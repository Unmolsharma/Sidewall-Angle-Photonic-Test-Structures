"""
Rib derivatives on a fixed mesh: the nominal mesh is reused for every point of the finite-difference
stencil, with its nodes moved onto the perturbed geometry, so there is no remeshing noise between
stencil points (the regular batch remeshes each perturbed geometry).

Runs in the femwell virtual environment, from src/:
    C:/Users/unmol/.venvs/femwell/Scripts/python.exe fixedmesh_rib.py [max_workers] [battery_floor]

Node map (nominal: bottom width w0, rib top t0 = 220 nm, slab s0, vertical walls). Every material boundary
is a horizontal line or a rib wall, so a continuous piecewise-linear map that sends the corners to their
new places keeps the mesh conforming:
  1. vertical: y -> piecewise linear through (0, 0), (s0, s), (t0, t), (t0 + 0.3, t0 + 0.3) um; identity
     below 0 and above t0 + 0.3 um (moves the slab top and rib top, i.e. the 's' and 't' steps);
  2. horizontal, at the mapped y: inside the rib |x| <= w0/2 is scaled onto the new wall
     x_w(y) = w/2 - tan(a)(y - s); between the wall and |x| = w0/2 + 0.4 um the gap is stretched
     linearly; the shift is tapered linearly to zero over 0.2 um below the slab top and above the rib top.
The element permittivities are those of the nominal mesh (each element keeps its material).
Steps and wavelengths as in femwell_batch.py (h = 5 nm w, 2 nm t, 1 deg a, 2 nm slab; h/2, h, 2h).

Writes results_montecarlo/fixedmesh/cache/solves_r12.jsonl (keys 'fixedmesh|' + femwell_batch key) and
logs/fixedmesh_r12.log.
"""
import dataclasses
import json
import os
import sys
import time

import numpy as np
from skfem import Basis, ElementTriP0

import femwell_batch as fb
from femwell.maxwell.waveguide import compute_modes

OUT = os.path.join(fb.REPO, 'results_montecarlo', 'fixedmesh')
CACHE = os.path.join(OUT, 'cache', 'solves_r12.jsonl')
UNITS = [('rib', 'TE', 300, 50), ('rib', 'TE', 300, 90)]
SRC_CACHE = os.path.join(fb.REPO, 'results_femwell', 'cache', 'solves_r12.jsonl')


def ramp(y, knots):
    xs, ys = zip(*knots)
    return np.interp(y, xs, ys)


def morph(p, w0, t0, s0, w, t, a, s):
    """Node coordinates (um) of the nominal mesh moved onto the geometry (w, t, a, s)."""
    x, y = p
    top = t0 + 0.3
    yn = np.where((y > 0) & (y < top), ramp(y, [(0, 0), (s0, s), (t0, t), (top, top)]), y)
    xw0, X0 = w0 / 2, w0 / 2 + 0.4
    xw = w / 2 - np.tan(np.radians(a)) * (np.clip(yn, s, t) - s)
    ax = np.abs(x)
    xn = np.where(ax <= xw0, ax * xw / xw0,
                  np.where(ax <= X0, xw + (ax - xw0) * (X0 - xw) / (X0 - xw0), ax))
    g = np.clip(np.minimum((yn - (s - 0.2)) / 0.2, ((t + 0.2) - yn) / 0.2), 0, 1)
    return np.vstack([np.sign(x) * (ax + g * (xn - ax)), yn])


def solve(task):
    geom, pol, W, slab_nm, clad, mesh, p, d, guess = task
    t0 = time.time()
    w0, tt0, s0 = W / 1000, 0.22, slab_nm / 1000
    w, t, a, s = w0, tt0, 0.0, s0
    if p == 'w':
        w += d / 1000
    elif p == 't':
        t += d / 1000
    elif p == 'a':
        a += d
    elif p == 's':
        s += d / 1000
    basis0, eps = fb.build(geom, w0, tt0, 0.0, s0, clad, mesh)
    m = basis0.mesh
    basis1 = Basis(dataclasses.replace(m, doflocs=morph(m.p, w0, tt0, s0, w, t, a, s)), ElementTriP0())
    ns, tes, notes = [], [], []
    for wl in fb.WLS:
        md = fb.pick(compute_modes(basis1, eps, wavelength=wl, num_modes=2, order=2, n_guess=guess), pol)
        if md is None or abs(md.n_eff.real - guess) > 0.03 * guess:
            notes.append(f'wl {wl}: fallback to 6 modes without guess')
            md = fb.pick(compute_modes(basis1, eps, wavelength=wl, num_modes=6, order=2), pol)
        ns.append(float(md.n_eff.real))
        tes.append(float(md.te_fraction))
    return dict(key=key_of(task), n=ns, te_fraction=tes, n_elements=int(m.nelements),
                seconds=time.time() - t0, notes=notes)


def key_of(task):
    return 'fixedmesh|' + fb.key_of(task)


def tasks():
    src = fb.load_cache(SRC_CACHE)
    out = []
    for geom, pol, W, S in UNITS:
        guess = src[f'{geom}|{pol}|{W}|{S}|oxide|r12|nom|+0.0000']['n'][0]
        out += fb.unit_tasks(geom, pol, W, S, 'oxide', 'r12', guess, params=['w', 't', 'a', 's'])
    return out


if __name__ == '__main__':
    os.makedirs(os.path.join(OUT, 'cache'), exist_ok=True)
    os.makedirs(os.path.join(OUT, 'logs'), exist_ok=True)
    cap = int(sys.argv[1]) if len(sys.argv) > 1 else 2
    if len(sys.argv) > 2:
        fb.BATTERY_FLOOR = int(sys.argv[2])
    logpath = os.path.join(OUT, 'logs', 'fixedmesh_r12.log')
    failed = fb.run(tasks(), CACHE, logpath, fb.choose_workers(2.0, max_workers=cap, logpath=logpath), 2.0,
                    solver=solve, keyer=key_of)
    with open(os.path.join(OUT, 'logs', 'fixedmesh_r12_failed.json'), 'w') as f:
        json.dump(failed, f, indent=1)
