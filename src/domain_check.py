"""
Domain-size convergence check: the structures of the two-width headline set (strip TE 400-600 nm, rib TE 300 and
500 nm on a 50 nm slab) solved on a larger computational domain, x in [-2.1, 2.1] um and y in [-1.8, 1.8] um
(1.5x the standard [-1.4, 1.4] x [-1.2, 1.2] um), with the same mesh settings (r12), wavelengths and mode selection.
Nominal + geometric steps at h (w, t, a; slab for ribs). Compared with the standard-domain solves in
validation_checks.py.

The geometry builder is femwell_batch.build with the domain half-sizes as parameters (copied here so that the
production code is unchanged).

Runs in the femwell virtual environment, from src/:
    C:/Users/unmol/.venvs/femwell/Scripts/python.exe domain_check.py [max_workers]
Writes results_final/domain_check/cache/solves_r12.jsonl (keys 'domain|1.5|' + femwell_batch key), logs/domain_check.log
"""
import json
import os
import sys
import time
from collections import OrderedDict

from shapely.geometry import box
from skfem import Basis, ElementTriP0
from skfem.io.meshio import from_meshio

import femwell_batch as fb
from femwell.maxwell.waveguide import compute_modes
from femwell.mesh import mesh_from_OrderedDict

OUT = os.path.join(fb.REPO, 'results_final', 'domain_check')
CACHE = os.path.join(OUT, 'cache', 'solves_r12.jsonl')
SCALE = 1.5
HX, HY = 1.4 * SCALE, 1.2 * SCALE
UNITS = [('strip', 'TE', w, None) for w in (400, 450, 500, 550, 600)] + [('rib', 'TE', w, 50) for w in (300, 500)]
SRC = [os.path.join(fb.REPO, 'results_femwell', 'cache', 'solves_r12.jsonl')]


def build(geom, w, t, a, slab, mesh):
    dom = box(-HX, -HY, HX, HY)
    core = fb.shapes(geom, w, t, a, slab)[0].intersection(dom)
    polys = OrderedDict(core0=core)
    polys['box'] = box(-HX, -HY, HX, 0).difference(core)
    polys['clad'] = box(-HX, 0, HX, HY).difference(core)
    res, dist, rmax = fb.MESH[mesh]
    m = from_meshio(mesh_from_OrderedDict(polys, {'core0': dict(resolution=res, distance=dist)}, default_resolution_max=rmax))
    basis0 = Basis(m, ElementTriP0())
    eps = basis0.zeros()
    eps[basis0.get_dofs(elements='core0')] = fb.N_SI ** 2
    eps[basis0.get_dofs(elements='box')] = fb.N_OX ** 2
    eps[basis0.get_dofs(elements='clad')] = fb.N_OX ** 2
    return basis0, eps


def solve(task):
    geom, pol, W, slab_nm, clad, mesh, p, d, guess = task
    t0 = time.time()
    w, t, a = W / 1000, 0.22, 0.0
    slab = (slab_nm if slab_nm is not None else fb.SLAB0 * 1000) / 1000
    if p == 'w':
        w += d / 1000
    elif p == 't':
        t += d / 1000
    elif p == 'a':
        a += d
    elif p == 's':
        slab += d / 1000
    basis0, eps = build(geom, w, t, a, slab, mesh)
    ns, tes, notes = [], [], []
    for wl in fb.WLS:
        m = fb.pick(compute_modes(basis0, eps, wavelength=wl, num_modes=2, order=2, n_guess=guess), pol)
        if m is None or abs(m.n_eff.real - guess) > 0.03 * guess:
            notes.append(f'wl {wl}: fallback to 6 modes without guess')
            m = fb.pick(compute_modes(basis0, eps, wavelength=wl, num_modes=6, order=2), pol)
        ns.append(float(m.n_eff.real))
        tes.append(float(m.te_fraction))
    return dict(key=key_of(task), n=ns, te_fraction=tes, n_elements=int(basis0.mesh.nelements),
                seconds=time.time() - t0, notes=notes, domain_um=[HX, HY])


def key_of(task):
    return f'domain|{SCALE}|' + fb.key_of(task)


def tasks():
    rows = fb.load_cache(SRC[0])
    out = []
    for geom, pol, W, S in UNITS:
        guess = rows[f'{geom}|{pol}|{W}|{S}|oxide|r12|nom|+0.0000']['n'][0]
        out.append((geom, pol, W, S, 'oxide', 'r12', 'nom', 0.0, guess))
        for p in ('w', 't', 'a') + (('s',) if geom == 'rib' else ()):
            for sgn in (1, -1):
                out.append((geom, pol, W, S, 'oxide', 'r12', p, sgn * fb.H[p], guess))
    return out


if __name__ == '__main__':
    os.makedirs(os.path.join(OUT, 'cache'), exist_ok=True)
    os.makedirs(os.path.join(OUT, 'logs'), exist_ok=True)
    cap = int(sys.argv[1]) if len(sys.argv) > 1 else 2
    logpath = os.path.join(OUT, 'logs', 'domain_check.log')
    failed = fb.run(tasks(), CACHE, logpath, max(cap, fb.choose_workers(2.5, max_workers=cap, logpath=logpath)), 2.5,
                    solver=solve, keyer=key_of)
    with open(os.path.join(OUT, 'logs', 'domain_check_failed.json'), 'w') as f:
        json.dump(failed, f, indent=1)
