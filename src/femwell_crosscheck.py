"""
Independent reference solver: femwell (finite elements on a conforming
triangular mesh; every material boundary is a mesh edge, so no sub-pixel
averaging is involved).

Runs in its own virtual environment (femwell pulls newer numpy/scipy than the
EMpy environment uses):
    C:/Users/unmol/.venvs/femwell/Scripts/python.exe src/femwell_crosscheck.py convergence
    C:/Users/unmol/.venvs/femwell/Scripts/python.exe src/femwell_crosscheck.py derivatives
It imports nothing from the EMpy code; the constants are copied from
fullvector_sensitivity.py / rib_waveguide.py.

Geometry conventions identical to the FD code: sidewall angle a = deviation
from vertical in degrees, top narrower, the tilt pivots about the BOTTOM edge
of the tilted wall (y = 0 for strip, y = slab top for rib), so the wall moves
inward by tan(a)*(y - y_pivot); width = bottom width of that wall; rib
thickness t perturbed with the slab held fixed. Optional air top cladding
(oxide below y = 0) for the cladding check.

Modes: one mode near n_guess (shift-invert), order-2 elements; the solve is
accepted only if it is TE-like (te_fraction > 0.5) and is checked against a
4-mode solve at the nominal geometry.

  convergence : n_eff and n_g vs mesh resolution and domain size
  derivatives : n_eff, n_g and their w/t/a derivatives at h/2, h, 2h
                (h = 5 nm, 2 nm, 1 deg; every geometry remeshed), step check

Writes results_anisotropic/data/femwell_{convergence,derivatives}.json
"""
import json
import os
import sys
import time
from collections import OrderedDict
from multiprocessing import Pool

import numpy as np
from shapely.geometry import Polygon, box
from skfem import Basis, ElementTriP0
from skfem.io.meshio import from_meshio

from femwell.maxwell.waveguide import compute_modes
from femwell.mesh import mesh_from_OrderedDict

N_SI, N_OX, N_AIR, WL0 = 3.48, 1.44, 1.0, 1.55
SLAB = 0.09
H = dict(w=5.0, t=2.0, a=1.0)
FACTORS = (0.5, 1.0, 2.0)
DLAM = 0.005
WLS = (WL0, WL0 + DLAM, WL0 - DLAM)
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, 'results_anisotropic', 'data')
# mesh: resolution at the silicon (um), distance over which it relaxes, max element size elsewhere
MESH = {'r20': (0.020, 0.3, 0.2), 'r12': (0.012, 0.3, 0.2), 'r8': (0.008, 0.3, 0.2),
        'r6': (0.006, 0.3, 0.2), 'r4': (0.004, 0.3, 0.2), 'r12_max0.1': (0.012, 0.3, 0.1)}
N_GUESS = dict(strip=2.45, rib=2.57)


def core_polygon(geom, w, t, a):
    tan_a = np.tan(np.radians(a))
    if geom == 'strip':
        hw0, hwt = w / 2, w / 2 - tan_a * t
        return Polygon([(-hw0, 0), (hw0, 0), (hwt, t), (-hwt, t)])
    hw0, hwt = w / 2, w / 2 - tan_a * (t - SLAB)
    return Polygon([(-50, 0), (50, 0), (50, SLAB), (hw0, SLAB), (hwt, t), (-hwt, t), (-hw0, SLAB), (-50, SLAB)])


def build(geom, w, t, a, mesh, half, clad):
    hx, hy = half
    core = core_polygon(geom, w, t, a).intersection(box(-hx, -hy, hx, hy))
    polys = OrderedDict(core=core, box=box(-hx, -hy, hx, 0).difference(core),
                        clad=box(-hx, 0, hx, hy).difference(core))
    res, dist, rmax = MESH[mesh]
    m = from_meshio(mesh_from_OrderedDict(polys, {'core': dict(resolution=res, distance=dist)},
                                          default_resolution_max=rmax))
    basis0 = Basis(m, ElementTriP0())
    eps = basis0.zeros()
    for name, n in (('core', N_SI), ('box', N_OX), ('clad', N_OX if clad == 'oxide' else N_AIR)):
        eps[basis0.get_dofs(elements=name)] = n ** 2
    return basis0, eps


def solve_nng(task):
    """task = (geom, w, t, a, mesh, half, clad) -> n_eff at WL0, n_g, n_elements."""
    geom, w, t, a, mesh, half, clad = task
    basis0, eps = build(geom, w, t, a, mesh, half, clad)
    ns = []
    for wl in WLS:
        m = compute_modes(basis0, eps, wavelength=wl, num_modes=1, order=2,
                          n_guess=N_GUESS[geom] if clad == 'oxide' else N_GUESS[geom] - 0.15)[0]
        if m.te_fraction <= 0.5:
            raise RuntimeError(f'non-TE mode picked for {task} at {wl}: te_fraction {m.te_fraction}')
        ns.append(float(m.n_eff.real))
    n0, npl, nmi = ns
    return task, (n0, n0 - WL0 * (npl - nmi) / (2 * DLAM), basis0.mesh.nelements)


def check_mode_order(geom, mesh, half, clad):
    basis0, eps = build(geom, 0.5, 0.22, 0.0, mesh, half, clad)
    ms = compute_modes(basis0, eps, wavelength=WL0, num_modes=4, order=2)
    te = [m for m in ms if m.te_fraction > 0.5]
    return float(max(m.n_eff.real for m in te))


def convergence():
    tasks = []
    for geom in ('strip', 'rib'):
        for mesh in ('r20', 'r12', 'r8', 'r6', 'r4', 'r12_max0.1'):
            tasks.append((geom, 0.5, 0.22, 0.0, mesh, (1.4, 1.2), 'oxide'))
        for half in ((1.0, 0.9), (2.0, 1.8)):
            tasks.append((geom, 0.5, 0.22, 0.0, 'r8', half, 'oxide'))
    t0 = time.time()
    with Pool(4) as p:
        res = dict(p.map(solve_nng, tasks, chunksize=1))
    out = {}
    for task, (n, g, nel) in res.items():
        geom, _, _, _, mesh, half, clad = task
        key = f'{geom} {mesh} domain {half[0]}x{half[1]}'
        out[key] = dict(neff=n, ng=g, n_elements=nel)
        print(f"  {key:34s} {nel:7d} el  n_eff {n:.7f}  n_g {g:.6f}", flush=True)
    for geom in ('strip', 'rib'):
        out[f'{geom} 4-mode check (r8)'] = check_mode_order(geom, 'r8', (1.4, 1.2), 'oxide')
        print(f"  {geom}: highest TE n_eff from a 4-mode solve (r8) {out[f'{geom} 4-mode check (r8)']:.7f}")
    print(f"  ({time.time()-t0:.0f}s)")
    with open(os.path.join(OUT, 'femwell_convergence.json'), 'w') as f:
        json.dump(out, f, indent=1)


def perturbed(W, p, d):
    w, t, a = W / 1000, 0.22, 0.0
    if p == 'w':
        w += d / 1000
    elif p == 't':
        t += d / 1000
    elif p == 'a':
        a += d
    return w, t, a


def derivatives(cases, mesh, half):
    tasks = []
    for geom, W, clad, params in cases:
        tasks.append((geom,) + perturbed(W, 'nom', 0) + (mesh, half, clad))
        for p in params:
            for f in FACTORS:
                for s in (1, -1):
                    tasks.append((geom,) + perturbed(W, p, s * f * H[p]) + (mesh, half, clad))
    t0 = time.time()
    res = {}
    with Pool(4) as pool:
        for i, (task, v) in enumerate(pool.imap_unordered(solve_nng, tasks, chunksize=1)):
            res[task] = v
            print(f"  {i+1}/{len(tasks)}  {time.time()-t0:.0f}s", flush=True)
    out = {}
    for geom, W, clad, params in cases:
        key = f'{geom}_TE_{W}_{clad}'
        nom = res[(geom,) + perturbed(W, 'nom', 0) + (mesh, half, clad)]
        D = {}
        for p in params:
            for f in FACTORS:
                h = f * H[p]
                (n1, g1, _), (n2, g2, _) = (res[(geom,) + perturbed(W, p, s * h) + (mesh, half, clad)] for s in (1, -1))
                D[f'{p}|neff|{f}'] = (n1 - n2) / (2 * h)
                D[f'{p}|ng|{f}'] = (g1 - g2) / (2 * h)
        steps = {}
        for obs in ('neff', 'ng'):
            big = max(abs(D[f'{p}|{obs}|1.0']) for p in params)
            for p in params:
                dh, dhalf, d2h = (D[f'{p}|{obs}|{f}'] for f in (1.0, 0.5, 2.0))
                dev = max(abs(dhalf - dh), abs(d2h - dh))
                near_zero = abs(dh) < 0.1 * big
                ok = dev <= 0.02 * abs(dh) or (near_zero and dev <= 0.02 * big)
                steps[f'd{obs}/d{p}'] = dict(D_half=dhalf, D_h=dh, D_2h=d2h, max_rel_dev=dev / abs(dh),
                                             near_zero=near_zero, passed=bool(ok))
        R = None
        if 'w' in params and 'a' in params:
            hw = (0.22 - SLAB) if geom == 'rib' else 0.22
            k = -hw * 1000 * np.tan(np.radians(1.0))
            R = dict(R_neff=D['a|neff|1.0'] / (k * D['w|neff|1.0']), R_ng=D['a|ng|1.0'] / (k * D['w|ng|1.0']),
                     ratio_neff=D['a|neff|1.0'] / D['w|neff|1.0'], ratio_ng=D['a|ng|1.0'] / D['w|ng|1.0'])
        out[key] = dict(neff=nom[0], ng=nom[1], n_elements=nom[2], mesh=mesh, domain=list(half), D=D, step_check=steps, R=R)
        print(f"\n{key}: n_eff {nom[0]:.6f}  n_g {nom[1]:.5f}  ({nom[2]} elements, mesh {mesh}, domain {half})")
        for k2, s in steps.items():
            print(f"  {'PASS' if s['passed'] else 'FAIL'}  {k2:10s} D(h/2,h,2h)=({s['D_half']:+.5e}, {s['D_h']:+.5e}, {s['D_2h']:+.5e})"
                  f"  rel dev {100*s['max_rel_dev']:.2f}%{'  (near-zero)' if s['near_zero'] else ''}")
        if R:
            print(f"  (dn/da)/(dn/dw): n_eff {R['ratio_neff']:+.3f} nm/deg, n_g {R['ratio_ng']:+.3f} nm/deg;  "
                  f"R = (dn/da)/(-h tan1deg dn/dw): n_eff {R['R_neff']:.4f}, n_g {R['R_ng']:.4f}")
    print(f"({time.time()-t0:.0f}s)")
    return out


if __name__ == '__main__':
    os.makedirs(OUT, exist_ok=True)
    what = sys.argv[1]
    if what == 'convergence':
        convergence()
    elif what == 'derivatives':
        mesh = sys.argv[2]
        half = tuple(float(x) for x in sys.argv[3].split('x'))
        cases = [('strip', 500, 'oxide', 'wta'), ('rib', 500, 'oxide', 'wta')]
        if len(sys.argv) > 4 and sys.argv[4] == '+air':
            cases.append(('strip', 500, 'air', 'wa'))
        out = derivatives(cases, mesh, half)
        path = os.path.join(OUT, 'femwell_derivatives.json')
        prev = json.load(open(path)) if os.path.exists(path) else {}
        prev.update(out)
        with open(path, 'w') as f:
            json.dump(prev, f, indent=1)
