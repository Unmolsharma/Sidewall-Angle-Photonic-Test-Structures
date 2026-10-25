"""
Derivative engine for the anisotropic-averaging runs (Phase 2 and its pilot).

Same derivative definitions and stability criterion as
recompute_sensitivities.py (central differences at h/2, h, 2h with
h = 5 nm / 2 nm / 1 deg; n_g = n_eff - lambda*dn_eff/dlambda, dlambda = 5 nm;
2% relative, or 2% of the largest |D| for near-zero entries), but every solve
can also take an explicit mesh, so production and finer meshes, and the three
averagings, run through one code path.

A unit is (averaging, geometry, pol, width_nm, slab_nm or None, mesh or None);
mesh is None for the production grid of that geometry, else a key of MESHES.
"""
import os
os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('MKL_NUM_THREADS', '1')
import time
from multiprocessing import Pool

import fullvector_sensitivity as fv
import slot_waveguide as sw
import anisotropic_averaging as an
from recompute_sensitivities import H, FACTORS, DLAM, step_rows, as_s

# (half_x, half_y, n_core, n_clad); strip/rib/swg use fv.make_grid, slot uses sw.make_slot_grid
MESHES = {'30/18': (1.4, 1.2, 30, 18), '45/18': (1.4, 1.2, 45, 18), '60/24': (1.4, 1.2, 60, 24),
          '90/36': (1.4, 1.2, 90, 36),
          '50/20': (1.6, 1.2, 50, 20), '70/28': (1.6, 1.2, 70, 28), '100/40': (1.6, 1.2, 100, 40)}
WLS = (fv.WL0, fv.WL0 + DLAM, fv.WL0 - DLAM)


def all_geometries():
    yield ('nom', 0.0)
    for p in ('w', 't', 'a'):
        for f in FACTORS:
            for sgn in (+1, -1):
                yield (p, sgn * f * H[p])


def h_geometries(params):
    """Only the +-h perturbations of the given parameters (mesh checks)."""
    return [(p, s * H[p]) for p in params for s in (1, -1)]


def tasks_for(unit, geoms):
    return [unit + (g, wl) for g in geoms for wl in WLS]


def solve(task):
    avg, geom, pol, W, slab, mesh, (p, d), wl = task
    w, t, a = W / 1000, 0.22, 0.0
    if p == 'w':
        w += d / 1000
    elif p == 't':
        t += d / 1000
    elif p == 'a':
        a += d
    kw = {}
    if mesh is not None:
        X, Y = (sw.make_slot_grid if geom == 'slot' else fv.make_grid)(*MESHES[mesh])
        kw.update(X=X, Y=Y)
        if avg == 'original':
            raise ValueError('explicit meshes are only wired for cell_exact / anisotropic')
    if slab is not None:
        kw.update(slab_h_um=slab / 1000)
    return task, an.SOLVERS[avg][geom](w, t, a, wl, pol, **kw)


def run(tasks, procs=4, label=''):  # Pool(4) was fastest on this laptop (results_anisotropic/data/timing.json)
    t0 = time.time()
    tasks = list(dict.fromkeys(tasks))
    res = {}
    with Pool(procs) as pool:
        for i, (task, n) in enumerate(pool.imap_unordered(solve, tasks, chunksize=1)):
            res[task] = n
            if (i + 1) % 200 == 0 or i + 1 == len(tasks):
                print(f"  {label} {i+1}/{len(tasks)} solves  {time.time()-t0:.0f}s", flush=True)
    return res


def nng(res, unit, g):
    n0, npl, nmi = (res[unit + (g, wl)] for wl in WLS)
    return n0, n0 - fv.WL0 * (npl - nmi) / (2 * DLAM)


def derivs(res, unit, factors=FACTORS, params=('w', 't', 'a')):
    D = {}
    for p in params:
        for f in factors:
            h = f * H[p]
            (n1, g1), (n2, g2) = nng(res, unit, (p, h)), nng(res, unit, (p, -h))
            D[(p, 'neff', f)] = (n1 - n2) / (2 * h)
            D[(p, 'ng', f)] = (g1 - g2) / (2 * h)
    return D


def full(res, unit):
    """n_eff, n_g, the step-check rows and the as_s() record at step h."""
    n0, g0 = nng(res, unit, ('nom', 0.0))
    D = derivs(res, unit)
    return n0, g0, D, as_s(n0, g0, D)


def unit_id(unit):
    avg, geom, pol, W, slab, mesh = unit
    s = f'{geom}_{pol}_{W}' if slab is None else f'rib_slab{slab}_TE_{W}'
    return s if mesh is None else f'{s}@{mesh}'


__all__ = ['MESHES', 'WLS', 'all_geometries', 'h_geometries', 'tasks_for', 'run', 'nng', 'derivs', 'full',
           'step_rows', 'unit_id']
