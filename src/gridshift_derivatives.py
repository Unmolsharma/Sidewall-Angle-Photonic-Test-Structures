"""
Option 2: grid-shift-averaged finite-difference derivatives.

The width scan (results_anisotropic/logs/diag_width_scan.log) showed that on
EMpy's grid the slope dn_eff/dw depends on where the wall sits inside its cell
(periodic in the wall position, period = one cell). A central difference
therefore carries an error that is periodic in the offset between the
geometry and the grid. Averaging the same derivative over N grid offsets
spread evenly over one cell cancels the first N-1 harmonics of that error.

Shifts: the whole grid (X, Y) is translated by (k/N)*(dx_core, dy_core),
k = 0..N-1 (diagonal shifts, so both the sidewall x-phase and the top/bottom
y-phase are sampled). The domain moves with it by at most one core cell,
which changes n_eff by far less than the derivatives resolved here (domain
check in validate_anisotropic.py).

Derivative definitions, steps (h/2, h, 2h), n_g definition and the step
criterion are those of recompute_sensitivities.py / aniso_derivatives.py.
A unit is (averaging, geometry, pol, width_nm, slab_nm or None, mesh, N).
"""
import os
os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('MKL_NUM_THREADS', '1')
import pickle
import time
from multiprocessing import Pool

import numpy as np

import fullvector_sensitivity as fv
import slot_waveguide as sw
import anisotropic_averaging as an
from aniso_derivatives import MESHES, WLS, all_geometries, h_geometries
from recompute_sensitivities import H, FACTORS, DLAM, step_rows, as_s

PROD = dict(strip='30/18', swg='30/18', rib='30/18', slot='50/20')


def core_steps(geom, mesh):
    hx, hy, nc, ncl = MESHES[mesh]
    span_x = 1.0 if geom == 'slot' else 0.8
    return span_x / nc, 0.6 / nc


def grid(geom, mesh, k, N):
    hx, hy, nc, ncl = MESHES[mesh]
    X, Y = (sw.make_slot_grid if geom == 'slot' else fv.make_grid)(hx, hy, nc, ncl)
    dx, dy = core_steps(geom, mesh)
    return X + dx * k / N, Y + dy * k / N


def tasks_for(unit, geoms):
    avg, geom, pol, W, slab, mesh, N = unit
    geoms = list(geoms)  # may be a generator; it is reused for every shift
    return [(avg, geom, pol, W, slab, mesh, k, N, g, wl) for k in range(N) for g in geoms for wl in WLS]


def solve(task):
    avg, geom, pol, W, slab, mesh, k, N, (p, d), wl = task
    w, t, a = W / 1000, 0.22, 0.0
    if p == 'w':
        w += d / 1000
    elif p == 't':
        t += d / 1000
    elif p == 'a':
        a += d
    X, Y = grid(geom, mesh, k, N)
    kw = dict(X=X, Y=Y)
    if slab is not None:
        kw.update(slab_h_um=slab / 1000)
    return task, an.SOLVERS[avg][geom](w, t, a, wl, pol, **kw)


def run(tasks, procs=4, label='', cache=None):
    """Solve all tasks; with cache (a .pkl path) every finished solve is kept on
    disk and reused, so a crash in post-processing never costs the solves."""
    t0 = time.time()
    tasks = list(dict.fromkeys(tasks))
    res = {}
    if cache and os.path.exists(cache):
        with open(cache, 'rb') as f:
            res = {k: v for k, v in pickle.load(f).items()}
    todo = [t for t in tasks if t not in res]
    print(f"  {label}: {len(tasks)} tasks, {len(tasks) - len(todo)} cached", flush=True)
    with Pool(procs) as pool:
        for i, (task, n) in enumerate(pool.imap_unordered(solve, todo, chunksize=1)):
            res[task] = n
            if (i + 1) % 200 == 0 or i + 1 == len(todo):
                print(f"  {label} {i+1}/{len(todo)} solves  {time.time()-t0:.0f}s", flush=True)
                if cache:
                    with open(cache, 'wb') as f:
                        pickle.dump(res, f)
    return res


def _nng(res, unit, k, g):
    n0, npl, nmi = (res[unit[:6] + (k, unit[6], g, wl)] for wl in WLS)
    return n0, n0 - fv.WL0 * (npl - nmi) / (2 * DLAM)


def derivs(res, unit, factors=FACTORS, params=('w', 't', 'a'), per_shift=False):
    """Shift-averaged D[(p, obs, f)]; with per_shift also the list of per-shift values."""
    N = unit[6]
    D, Dk = {}, {}
    for p in params:
        for f in factors:
            h = f * H[p]
            vals_n, vals_g = [], []
            for k in range(N):
                (n1, g1), (n2, g2) = _nng(res, unit, k, (p, h)), _nng(res, unit, k, (p, -h))
                vals_n.append((n1 - n2) / (2 * h))
                vals_g.append((g1 - g2) / (2 * h))
            D[(p, 'neff', f)], D[(p, 'ng', f)] = float(np.mean(vals_n)), float(np.mean(vals_g))
            Dk[(p, 'neff', f)], Dk[(p, 'ng', f)] = vals_n, vals_g
    return (D, Dk) if per_shift else D


def full(res, unit):
    N = unit[6]
    nom = [_nng(res, unit, k, ('nom', 0.0)) for k in range(N)]
    n0, g0 = float(np.mean([x[0] for x in nom])), float(np.mean([x[1] for x in nom]))
    D = derivs(res, unit)
    return n0, g0, D, as_s(n0, g0, D)


def unit_id(unit):
    avg, geom, pol, W, slab, mesh, N = unit
    s = f'{geom}_{pol}_{W}' if slab is None else f'rib_slab{slab}_TE_{W}'
    return s
