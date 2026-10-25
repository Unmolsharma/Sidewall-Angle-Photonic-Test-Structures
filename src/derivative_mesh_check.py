"""
Mesh convergence of the DERIVATIVES (not just n_eff) with cell-exact
averaging, at step h. The step-size check cannot see a symmetric kink in
n_eff(param) at a cell edge (a central difference across a kink is constant
in h), so this checks whether D(h) itself is mesh-converged.

Structures: strip TE 500, strip TM 500, rib TE 400/500 (90 nm slab),
slot TE 180, slot TM 180. n_eff derivatives only (w, t, a at h = 5 nm, 2 nm, 1 deg).
Strip/rib meshes: production 30/18, 45/18, 60/24, 90/36 (same domain).
Slot meshes: production 50/20, 70/28, 100/40 (same domain).

Run from src/:  python derivative_mesh_check.py
"""
import os
os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
import json
from multiprocessing import Pool

import fullvector_sensitivity as fv
import slot_waveguide as sw
import cellexact_averaging as ce

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
H = dict(w=5.0, t=2.0, a=1.0)
STRIP_MESHES = {'30/18': (1.4, 1.2, 30, 18), '45/18': (1.4, 1.2, 45, 18), '60/24': (1.4, 1.2, 60, 24),
                '90/36': (1.4, 1.2, 90, 36)}
SLOT_MESHES = {'50/20': (1.6, 1.2, 50, 20), '70/28': (1.6, 1.2, 70, 28), '100/40': (1.6, 1.2, 100, 40)}
STRUCTS = [('strip', 'TE', 500), ('strip', 'TM', 500), ('rib', 'TE', 400), ('rib', 'TE', 500),
           ('slot', 'TE', 180), ('slot', 'TM', 180)]


def solve(task):
    geom, pol, W, mesh, p, d = task
    if geom == 'slot':
        X, Y = sw.make_slot_grid(*SLOT_MESHES[mesh])
    else:
        X, Y = fv.make_grid(*STRIP_MESHES[mesh])
    w, t, a = W / 1000 + (d / 1000 if p == 'w' else 0), 0.22 + (d / 1000 if p == 't' else 0), (d if p == 'a' else 0.0)
    return task, ce.SOLVERS['cell_exact'][geom](w, t, a, fv.WL0, pol, X=X, Y=Y)


if __name__ == '__main__':
    tasks = [(g, pol, W, m, p, s * H[p]) for g, pol, W in STRUCTS
             for m in (SLOT_MESHES if g == 'slot' else STRIP_MESHES)
             for p in ('w', 't', 'a') for s in (1, -1)]
    with Pool(10) as pool:
        res = dict(pool.map(solve, tasks, chunksize=1))
    out = {}
    for g, pol, W in STRUCTS:
        meshes = SLOT_MESHES if g == 'slot' else STRIP_MESHES
        for p in ('w', 't', 'a'):
            D = {m: (res[(g, pol, W, m, p, H[p])] - res[(g, pol, W, m, p, -H[p])]) / (2 * H[p]) for m in meshes}
            vals = list(D.values())
            finest = vals[-1]
            out[f'{g}_{pol}_{W}|dneff/d{p}'] = dict(D=D, production_vs_finest_pct=100 * (vals[0] - finest) / abs(finest))
            print(f"{g}_{pol}_{W:<4} dneff/d{p}: " + '  '.join(f"{m}: {v:+.4e}" for m, v in D.items()) +
                  f"   production vs finest {100*(vals[0]-finest)/abs(finest):+.1f}%", flush=True)
    with open(os.path.join(REPO, 'results_corrected', 'data', 'derivative_mesh_check.json'), 'w') as f:
        json.dump(out, f, indent=1)
