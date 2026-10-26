"""
Validate option 2 (grid-shift-averaged FD derivatives) for strip TE 500 and
rib TE 500 (90 nm slab), cell-exact and anisotropic averaging:
  - production mesh, N = 8 shifts (k/8 of a core cell); the N = 4, 2, 1
    averages are the subsets k even, k in {0, 4}, k = 0
  - all 12 derivatives at h/2, h, 2h -> step-size check
  - mesh check: the same shift-averaged derivatives (step h) on the 45/18 mesh
Compared with femwell in compare_derivatives.py.

Run from src/:  python gridshift_validate.py
Writes results_anisotropic/data/gridshift_validation.json
"""
import json
import os
import time

import numpy as np

import gridshift_derivatives as gd
from aniso_derivatives import all_geometries, h_geometries
from recompute_sensitivities import H, FACTORS

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, 'results_anisotropic', 'data')
STRUCTS = [('strip', 'TE', 500, None), ('rib', 'TE', 500, None)]
AVGS = ('cell_exact', 'anisotropic')
NMAX = 8
SUBSETS = {8: range(8), 4: range(0, 8, 2), 2: (0, 4), 1: (0,)}

if __name__ == '__main__':
    t0 = time.time()
    tasks = []
    for avg in AVGS:
        for g, pol, W, slab in STRUCTS:
            tasks += gd.tasks_for((avg, g, pol, W, slab, '30/18', NMAX), all_geometries())
            tasks += gd.tasks_for((avg, g, pol, W, slab, '45/18', NMAX), [('nom', 0.0)] + h_geometries(('w', 't', 'a')))
    res = gd.run(tasks, label='gridshift', cache=os.path.join(OUT, 'gridshift_validation_solves.pkl'))
    out = {}
    for avg in AVGS:
        for g, pol, W, slab in STRUCTS:
            key = f'{avg}|{g}_{pol}_{W}'
            o = out[key] = {}
            for mesh, factors in (('30/18', FACTORS), ('45/18', (1.0,))):
                u = (avg, g, pol, W, slab, mesh, NMAX)
                _, Dk = gd.derivs(res, u, factors=factors, per_shift=True)
                o[mesh] = {}
                for n, ks in SUBSETS.items():
                    o[mesh][f'N{n}'] = {f'{p}|{obs}|{f}': float(np.mean([v[k] for k in ks]))
                                        for (p, obs, f), v in Dk.items()}
                o[mesh]['per_shift'] = {f'{p}|{obs}|{f}': v for (p, obs, f), v in Dk.items()}
                nom = [gd._nng(res, u, k, ('nom', 0.0)) for k in range(NMAX)]
                o[mesh]['neff_per_shift'] = [x[0] for x in nom]
                o[mesh]['ng_per_shift'] = [x[1] for x in nom]
            print(f"\n{key}: n_eff over shifts {min(o['30/18']['neff_per_shift']):.5f}..{max(o['30/18']['neff_per_shift']):.5f}")
            for p in ('w', 't', 'a'):
                for obs in ('neff', 'ng'):
                    k1 = f'{p}|{obs}|1.0'
                    per = o['30/18']['per_shift'][k1]
                    print(f"  d{obs:4s}/d{p}  per-shift spread {100*(max(per)-min(per))/abs(np.mean(per)):5.1f}%   "
                          + '  '.join(f"N{n}: {o['30/18'][f'N{n}'][k1]:+.4e}" for n in (1, 2, 4, 8))
                          + f"   45/18 N8: {o['45/18']['N8'][k1]:+.4e}")
    with open(os.path.join(OUT, 'gridshift_validation.json'), 'w') as f:
        json.dump(out, f, indent=1)
    print(f"\nElapsed {time.time()-t0:.0f}s")
