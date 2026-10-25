"""
Phase 2 pilot: anisotropic vs cell-exact derivatives for
  strip TE 500, rib TE 500 (90 nm slab), rib TE 300 (90 nm slab), slot TE 180 (middle width)
  (1) step-size check (h/2, h, 2h) on the production mesh, all 12 derivatives per structure
  (2) mesh check of the width and angle derivatives (n_eff and n_g, step h):
      strip/rib meshes 30/18 (production), 45/18, 60/24, 90/36; slot 50/20 (production), 70/28, 100/40
Both averagings are computed here through the same code path.

Run from src/:  python pilot_anisotropic.py
Writes results_anisotropic/data/pilot.json
"""
import json
import os
import time

import aniso_derivatives as ad

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, 'results_anisotropic', 'data')
STRUCTS = [('strip', 'TE', 500, None), ('rib', 'TE', 500, None), ('rib', 'TE', 300, 90), ('slot', 'TE', 180, None)]
FINE = dict(strip=['45/18', '60/24', '90/36'], rib=['45/18', '60/24', '90/36'], slot=['70/28', '100/40'])
PROD = dict(strip='30/18', rib='30/18', slot='50/20')
AVGS = ('cell_exact', 'anisotropic')

if __name__ == '__main__':
    t0 = time.time()
    tasks = []
    for avg in AVGS:
        for g, pol, W, slab in STRUCTS:
            tasks += ad.tasks_for((avg, g, pol, W, slab, None), ad.all_geometries())
            for m in FINE[g]:
                tasks += ad.tasks_for((avg, g, pol, W, slab, m), ad.h_geometries(('w', 'a')))
    res = ad.run(tasks, label='pilot')

    out = dict(step_check={}, mesh_check={}, summary={})
    for avg in AVGS:
        rows, mesh = [], {}
        for g, pol, W, slab in STRUCTS:
            u = (avg, g, pol, W, slab, None)
            _, _, D, _ = ad.full(res, u)
            rows += ad.step_rows(avg, ad.unit_id(u), D)
            for p in ('w', 'a'):
                for obs in ('neff', 'ng'):
                    vals = {PROD[g]: D[(p, obs, 1.0)]}
                    for m in FINE[g]:
                        vals[m] = ad.derivs(res, (avg, g, pol, W, slab, m), factors=(1.0,), params=(p,))[(p, obs, 1.0)]
                    finest = list(vals.values())[-1]
                    mesh[f'{ad.unit_id(u)}|d{obs}/d{p}'] = dict(
                        D=vals, production_vs_finest_pct=100 * (vals[PROD[g]] - finest) / abs(finest),
                        max_vs_finest_pct=max(abs(100 * (v - finest) / finest) for v in list(vals.values())[:-1]))
        out['step_check'][avg] = rows
        out['mesh_check'][avg] = mesh
        out['summary'][avg] = dict(
            step_fail=[f"{r['structure']} d{r['output']}/d{r['param']}" for r in rows if r['verdict'] == 'FAIL'],
            n_pass=sum(r['verdict'] == 'PASS' for r in rows),
            n_near_zero=sum(r['verdict'] == 'PASS_NEAR_ZERO_ABS' for r in rows),
            n_fail=sum(r['verdict'] == 'FAIL' for r in rows))

    for avg in AVGS:
        s = out['summary'][avg]
        print(f"\n=== {avg}: step check {s['n_pass']} pass, {s['n_near_zero']} near-zero pass, {s['n_fail']} fail / {len(out['step_check'][avg])}")
        for r in out['step_check'][avg]:
            print(f"  {r['verdict']:18s} {r['structure']:20s} d{r['output']:4s}/d{r['param']}  "
                  f"D(h/2,h,2h)=({r['D_half']:+.4e}, {r['D_h']:+.4e}, {r['D_2h']:+.4e})  rel dev {r['max_rel_dev']:6.1%}")
    print("\n=== mesh check (step h): production vs finest, % ===")
    print(f"  {'derivative':34s} {'cell-exact':>11s} {'anisotropic':>12s}   anisotropic values per mesh")
    for k in out['mesh_check']['anisotropic']:
        ce_, an_ = out['mesh_check']['cell_exact'][k], out['mesh_check']['anisotropic'][k]
        print(f"  {k:34s} {ce_['production_vs_finest_pct']:+10.1f}% {an_['production_vs_finest_pct']:+11.1f}%   "
              + '  '.join(f"{m}: {v:+.4e}" for m, v in an_['D'].items()))
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, 'pilot.json'), 'w') as f:
        json.dump(out, f, indent=1, default=float)
    print(f"\nElapsed {time.time()-t0:.0f}s")
