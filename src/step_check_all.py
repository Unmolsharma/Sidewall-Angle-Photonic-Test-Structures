"""
Step check (h/2, h, 2h; same rule as etchdepth_femwell_analysis.step_check) over every structure in the
femwell caches, including the strip TM index solves and the rib 300/500 nm solves at 50/110/130/150 nm
slabs added on 2026-10-28. Raw derivatives (no strip symmetry), since the check is about the solves.

Run from src/ (main environment):  python step_check_all.py
Writes  results_montecarlo/sym_search/data/step_check_all.csv
"""
import collections
import csv
import os

import etchdepth_femwell_analysis as ed
import priors_search as ps

OUT = os.path.join(ed.REPO, 'results_montecarlo', 'sym_search', 'data', 'step_check_all.csv')

if __name__ == '__main__':
    steps = ed.step_check(ps.derivatives(ed.load(), strip_symmetry=False))   # every structure with solves
    with open(OUT, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(steps[0].keys()))
        w.writeheader()
        w.writerows(steps)
    print(f'{len(steps)} derivatives:', dict(collections.Counter(r['verdict'] for r in steps)))
    for r in steps:
        if r['verdict'] == 'FAIL':
            print(f"FAIL {r['structure']} {r['param']} {r['output']}: max rel dev {r['max_rel_dev']:.4f}, "
                  f"{r['frac_of_largest']:.3f} of largest")
