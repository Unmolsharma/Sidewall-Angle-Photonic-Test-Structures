"""
Numbers behind qualitative claims that no other script computes:
  1. MZI redundancy: for the headline sets at the primary noise (1 pm), the CRB with ring + MZI equals the CRB
     with the ring alone when the ring's n_eff noise is replaced by the combination of the ring and MZI n_eff
     noises (the MZI phase row is parallel to the ring n_eff row).
  2. Greedy vs exhaustive: final greedy set (8 structures) vs the exhaustive best of 8, every stage / prior /
     slab of the 1 pm search.
  3. EMpy vs femwell derivatives: relative differences per structure and derivative between the committed EMpy
     (cell-exact) database results_corrected/data/sensitivity_db.json and the femwell database
     results_femwell/data/sensitivity_db.json, and the EMpy step-check verdicts
     (results_corrected/data/step_check.csv). These are committed outputs of earlier sessions, read as is.
  4. Slab trend of the two-width set (strip TE 400-600 + rib 300 + 500 nm) in the 1 pm search.

Run from src/ (main environment):  python claims_checks.py
Writes results_final/data/claims_checks.json
"""
import collections
import csv
import json
import os

import numpy as np

import etchdepth_femwell_analysis as ed
import fullvector_sensitivity as fv
import montecarlo
import priors_search as ps
from assumption_sensitivity import headline_sets

REPO = ed.REPO
OUT = os.path.join(REPO, 'results_final', 'data')
S1PM = os.path.join(REPO, 'results_final', 'sym_search_1pm', 'data')


def crb(units, Dall, readouts):
    return ps.evaluate([(u, readouts) for u in units], Dall, 'D', ps.PRIORS['realistic'], factors=(1.0,))[1.0]['a']


def main():
    os.makedirs(OUT, exist_ok=True)
    res = {}
    Dall = ps.derivatives(ed.load())
    sets, _ = headline_sets()

    # 1. MZI redundancy
    sig = montecarlo.set_noise('1pm')
    s_ph_neff = fv.SIGMA_PHASE / montecarlo.KPH
    s_comb = (fv.SIGMA_NEFF ** -2 + s_ph_neff ** -2) ** -0.5
    mzi = {}
    for name, units in sets.items():
        both = crb(units, Dall, ('ring', 'mzi'))
        ring = crb(units, Dall, ('ring',))
        fv.SIGMA_NEFF = s_comb
        ring_comb = crb(units, Dall, ('ring',))
        fv.SIGMA_NEFF = sig['SIGMA_NEFF']
        mzi[name] = dict(ring_and_mzi=both, ring_only=ring, ring_with_combined_neff_noise=ring_comb,
                         rel_diff_combined_vs_both=ring_comb / both - 1)
    res['mzi_redundancy'] = dict(noise='1pm', sigma_neff=sig['SIGMA_NEFF'], sigma_phase_as_neff=s_ph_neff,
                                 combined_sigma_neff=s_comb, sets=mzi)

    # 2. greedy vs exhaustive (1 pm search)
    S = json.load(open(os.path.join(S1PM, 'search.json')))
    gv = []
    for tag, e in S.items():
        for slab, d in e['one_etch_per_slab'].items():
            if '8' in d and d.get('greedy'):
                g = d['greedy'][-1]['a_deg']
                x = d['8']['a_deg_h']
                gv.append(dict(search=tag, slab=int(slab), greedy=g, exhaustive=x,
                               ratio=(g / x) if np.isfinite(g) and np.isfinite(x) and x > 0 else None))
    fin = [r for r in gv if r['ratio'] is not None]
    res['greedy_vs_exhaustive'] = dict(cases=gv, n_cases=len(fin), n_greedy_worse=sum(r['ratio'] > 1 + 1e-9 for r in fin),
                                       n_equal=sum(abs(r['ratio'] - 1) <= 1e-9 for r in fin),
                                       max_ratio=max(r['ratio'] for r in fin), median_ratio=float(np.median([r['ratio'] for r in fin])))

    # 3. EMpy vs femwell
    E = {r['id']: r for r in json.load(open(os.path.join(REPO, 'results_corrected', 'data', 'sensitivity_db.json')))}
    F = {r['id']: r for r in json.load(open(os.path.join(REPO, 'results_femwell', 'data', 'sensitivity_db.json')))}
    keys = ['dneff_dw', 'dneff_dt', 'dneff_da', 'dng_dw', 'dng_dt', 'dng_da']
    rel = {k: {} for k in keys}
    for i in sorted(set(E) & set(F)):
        for k in keys:
            if F[i][k]:
                rel[k][i] = E[i][k] / F[i][k] - 1
    summ = {k: dict(median_abs=float(np.median(np.abs(list(v.values())))), max_abs=float(np.max(np.abs(list(v.values())))),
                    n=len(v)) for k, v in rel.items()}
    with open(os.path.join(REPO, 'results_corrected', 'data', 'step_check.csv'), newline='') as f:
        verdicts = collections.Counter(r['verdict'] for r in csv.DictReader(f))
    res['empy_vs_femwell'] = dict(n_structures=len(set(E) & set(F)), per_derivative=summ, rel_diff=rel,
                                  empy_step_check_verdicts=dict(verdicts))

    # 4. slab trend of the two-width set
    P = json.load(open(os.path.join(S1PM, 'slab_pairs.json')))['D|realistic']
    two = {v['slab']: v['a_deg_h'] for k, v in P.items() if 'rib_TE_300' in k and 'rib_TE_500' in k}
    sl = sorted(two)
    res['two_width_slab_trend'] = dict(noise='1pm', stage='D realistic', a_deg_by_slab={str(s): two[s] for s in sl},
                                       best_slab=min(two, key=two.get),
                                       monotonic_increasing_with_slab=bool(all(two[a] < two[b] for a, b in zip(sl, sl[1:]))))
    with open(os.path.join(OUT, 'claims_checks.json'), 'w') as f:
        json.dump(res, f, indent=1)
    print(json.dumps({k: v for k, v in res['mzi_redundancy'].items()}, indent=1))
    g = res['greedy_vs_exhaustive']
    print(f"greedy: {g['n_cases']} cases, worse {g['n_greedy_worse']}, equal {g['n_equal']}, max ratio {g['max_ratio']:.3g}")
    print('EMpy vs femwell:', {k: round(v['median_abs'], 4) for k, v in summ.items()}, dict(verdicts))
    print('two-width trend:', res['two_width_slab_trend'])


if __name__ == '__main__':
    main()
