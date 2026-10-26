"""
Cross-check that the analysis scripts agree with each other on the femwell
database: the same quantity read from every output that reports it, plus an
independent Fisher computation straight from probe_db.json.

  writeup_numbers_femwell.json   (writeup_numbers.main via downstream_analysis.py)
  fisher_reanalysis.json         (original fisher_reanalysis.py)
  angle_crb.log                  (original angle_crb.py; agrees if equal to the digits printed)
  robustness_checks.json         (degeneracy R, downstream_analysis.py)
  crosscheck_and_ratios.json     (R from femwell_postprocess / comparison table)
  impact_sidewall_angle.json     (femwell_impact.py)
  sensitivity_db.json, rib_slab_sweep_db_femwell.json

Run from src/:  python femwell_script_consistency.py
Writes results_femwell/data/script_consistency.json
"""
import json
import os
import re

import numpy as np

from fisher_reanalysis import fisher_from_probes

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D = os.path.join(REPO, 'results_femwell', 'data')
LOGS = os.path.join(REPO, 'results_femwell', 'logs', 'original_scripts_on_femwell_db')
TOL_JSON = 1e-9


def angle_std(ids, by_id):
    F, _ = fisher_from_probes(ids, by_id)
    return float(np.sqrt(np.linalg.inv(F)[2, 2]))


def main():
    N = json.load(open(os.path.join(D, 'writeup_numbers_femwell.json')))
    FR = json.load(open(os.path.join(D, 'fisher_reanalysis.json')))
    rob = json.load(open(os.path.join(D, 'robustness_checks.json')))['femwell']['degeneracy']
    xr = json.load(open(os.path.join(D, 'crosscheck_and_ratios.json')))['ratios']
    imp = json.load(open(os.path.join(D, 'impact_sidewall_angle.json')))
    recs = {r['id']: r for r in json.load(open(os.path.join(D, 'sensitivity_db.json')))}
    probes = json.load(open(os.path.join(D, 'probe_db.json')))
    by_id = {p['probe_id']: p for p in probes}
    sweep = json.load(open(os.path.join(D, 'rib_slab_sweep_db_femwell.json')))
    log = open(os.path.join(LOGS, 'angle_crb.log'), encoding='utf-8').read()

    checks = []

    def chk(name, a, b, tol, src):
        rel = abs(a - b) / max(abs(b), 1e-300)
        checks.append(dict(quantity=name, sources=src, a=a, b=b, rel_diff=rel, agree=bool(rel <= tol)))

    def chk_log(name, a, printed, src):
        # a printed log value agrees if it is a rounded to the digits shown
        b = float(printed)
        dec = len(printed.split('.')[1]) if '.' in printed else 0
        ok = abs(a - b) <= 0.5 * 10 ** -dec * (1 + 1e-9)
        checks.append(dict(quantity=name, sources=src, a=a, b=b, printed=printed,
                           rel_diff=abs(a - b) / max(abs(b), 1e-300), agree=bool(ok)))

    # strip baseline (TE+TM, ring+mzi, all widths) and the two-ring minimal budget
    strip_ids = [p['probe_id'] for p in probes if p['geometry'] == 'strip']
    base = angle_std(strip_ids, by_id)
    two = angle_std(['strip_TE_550_ring', 'strip_TE_600_ring'], by_id)
    m = re.search(r'strip-only \(TE\+TM, ring\+mzi\) baseline: .*angle std=([\d.]+) deg', log)
    chk('strip baseline angle CRB (deg)', N['crb_large_budget']['strip_baseline'], base, TOL_JSON,
        'writeup_numbers vs probe_db recompute')
    chk_log('strip baseline angle CRB (deg)', base, m.group(1), 'probe_db recompute vs angle_crb.log')
    m = re.search(r'strip TE@550\+TE@600 alone: .*angle std=([\d.]+) deg', log)
    chk('two strip TE rings 550+600 angle CRB (deg)', N['crb_two_rings']['angle_std'], two, TOL_JSON,
        'writeup_numbers vs probe_db recompute')
    chk_log('two strip TE rings 550+600 angle CRB (deg)', two, m.group(1), 'probe_db recompute vs angle_crb.log')
    # additions with noise multipliers: writeup_numbers vs angle_crb (independent code paths)
    for key, pat in (('+rib x1.0', r'rib \[conservative 1\.0x\]'), ('+rib x0.2', r'rib \[optimistic 0\.2x\]'),
                     ('+slot x4.0', r'slot \[sourced 4\.0x\]')):
        m = re.search(pat + r'\s+[\d.e+-]+\s+([\d.]+)\s+([\d.]+)x', log)
        chk_log(f'strip baseline {key}: improvement factor', N['crb_large_budget'][key]['improvement'],
                m.group(2), 'writeup_numbers vs angle_crb.log')
    for key, pat in (('strip x1.0', r'strip \[reference, unchanged\]'), ('rib x1.0', r'rib \[conservative 1\.0x, no benefit assumed\]'),
                     ('rib x0.2', r'rib \[optimistic 0\.2x, sourced wide-rib lit\.\]'), ('slot x4.0', r'slot \[sourced 4\.0x\]')):
        m = re.search(pat + r'\s+[\d.e+-]+\s+[\d.e+-]+\s+([\d.]+)', log)
        chk_log(f'{key} alone angle CRB (deg)', N['crb_alone'][key], m.group(1), 'writeup_numbers vs angle_crb.log')
    m = re.search(r'rib TE@400 \[cons\. 1\.0x\]\s+[\d.e+-]+\s+[\d.e+-]+\s+([\d.]+)\s+([\d.]+)x', log)
    chk_log('two rings + rib 400 (cons.): improvement factor', N['crb_two_rings']['+rib400 x1.0']['improvement'],
            m.group(2), 'writeup_numbers vs angle_crb.log')

    # ring-pair landscape: writeup_numbers vs fisher_reanalysis (same combos, all included)
    conds = [mm['cond'] for _, mm in FR['(a) TE ring + TE ring']]
    chk('TE+TE ring pairs: best cond', N['ring_pairs']['cond_best'], min(conds), TOL_JSON,
        'writeup_numbers vs fisher_reanalysis.json')
    chk('TE+TE ring pairs: worst cond', N['ring_pairs']['cond_worst'], max(conds), TOL_JSON,
        'writeup_numbers vs fisher_reanalysis.json')
    F, _ = fisher_from_probes(['strip_TE_550_ring', 'strip_TE_600_ring'], by_id)
    chk('TE 550+600 pair lambda_min', N['additions_550_600']['baseline']['lambda_min'],
        float(np.linalg.eigvalsh(F)[0]), TOL_JSON, 'writeup_numbers vs probe_db recompute')

    # degeneracy R: robustness_checks vs crosscheck_and_ratios vs direct from sensitivity_db
    wall = dict(strip=0.22, swg=0.22, slot=0.22, rib=0.13)
    for sid, r in rob.items():
        rec = recs[sid]
        direct = rec['dneff_da'] / rec['dneff_dw'] / (-wall[rec['geometry']] * 1000 * np.tan(np.radians(1.0)))
        chk(f'R_neff {sid}', r['R_neff'], direct, TOL_JSON, 'robustness_checks vs sensitivity_db')
        if sid in xr:
            chk(f'R_neff {sid}', xr[sid]['R_neff'], direct, TOL_JSON, 'crosscheck_and_ratios vs sensitivity_db')

    # 6b impact: the "true" rib change must be the database derivative times the angle step
    rib = recs['rib_TE_500']
    for da in (1.0, 2.0):
        for fit in ('width_only', 'width_and_thickness'):
            for obs in ('neff', 'ng'):
                chk(f'impact da={da:g} {fit} true d{obs}', imp[f'da={da:g}deg'][fit]['rib'][obs]['true'],
                    rib[f'd{obs}_da'] * da, TOL_JSON, 'femwell_impact vs sensitivity_db')

    # slab sweep: same slabs in the sweep file and in writeup_numbers
    slabs_db = sorted({e['slab_nm'] for e in sweep})
    slabs_wn = sorted(int(s) for s in N['slab_sweep'])
    checks.append(dict(quantity='slab sweep slabs', sources='writeup_numbers vs rib_slab_sweep_db_femwell.json',
                       a=slabs_wn, b=slabs_db, agree=slabs_wn == slabs_db))

    bad = [c for c in checks if not c['agree']]
    for c in checks:
        if 'rel_diff' in c:
            print(f"{'OK ' if c['agree'] else 'BAD'} {c['quantity']:48s} {c['sources']:42s} rel diff {c['rel_diff']:.1e}")
        else:
            print(f"{'OK ' if c['agree'] else 'BAD'} {c['quantity']:48s} {c['a']} vs {c['b']}")
    print(f'{len(checks)} checks, {len(bad)} disagree')
    with open(os.path.join(D, 'script_consistency.json'), 'w') as f:
        json.dump(dict(n_checks=len(checks), n_disagree=len(bad), checks=checks), f, indent=1)


if __name__ == '__main__':
    main()
