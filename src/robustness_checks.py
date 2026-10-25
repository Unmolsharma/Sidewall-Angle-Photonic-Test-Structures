"""
Robustness of the downstream numbers to the derivative values themselves.

(1) Step-choice variants: rebuild the sensitivity database with every
    derivative taken at h/2 or at 2h instead of h (from step_check.csv; all
    values computed in recompute_sensitivities.py), and recompute every
    writeup number on each. A conclusion that only holds for one step choice
    is not supported by the derivatives.

(2) Sidewall-angle / width degeneracy. To first order in boundary
    perturbation theory, tilting a sidewall by a (top narrower, bottom fixed)
    moves the wall inward by tan(a)*y at height y. If the structure is
    mirror-symmetric top-to-bottom (strip, slot and SWG in this model: oxide
    above and below), the sidewall field weighting is symmetric about
    mid-height, so the tilt is exactly equivalent to a uniform width change
    dw = -h*tan(a) (h = height of the tilted wall). Then
        R = dn/da / (-h*tan(1 deg) * dn/dw) = 1   (for n_eff and n_g alike)
    and sidewall angle carries no first-order information beyond width.
    The rib's slab breaks the symmetry (h = t - slab = 130 nm). R - 1 is
    reported for every structure and both averagings.

Run from src/:  python robustness_checks.py
"""
import csv
import json
import os

import numpy as np

import writeup_numbers as wn

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D = os.path.join(REPO, 'results_corrected', 'data')
KEY = {('w', 'neff'): 'dneff_dw', ('t', 'neff'): 'dneff_dt', ('a', 'neff'): 'dneff_da',
       ('w', 'ng'): 'dng_dw', ('t', 'ng'): 'dng_dt', ('a', 'ng'): 'dng_da'}


def variant(records, steprows, col):
    out = []
    for r in records:
        r = json.loads(json.dumps(r))
        for (p, o), k in KEY.items():
            r[k] = steprows[(r['id'], p, o)][col]
        raw = np.array([r['dneff_dw'], r['dneff_dt'], r['dneff_da']])
        scaled = raw * wn.FAB_SCALE
        r['sens_vector_raw'], r['sens_vector_scaled'] = raw.tolist(), scaled.tolist()
        r['sens_vector_unit'] = (scaled / np.linalg.norm(scaled)).tolist()
        out.append(r)
    return out


def slab_variant(sweep, steprows, col):
    out = []
    for e in sweep:
        e = dict(e)
        sid = f"rib_slab{e['slab_nm']}_TE_{e['width_nm']}"
        for (p, o), k in KEY.items():
            e[k] = steprows[(sid, p, o)][col]
        out.append(e)
    return out


def degeneracy(records):
    rows = {}
    for r in records:
        h = 0.130 if r['geometry'] == 'rib' else 0.220
        k = -h * 1000 * np.tan(np.radians(1.0))
        rows[r['id']] = dict(R_neff=r['dneff_da'] / (k * r['dneff_dw']), R_ng=r['dng_da'] / (k * r['dng_dw']))
    return rows


if __name__ == '__main__':
    out = {}
    for avg, dbfile, slabfile in (('cell_exact', 'sensitivity_db.json', 'rib_slab_sweep_db_cell_exact.json'),
                                  ('original', 'sensitivity_db_original_rerun.json', 'rib_slab_sweep_db_original.json')):
        with open(os.path.join(D, dbfile)) as f:
            records = json.load(f)
        with open(os.path.join(D, slabfile)) as f:
            sweep = json.load(f)
        with open(os.path.join(D, 'step_check.csv')) as f:
            steprows = {(r['structure'], r['param'], r['output']): {c: float(r[c]) for c in ('D_half', 'D_h', 'D_2h')}
                        for r in csv.DictReader(f) if r['averaging'] == avg}
        out[avg] = {}
        for col in ('D_half', 'D_2h'):
            recs = variant(records, steprows, col)
            sv = slab_variant(sweep, steprows, col)
            p1, p2 = os.path.join(D, f'_tmp_db.json'), os.path.join(D, f'_tmp_slab.json')
            json.dump(recs, open(p1, 'w'))
            json.dump(sv, open(p2, 'w'))
            out[avg][col] = wn.main(p1, p2)
            os.remove(p1); os.remove(p2)
        out[avg]['degeneracy'] = degeneracy(records)
        print(f"\n=== {avg}: R = dn/da / (-h tan1deg dn/dw), 1.000 = exactly degenerate with width ===")
        for fam in ('strip_TE', 'strip_TM', 'slot_TE', 'slot_TM', 'swg_TE', 'rib_TE'):
            v = [(k, x) for k, x in out[avg]['degeneracy'].items() if k.startswith(fam)]
            print(f"  {fam:9s} R_neff " + ' '.join(f"{x['R_neff']:.4f}" for _, x in v) +
                  "   R_ng " + ' '.join(f"{x['R_ng']:.4f}" for _, x in v))
    with open(os.path.join(D, 'robustness_checks.json'), 'w') as f:
        json.dump(out, f, indent=1, default=float)
    print("\nSaved robustness_checks.json")
