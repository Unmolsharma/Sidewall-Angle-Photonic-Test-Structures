"""
Step 6c cross-check and angle/width ratios on the femwell database.

  - femwell batch n_eff derivatives (strip/rib TE 500) vs the earlier femwell_crosscheck r12 run and
    vs EMpy grid-shift averaging (results_anisotropic/data/gridshift_validation.json; no new EMpy runs)
  - (dn/da)/(dn/dw) and R = that / (-h tan 1 deg) for every structure and every slab-sweep rib,
    h = height of the tilted wall (220 nm; rib 220 nm - slab)

Run from src/:  python femwell_crosscheck_ratios.py
Writes results_femwell/data/crosscheck_and_ratios.json
"""
import json
import os

import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D = os.path.join(REPO, 'results_femwell', 'data')
A = os.path.join(REPO, 'results_anisotropic', 'data')


def main():
    F = {r['id']: r for r in json.load(open(os.path.join(D, 'sensitivity_db.json')))}
    ref = json.load(open(os.path.join(A, 'femwell_derivatives_r12.json')))
    gs = json.load(open(os.path.join(A, 'gridshift_validation.json')))
    out = {'crosscheck': {}, 'ratios': {}}
    print("6c: femwell batch vs EMpy grid-shift (production mesh, N=4 / N=8), n_eff derivatives at h")
    for sid, g in (('strip_TE_500', 'strip'), ('rib_TE_500', 'rib')):
        for p in 'wta':
            f = F[sid][f'dneff_d{p}']
            r0 = ref[f'{g}_TE_500_oxide']['D'][f'{p}|neff|1.0']
            row = {'femwell_batch': f, 'femwell_crosscheck_r12_pct': 100 * (r0 - f) / abs(f)}
            for avg in ('anisotropic', 'cell_exact'):
                for n in (4, 8):
                    row[f'{avg}_N{n}_pct'] = 100 * (gs[f'{avg}|{g}_TE_500']['30/18'][f'N{n}'][f'{p}|neff|1.0'] - f) / abs(f)
            out['crosscheck'][f'{sid}|dneff/d{p}'] = row
            print(f"  {sid} dneff/d{p}: femwell {f:+.4e} (earlier r12 run {row['femwell_crosscheck_r12_pct']:+.3f}%)  "
                  f"aniso N4 {row['anisotropic_N4_pct']:+.1f}% N8 {row['anisotropic_N8_pct']:+.1f}%  "
                  f"cell-exact N4 {row['cell_exact_N4_pct']:+.1f}% N8 {row['cell_exact_N8_pct']:+.1f}%")
    print("\nangle/width ratio (dn_eff/da)/(dn_eff/dw) [nm/deg] and R, femwell r12, oxide clad")
    for sid, r in F.items():
        h = 0.13 if r['geometry'] == 'rib' else 0.22
        k = -h * 1000 * np.tan(np.radians(1))
        ra, rg = r['dneff_da'] / r['dneff_dw'], r['dng_da'] / r['dng_dw']
        out['ratios'][sid] = dict(ratio_neff=ra, R_neff=ra / k, ratio_ng=rg, R_ng=rg / k)
        print(f"  {sid:14s} n_eff {ra:+.4f} (R {ra / k:.4f})   n_g {rg:+.4f} (R {rg / k:.4f})")
    for e in json.load(open(os.path.join(D, 'rib_slab_sweep_db_femwell.json'))):
        k = -(0.22 - e['slab_nm'] / 1000) * 1000 * np.tan(np.radians(1))
        ra = e['dneff_da'] / e['dneff_dw']
        out['ratios'][f"rib_slab{e['slab_nm']}_TE_{e['width_nm']}"] = dict(ratio_neff=ra, R_neff=ra / k)
        print(f"  rib slab {e['slab_nm']:3d} w {e['width_nm']}: n_eff {ra:+.4f} (R {ra / k:.4f})")
    with open(os.path.join(D, 'crosscheck_and_ratios.json'), 'w') as f:
        json.dump(out, f, indent=1)


if __name__ == '__main__':
    main()
