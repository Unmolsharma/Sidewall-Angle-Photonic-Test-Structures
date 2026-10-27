"""
Where does the strips-only sidewall-angle information come from, and how robust is it?

Run from src/ (main environment):  python strips_angle_info.py
Reads   femwell caches (via priors_search / etchdepth_femwell_analysis), results_priors_search/data/search.json
Writes  results_montecarlo/data/strips_info.json, strips_info_perturbation.csv

Symmetry argument: an oxide-clad strip is mirror-symmetric about mid-height. A wall tilt about the bottom
edge is a uniform wall shift of (t/2) tan(a) per side (identical to a width change of -t tan(a)) plus a
rotation about mid-height, which is odd under the mirror and changes n_eff, n_g only at second order. So
at first order dn/da = -t (pi/180) dn/dw [a in degrees] exactly, for TE and TM and for n_eff and n_g:
R = (pi/180) / tan(1 deg) = 0.99990. Any strip angle information beyond that comes from the residual
R - 0.99990, i.e. finite-difference truncation, mesh error, and the slightly off-centre domain boundary.

Checks:
  1. R - R_exact per strip at h/2, h, 2h (truncation scales as h^2)
  2. information decomposition: the angle's marginal Fisher information is |P J_a|^2, with P the projector
     orthogonal to the other columns; split by row (TE/TM, n_eff/n_g/MZI phase)
  3. exact-symmetry model: strip a-columns replaced by -t (pi/180) x w-column
  4. random perturbation of the strip derivatives (each entry x (1 + u), u uniform in +-p), many draws
Sets: strip TE 5 widths; best strips-only set (stage C, 8); strip + rib70 500; best one-etch multi-width
rib sets (stage C none, stage D realistic, 8 structures, from search.json). Rib derivatives untouched.
"""
import csv
import json
import os

import numpy as np

import etchdepth_femwell_analysis as ed
import priors_search as ps

REPO = ed.REPO
OUT = os.path.join(REPO, 'results_montecarlo', 'data')
R_EXACT = np.radians(1.0) / np.tan(np.radians(1.0))
T_UM = 220.0  # nm
NDRAW = 2000


def sets_from_search():
    S = json.load(open(os.path.join(REPO, 'results_priors_search', 'data', 'search.json')))
    rm = lambda combo: [(u, ('ring', 'mzi')) for u in combo]
    strip5 = [f'strip_TE_{w}' for w in (400, 450, 500, 550, 600)]
    return {
        'strip TE 5 widths': ('C', 'none', rm(strip5)),
        'best strips-only (C, 8)': ('C', 'none', rm(S['C|none|strips only']['best']['8']['set'])),
        'strip TE 5 + rib70 500 (C)': ('C', 'none', rm(strip5 + ['rib_TE_500_slab70'])),
        'strip TE 5 + rib70 500 + 300 (C)': ('C', 'none', rm(strip5 + ['rib_TE_500_slab70', 'rib_TE_300_slab70'])),
        'strip TE 5 + rib70 500 + 300 (E, realistic)': ('E', 'realistic', rm(strip5 + ['rib_TE_500_slab70', 'rib_TE_300_slab70'])),
        'best one-etch set (C, none, 8)': ('C', 'none', rm(S['C|none|one partial etch']['best']['8']['set'])),
        'best one-etch set (D, realistic, 8)': ('D', 'realistic', rm(S['D|realistic|one partial etch']['best']['8']['set'])),
    }


def crb(items, Dall, stage, prior):
    params, J = ps.jacobian(items, Dall, 1.0, stage, ps.PRIORS[prior])
    return ps.strip_angle(J, params)


def r_check(Dall):
    out = {}
    for uid, Du in Dall.items():
        if not uid.startswith('strip'):
            continue
        out[uid] = {}
        for f in ed.FACTORS:
            D = Du[f]
            k = -T_UM * np.tan(np.radians(1.0))
            out[uid][str(f)] = {obs: float(D['a'][i] / D['w'][i] / k - R_EXACT) for i, obs in enumerate(('neff', 'ng'))}
    return out


def decomposition(items, Dall, stage, prior):
    """Per-row share of the angle's marginal information (no rows for priors are split out)."""
    params, J = ps.jacobian(items, Dall, 1.0, stage, ps.PRIORS[prior])
    labels = []
    for uid, readouts in items:
        if 'ring' in readouts:
            labels += [(uid, 'neff'), (uid, 'ng')]
        if 'mzi' in readouts:
            labels += [(uid, 'phase')]
    labels += [('prior', 'prior')] * (len(J) - len(labels))
    i = params.index('a')
    others = np.delete(J, i, axis=1)
    Q, s, _ = np.linalg.svd(others, full_matrices=False)
    Q = Q[:, s > ed.RANK_TOL * s[0]]
    resid = J[:, i] - Q @ (Q.T @ J[:, i])          # component of the angle column no other parameter can mimic
    info = resid ** 2
    tot = info.sum()
    groups = {}
    for (uid, obs), x in zip(labels, info):
        pol = uid.split('_')[1] if uid.startswith('strip') else 'rib' if uid.startswith('rib') else 'prior'
        key = f"{'strip ' + pol if pol in ('TE', 'TM') else pol} {obs}"
        groups[key] = groups.get(key, 0.0) + float(x)
    return dict(total_info=float(tot), crb_deg=float(1 / np.sqrt(tot)) if tot > 0 else float('inf'),
                share={k: v / tot for k, v in sorted(groups.items())})


def subset_crbs(Dall):
    """Strips-only precision from TE-only / TM-only structures and from n_eff-type or n_g rows only."""
    base = sets_from_search()['best strips-only (C, 8)'][2]
    out = {}
    for name, keep in (('all', lambda u, o: True), ('TE structures only', lambda u, o: '_TE_' in u),
                       ('TM structures only', lambda u, o: '_TM_' in u),
                       ('n_eff + MZI phase only (no n_g)', lambda u, o: o != 'ng'),
                       ('n_g only', lambda u, o: o == 'ng')):
        params, J = ps.jacobian(base, Dall, 1.0, 'C', {})
        labels = [(u, o) for u, r in base for o in ('neff', 'ng', 'phase')]
        rows = [j for j, (u, o) in enumerate(labels) if keep(u, o)]
        Js = J[rows]
        Js = Js[:, np.any(Js != 0, axis=0)]
        out[name] = ps.strip_angle(Js, [p for p, c in zip(params, np.any(J[rows] != 0, axis=0)) if c])
    return out


def symmetrized(Dall):
    """Copy of the derivatives with every strip a-column set to the exact first-order value."""
    S = {}
    for uid, Du in Dall.items():
        S[uid] = {}
        for f, D in Du.items():
            D2 = dict(D)
            if uid.startswith('strip'):
                D2['a'] = -T_UM * np.radians(1.0) * np.asarray(D['w'])
            S[uid][f] = D2
    return S


def perturbed(Dall, p, rng):
    S = {}
    for uid, Du in Dall.items():
        S[uid] = {}
        for f, D in Du.items():
            if uid.startswith('strip'):
                S[uid][f] = {k: np.asarray(v) * (1 + rng.uniform(-p, p, size=np.shape(v))) for k, v in D.items()}
            else:
                S[uid][f] = D
    return S


def main():
    os.makedirs(OUT, exist_ok=True)
    Dall = ps.derivatives(ed.load())
    sets = sets_from_search()
    res = dict(R_exact=R_EXACT, R_minus_exact=r_check(Dall))
    res['nominal'] = {n: crb(it, Dall, st, pr) for n, (st, pr, it) in sets.items()}
    res['decomposition'] = {n: decomposition(it, Dall, st, pr) for n, (st, pr, it) in sets.items()}
    res['strips_only_subsets'] = subset_crbs(Dall)
    Ds = symmetrized(Dall)
    res['exact_symmetry'] = {n: crb(it, Ds, st, pr) for n, (st, pr, it) in sets.items()}
    rng = np.random.default_rng(20261027)
    rows, summ = [], {}
    for p in (0.001, 0.003, 0.01, 0.03):
        vals = {n: [] for n in sets}
        for k in range(NDRAW):
            Dp = perturbed(Dall, p, rng)
            for n, (st, pr, it) in sets.items():
                v = crb(it, Dp, st, pr)
                vals[n].append(v)
                rows.append(dict(p=p, draw=k, set=n, a_deg=v))
        summ[str(p)] = {n: dict(nominal=res['nominal'][n], median=float(np.median(v)), p05=float(np.percentile(v, 5)),
                                p95=float(np.percentile(v, 95)), min=float(np.min(v)), max=float(np.max(v)),
                                frac_better_than_nominal_by_2x=float(np.mean(np.array(v) < res['nominal'][n] / 2)))
                        for n, v in vals.items()}
        print(f'perturbation +-{p:.1%}:')
        for n, s in summ[str(p)].items():
            print(f"  {n:36s} nominal {s['nominal']:.3g}  median {s['median']:.3g}  5-95% {s['p05']:.3g}..{s['p95']:.3g}")
    res['perturbation'] = summ
    with open(os.path.join(OUT, 'strips_info.json'), 'w') as f:
        json.dump(res, f, indent=1)
    with open(os.path.join(OUT, 'strips_info_perturbation.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print('R - R_exact vs step (h/2, h, 2h):')
    for u, d in res['R_minus_exact'].items():
        print(f'  {u:14s} ' + '  '.join(f"{o}: " + ' '.join(f"{d[f][o]:+.1e}" for f in ('0.5', '1.0', '2.0')) for o in ('neff', 'ng')))
    print('R - R_exact (strips, h):', {u: {k: f'{v:+.1e}' for k, v in d['1.0'].items()} for u, d in res['R_minus_exact'].items()})
    print('nominal:', res['nominal'])
    print('exact symmetry:', res['exact_symmetry'])
    print('strips-only subsets:', res['strips_only_subsets'])
    for n, d in res['decomposition'].items():
        print(n, {k: f'{v:.2f}' for k, v in d['share'].items()})


if __name__ == '__main__':
    main()
