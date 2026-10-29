"""
Robustness of the two headline CRBs to the rib 300 nm derivatives (whose femwell solves carry remeshing
noise; see results_montecarlo/SUMMARY.md, step 3):
  (a) best one-etch set, stage D, realistic priors: strip TE 400, TE 600, TM 400, TM 600 + rib
      300/400/500/600 nm, 90 nm slab (0.162 deg);
  (b) two-width set at the 50 nm slab, stage D, realistic priors: strip TE 400-600 + rib 300 + 500 nm
      (0.367 deg).

1. Random perturbation: every derivative entry of the rib 300 nm structure in the set (all parameters,
   n_eff and n_g, the same multiplier at h/2, h, 2h) is multiplied by (1 + u), u uniform in +-p,
   p = 3% and 5% (as asked) plus 0.3% and 1% (the size of the mesh and remeshing errors), 2,000 draws each; the strip-angle CRB at h is reported as median and 5-95%.
2. Fixed mesh: the rib 300 nm geometric derivatives (w, t, a, slab) are replaced by those from
   src/fixedmesh_rib.py (nominal mesh reused across the stencil, no remeshing), if those solves exist.

Run from src/ (main environment):  python robustness_rib300.py
Writes results_montecarlo/robustness/data/robustness.json (log logs/robustness.log)
"""
import copy
import json
import os

import numpy as np

import etchdepth_femwell_analysis as ed
import priors_search as ps

OUT = os.path.join(ed.REPO, 'results_montecarlo', 'robustness', 'data')
FIXED = os.path.join(ed.REPO, 'results_montecarlo', 'fixedmesh', 'cache', 'solves_r12.jsonl')
RM = ('ring', 'mzi')
SETS = {
    'best_one_etch_D_realistic': ['strip_TE_400', 'strip_TE_600', 'strip_TM_400', 'strip_TM_600'] +
                                 [f'rib_TE_{w}_slab90' for w in (300, 400, 500, 600)],
    'two_width_slab50_D_realistic': [f'strip_TE_{w}' for w in (400, 450, 500, 550, 600)] +
                                    ['rib_TE_300_slab50', 'rib_TE_500_slab50'],
}
LEVELS, NDRAW, SEED = (0.003, 0.01, 0.03, 0.05), 2000, 20261029


def crb(units, Dall):
    return ps.rng(ps.evaluate([(u, RM) for u in units], Dall, 'D', ps.PRIORS['realistic']))


def perturbed(Dall, uid, p, rng):
    D = copy.deepcopy(Dall)
    keys = list(D[uid][1.0])
    mult = {k: 1 + rng.uniform(-p, p, 2) for k in keys}
    for f in D[uid]:
        for k in keys:
            if k in D[uid][f]:
                D[uid][f][k] = np.asarray(D[uid][f][k]) * mult[k]
    if 't_d' in keys:                          # keep the derived column consistent with its parts
        for f in D[uid]:
            D[uid][f]['t_d'] = D[uid][f]['t'] + D[uid][f]['s']
    return D


def fixed_rows(rows):
    """The regular solve rows with the rib 300 nm geometric steps replaced by the fixed-mesh solves, for
    each structure whose 25 fixed-mesh solves (nominal + w, t, a, slab at +-h/2, h, 2h) are all present.
    Returns (rows, {uid: number of fixed-mesh solves present})."""
    if not os.path.exists(FIXED):
        return None, {}
    got = {}
    with open(FIXED, encoding='utf-8') as f:
        for line in f:
            if line.strip():
                r = json.loads(line)
                k = r['key'].removeprefix('fixedmesh|')
                g, pol, W, S = k.split('|')[:4]
                got.setdefault(f'{g}_{pol}_{W}_slab{S}', {})[k] = r
    out = dict(rows)
    for uid, rs in got.items():
        if len(rs) == 25:
            for k, r in rs.items():
                out[k] = dict(rows.get(k, {}), **{kk: v for kk, v in r.items() if kk != 'key'}, key=k)
    return out, {uid: len(rs) for uid, rs in got.items()}


def main():
    os.makedirs(OUT, exist_ok=True)
    rows = ed.load()
    Dall = ps.derivatives(rows)                # strip symmetry imposed (default)
    res = {}
    for name, units in SETS.items():
        v, lo, hi = crb(units, Dall)
        rib300 = [u for u in units if u.startswith('rib_TE_300')][0]
        e = dict(set=units, rib300=rib300, nominal=dict(a_deg_h=v, a_deg_min=lo, a_deg_max=hi), perturbation={})
        print(f'{name}: nominal {v:.4g} [{lo:.4g}-{hi:.4g}]')
        for p in LEVELS:
            rng = np.random.default_rng([SEED, int(p * 1000), len(name)])
            a = np.array([crb(units, perturbed(Dall, rib300, p, rng))[0] for _ in range(NDRAW)])
            q = np.percentile(a, [5, 50, 95])
            e['perturbation'][f'{p:.1%}'] = dict(median=float(q[1]), p05=float(q[0]), p95=float(q[2]),
                                                 min=float(a.min()), max=float(a.max()),
                                                 median_rel_change=float(q[1] / v - 1),
                                                 p05_rel_change=float(q[0] / v - 1), p95_rel_change=float(q[2] / v - 1))
            print(f'  +-{p:.1%}: median {q[1]:.4g}, 5-95% {q[0]:.4g}-{q[2]:.4g} '
                  f'({q[0] / v - 1:+.1%} to {q[2] / v - 1:+.1%})')
        res[name] = e
    frows, counts = fixed_rows(rows)
    if frows is not None:
        Dfix = ps.derivatives(frows)
        for name, units in SETS.items():
            rib300 = res[name]['rib300']
            if counts.get(rib300, 0) < 25:
                res[name]['fixed_mesh'] = dict(status='incomplete', n_fixed_solves=counts.get(rib300, 0))
                print(f'{name}: fixed mesh incomplete ({counts.get(rib300, 0)} of 25 solves for {rib300})')
                continue
            v, lo, hi = crb(units, Dfix)
            dev = {k: {str(f): (np.asarray(Dfix[rib300][f][k]) / np.asarray(Dall[rib300][f][k]) - 1).tolist()
                       for f in ed.FACTORS} for k in ('w', 't', 'a', 's')}
            res[name]['fixed_mesh'] = dict(status='complete', n_fixed_solves=25, a_deg_h=v, a_deg_min=lo, a_deg_max=hi,
                                           rel_change=v / res[name]['nominal']['a_deg_h'] - 1,
                                           rib300_derivative_rel_diff_fixed_vs_remeshed=dev)
            print(f'{name}: fixed mesh {v:.4g} [{lo:.4g}-{hi:.4g}] ({v / res[name]["nominal"]["a_deg_h"] - 1:+.1%})')
            for k, d in dev.items():
                print(f'    d/d{k} fixed vs remeshed (n_eff, n_g) at h/2, h, 2h: ' +
                      '; '.join(f'{x[0]:+.2%}, {x[1]:+.2%}' for x in d.values()))
    with open(os.path.join(OUT, 'robustness.json'), 'w') as f:
        json.dump(res, f, indent=1)


if __name__ == '__main__':
    main()
