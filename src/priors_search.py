"""
Priors study and structure search on the existing femwell r12 derivatives (no new solves).

Run from src/ (main environment):  python priors_search.py
Reads   results_femwell/cache/solves_r12.jsonl, results_etchdepth_femwell/cache/solves_r12.jsonl
Writes  results_priors_search/data/priors.json, priors_table.csv, search.json, search_table.csv

Model, parameters, units and stages as in etchdepth_femwell_analysis.py (w, t, d [nm]; a, ar [deg];
n_Si, n_ox [1e-3]; T [K]). Gaussian priors enter as pseudo-observations: one extra Jacobian row
e_p / sigma_p per prior, so F_total = F + diag(1/sigma^2); a prior on the rib-angle offset adds the row
(e_ar - e_a) / sigma. The identifiability test and CRB are the same as in that script.

1. Priors: strip-angle CRB for strip alone, strip + rib70 (500 nm), strip + rib70 (500 + 300 nm), stages
   B-E, for each prior on its own and for combinations.
2. Search, stages C and D, no priors and realistic priors: exhaustive best set of each size 1-8 and
   greedy selection, over every structure that has the derivatives the stage needs. Each candidate is
   one structure read out by a ring (n_eff, n_g) and an MZI (phase). Also restricted to ribs from one
   partial etch, since two slabs on one chip need two separate partial etches.
Noise unchanged (placeholders, see fullvector_sensitivity.py); rib multiplier 1.0 unless stated.
"""
import csv
import itertools
import json
import os

import numpy as np

import etchdepth_femwell_analysis as ed

OUT = os.path.join(ed.REPO, 'results_priors_search')
FACTORS = ed.FACTORS
STRIP_UNITS = [('strip', pol, w, None) for pol in ('TE', 'TM') for w in (400, 450, 500, 550, 600)]
RIB_UNITS = [('rib', 'TE', w, s) for s in (50, 70, 90, 110, 130, 150) for w in (300, 400, 500, 600)]

# prior options; sigma in parameter units (n in 1e-3, so 1e-3 -> 1.0 and 1e-4 -> 0.1)
PRIORS = {
    'none': {},
    'd 1 nm': {'d': 1.0}, 'd 2 nm': {'d': 2.0}, 'd 5 nm': {'d': 5.0},
    'n_Si 1e-3': {'nsi': 1.0}, 'n_Si 1e-4': {'nsi': 0.1},
    'n_ox 1e-3': {'nox': 1.0}, 'n_ox 1e-4': {'nox': 0.1},
    'n_Si, n_ox 1e-3': {'nsi': 1.0, 'nox': 1.0}, 'n_Si, n_ox 1e-4': {'nsi': 0.1, 'nox': 0.1},
    'ar-a 1 deg': {'dar': 1.0}, 'ar-a 2 deg': {'dar': 2.0}, 'ar-a 5 deg': {'dar': 5.0},
    'T 1 K': {'T': 1.0},
    'realistic': {'d': 2.0, 'nsi': 0.1, 'nox': 1.0, 'T': 1.0},
    'realistic + ar-a 2 deg': {'d': 2.0, 'nsi': 0.1, 'nox': 1.0, 'T': 1.0, 'dar': 2.0},
    'all optimistic': {'d': 1.0, 'nsi': 0.1, 'nox': 0.1, 'T': 0.1, 'dar': 1.0},
}
REALISTIC = 'realistic'


def derivatives(rows):
    """Like ed.derivatives, but keeps whatever parameters a structure has solves for."""
    out = {}
    for u in STRIP_UNITS + RIB_UNITS:
        uid = f'{u[0]}_{u[1]}_{u[2]}' + (f'_slab{u[3]}' if u[3] else '')
        try:
            ed.nng(rows, u, 'nom', 0.0)
        except KeyError:
            continue
        out[uid] = {}
        for f in FACTORS:
            D = {}
            for p in ('w', 't', 'a', 's', 'nsi', 'nox'):
                h = f * ed.H[p]
                try:
                    D[p] = (ed.nng(rows, u, p, h) - ed.nng(rows, u, p, -h)) / (2 * h)
                except KeyError:
                    continue
                if p in ed.MAT:
                    D[p] = D[p] * 1e-3
            if 's' in D:
                D['t_d'] = D['t'] + D['s']
            out[uid][f] = D
    return out


def eligible(Dall, stage):
    need = {'strip': {'w', 't', 'a'}, 'rib': {'w', 't', 'a', 's'}}
    mat = {'nsi', 'nox'} if stage in ('D', 'E') else set()
    return [u for u, Du in Dall.items() if (need[u.split('_')[0]] | mat) <= set(Du[1.0])]


def jacobian(items, Dall, f, stage, prior, rib_mult=1.0):
    params, J = ed.rows_for(items, Dall, f, stage, rib_mult)
    extra = []
    for key, sig in prior.items():
        if key == 'dar':
            for p in params:
                if p.startswith('ar'):
                    r = np.zeros(len(params))
                    r[params.index(p)], r[params.index('a')] = 1 / sig, -1 / sig
                    extra.append(r)
        else:
            for p in params:
                if p == key or (key == 'd' and p.startswith('d')):
                    r = np.zeros(len(params))
                    r[params.index(p)] = 1 / sig
                    extra.append(r)
    return params, (np.vstack([J] + extra) if extra else J)


def evaluate(items, Dall, stage, prior, rib_mult=1.0, factors=FACTORS):
    out = {}
    for f in factors:
        params, J = jacobian(items, Dall, f, stage, prior, rib_mult)
        ev, evec = np.linalg.eigh(J.T @ J)
        out[f] = dict(a=strip_angle(J, params), eig_min=float(ev[0]),
                      eigvec_min={p: float(x) for p, x in zip(params, evec[:, 0])},
                      rib_angles={p: ed.marginal_std(J, i) for i, p in enumerate(params) if p.startswith('ar')})
    return out


def strip_angle(J, params):
    # a set with no strip has no strip-angle column: the strip angle is not determined at all
    return ed.marginal_std(J, params.index('a')) if 'a' in params else float('inf')


def rng(e):
    v = [e[f]['a'] for f in FACTORS]
    return e[1.0]['a'], min(v), max(v)


def priors_study(Dall):
    strip = [(f'strip_TE_{w}', ('ring', 'mzi')) for w in (400, 450, 500, 550, 600)]
    sets = {'strip alone': strip,
            'strip + rib70 (500 nm)': strip + [('rib_TE_500_slab70', ('ring', 'mzi'))],
            'strip + rib70 (500 + 300 nm)': strip + [('rib_TE_500_slab70', ('ring', 'mzi')), ('rib_TE_300_slab70', ('ring', 'mzi'))]}
    res, table = {}, []
    for stage in ('B', 'C', 'D', 'E'):
        for pname, prior in PRIORS.items():
            for sname, items in sets.items():
                for mult in (1.0, 0.2):
                    e = evaluate(items, Dall, stage, prior, mult)
                    res.setdefault(stage, {}).setdefault(pname, {}).setdefault(sname, {})[mult] = e
                    v, lo, hi = rng(e)
                    table.append(dict(stage=stage, prior=pname, set=sname, rib_mult=mult, a_deg_h=v, a_deg_min=lo,
                                      a_deg_max=hi, eig_min_h=e[1.0]['eig_min'],
                                      rib_angle_deg_h=json.dumps(e[1.0]['rib_angles'])))
    return res, table


def crb_h(items, Dall, stage, prior):
    params, J = jacobian(items, Dall, 1.0, stage, prior)
    norms = np.linalg.norm(J, axis=0)
    s = np.linalg.svd(J / norms, compute_uv=False)
    return strip_angle(J, params), int(np.sum(s > ed.RANK_TOL * s[0]))


def search(Dall, stage, prior, mode, kmax=8):
    """mode: 'any slabs', 'one partial etch' (ribs from at most one slab), or 'strips only' (the fair
    baseline for the rib's contribution)."""
    pool = eligible(Dall, stage)
    if mode == 'strips only':
        pool = [u for u in pool if u.startswith('strip')]
    cand = {u: (u, ('ring', 'mzi')) for u in pool}

    def ok(combo):
        slabs = {u.split('slab')[1] for u in combo if 'slab' in u}
        return mode != 'one partial etch' or len(slabs) <= 1

    best = {}
    for k in range(1, min(kmax, len(pool)) + 1):
        top = None
        for combo in itertools.combinations(pool, k):
            if not ok(combo):
                continue
            a, _ = crb_h([cand[u] for u in combo], Dall, stage, prior)
            if top is None or a < top[0]:
                top = (a, combo)
        best[k] = top
    # greedy: smallest strip-angle CRB; while every candidate is unidentifiable, prefer the largest rank
    chosen, steps = [], []
    for k in range(1, kmax + 1):
        scored = []
        for u in pool:
            if u in chosen or not ok(chosen + [u]):
                continue
            a, r = crb_h([cand[x] for x in chosen + [u]], Dall, stage, prior)
            scored.append((a, -r, u))
        if not scored:
            break
        a, negr, u = min(scored)
        chosen.append(u)
        steps.append(dict(added=u, a_deg=a, rank=-negr))
    return pool, best, steps, cand


def main(parts=('priors', 'dsweep', 'search')):
    Dall = derivatives(ed.load())
    os.makedirs(os.path.join(OUT, 'data'), exist_ok=True)
    if 'priors' in parts:
        run_priors(Dall)
    if 'dsweep' in parts:
        run_dsweep(Dall)
    if 'search' in parts:
        run_search(Dall)


def run_priors(Dall):
    res, table = priors_study(Dall)
    with open(os.path.join(OUT, 'data', 'priors.json'), 'w') as f:
        json.dump({s: {p: {n: {str(m): {str(k): v for k, v in e.items()} for m, e in d.items()} for n, d in x.items()}
                       for p, x in y.items()} for s, y in res.items()}, f, indent=1)
    with open(os.path.join(OUT, 'data', 'priors_table.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(table[0].keys()))
        w.writeheader()
        w.writerows(table)
    print('priors: strip-angle CRB at h [min..max], rib x1.0')
    for r in table:
        if r['rib_mult'] == 1.0:
            print(f"  {r['stage']} {r['prior']:24s} {r['set']:30s} {r['a_deg_h']:.4g} [{r['a_deg_min']:.4g}..{r['a_deg_max']:.4g}]")


def run_dsweep(Dall):
    """How well must the etch depth be known to recover the depth-known (A*) precision? Stage B and C,
    strip + rib70 500 nm and strip + rib70 500 + 300 nm, depth prior from 0.01 to 5 nm."""
    strip = [(f'strip_TE_{w}', ('ring', 'mzi')) for w in (400, 450, 500, 550, 600)]
    sets = {'strip + rib70 (500 nm)': strip + [('rib_TE_500_slab70', ('ring', 'mzi'))],
            'strip + rib70 (500 + 300 nm)': strip + [('rib_TE_500_slab70', ('ring', 'mzi')), ('rib_TE_300_slab70', ('ring', 'mzi'))]}
    out = {}
    for stage in ('B', 'C'):
        for name, items in sets.items():
            for sig in (0.01, 0.03, 0.1, 0.3, 1.0, 2.0, 5.0, None):
                e = evaluate(items, Dall, stage, {'d': sig} if sig else {})
                v, lo, hi = rng(e)
                out.setdefault(stage, {}).setdefault(name, {})['free' if sig is None else str(sig)] = dict(a_deg_h=v, a_deg_min=lo, a_deg_max=hi)
                print(f"  {stage} {name:30s} sigma_d {sig}: {v:.4g} [{lo:.4g}..{hi:.4g}]")
        # reference: depth known exactly (stage A*)
        for name, items in sets.items():
            e = evaluate(items, Dall, 'A*', {})
            out.setdefault('A* (depth known)', {})[name] = dict(zip(('a_deg_h', 'a_deg_min', 'a_deg_max'), rng(e)))
    with open(os.path.join(OUT, 'data', 'd_prior_sweep.json'), 'w') as f:
        json.dump(out, f, indent=1)


def run_search(Dall):
    out, rows = {}, []
    for stage in ('C', 'D'):
        for pname in ('none', REALISTIC):
            for mode in ('any slabs', 'one partial etch', 'strips only'):
                pool, best, steps, cand = search(Dall, stage, PRIORS[pname], mode)
                tag = f'{stage}|{pname}|{mode}'
                ent = dict(pool=pool, best={}, greedy=[])
                for k, (a, combo) in best.items():
                    e = evaluate([cand[u] for u in combo], Dall, stage, PRIORS[pname])
                    v, lo, hi = rng(e)
                    ent['best'][k] = dict(set=list(combo), a_deg_h=v, a_deg_min=lo, a_deg_max=hi,
                                          eig_min_h=e[1.0]['eig_min'], rib_angles_h=e[1.0]['rib_angles'])
                    rows.append(dict(search=tag, method='exhaustive', k=k, set=' + '.join(combo), a_deg_h=v,
                                     a_deg_min=lo, a_deg_max=hi, eig_min_h=e[1.0]['eig_min']))
                for i, st in enumerate(steps):
                    e = evaluate([cand[u] for u in [s['added'] for s in steps[:i + 1]]], Dall, stage, PRIORS[pname])
                    v, lo, hi = rng(e)
                    ent['greedy'].append(dict(st, a_deg_h=v, a_deg_min=lo, a_deg_max=hi, eig_min_h=e[1.0]['eig_min']))
                    rows.append(dict(search=tag, method='greedy', k=i + 1, set=st['added'], a_deg_h=v,
                                     a_deg_min=lo, a_deg_max=hi, eig_min_h=e[1.0]['eig_min']))
                out[tag] = ent
                print(f'\nsearch {tag}: pool {len(pool)}')
                for k, b in ent['best'].items():
                    print(f"  best {k}: {b['a_deg_h']:.4g} [{b['a_deg_min']:.4g}..{b['a_deg_max']:.4g}]  {' + '.join(b['set'])}")
                for g in ent['greedy']:
                    print(f"  greedy +{g['added']:20s} -> {g['a_deg_h']:.4g} [{g['a_deg_min']:.4g}..{g['a_deg_max']:.4g}] eig_min {g['eig_min_h']:.3g}")
    with open(os.path.join(OUT, 'data', 'search.json'), 'w') as f:
        json.dump(out, f, indent=1)
    with open(os.path.join(OUT, 'data', 'search_table.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


if __name__ == '__main__':
    import sys
    main(tuple(sys.argv[1:]) or ('priors', 'dsweep', 'search'))
