"""
Review analyses of the two headline sets (best one-etch set, 90 nm slab; two-width set, 50 nm slab), stage D,
realistic priors, femwell derivatives at step h with strip symmetry imposed, sigma_lambda = 1 pm unless stated.

  noise       correlated n_eff/n_g errors of one ring spectrum (rho sweep and the derived value), a common
              wavelength-calibration offset (unknown or with a prior), temperature as a shared nuisance (prior 0.1,
              0.5, 1 K) and, for contrast, an independent temperature drift per structure
  mzi         MZI fringe precision x1, x3, x10 of the ring precision (independent), and ring-only; plus an exhaustive
              ring-only search of the one-etch sets
  variation   per-structure departures from the shared parameters: bias of the shared fit for random angle offsets
              (0.25-2 deg), random width offsets (0.5-2 nm) and a linear width trend of the angle; and the CRLB of
              the shared angle when per-structure angle, width or thickness deviations are unknowns with a Gaussian
              prior (random-effects model), with the spread at which the CRLB doubles
  practical   resonance prediction error of a strip and a rib ring from the fitted parameters (full covariance)

Noise conversion as in noise_sourced.py (n_g 4.26, L 100 um, lambda 1.55 um).
Derived n_eff/n_g correlation: n_eff from one resonance lambda_1 and n_g from one FSR, lambda_1 - lambda_2, give
corr = -1/sqrt(2) for independent resonance errors (n_g = lambda^2/(L FSR) falls when lambda_1 rises); n_eff from
the mean of the two resonances gives 0.

Run from src/ (main environment):  python review_analyses.py
Writes results_review/data/review_analyses.json (log results_review/logs/review_analyses.log)
"""
import itertools
import json
import os

import numpy as np

import etchdepth_femwell_analysis as ed
import noise_sourced
import priors_search as ps

REPO = ed.REPO
OUT = os.path.join(REPO, 'results_review', 'data')
LAM = 1.55e-6
L = noise_sourced.L
SETS = {'best': ['strip_TE_400', 'strip_TE_600', 'strip_TM_400', 'strip_TM_600'] +
                 [f'rib_TE_{w}_slab90' for w in (300, 400, 500, 600)],
        'two': [f'strip_TE_{w}' for w in (400, 450, 500, 550, 600)] + ['rib_TE_300_slab50', 'rib_TE_500_slab50']}
STAGE, PRIOR = 'D', ps.PRIORS['realistic']
ROWS = ed.load()
D = ps.derivatives(ROWS)                       # strip symmetry imposed


def unit(uid):
    g, pol, W = uid.split('_')[:3]
    return (g, pol, int(W), int(uid.split('slab')[1]) if 'slab' in uid else None)


def ng0(uid):
    return float(ed.nng(ROWS, unit(uid), 'nom', 0.0)[1])


def raw_rows(units, readouts=('ring', 'mzi'), stage=STAGE, f=1.0):
    """Unwhitened rows: list of (uid, kind, {param: derivative}) with kind in neff, ng, phase."""
    k = 2 * np.pi / LAM * L
    out = []
    for u in units:
        c = ed.columns(u, D[u][f], stage)
        if 'ring' in readouts:
            out.append((u, 'neff', {p: v[0] for p, v in c.items()}))
            out.append((u, 'ng', {p: v[1] for p, v in c.items()}))
        if 'mzi' in readouts:
            out.append((u, 'phase', {p: k * v[0] for p, v in c.items()}))
    return out


def sig(sl=1e-12, phase_mult=1.0):
    s = noise_sourced.sigmas(sl)
    return dict(neff=s['SIGMA_NEFF'], ng=s['SIGMA_NG'], phase=s['SIGMA_PHASE'] * phase_mult)


def crb(units, readouts=('ring', 'mzi'), sl=1e-12, phase_mult=1.0, rho=0.0, extra_cols=None, extra_prior=None,
        prior=None, stage=STAGE, drift_K=0.0, f=1.0, return_cov=False):
    """Strip-angle CRLB with optional n_eff/n_g correlation rho per ring, extra nuisance columns
    (extra_cols: {name: fn(uid, kind) -> derivative}), priors on them (extra_prior: {name: sigma or None}), and an
    independent temperature drift per structure (adds (dn/dT sigma_T)^2 to that structure's noise)."""
    rows = raw_rows(units, readouts, stage, f)
    s = sig(sl, phase_mult)
    base = {p for _, _, r in rows for p in r}
    order = ['w', 't', 'a'] + sorted(base - {'w', 't', 'a'})
    params = order + list((extra_cols or {}).keys())
    J = np.array([[r.get(p, 0.0) for p in order] + [fn(u, kind) for fn in (extra_cols or {}).values()] for u, kind, r in rows])
    n = len(rows)
    S = np.diag([s[kind] ** 2 for _, kind, _ in rows])
    if rho:
        for i, (u, ki, _) in enumerate(rows):
            for j, (v, kj, _) in enumerate(rows):
                if u == v and ki == 'neff' and kj == 'ng':
                    S[i, j] = S[j, i] = rho * np.sqrt(S[i, i] * S[j, j])
    if drift_K:   # independent temperature excursion per structure, shared by that structure's rows
        dT = np.array([(ed.DNDT_SI * r['nsi'] + ed.DNDT_OX * r['nox']) / 1e-3 for _, _, r in rows])   # per K
        same = np.array([[u == v for v, _, _ in rows] for u, _, _ in rows])
        S = S + same * np.outer(dT, dT) * drift_K ** 2
    Lc = np.linalg.cholesky(S)
    Jw = np.linalg.solve(Lc, J)
    pri = []
    prior = PRIOR if prior is None else prior
    for key_, sp in prior.items():
        for i, p in enumerate(params):
            if p == key_ or (key_ == 'd' and p.startswith('d') and p[1:].isdigit()):
                r = np.zeros(len(params))
                r[i] = 1 / sp
                pri.append(r)
    for name, sp in (extra_prior or {}).items():
        if sp:
            r = np.zeros(len(params))
            r[params.index(name)] = 1 / sp
            pri.append(r)
    Jt = np.vstack([Jw] + pri) if pri else Jw
    try:
        a = ed.marginal_std(Jt, params.index('a'))
    except (np.linalg.LinAlgError, ValueError):
        a = float('inf')
    if not np.isfinite(np.linalg.norm(Jt[:, params.index('a')])) or np.linalg.norm(Jt[:, params.index('a')]) == 0:
        a = float('inf')
    if return_cov:
        return a, params, Jt, Jw
    return a


def noise_analysis():
    res = {}
    for name, units in SETS.items():
        base = crb(units)
        r = dict(base=base, rho={}, calibration={}, temperature_shared={}, temperature_drift_independent={})
        for rho in (-0.9, -1 / np.sqrt(2), -0.6, -0.3, 0.0, 0.3, 0.6, 0.9):
            r['rho'][f'{rho:.3f}'] = crb(units, rho=rho)
        cal = {'cal': lambda u, kind: (ng0(u) / LAM * 1e-12) if kind == 'neff' else
               (2 * np.pi * ng0(u) * L / LAM ** 2 * 1e-12 if kind == 'phase' else 0.0)}     # per pm of offset
        for label, sp in (('unknown', None), ('10 pm', 10.0), ('1 pm', 1.0)):
            r['calibration'][label] = crb(units, extra_cols=cal, extra_prior={'cal': sp})
        for TK in (0.1, 0.5, 1.0, None):
            pr = dict(PRIOR)
            if TK is None:
                pr.pop('T', None)
            else:
                pr['T'] = TK
            r['temperature_shared']['unknown' if TK is None else f'{TK} K'] = crb(units, prior=pr, stage='E')
        for dK in (0.001, 0.01, 0.1):
            r['temperature_drift_independent'][f'{dK * 1000:g} mK'] = crb(units, drift_K=dK)
        res[name] = r
        print(name, json.dumps(r, indent=None, default=float)[:900])
    return res


def mzi_analysis():
    res = {}
    for name, units in SETS.items():
        r = dict(ring_only=crb(units, readouts=('ring',)))
        for m in (1, 3, 10):
            r[f'mzi_x{m}'] = crb(units, phase_mult=m)
        res[name] = r
        print(name, r)
    # exhaustive ring-only search, one partial etch, stage D realistic, 1 pm
    strips = [u for u in D if u.startswith('strip') and all(p in D[u][1.0] for p in ('w', 't', 'a', 'nsi', 'nox'))]
    ribs = [u for u in D if u.startswith('rib') and all(p in D[u][1.0] for p in ('w', 't', 'a', 's', 'nsi', 'nox'))]
    slabs = sorted({int(u.split('slab')[1]) for u in ribs})
    best = {}
    for k in range(4, 9):
        bk = (np.inf, None)
        for s in slabs:
            pool = strips + [u for u in ribs if u.endswith(f'slab{s}')]
            for combo in itertools.combinations(pool, k):
                if not any(u.startswith('rib') for u in combo) or not any(u.startswith('strip') for u in combo):
                    continue
                v = crb(list(combo), readouts=('ring',))
                if v < bk[0]:
                    bk = (v, combo)
        best[str(k)] = dict(a_deg=bk[0], set=list(bk[1]), slab=int(bk[1][-1].split('slab')[1]))
        print('ring-only best', k, bk[0], bk[1])
    res['ring_only_search'] = dict(stage='D', priors='realistic', noise='1pm', best=best)
    two = SETS['two']
    res['ring_only_two_width_by_slab'] = {str(s): crb([u for u in two if u.startswith('strip')] + [f'rib_TE_300_slab{s}', f'rib_TE_500_slab{s}'],
                                                     readouts=('ring',)) for s in slabs}
    return res


def angle_col(uid, cols_key):
    """Per-structure angle column: the strip angle for a strip (symmetric column), the rib angle for a rib."""
    return 'a' if uid.startswith('strip') else cols_key


def variation_analysis():
    res = {}
    for name, units in SETS.items():
        a, params, Jt, Jw = crb(units, return_cov=True)
        rows = raw_rows(units)
        F = Jt.T @ Jt
        Finv = np.linalg.pinv(F)
        G = Finv @ Jw.T                                     # estimator gain: dp = G (whitened data error)
        ia = params.index('a')
        s = sig()
        # whitened data error of a unit departure of structure i's angle (deg) or width (nm)
        def dy(kind_param, i_unit):
            v = np.zeros(len(rows))
            for k, (u, kind, r) in enumerate(rows):
                if u == units[i_unit]:
                    p = 'a' if (kind_param == 'angle' and u.startswith('strip')) else (
                        [q for q in r if q.startswith('ar')][0] if kind_param == 'angle' else 'w')
                    v[k] = r[p] / s[kind]
            return v
        ga = np.array([G[ia] @ dy('angle', i) for i in range(len(units))])     # bias of theta per deg at structure i
        gw = np.array([G[ia] @ dy('width', i) for i in range(len(units))])     # per nm at structure i
        widths = np.array([int(u.split('_')[2]) for u in units], float)
        r = dict(crb=a, units=units, bias_per_deg_by_structure=ga.tolist(), bias_per_nm_by_structure=gw.tolist(),
                 bias_gain_angle_deg_per_deg=float(np.linalg.norm(ga)), bias_gain_angle_deg_per_nm=float(np.linalg.norm(gw)),
                 random_angle={}, random_width={}, linear_angle_trend={})
        for sd in (0.25, 0.5, 1.0, 2.0):
            b = sd * np.linalg.norm(ga)
            r['random_angle'][f'{sd} deg'] = dict(bias_std=b, total_rms=float(np.hypot(b, a)))
        for sw in (0.5, 1.0, 2.0):
            b = sw * np.linalg.norm(gw)
            r['random_width'][f'{sw} nm'] = dict(bias_std=b, total_rms=float(np.hypot(b, a)))
        for k in (0.25, 0.5, 1.0, 2.0):                    # deg per 100 nm of width, about the mean width of the set
            d = k * (widths - widths.mean()) / 100.0
            r['linear_angle_trend'][f'{k} deg per 100 nm'] = dict(bias=float(ga @ d))
        # random-effects model: per-structure deviations of the angle (strip angle or rib angle, deg), the width (nm) or
        # the thickness (nm) as unknowns with a Gaussian prior sigma_dev; CRLB of the shared strip angle
        def dev_cols(kind_param):
            cols = []
            for i in range(len(units)):
                v = np.zeros(Jt.shape[0])
                for k, (u, kind, rr) in enumerate(rows):
                    if u == units[i]:
                        q = ('a' if u.startswith('strip') else [x for x in rr if x.startswith('ar')][0]) if kind_param == 'angle' else kind_param
                        v[k] = rr[q] / s[kind]
                cols.append(v)
            return np.array(cols).T
        for kp, levels in (('angle', (0.0001, 0.001, 0.01, 0.1, 0.25, 0.5, 1.0, 2.0)), ('w', (0.001, 0.01, 0.1, 0.5, 1.0, 2.0)),
                           ('t', (0.001, 0.01, 0.1, 0.5, 1.0))):
            cols = dev_cols(kp)
            out = {}
            for sd in levels:
                Jre = np.hstack([Jt, cols])
                Jre = np.vstack([Jre, np.hstack([np.zeros((len(units), Jt.shape[1])), np.eye(len(units)) / sd])])
                out[f'{sd:g}'] = ed.marginal_std(Jre, ia)
            r['random_effects_crb' if kp == 'angle' else f'random_effects_crb_{kp}'] = out
            # spread at which the CRLB doubles (log interpolation over the levels)
            xs, ys = np.log(np.array(levels)), np.log(np.array([out[f'{x:g}'] for x in levels]))
            r[f'doubling_spread_{kp}'] = float(np.exp(np.interp(np.log(2 * a), ys, xs))) if ys[-1] > np.log(2 * a) > ys[0] else None
        res[name] = r
        print(name, {k: v for k, v in r.items() if k not in ('units',)})
    return res


def practical_analysis():
    """Resonance prediction error d lambda = lambda sigma_neff_pred / n_g of a strip ring (TE 500) and a rib ring (TE
    500 at the set's slab) from the fitted parameters, at 1 pm and in the correlated-noise variants."""
    res = {}
    devices = {'best': ['strip_TE_500', 'rib_TE_500_slab90'], 'two': ['strip_TE_500', 'rib_TE_500_slab50']}
    variants = {'1 pm': {}, '1 pm, rho = -0.707': dict(rho=-1 / np.sqrt(2)),
                '1 pm, unknown calibration offset': 'cal', '0.1 pm': dict(sl=1e-13)}
    for name, units in SETS.items():
        r = {}
        for vname, kw in variants.items():
            if kw == 'cal':
                cal = {'cal': lambda u, kind: (ng0(u) / LAM * 1e-12) if kind == 'neff' else
                       (2 * np.pi * ng0(u) * L / LAM ** 2 * 1e-12 if kind == 'phase' else 0.0)}
                a, params, Jt, _ = crb(units, extra_cols=cal, extra_prior={'cal': None}, return_cov=True)
            else:
                a, params, Jt, _ = crb(units, return_cov=True, **kw)
            C = np.linalg.pinv(Jt.T @ Jt)
            dev = {}
            for d in devices[name]:
                c = ed.columns(d, D[d][1.0], STAGE)
                g = np.array([c.get(p, [0.0])[0] if p in c else 0.0 for p in params])
                sn = float(np.sqrt(g @ C @ g))
                ga = np.zeros(len(params))
                ga[params.index('a')] = g[params.index('a')]
                arp = [p for p in params if p.startswith('ar')]
                if d.startswith('rib') and arp:
                    ga[params.index(arp[0])] = g[params.index(arp[0])]
                dev[d] = dict(sigma_resonance_pm=float(LAM * sn / ng0(d) * 1e12),
                              angle_terms_only_pm=float(LAM * np.sqrt(ga @ C @ ga) / ng0(d) * 1e12))
            r[vname] = dict(crb_angle=a, devices=dev)
        res[name] = r
        print(name, json.dumps(r, default=float)[:700])
    return res


def main():
    os.makedirs(OUT, exist_ok=True)
    base_check = {k: crb(v) for k, v in SETS.items()}
    ref = {k: ps.evaluate([(u, ('ring', 'mzi')) for u in v], D, 'D', PRIOR, factors=(1.0,))[1.0]['a'] for k, v in SETS.items()}
    print('base check (this module / priors_search at 1 pm):', base_check, ref)
    res = dict(base_check=dict(this_module=base_check), noise=noise_analysis(), mzi=mzi_analysis(),
               variation=variation_analysis(), practical=practical_analysis())
    with open(os.path.join(OUT, 'review_analyses.json'), 'w') as f:
        json.dump(res, f, indent=1, default=float)


if __name__ == '__main__':
    import fullvector_sensitivity as fv
    s = noise_sourced.sigmas(1e-12)
    fv.SIGMA_NEFF, fv.SIGMA_NG, fv.SIGMA_PHASE = s['SIGMA_NEFF'], s['SIGMA_NG'], s['SIGMA_PHASE']   # for the reference
    main()
