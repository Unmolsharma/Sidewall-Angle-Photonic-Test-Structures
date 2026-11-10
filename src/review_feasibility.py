"""
Practical-feasibility analyses after the review of the 2026-11-08 draft. Two headline sets of review_analyses.py
(two-width, 50 nm slab; best set of eight, 90 nm slab), stage D, realistic priors, femwell derivatives at step h with
the strip symmetry imposed, sigma_lambda = 1 pm unless stated.

Model. Structure i is repeated N times. Copy j reads y_ij = J_i theta + D_i (s_i + r_ij) + e_ij, with theta the shared
parameters, D_i the columns of the structure's own width, thickness and angle (strip angle for a strip, rib angle for a
rib), r_ij ~ N(0, diag(sw, st, sa)^2) independent between copies, s_i a per-design offset common to all copies of design
i (zero unless stated) and e_ij the readout noise S_i. The copies of one design enter only through their mean, whose
covariance is V_i(N) = (S_i + D_i R D_i^T) / N + D_i Q D_i^T, so that
    F = sum_i J_i^T V_i(N)^-1 J_i + priors.
This is the exact Schur complement of the Fisher matrix with every offset an explicit unknown with a Gaussian prior
(checked against the explicit augmented matrix in main()).

  variation    bound at every level of the requested grid (width 0.1, 0.3, 0.6, 1 nm; angle 0.1, 0.3, 0.5, 1 deg;
               thickness 0.1, 0.3 nm)
  replication  bound vs N (1, 4, 16, 64, 256) and the N that returns to 1.5x and 2x of the shared-parameter bound (N = 1)
  layout       linear gradient along the array for blocked, interleaved and mirrored placement: bias of the shared fit
               per unit gradient and bound with the gradient fitted. Ring pitch 50 um (ASSUMED).
  calibration  common wavelength offset: differential formulation (explicit contrasts) and absolute reference of given
               accuracy; reference error independent per structure; resonance prediction in the true and instrument frame
  ablation     (separate command) femwell-data Monte Carlo of review_surrogate_check.py with the data of the 300 nm rib
               or of all other structures replaced by the surrogate

Run from src/ (main environment):  python review_feasibility.py [ablation]
Writes results_review/data/feasibility.json and feasibility_ablation.json
"""
import itertools
import json
import os
import sys

import numpy as np

import review_analyses as ra

ed = ra.ed
OUT = ra.OUT
LAM, L = ra.LAM, ra.L
WIDTHS = (0.1, 0.3, 0.6, 1.0)
ANGLES = (0.1, 0.3, 0.5, 1.0)
THICK = (0.1, 0.3)
COPIES = (1, 4, 16, 64, 256)
SCEN = {'low': (0.1, 0.1, 0.1), 'moderate': (0.3, 0.3, 0.1), 'high': (0.6, 0.5, 0.3), 'worst': (1.0, 1.0, 0.3)}   # w nm, angle deg, t nm
PITCH_MM = 0.05            # ASSUMED ring pitch
READ = {'ring_mzi': ('ring', 'mzi'), 'ring': ('ring',)}


def blocks(units, readouts=('ring', 'mzi'), stage=ra.STAGE, sl=1e-12):
    rows = ra.raw_rows(units, readouts, stage)
    s = ra.sig(sl)
    base = {p for _, _, r in rows for p in r}
    params = ['w', 't', 'a'] + sorted(base - {'w', 't', 'a'})
    out = []
    for u in units:
        rr = [r for r in rows if r[0] == u]
        akey = 'a' if u.startswith('strip') else [q for q in rr[0][2] if q.startswith('ar')][0]
        out.append(dict(uid=u, kinds=[r[1] for r in rr],
                        J=np.array([[r[2].get(p, 0.0) for p in params] for r in rr]),
                        S=np.diag([s[r[1]] ** 2 for r in rr]),
                        D=np.array([[r[2]['w'], r[2]['t'], r[2][akey]] for r in rr])))
    return params, out


def prior_rows(params, prior=None):
    prior = ra.PRIOR if prior is None else prior
    pri = []
    for key_, sp in prior.items():
        for i, p in enumerate(params):
            if p == key_ or (key_ == 'd' and p.startswith('d') and p[1:].isdigit()):
                r = np.zeros(len(params))
                r[i] = 1 / sp
                pri.append(r)
    return pri


def whiten(J, V):
    d = np.sqrt(np.diag(V))                    # scaled first: V is nearly singular at very large N
    return np.linalg.solve(np.linalg.cholesky(V / np.outer(d, d)), J / d[:, None])


def V_of(b, N, sd, sys_sd=(0.0, 0.0, 0.0)):
    R, Q = np.diag(np.array(sd, float) ** 2), np.diag(np.array(sys_sd, float) ** 2)
    return (b['S'] + b['D'] @ R @ b['D'].T) / N + b['D'] @ Q @ b['D'].T


def bound(params, bl, N=1, sd=(0.0, 0.0, 0.0), sys_sd=(0.0, 0.0, 0.0)):
    Jt = np.vstack([whiten(b['J'], V_of(b, N, sd, sys_sd)) for b in bl] + prior_rows(params))
    return ed.marginal_std(Jt, params.index('a'))


def explicit_bound(params, bl, N, sd):
    """Same bound with every offset of every copy an explicit unknown (check of the marginal form)."""
    n = len(bl) * N
    ndev = 3 * n
    rows = []
    pos = 0
    for b in bl:
        Sw = np.linalg.cholesky(b['S'])
        for _ in range(N):
            A = np.linalg.solve(Sw, b['J'])
            Dd = np.zeros((A.shape[0], ndev))
            Dd[:, pos:pos + 3] = np.linalg.solve(Sw, b['D'])
            rows.append(np.hstack([A, Dd]))
            pos += 3
    J = np.vstack(rows)
    pri = [np.hstack([r, np.zeros(ndev)]) for r in prior_rows(params)]
    dev = np.hstack([np.zeros((ndev, len(params))), np.diag(np.tile(1 / np.array(sd, float), n))])
    return ed.marginal_std(np.vstack([J] + pri + [dev]), params.index('a'))


def n_needed(params, bl, target, sd, sys_sd=(0.0, 0.0, 0.0), nmax=1e7):
    """Smallest integer N with bound <= target, or None below nmax."""
    if bound(params, bl, 1, sd, sys_sd) <= target:
        return 1
    if bound(params, bl, nmax, sd, sys_sd) > target:
        return None
    lo, hi = 0.0, np.log(nmax)
    for _ in range(70):
        mid = 0.5 * (lo + hi)
        if bound(params, bl, np.exp(mid), sd, sys_sd) > target:
            lo = mid
        else:
            hi = mid
    return int(np.ceil(np.exp(hi) - 1e-6))


def variation():
    res = {}
    for name, units in ra.SETS.items():
        res[name] = {}
        for rname, ro in READ.items():
            params, bl = blocks(units, ro)
            base = bound(params, bl)
            single = dict(width={f'{x:g}': bound(params, bl, 1, (x, 0, 0)) for x in WIDTHS},
                          angle={f'{x:g}': bound(params, bl, 1, (0, 0, x)) for x in ANGLES},
                          thickness={f'{x:g}': bound(params, bl, 1, (0, x, 0)) for x in THICK})
            grid = []
            for t, w, a in itertools.product(THICK, WIDTHS, ANGLES):
                b = bound(params, bl, 1, (w, t, a))
                grid.append(dict(width_nm=w, angle_deg=a, thickness_nm=t, bound_deg=b, ratio_to_shared=b / base))
            res[name][rname] = dict(shared=base, single=single, grid=grid)
    return res


def replication():
    res = {}
    for name, units in ra.SETS.items():
        res[name] = {}
        for rname, ro in READ.items():
            params, bl = blocks(units, ro)
            base = bound(params, bl)
            shared_N = {str(N): bound(params, bl, N) for N in COPIES}
            r = dict(shared_N1=base, target_1p5=1.5 * base, target_2=2 * base, shared_same_N=shared_N, scenarios={}, grid=[])
            for sname, (w, a, t) in SCEN.items():
                sd = (w, t, a)
                r['scenarios'][sname] = dict(
                    width_nm=w, angle_deg=a, thickness_nm=t, bound={str(N): bound(params, bl, N, sd) for N in COPIES},
                    N_1p5=n_needed(params, bl, 1.5 * base, sd), N_2=n_needed(params, bl, 2 * base, sd),
                    # part of each offset that is common to all copies of a design and cannot average
                    per_design={f'{f:g}': dict(bound_N256=bound(params, bl, 256, sd, tuple(f * x for x in sd)),
                                               floor=bound(params, bl, 1e7, sd, tuple(f * x for x in sd)),
                                               N_2=n_needed(params, bl, 2 * base, sd, tuple(f * x for x in sd)))
                                for f in (0.01, 0.03, 0.1, 0.3)})
            for t, w, a in itertools.product(THICK, WIDTHS, ANGLES):
                sd = (w, t, a)
                r['grid'].append(dict(width_nm=w, angle_deg=a, thickness_nm=t, N_1p5=n_needed(params, bl, 1.5 * base, sd),
                                      N_2=n_needed(params, bl, 2 * base, sd),
                                      **{f'bound_N{N}': bound(params, bl, N, sd) for N in COPIES}))
            res[name][rname] = r
    return res


def positions(ns, N, layout):
    """x position (mm, centred) of each copy of each design: n_s sqrt(N) columns by sqrt(N) rows, gradient along x."""
    q = int(round(np.sqrt(N)))
    columns = []                                   # design index of each column
    if layout == 'blocked':
        for i in range(ns):
            columns += [i] * q
    elif layout == 'interleaved':
        columns = list(range(ns)) * q
    elif layout == 'mirrored':
        for k in range(q):
            columns += list(range(ns)) if k % 2 == 0 else list(range(ns))[::-1]
    columns = np.array(columns)
    x = (np.arange(len(columns)) - (len(columns) - 1) / 2) * PITCH_MM
    return [np.repeat(x[columns == i], q) for i in range(ns)], float(ns * q * PITCH_MM), float(q * PITCH_MM)


def layout():
    res = dict(pitch_mm=PITCH_MM)
    for name, units in ra.SETS.items():
        params, bl = blocks(units)
        ia = params.index('a')
        res[name] = {}
        for sname, sd in (('shared', (0.0, 0.0, 0.0)), ('moderate', SCEN['moderate']), ('high', SCEN['high'])):
            sd3 = (sd[0], sd[2], sd[1])                 # (w, t, a)
            res[name][sname] = {}
            for lay in ('blocked', 'interleaved', 'mirrored'):
                res[name][sname][lay] = {}
                for N in COPIES:
                    if N == 1 and lay == 'mirrored':
                        continue
                    xs, fx, fy = positions(len(units), N, lay)
                    # bias of the shared fit per unit linear gradient
                    A = [np.linalg.inv(V_of(b, N, sd3)) for b in bl]
                    F = sum(b['J'].T @ a @ b['J'] for b, a in zip(bl, A)) + sum(np.outer(p, p) for p in prior_rows(params))
                    C = np.linalg.inv(F)
                    bias = {}
                    for k, lab in enumerate(('width_deg_per_nm_per_mm', 'thickness_deg_per_nm_per_mm', 'angle_deg_per_deg_per_mm')):
                        g = sum(b['J'].T @ a @ b['D'][:, k] * float(np.mean(x)) for b, a, x in zip(bl, A, xs))
                        bias[lab] = float((C @ g)[ia])
                    # bound with the three gradients fitted (explicit sum over copies, offsets marginalized per copy)
                    nparam = len(params)
                    Fg = np.zeros((nparam + 3, nparam + 3))
                    for b, x in zip(bl, xs):
                        W = np.linalg.inv(b['S'] + b['D'] @ np.diag(np.array(sd3, float) ** 2) @ b['D'].T)
                        n1, s1, s2 = len(x), float(np.sum(x)), float(np.sum(x ** 2))
                        Fg[:nparam, :nparam] += n1 * b['J'].T @ W @ b['J']
                        c = s1 * b['J'].T @ W @ b['D']
                        Fg[:nparam, nparam:] += c
                        Fg[nparam:, :nparam] += c.T
                        Fg[nparam:, nparam:] += s2 * b['D'].T @ W @ b['D']
                    for p in prior_rows(params):
                        pp = np.r_[p, np.zeros(3)]
                        Fg += np.outer(pp, pp)
                    scale = np.sqrt(np.diag(Fg))
                    try:
                        Cg = np.linalg.inv(Fg / np.outer(scale, scale)) / np.outer(scale, scale)
                        fitted = float(np.sqrt(Cg[ia, ia]))
                    except np.linalg.LinAlgError:
                        fitted = float('inf')
                    res[name][sname][lay][str(N)] = dict(footprint_x_mm=fx, footprint_y_mm=fy,
                                                         centroid_span_mm=float(max(np.mean(x) for x in xs) - min(np.mean(x) for x in xs)),
                                                         bound_no_gradient=bound(params, bl, N, sd3),
                                                         bound_gradient_fitted=fitted, **bias)
    return res


def cal_factor(b, kind, u):
    """Change of the observable per pm of common wavelength offset."""
    ng = ra.ng0(u)
    return {'neff': ng / LAM * 1e-12, 'ng': 0.0, 'phase': 2 * np.pi * ng * L / LAM ** 2 * 1e-12}[kind]


def differential(params, bl):
    """Explicit differential formulation: wavelength-like rows (resonance, MZI fringe) expressed in pm and differenced
    against the first structure; FSR-based n_g rows unchanged. Returns the bound."""
    J, sig2, gJ, gS = [], [], [], []
    for b in bl:
        for k, kind in enumerate(b['kinds']):
            c = cal_factor(b, kind, b['uid'])
            if kind == 'ng':
                gJ.append(b['J'][k])
                gS.append(b['S'][k, k])
            else:
                J.append(b['J'][k] / c)
                sig2.append(b['S'][k, k] / c ** 2)
    J, sig2 = np.array(J), np.array(sig2)
    n = len(sig2)
    C = np.hstack([-np.ones((n - 1, 1)), np.eye(n - 1)])
    Jd = C @ J
    Vd = C @ np.diag(sig2) @ C.T
    rows = [whiten(Jd, Vd), np.array(gJ) / np.sqrt(np.array(gS))[:, None]] + [np.array(prior_rows(params))]
    return ed.marginal_std(np.vstack(rows), params.index('a'))


def offset_bound(params, bl, sigma_c=None, stage_T=False):
    """Bound with one common offset unknown (prior sigma_c in pm, None = no prior). Returns (bound, params+cal, Jt)."""
    rows, pri = [], prior_rows(params)
    for b in bl:
        c = np.array([[cal_factor(b, kind, b['uid'])] for kind in b['kinds']])
        rows.append(whiten(np.hstack([b['J'], c]), b['S']))
    pri = [np.r_[p, 0.0] for p in pri]
    if sigma_c:
        pri.append(np.r_[np.zeros(len(params)), 1 / sigma_c])
    Jt = np.vstack(rows + pri)
    return ed.marginal_std(Jt, params.index('a')), params + ['cal'], Jt


def predict(par, Jt, device, frame):
    """1-sigma resonance prediction (pm) of a ring. frame 'true' = actual wavelength; 'instrument' = what the same
    instrument (same offset) reads."""
    C = np.linalg.pinv(Jt.T @ Jt, rcond=1e-14)
    c = ed.columns(device, ra.D[device][1.0], ra.STAGE)
    g = np.array([c[p][0] if p in c else 0.0 for p in par])
    ng = ra.ng0(device)
    res = LAM / ng * 1e12                         # pm per unit n_eff
    g = g * res
    if frame == 'instrument' and 'cal' in par:
        g[par.index('cal')] = 1.0
    return float(np.sqrt(g @ C @ g))


def calibration():
    res = {}
    devices = {'best': ['strip_TE_500', 'rib_TE_500_slab90'], 'two': ['strip_TE_500', 'rib_TE_500_slab50']}
    for name, units in ra.SETS.items():
        res[name] = {}
        for rname, ro in READ.items():
            params, bl = blocks(units, ro)
            r = dict(no_offset=bound(params, bl))
            r['offset_unknown'] = offset_bound(params, bl)[0]
            r['differential'] = differential(params, bl)
            r['reference'] = {}
            for sc in (0.1, 1, 5, 10, 20, 50, 100, 200, 500):
                v, par, Jt = offset_bound(params, bl, sc)
                r['reference'][f'{sc:g}'] = dict(bound=v, pred={d: dict(true=predict(par, Jt, d, 'true'),
                                                                         instrument=predict(par, Jt, d, 'instrument'))
                                                                 for d in devices[name]})
            _, par, Jt = offset_bound(params, bl)
            r['unknown_pred'] = {d: dict(true=predict(par, Jt, d, 'true'), instrument=predict(par, Jt, d, 'instrument')) for d in devices[name]}
            # reference accuracy at which the bound is 1.1x, 1.5x, 2x of no-offset (log interpolation on a fine grid)
            xs = np.logspace(-1, 4, 501)
            ys = np.array([offset_bound(params, bl, x)[0] for x in xs])
            r['reference_for_ratio'] = {f'{f:g}': (float(np.interp(np.log(f), np.log(ys / r['no_offset']), np.log(xs))) if (ys / r['no_offset']).max() > f else None)
                                        for f in (1.1, 1.5, 2.0)}
            r['reference_for_ratio'] = {k: (float(np.exp(v)) if v is not None else None) for k, v in r['reference_for_ratio'].items()}
            # reference error independent per structure (separate sweeps): adds to the wavelength-like rows of a structure
            r['independent_reference'] = {}
            for sc in (0.1, 1, 5):
                rows = []
                for b in bl:
                    c = np.array([cal_factor(b, kind, b['uid']) for kind in b['kinds']])
                    rows.append(whiten(b['J'], b['S'] + np.outer(c, c) * sc ** 2))
                r['independent_reference'][f'{sc:g}'] = ed.marginal_std(np.vstack(rows + prior_rows(params)), params.index('a'))
            res[name][rname] = r
    return res


def ablation():
    import montecarlo as mc
    import review_surrogate_check as sc
    D = json.load(open(os.path.join(OUT, 'surrogate_check_design.json')))
    sol = {}
    for line in open(sc.CACHE, encoding='utf-8'):
        if line.strip():
            r = json.loads(line)
            sol[r['key']] = r
    res = {}
    for name in sc.CASES:
        d = D['cases'][name]
        c = sc.case(name)
        res[name] = dict(units=c.units, crb_a_fd=c.crb_fd, crb_a_surrogate=c.crb_sur)
        rib = [k for k, u in enumerate(c.units) if u.startswith('rib_TE_300')][0]
        variants = {'all_femwell': lambda k: True, 'rib300_surrogate': lambda k: k != rib,
                    'only_rib300_femwell': lambda k: k == rib, 'all_surrogate': lambda k: False}
        for vname, use_fw in variants.items():
            rng = np.random.default_rng([sc.SEED, ord(name), 1])            # same noise draws as review_surrogate_check
            errs, status = [], {'ok': 0, 'at_bound': 0, 'failed': 0}
            for j, t in enumerate(d['truths']):
                t = np.array(t)
                pred = c.M.predict(t)
                Y = []
                for k, u in enumerate(c.units):
                    p = [q for q in d['points'] if q['truth'] == j and q['unit'] == u][0]
                    Y.append(mc.nng_of(sol[sc.key(p)]['n']) if use_fw(k) else pred[k])
                y0 = c.M.whiten(np.array(Y))
                for _ in range(sc.N_NOISE):
                    ymeas = y0 + rng.standard_normal(len(y0))
                    yprior = t[c.idx_prior] + rng.standard_normal(len(c.idx_prior)) * c.sig_prior
                    starts = [t * 0, t] + [rng.uniform(c.lo, c.hi) * 0.5 for _ in range(mc.NSTART - 2)]
                    ok = [s for s in c.fit(ymeas, yprior, starts) if s[2] > 0]
                    if not ok:
                        status['failed'] += 1
                        continue
                    chi2, x, _ = min(ok, key=lambda s: s[0])
                    status['at_bound' if np.any(np.abs(x) >= 0.999 * c.hi) else 'ok'] += 1
                    errs.append(x[c.ia] - t[c.ia])
            errs = np.array(errs)
            boot = np.random.default_rng(sc.SEED).choice(errs, size=(1000, len(errs))).std(axis=1, ddof=1)
            res[name][vname] = dict(angle_error_std=float(errs.std(ddof=1)), std_ci95=[float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))],
                                    std_over_crb_fd=float(errs.std(ddof=1) / c.crb_fd), angle_error_mean=float(errs.mean()), status=status)
            print(name, vname, res[name][vname], flush=True)
    json.dump(res, open(os.path.join(OUT, 'feasibility_ablation.json'), 'w'), indent=1)


def checks():
    """Marginal form against review_analyses.py (N = 1, one factor) and against explicit unknowns (N = 4)."""
    R = json.load(open(os.path.join(OUT, 'review_analyses.json')))['variation']
    out = {}
    for name, units in ra.SETS.items():
        params, bl = blocks(units)
        out[name] = dict(
            shared=bound(params, bl), shared_ref=R[name]['crb'],
            width_0p1=bound(params, bl, 1, (0.1, 0, 0)), width_0p1_ref=R[name]['random_effects_crb_w']['0.1'],
            thickness_0p1=bound(params, bl, 1, (0, 0.1, 0)), thickness_0p1_ref=R[name]['random_effects_crb_t']['0.1'],
            angle_0p1=bound(params, bl, 1, (0, 0, 0.1)), angle_0p1_ref=R[name]['random_effects_crb']['0.1'],
            combined_N4=bound(params, bl, 4, (0.3, 0.1, 0.3)), combined_N4_explicit=explicit_bound(params, bl, 4, (0.3, 0.1, 0.3)))
    return out


def main():
    os.makedirs(OUT, exist_ok=True)
    res = dict(checks=checks())
    print(json.dumps(res['checks'], indent=0))
    res['scenarios'] = {k: dict(width_nm=v[0], angle_deg=v[1], thickness_nm=v[2]) for k, v in SCEN.items()}
    for key, fn in (('variation', variation), ('replication', replication), ('layout', layout), ('calibration', calibration)):
        res[key] = fn()
        print('done', key, flush=True)
    json.dump(res, open(os.path.join(OUT, 'feasibility.json'), 'w'), indent=1, default=float)


if __name__ == '__main__':
    if 'ablation' in sys.argv:
        ablation()
    else:
        main()
