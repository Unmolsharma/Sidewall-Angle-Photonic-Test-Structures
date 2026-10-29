"""
Monte Carlo fitting check of the strip sidewall-angle Cramér-Rao bound, with strip symmetry imposed.

Run from src/ (main environment):
    python montecarlo.py plan     cases, surrogate domains and design points -> data/mc_design.json, estimate
    (femwell venv) montecarlo_batch.py   solves the design points
    python montecarlo.py fit      surrogates + held-out accuracy -> data/surrogate_accuracy.{json,csv}
    python montecarlo.py mc       noisy measurement sets and fits -> data/mc_results.json, mc_fits.csv
Reads   results_montecarlo/sym_search/data/{search.json, slab_pairs.json} (sym_search.py), femwell caches,
        results_montecarlo/cache/solves_r12.jsonl
Writes  results_montecarlo/data/...

Cases (stage D: w, t, a, d, ar, n_Si, n_ox; realistic priors: d 2 nm, n_Si 1e-4, n_ox 1e-3):
  a  best one-etch set of 8 structures (sym_search, D realistic)
  b  strip TE (5 widths) + rib 300 + 500 nm at the best slab for a two-width set (D realistic)
  c  strip TE (5 widths) + rib 500 nm at the same slab (expected: width/angle valley)
Each structure is read out by a ring (n_eff, n_g) and an MZI (phase); noise placeholders as before.

Forward model: per structure a smooth surrogate (full quadratic + pure cubes) fitted to femwell r12 solves
(Latin hypercube in the structure's own inputs, plus the existing derivative solves that fall inside the box).
Strip symmetry: a strip sees the angle only through its mid-height width, w_mid = W + w - t tan(a), so strip
surrogates are in (w_mid, t, n_Si, n_ox) and are solved with vertical walls. Ribs are in (w, t, ar, slab,
n_Si, n_ox) with slab = S + t - d. Surrogate boxes are +-K_BOX marginal CRB standard deviations of the cases
(a, b) mapped onto each structure's inputs, with floors; case c gets a wide rib box (valley).

Monte Carlo: truth = nominal; each set draws measurement noise and prior centres (prior pseudo-
observations = truth + sigma_p N(0,1)); bounded nonlinear least squares from NSTART starts (nominal + random
in the box). Per case: spread and bias of the fitted strip angle vs the CRB (finite-difference Jacobian and
surrogate Jacobian), failed fits, fits at the box edge, and ambiguous fits (another start within
delta-chi^2 < 1 with a strip angle differing by more than min(CRB, 1 deg)).
"""
import csv
import json
import os
import sys

import numpy as np
from scipy.optimize import least_squares
from scipy.stats import qmc

import etchdepth_femwell_analysis as ed
import fullvector_sensitivity as fv
import priors_search as ps

REPO = ed.REPO
OUT = os.path.join(REPO, 'results_montecarlo', 'data')
SYM = os.path.join(REPO, 'results_montecarlo', 'sym_search', 'data')
MC_CACHE = os.path.join(REPO, 'results_montecarlo', 'cache', 'solves_r12.jsonl')
T0, N_SI, N_OX = 220.0, 3.48, 1.44
STAGE, PRIOR = 'D', 'realistic'
STRIP5 = [f'strip_TE_{w}' for w in (400, 450, 500, 550, 600)]
K_BOX = 6.0
FLOOR = dict(w=1.0, t=1.0, ar=0.25, slab=1.0, nsi=1e-4, nox=1e-4)     # half-widths: nm, deg, absolute index
WIDE = dict(w=12.0, a=3.0, ar=5.0)                                     # case c global box (nm, deg, deg)
N_TRAIN = dict(strip=20, rib=30)
N_TEST = dict(strip=6, rib=8)
NSETS, NSTART, SEED = 2000, 5, 20261027
KPH = 2 * np.pi / fv.WL0 * fv.L_MZI
SIG_PH_NEFF = fv.SIGMA_PHASE / KPH        # phase noise as an n_eff noise


_ROWS = None


def rows():
    global _ROWS
    if _ROWS is None:
        _ROWS = ed.load()
    return _ROWS


def mc_key(d, p):
    # same format as montecarlo_batch.key_of (that module needs the femwell environment)
    return (f"mc|{d['geom']}|{d['pol']}|{d['W']}|{d['S']}|{p['w']:.4f}|{p['t']:.4f}|{p['a']:.5f}|"
            f"{p['slab'] if p['slab'] is not None else 'None'}|{p['nsi']:.7f}|{p['nox']:.7f}")


def unit_of(uid):
    p = uid.split('_')
    return (p[0], p[1], int(p[2]), int(p[3].replace('slab', '')) if len(p) > 3 else None)


def inputs_names(geom):
    return ('w_mid', 't', 'nsi', 'nox') if geom == 'strip' else ('w', 't', 'ar', 'slab', 'nsi', 'nox')


# ---------------------------------------------------------------------------------------------- cases
def cases():
    S = json.load(open(os.path.join(SYM, 'search.json')))
    P = json.load(open(os.path.join(SYM, 'slab_pairs.json')))[f'{STAGE}|{PRIOR}']
    pairs = {v['slab']: v['a_deg_h'] for k, v in P.items() if k.startswith('rib_TE_300') and 'rib_TE_500' in k}
    slab = min(pairs, key=pairs.get)
    return {'a': S[f'{STAGE}|{PRIOR}']['one_etch_overall']['8']['set'],
            'b': STRIP5 + [f'rib_TE_300_slab{slab}', f'rib_TE_500_slab{slab}'],
            'c': STRIP5 + [f'rib_TE_500_slab{slab}']}, slab, pairs


def fd_model(units, Dall):
    items = [(u, ('ring', 'mzi')) for u in units]
    return ps.jacobian(items, Dall, 1.0, STAGE, ps.PRIORS[PRIOR])


def input_map(uid, params):
    """Linearised map (rows: structure inputs, columns: global params) at the nominal point."""
    geom, pol, W, S = unit_of(uid)
    names = inputs_names(geom)
    M = np.zeros((len(names), len(params)))
    col = {p: i for i, p in enumerate(params)}
    if geom == 'strip':
        M[0, col['w']], M[0, col['a']], M[1, col['t']] = 1.0, -T0 * np.pi / 180, 1.0
        M[2, col['nsi']], M[3, col['nox']] = 1e-3, 1e-3
    else:
        M[0, col['w']], M[1, col['t']], M[2, col[f'ar{S}']] = 1.0, 1.0, 1.0
        M[3, col['t']], M[3, col[f'd{S}']] = 1.0, -1.0
        M[4, col['nsi']], M[5, col['nox']] = 1e-3, 1e-3
    return names, M


def floor_of(name):
    return FLOOR['w'] if name in ('w', 'w_mid') else FLOOR[name]


# ---------------------------------------------------------------------------------------------- plan
def plan():
    Dall = ps.derivatives(rows())
    C, slab, pairs = cases()
    boxes, theta_box, crb_fd = {}, {}, {}
    for name in ('a', 'b'):
        params, J = fd_model(C[name], Dall)
        cov = np.linalg.pinv(J.T @ J)
        sd = np.sqrt(np.diag(cov))
        theta_box[name] = {p: float(K_BOX * s) for p, s in zip(params, sd)}
        crb_fd[name] = dict(params=params, sd=dict(zip(params, sd.tolist())),
                            a_marginal=ed.marginal_std(J, params.index('a')))
        for u in C[name]:
            names, M = input_map(u, params)
            hw = K_BOX * np.sqrt(np.einsum('ij,jk,ik->i', M, cov, M))
            hw = np.maximum(hw, [floor_of(n) for n in names])
            old = boxes.get(u, np.zeros(len(names)))
            boxes[u] = np.maximum(old, hw)
    # case c: wide box on w, a, ar; the other parameters as case b
    params, J = fd_model(C['c'], Dall)
    tb = dict(theta_box['b'])
    tb.update(w=WIDE['w'], a=WIDE['a'], **{f'ar{slab}': WIDE['ar']})
    theta_box['c'] = {p: tb[p] for p in params}
    crb_fd['c'] = dict(params=params, a_marginal=ed.marginal_std(J, params.index('a')))
    wide = {}
    for u in C['c']:
        if u.startswith('rib'):
            names, M = input_map(u, params)
            wide[u] = np.maximum(np.abs(M) @ np.array([theta_box['c'][p] for p in params]), [floor_of(n) for n in names])
    designs = []
    rng_seed = SEED
    for tag, bx in (('narrow', boxes), ('wide', wide)):
        for u, hw in bx.items():
            geom, pol, W, S = unit_of(u)
            n_tr, n_te = N_TRAIN[geom], N_TEST[geom]
            z = qmc.LatinHypercube(d=len(hw), seed=rng_seed).random(n_tr + n_te) * 2 - 1
            rng_seed += 1
            pts = [absolute(geom, W, S, zi * hw) | dict(role='train' if i < n_tr else 'test', z=zi.tolist())
                   for i, zi in enumerate(z)]
            designs.append(dict(uid=u, box=tag, geom=geom, pol=pol, W=W, S=S, names=list(inputs_names(geom)),
                                halfwidth=hw.tolist(), guess=nominal_guess(u), points=pts))
    os.makedirs(OUT, exist_ok=True)
    out = dict(cases=C, two_width_slab=slab, two_width_by_slab=pairs, theta_box=theta_box, crb_fd=crb_fd,
               K_BOX=K_BOX, FLOOR=FLOOR, WIDE=WIDE, designs=designs)
    with open(os.path.join(OUT, 'mc_design.json'), 'w') as f:
        json.dump(out, f, indent=1)
    estimate(designs)


def absolute(geom, W, S, x):
    """Structure inputs (offsets) -> absolute solve values."""
    if geom == 'strip':
        return dict(w=W + x[0], t=T0 + x[1], a=0.0, slab=None, nsi=N_SI + x[2], nox=N_OX + x[3])
    return dict(w=W + x[0], t=T0 + x[1], a=x[2], slab=S + x[3], nsi=N_SI + x[4], nox=N_OX + x[5])


def nominal_guess(uid):
    geom, pol, W, S = unit_of(uid)
    return rows()[f'{geom}|{pol}|{W}|{S}|oxide|r12|nom|+0.0000']['n'][0]


def estimate(designs):
    """Solve count and runtime from the per-solve times already in the caches (median and mean)."""
    sec = {}
    for r in rows().values():
        g = r['key'].split('|')[0]
        sec.setdefault(g, []).append(r['seconds'])
    n = {}
    for d in designs:
        n[d['geom']] = n.get(d['geom'], 0) + len(d['points'])
    print('designs:', [(d['uid'], d['box'], [round(h, 4) for h in d['halfwidth']]) for d in designs])
    tot_med = tot_mean = 0.0
    for g, k in n.items():
        med, mean = float(np.median(sec[g])), float(np.mean(sec[g]))
        tot_med += k * med
        tot_mean += k * mean
        print(f'{g}: {k} solves, per-solve median {med:.1f} s, mean {mean:.1f} s (one worker)')
    print(f'total {sum(n.values())} solves; worker-hours {tot_med/3600:.2f} (median) .. {tot_mean/3600:.2f} (mean); '
          f'with 4 workers ~{tot_med/4/3600:.2f} .. {tot_mean/4/3600:.2f} h if per-solve times hold under load')


# ---------------------------------------------------------------------------------------------- surrogates
_PAIRS = {}


def _pairs(d):
    if d not in _PAIRS:
        _PAIRS[d] = np.triu_indices(d)
    return _PAIRS[d]


def basis(z):
    """z: (..., d) -> (..., nb): 1, z_i, z_i z_j (i <= j), z_i^3."""
    I, J = _pairs(z.shape[-1])
    return np.concatenate([np.ones(z.shape[:-1] + (1,)), z, z[..., I] * z[..., J], z ** 3], axis=-1)


def dbasis(z):
    """d basis / d z: (d,) -> (nb, d)."""
    d = len(z)
    I, J = _pairs(d)
    k = np.arange(len(I))
    Dp = np.zeros((len(I), d))
    Dp[k, I] = z[J]
    Dp[k, J] += z[I]          # (k, J) pairs are unique, so no repeated-index problem; i == j gives 2 z_i
    return np.vstack([np.zeros((1, d)), np.eye(d), Dp, np.diag(3 * z ** 2)])


def basis_masks(d):
    """Candidate term subsets of basis(); the fit keeps the one with the smallest leave-one-out error on the
    training points. 'full' = all terms; 'quad' = no cubes; 'reduced' = linear in every input, quadratic only
    among the geometry inputs (w, t, ar, slab for ribs; w_mid, t for strips), cubes only in w and ar (w_mid)."""
    I, J = _pairs(d)
    geo = range(4) if d == 6 else range(2)
    cub = (0, 2) if d == 6 else (0,)
    return {'full': np.ones(1 + d + len(I) + d, bool),
            'quad': np.r_[np.ones(1 + d + len(I), bool), np.zeros(d, bool)],
            'reduced': np.r_[np.ones(1 + d, bool), [i in geo and j in geo for i, j in zip(I, J)],
                             [i in cub for i in range(d)]]}


def nng_of(n):
    return np.array([n[0], n[0] - fv.WL0 * (n[1] - n[2]) / (2 * fb_DLAM)])


fb_DLAM = 0.005


def axis_points(uid, hw):
    """Existing derivative solves (nominal and single-parameter steps) that fall inside the box, as offsets."""
    geom, pol, W, S = unit_of(uid)
    names = inputs_names(geom)
    which = {'w': 'w_mid' if geom == 'strip' else 'w', 't': 't', 'a': None if geom == 'strip' else 'ar',
             's': 'slab', 'nsi': 'nsi', 'nox': 'nox', 'nom': None}
    out = []
    for p in ('nom', 'w', 't', 'a', 's', 'nsi', 'nox'):
        steps = [0.0] if p == 'nom' else [s * f * ed.H[p] for f in ed.FACTORS for s in (1, -1)]
        for d in steps:
            r = rows().get(f'{geom}|{pol}|{W}|{S}|oxide|r12|{p}|{d:+.4f}')
            if r is None or (p != 'nom' and which[p] is None):
                continue
            x = np.zeros(len(names))
            if p != 'nom':
                x[names.index(which[p])] = d
            if np.all(np.abs(x) <= np.array(hw) * (1 + 1e-9)):
                out.append((x, nng_of(r['n'])))
    return out


def fit_surrogates():
    design = json.load(open(os.path.join(OUT, 'mc_design.json')))
    mc = {}
    with open(MC_CACHE, encoding='utf-8') as f:
        for line in f:
            if line.strip():
                r = json.loads(line)
                mc[r['key']] = r
    sur, acc = {}, []
    for d in design['designs']:
        hw = np.array(d['halfwidth'])
        X, Y, Xt, Yt = [], [], [], []
        for p in d['points']:
            r = mc.get(mc_key(d, p))
            if r is None:
                continue
            x = np.array(p['z']) * hw
            (X if p['role'] == 'train' else Xt).append(x)
            (Y if p['role'] == 'train' else Yt).append(nng_of(r['n']))
        ax = axis_points(d['uid'], hw)
        X += [x for x, _ in ax]
        Y += [y for _, y in ax]
        X, Y, Xt, Yt = map(np.array, (X, Y, Xt, Yt))
        Bf = basis(X / hw)
        wsig = np.array([SIG_PH_NEFF, fv.SIGMA_NG])
        loo = {}
        for bname, mk in basis_masks(X.shape[1]).items():     # choose the basis on training points only
            B = Bf[:, mk]
            c, *_ = np.linalg.lstsq(B, Y, rcond=None)
            h = np.einsum('ij,ji->i', B, np.linalg.pinv(B))
            e = (B @ c - Y) / (1 - h)[:, None] / wsig
            loo[bname] = float(np.sqrt(np.mean(e ** 2)))
        bname = min(loo, key=loo.get)
        mk = basis_masks(X.shape[1])[bname]
        coef = np.zeros((len(mk), Y.shape[1]))
        coef[mk], *_ = np.linalg.lstsq(Bf[:, mk], Y, rcond=None)
        name = f"{d['uid']}|{d['box']}"
        sur[name] = dict(uid=d['uid'], box=d['box'], halfwidth=hw.tolist(), coef=coef.tolist(), names=d['names'],
                         basis=bname)
        e_tr = Bf @ coef - Y
        e_te = basis(Xt / hw) @ coef - Yt
        row = dict(surrogate=name, n_train=len(X), n_train_axis=len(ax), n_test=len(Xt), basis=bname,
                   n_terms=int(mk.sum()), **{f'loo_{k}': v for k, v in loo.items()},
                   neff_rms_train=float(np.sqrt(np.mean(e_tr[:, 0] ** 2))),
                   neff_rms_test=float(np.sqrt(np.mean(e_te[:, 0] ** 2))), neff_max_test=float(np.max(np.abs(e_te[:, 0]))),
                   ng_rms_test=float(np.sqrt(np.mean(e_te[:, 1] ** 2))), ng_max_test=float(np.max(np.abs(e_te[:, 1]))),
                   neff_range=float(np.ptp(np.r_[Y[:, 0], Yt[:, 0]])), ng_range=float(np.ptp(np.r_[Y[:, 1], Yt[:, 1]])))
        row['neff_rms_test_over_sigma_phase'] = row['neff_rms_test'] / SIG_PH_NEFF
        row['neff_max_test_over_sigma_phase'] = row['neff_max_test'] / SIG_PH_NEFF
        row['neff_rms_test_over_sigma_neff'] = row['neff_rms_test'] / fv.SIGMA_NEFF
        row['ng_rms_test_over_sigma_ng'] = row['ng_rms_test'] / fv.SIGMA_NG
        acc.append(row)
        print(f"{name:28s} {bname:7s} train {len(X):3d} (axis {len(ax)}) test {len(Xt):2d}: n_eff rms {row['neff_rms_test']:.2e} "
              f"(= {row['neff_rms_test_over_sigma_phase']:.2f} sigma_phase), max {row['neff_max_test']:.2e}; "
              f"n_g rms {row['ng_rms_test']:.2e} (= {row['ng_rms_test_over_sigma_ng']:.3f} sigma_ng)")
    with open(os.path.join(OUT, 'surrogates.json'), 'w') as f:
        json.dump(sur, f, indent=1)
    with open(os.path.join(OUT, 'surrogate_accuracy.json'), 'w') as f:
        json.dump(dict(sigma_phase_as_neff=SIG_PH_NEFF, sigma_neff=fv.SIGMA_NEFF, sigma_ng=fv.SIGMA_NG, rows=acc), f, indent=1)
    with open(os.path.join(OUT, 'surrogate_accuracy.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(acc[0].keys()))
        w.writeheader()
        w.writerows(acc)


# ---------------------------------------------------------------------------------------------- forward model
class Model:
    def __init__(self, units, params, sur, box):
        self.units, self.params = units, params
        self.col = {p: i for i, p in enumerate(params)}
        self.S = [sur[f'{u}|{box if u.startswith("rib") else "narrow"}'] for u in units]
        self.coef = [np.array(s['coef']) for s in self.S]
        self.hw = [np.array(s['halfwidth']) for s in self.S]
        prior = ps.PRIORS[PRIOR]
        self.prior = [(self.col[p], prior['d'] if p.startswith('d') else prior[p]) for p in params
                      if p.startswith('d') or p in ('nsi', 'nox')]

    def x_and_dx(self, uid, th):
        geom, pol, W, S = unit_of(uid)
        c, P = self.col, len(self.params)
        if geom == 'strip':
            ta = np.tan(np.radians(th[c['a']]))
            x = np.array([th[c['w']] - (T0 + th[c['t']]) * ta, th[c['t']], 1e-3 * th[c['nsi']], 1e-3 * th[c['nox']]])
            dx = np.zeros((4, P))
            dx[0, c['w']], dx[0, c['t']] = 1.0, -ta
            dx[0, c['a']] = -(T0 + th[c['t']]) * np.radians(1.0) / np.cos(np.radians(th[c['a']])) ** 2
            dx[1, c['t']], dx[2, c['nsi']], dx[3, c['nox']] = 1.0, 1e-3, 1e-3
        else:
            x = np.array([th[c['w']], th[c['t']], th[c[f'ar{S}']], th[c['t']] - th[c[f'd{S}']],
                          1e-3 * th[c['nsi']], 1e-3 * th[c['nox']]])
            dx = np.zeros((6, P))
            dx[0, c['w']], dx[1, c['t']], dx[2, c[f'ar{S}']] = 1.0, 1.0, 1.0
            dx[3, c['t']], dx[3, c[f'd{S}']] = 1.0, -1.0
            dx[4, c['nsi']], dx[5, c['nox']] = 1e-3, 1e-3
        return x, dx

    def predict(self, th, jac=False):
        """-> observables per structure [n_eff, n_g] (k, 2) and optionally d/dtheta (k, 2, P)."""
        Y, G = [], []
        for u, coef, hw in zip(self.units, self.coef, self.hw):
            x, dx = self.x_and_dx(u, th)
            z = x / hw
            Y.append(basis(z) @ coef)
            if jac:
                G.append((dbasis(z).T @ coef).T / hw @ dx)   # (2, d) / hw -> (2, d) @ (d, P)
        return (np.array(Y), np.array(G)) if jac else np.array(Y)

    def whiten(self, Y):
        """observables -> measurement vector / sigma: per structure n_eff (ring), n_g (ring), phase (MZI)."""
        return np.concatenate([Y[..., 0] / fv.SIGMA_NEFF, Y[..., 1] / fv.SIGMA_NG, KPH * Y[..., 0] / fv.SIGMA_PHASE], axis=-1)

    def jac_rows(self, th):
        _, G = self.predict(th, jac=True)
        Jm = np.concatenate([G[:, 0, :] / fv.SIGMA_NEFF, G[:, 1, :] / fv.SIGMA_NG, KPH * G[:, 0, :] / fv.SIGMA_PHASE])
        Jp = np.zeros((len(self.prior), len(self.params)))
        for k, (i, s) in enumerate(self.prior):
            Jp[k, i] = 1 / s
        return np.vstack([Jm, Jp])


NPROC, CHUNK = 2, 100


class Case:
    """Everything one Monte Carlo case needs; rebuilt cheaply inside each worker process."""

    def __init__(self, name, out=None):
        out = out or OUT
        design = json.load(open(os.path.join(out, 'mc_design.json')))
        sur = json.load(open(os.path.join(out, 'surrogates.json')))
        self.name, self.units = name, design['cases'][name]
        self.params, self.Jfd = fd_model(self.units, ps.derivatives(rows()))
        self.M = Model(self.units, self.params, sur, 'wide' if name == 'c' else 'narrow')
        self.hi = np.array([design['theta_box'][name][p] for p in self.params])
        self.lo = -self.hi
        self.th0 = np.zeros(len(self.params))
        self.ia = self.params.index('a')
        self.J0 = self.M.jac_rows(self.th0)
        self.crb_fd, self.crb_sur = ed.marginal_std(self.Jfd, self.ia), ed.marginal_std(self.J0, self.ia)
        self.Y0 = self.M.whiten(self.M.predict(self.th0))
        self.sig_prior = np.array([s for _, s in self.M.prior])
        self.idx_prior = [i for i, _ in self.M.prior]
        self.amb_thr = min(self.crb_sur, 1.0)

    def resid(self, th, ymeas, yprior):
        return np.concatenate([self.M.whiten(self.M.predict(th)) - ymeas, (th[self.idx_prior] - yprior) / self.sig_prior])

    def jac(self, th, ymeas, yprior):
        return self.M.jac_rows(th)

    def fit(self, ymeas, yprior, starts):
        sols = []
        for s0 in starts:
            r = least_squares(self.resid, s0, jac=self.jac, bounds=(self.lo, self.hi), args=(ymeas, yprior),
                              method='trf', x_scale=self.hi, xtol=1e-10, ftol=1e-10, gtol=1e-10, max_nfev=200)
            sols.append((2 * r.cost, r.x, r.status))
        return sols

    def one_set(self, rng):
        ymeas = self.Y0 + rng.standard_normal(len(self.Y0))
        yprior = rng.standard_normal(len(self.idx_prior)) * self.sig_prior
        shrink = 1.0 if self.name == 'c' else 0.5
        starts = [self.th0] + [rng.uniform(self.lo, self.hi) * shrink for _ in range(NSTART - 1)]
        ok = [s for s in self.fit(ymeas, yprior, starts) if s[2] > 0]
        if not ok:
            return dict(status='failed')
        chi2, x, _ = min(ok, key=lambda s: s[0])
        near = [s for s in ok if s[0] <= chi2 + 1.0 and abs(s[1][self.ia] - x[self.ia]) > self.amb_thr]
        at_bound = bool(np.any(np.abs(x) >= 0.999 * self.hi))
        return dict(status='ambiguous' if near else ('at_bound' if at_bound else 'ok'), chi2=float(chi2),
                    a_hat=float(x[self.ia]), at_bound=at_bound, n_converged=len(ok),
                    **{f'hat_{p}': float(v) for p, v in zip(self.params, x)})


def _chunk(args):
    name, ci, k0, k1, out = args
    c = Case(name, out)
    rng = np.random.default_rng([SEED, ord(name), ci])     # reproducible, independent per chunk
    return [dict(set=k, **c.one_set(rng)) for k in range(k0, k1)]


def run_mc():
    from multiprocessing import Pool
    design = json.load(open(os.path.join(OUT, 'mc_design.json')))
    jobs = [(name, ci, k0, min(k0 + CHUNK, NSETS), OUT) for name in design['cases']
            for ci, k0 in enumerate(range(0, NSETS, CHUNK))]
    with Pool(NPROC) as pool:
        parts = pool.map(_chunk, jobs)
    by_case = {}
    for (name, *_), part in zip(jobs, parts):
        by_case.setdefault(name, []).extend(part)
    results, fits = {}, []
    rng = np.random.default_rng(SEED)
    for name, out in by_case.items():
        c = Case(name)
        params = c.params
        a_all = np.array([o['a_hat'] for o in out if 'a_hat' in o])
        a_ok = np.array([o['a_hat'] for o in out if o['status'] == 'ok'])
        # model-error check: noise-free data = the actual femwell nominal solves, priors at truth
        y_true = c.M.whiten(np.array([ed.nng(rows(), unit_of(u), 'nom', 0.0) for u in c.units]))
        a_modelerr = float(c.fit(y_true, np.zeros(len(c.idx_prior)), [c.th0])[0][1][c.ia])
        boot = rng.choice(a_all, size=(1000, len(a_all)), replace=True).std(axis=1, ddof=1)
        res = dict(units=c.units, params=params, theta_box=design['theta_box'][name],
                   surrogate_box='wide' if name == 'c' else 'narrow',
                   crb_a_fd=c.crb_fd, crb_a_surrogate=c.crb_sur, n_sets=len(out), n_starts=NSTART,
                   n_ok=int(sum(o['status'] == 'ok' for o in out)),
                   n_failed=int(sum(o['status'] == 'failed' for o in out)),
                   n_at_bound=int(sum(o['status'] == 'at_bound' for o in out)),
                   n_ambiguous=int(sum(o['status'] == 'ambiguous' for o in out)),
                   ambiguity_threshold_deg=c.amb_thr,
                   a_std_all=float(a_all.std(ddof=1)),
                   a_std_all_ci95=[float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))],
                   a_robust_std_all=float(1.4826 * np.median(np.abs(a_all - np.median(a_all)))),
                   a_bias_all=float(a_all.mean()), a_bias_se=float(a_all.std(ddof=1) / np.sqrt(len(a_all))),
                   a_std_ok=float(a_ok.std(ddof=1)) if len(a_ok) > 1 else None,
                   a_bias_ok=float(a_ok.mean()) if len(a_ok) else None,
                   a_p05_p95=[float(np.percentile(a_all, 5)), float(np.percentile(a_all, 95))],
                   frac_abs_a_over_3crb=float(np.mean(np.abs(a_all) > 3 * c.crb_sur)) if np.isfinite(c.crb_sur) else None,
                   chi2_mean=float(np.mean([o['chi2'] for o in out if 'chi2' in o])),
                   dof=int(len(c.Y0) + len(c.idx_prior) - len(params)),
                   a_fit_to_femwell_nominal=a_modelerr,
                   params_std_mc={p: float(np.std([o[f'hat_{p}'] for o in out if f'hat_{p}' in o], ddof=1)) for p in params},
                   params_crb_surrogate={p: ed.marginal_std(c.J0, i) for i, p in enumerate(params)},
                   params_crb_fd={p: ed.marginal_std(c.Jfd, i) for i, p in enumerate(params)})
        res['std_over_crb_surrogate'] = res['a_std_all'] / c.crb_sur if np.isfinite(c.crb_sur) else None
        res['std_over_crb_fd'] = res['a_std_all'] / c.crb_fd if np.isfinite(c.crb_fd) else None
        results[name] = res
        for o in out:
            fits.append(dict(case=name, **{k: o.get(k) for k in ['set', 'status', 'chi2', 'a_hat', 'n_converged'] +
                                           [f'hat_{p}' for p in params]}))
        print(f"case {name}: CRB fd {c.crb_fd:.4g}, surrogate {c.crb_sur:.4g}; MC std {res['a_std_all']:.4g} "
              f"[{res['a_std_all_ci95'][0]:.4g}, {res['a_std_all_ci95'][1]:.4g}], bias {res['a_bias_all']:+.4g} "
              f"+- {res['a_bias_se']:.2g}; ok {res['n_ok']}, failed {res['n_failed']}, at bound {res['n_at_bound']}, "
              f"ambiguous {res['n_ambiguous']}; chi2 mean {res['chi2_mean']:.2f} (dof {res['dof']}); "
              f"fit to femwell nominal: a = {a_modelerr:+.4g}")
    with open(os.path.join(OUT, 'mc_results.json'), 'w') as f:
        json.dump(results, f, indent=1)
    head = ['case', 'set', 'status', 'chi2', 'a_hat', 'n_converged']
    keys = head + sorted({k for r in fits for k in r} - set(head))
    with open(os.path.join(OUT, 'mc_fits.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(fits)


def valley():
    """Width/angle relation of the fitted values per case (from data/mc_fits.csv): correlation and slope of
    fitted width against fitted strip angle. Writes data/mc_valley.json."""
    out = {}
    with open(os.path.join(OUT, 'mc_fits.csv'), newline='') as f:
        fits = list(csv.DictReader(f))
    for name in sorted({r['case'] for r in fits}):
        r = [x for x in fits if x['case'] == name and x['a_hat'] and x['hat_w']]
        a, w = np.array([float(x['a_hat']) for x in r]), np.array([float(x['hat_w']) for x in r])
        out[name] = dict(n=len(r), corr_w_a=float(np.corrcoef(w, a)[0, 1]), slope_w_per_a_nm_per_deg=float(np.polyfit(a, w, 1)[0]),
                         expected_slope_strip_symmetry=float(T0 * np.pi / 180))
        print(name, out[name])
    with open(os.path.join(OUT, 'mc_valley.json'), 'w') as f:
        json.dump(out, f, indent=1)


if __name__ == '__main__':
    {'plan': plan, 'fit': fit_surrogates, 'mc': run_mc, 'valley': valley}[sys.argv[1]]()
