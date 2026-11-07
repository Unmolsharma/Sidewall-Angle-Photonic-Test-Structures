"""
Independent check of the Monte Carlo surrogate against direct femwell solves (review, step 7).

For each headline set (case a: best one-etch set, 90 nm slab; case b: two-width set, 50 nm slab), N_TRUTH truth
parameter vectors are drawn from the CRLB covariance at 1 pm (surrogate Jacobian at nominal, realistic priors),
clipped to the Monte Carlo box. Every structure of the set is then solved directly with femwell at its geometry
(same mesh, mode selection and wavelengths as the surrogate training solves; strips with vertical walls at the
mid-height width, which encodes the strip symmetry), and:
  1. the surrogate prediction is compared with the femwell values at each truth (n_eff in units of the MZI phase noise
     expressed as an n_eff noise, n_g in units of sigma_ng, both at 1 pm);
  2. a direct-data Monte Carlo is run: data = femwell values at the truth + noise (N_NOISE draws per truth), prior
     centres = truth + sigma_p N(0, 1), fitted with the surrogate model exactly as in montecarlo.py. The fitted-angle
     error (theta_hat - theta_true) is compared with the CRLB.

Run from src/ (main environment):
    python review_surrogate_check.py design     -> results_review/data/surrogate_check_design.json
    (femwell venv) python review_surrogate_batch.py   solves the design
    python review_surrogate_check.py analyze    -> results_review/data/surrogate_check.json
"""
import json
import os
import sys

import numpy as np

import etchdepth_femwell_analysis as ed
import montecarlo as mc

REPO = ed.REPO
OUT = os.path.join(REPO, 'results_review', 'data')
CACHE = os.path.join(REPO, 'results_review', 'cache', 'solves_r12.jsonl')
N_TRUTH, N_NOISE, SEED = 20, 25, 20261103
CASES = ('a', 'b')


def case(name):
    mc.set_noise('1pm')
    return mc.Case(name, mc.OUT)


def design():
    os.makedirs(OUT, exist_ok=True)
    out = {}
    for name in CASES:
        c = case(name)
        J = c.J0
        C = np.linalg.inv(J.T @ J)
        rng = np.random.default_rng([SEED, ord(name)])
        th = rng.multivariate_normal(np.zeros(len(c.params)), C, size=N_TRUTH)
        th = np.clip(th, 0.9 * c.lo, 0.9 * c.hi)
        pts = []
        for j, t in enumerate(th):
            for u in c.units:
                geom, pol, W, S = mc.unit_of(u)
                x, _ = c.M.x_and_dx(u, t)
                a = mc.absolute(geom, W, S, x)
                pts.append(dict(truth=j, unit=u, geom=geom, pol=pol, W=W, S=S, guess=mc.nominal_guess(u), **a))
        out[name] = dict(units=c.units, params=c.params, truths=th.tolist(), points=pts,
                         crb_a_surrogate=c.crb_sur, crb_a_fd=c.crb_fd)
        print(name, len(pts), 'solves')
    json.dump(dict(seed=SEED, n_truth=N_TRUTH, cases=out), open(os.path.join(OUT, 'surrogate_check_design.json'), 'w'), indent=1)


def key(p):
    return (f"review|{p['geom']}|{p['pol']}|{p['W']}|{p['S']}|{p['w']:.4f}|{p['t']:.4f}|{p['a']:.5f}|{p['slab']}|"
            f"{p['nsi']:.7f}|{p['nox']:.7f}")


def analyze():
    D = json.load(open(os.path.join(OUT, 'surrogate_check_design.json')))
    sol = {}
    for line in open(CACHE, encoding='utf-8'):
        if line.strip():
            r = json.loads(line)
            sol[r['key']] = r
    res = {}
    for name in CASES:
        d = D['cases'][name]
        c = case(name)
        s_ph = mc.fv.SIGMA_PHASE / mc.KPH
        err_neff, err_ng, femwell_y, unit_of = [], [], {}, []
        for j, t in enumerate(d['truths']):
            pred = c.M.predict(np.array(t))
            Y = []
            for k, u in enumerate(c.units):
                p = [q for q in d['points'] if q['truth'] == j and q['unit'] == u][0]
                y = mc.nng_of(sol[key(p)]['n'])
                Y.append(y)
                err_neff.append((pred[k][0] - y[0]) / s_ph)
                err_ng.append((pred[k][1] - y[1]) / mc.fv.SIGMA_NG)
                unit_of.append(str(u))
            femwell_y[j] = np.array(Y)
        e1, e2 = np.array(err_neff), np.array(err_ng)
        uo = np.array(unit_of)
        per_unit = {u: dict(neff_err_rms_sigma_phase=float(np.sqrt(np.mean(e1[uo == u] ** 2))),
                            neff_err_max_sigma_phase=float(np.max(np.abs(e1[uo == u]))),
                            ng_err_rms_sigma_ng=float(np.sqrt(np.mean(e2[uo == u] ** 2))))
                    for u in dict.fromkeys(unit_of)}
        rng = np.random.default_rng([SEED, ord(name), 1])
        errs, status = [], {'ok': 0, 'at_bound': 0, 'failed': 0}
        for j, t in enumerate(d['truths']):
            t = np.array(t)
            y0 = c.M.whiten(femwell_y[j])
            for _ in range(N_NOISE):
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
        boot = np.random.default_rng(SEED).choice(errs, size=(1000, len(errs))).std(axis=1, ddof=1)
        res[name] = dict(units=c.units, n_truths=len(d['truths']), n_solves=len(d['points']),
                         neff_err_rms_sigma_phase=float(np.sqrt(np.mean(e1 ** 2))), neff_err_max_sigma_phase=float(np.max(np.abs(e1))),
                         ng_err_rms_sigma_ng=float(np.sqrt(np.mean(e2 ** 2))), ng_err_max_sigma_ng=float(np.max(np.abs(e2))),
                         per_unit=per_unit, n_fits=int(len(errs) + status['failed']), status=status,
                         angle_error_std=float(errs.std(ddof=1)), angle_error_std_ci95=[float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))],
                         angle_error_mean=float(errs.mean()), angle_error_mean_se=float(errs.std(ddof=1) / np.sqrt(len(errs))),
                         crb_a_fd=c.crb_fd, crb_a_surrogate=c.crb_sur, std_over_crb_fd=float(errs.std(ddof=1) / c.crb_fd))
        print(name, {k: v for k, v in res[name].items() if k != 'units'})
    json.dump(res, open(os.path.join(OUT, 'surrogate_check.json'), 'w'), indent=1)


if __name__ == '__main__':
    {'design': design, 'analyze': analyze}[sys.argv[1]]()
