"""
Third-review analyses of the two-width set (strip TE 400-600 nm + ribs 300 and 500 nm on a 50 nm slab), stage D,
realistic priors, 1 pm, ring and MZI unless stated. Builds on review_feasibility.py (same blocks, priors and noise).

  systematic   width offsets that depend on the drawn width (an etch bias), modeled as fitted nuisance parameters with
               and without Gaussian priors: x = (W - mean W)/100 nm, offset = c1 x (+ c2 x^2), applied to all structures,
               to the strips only, or separately to strips and ribs; a constant width offset of the ribs relative to the
               strips; a free offset per design. Also combined with random local offsets and N copies.
  correlation  replication with spatially correlated local offsets, rho(d) = exp(-d/L) (L ASSUMED), copies laid out in
               cells of pitch s (one copy of each design per cell, 3 x 3 sub-grid at s/3); either one random field
               shared by all structures, or correlation only among copies of the same design.
  prediction   1-sigma resonance prediction of a new strip ring and rib ring (TE 500 nm) from the fit of N copies with
               random local offsets, with and without the new device's own local offsets.

Run from src/ (main environment):  python review_systematic.py
Writes results_review/data/systematic.json (log results_review/logs/systematic.log)
"""
import json
import os

import numpy as np

import review_analyses as ra
import review_feasibility as rf

ed = ra.ed
OUT = ra.OUT
UNITS = ra.SETS['two']
MOD = (0.3, 0.1, 0.3)                       # (w nm, t nm, angle deg): the 'moderate' scenario
PRIORS = (None, 0.01, 0.03, 0.1, 0.3, 1.0)  # nm (per 100 nm for c1, per (100 nm)^2 for c2)
L_MM = (0.05, 0.2, 1.0, 5.0)                # ASSUMED correlation lengths
S_MM = (0.15, 0.5, 1.5)                     # copy spacing (cell pitch)
NMAX = 256


def widths():
    return np.array([float(u.split('_')[2]) for u in UNITS])


def basis(model):
    """Columns (one per nuisance coefficient) of the width offset of each structure, per unit coefficient."""
    W = widths()
    x = (W - W.mean()) / 100.0
    strip = np.array([u.startswith('strip') for u in UNITS], float)
    rib = 1 - strip
    one = np.ones_like(x)
    cols = {'lin_all': [x], 'quad_all': [x, x ** 2], 'lin_strip': [x * strip], 'quad_strip': [x * strip, x ** 2 * strip],
            'rib_offset': [rib], 'rib_offset_lin': [rib, x * rib], 'separate_lin': [x * strip, rib, x * rib],
            'per_design': [np.eye(len(UNITS))[i] for i in range(len(UNITS))]}[model]
    return np.array(cols).T                  # structures x coefficients


def bound_sys(params, bl, model=None, prior=None, N=1, sd=(0.0, 0.0, 0.0), return_cov=False):
    """CRLB of the strip angle with systematic width-offset coefficients fitted (Gaussian prior or none)."""
    B = basis(model) if model else np.zeros((len(bl), 0))
    rows = []
    for b, coef in zip(bl, B):
        Jx = np.hstack([b['J'], np.outer(b['D'][:, 0], coef)]) if B.shape[1] else b['J']
        rows.append(rf.whiten(Jx, rf.V_of(b, N, sd)))
    k = B.shape[1]
    pri = [np.r_[p, np.zeros(k)] for p in rf.prior_rows(params)]
    if prior:
        pri += [np.r_[np.zeros(len(params)), np.eye(k)[j] / prior] for j in range(k)]
    Jt = np.vstack(rows + pri)
    a = ed.marginal_std(Jt, params.index('a'))
    return (a, Jt) if return_cov else a


def systematic():
    res = {}
    for ro in ('ring_mzi', 'ring'):
        params, bl = rf.blocks(UNITS, rf.READ[ro])
        base = rf.bound(params, bl)
        r = dict(shared=base, centre_nm=float(widths().mean()), models={})
        for m in ('lin_all', 'quad_all', 'lin_strip', 'quad_strip', 'rib_offset', 'rib_offset_lin', 'separate_lin', 'per_design'):
            r['models'][m] = {('none' if p is None else f'{p:g}'): bound_sys(params, bl, m, p) for p in PRIORS}
        # fitted systematic terms (no prior) together with random local offsets and N copies
        r['with_random'] = {}
        for m in ('none', 'lin_all', 'quad_all', 'quad_strip'):
            r['with_random'][m] = {str(N): bound_sys(params, bl, None if m == 'none' else m, None, N, MOD) for N in (1, 4, 16, 34, 61, 64, 256)}
            r['with_random'][m]['N_2'] = n_needed_sys(params, bl, None if m == 'none' else m, 2 * base)
            r['with_random'][m]['N_1p5'] = n_needed_sys(params, bl, None if m == 'none' else m, 1.5 * base)
        res[ro] = r
        print(ro, json.dumps(r, default=float)[:1500], flush=True)
    return res


def n_needed_sys(params, bl, model, target, nmax=1e7):
    f = lambda N: bound_sys(params, bl, model, None, N, MOD)
    if f(1) <= target:
        return 1
    if f(nmax) > target:
        return None
    lo, hi = 0.0, np.log(nmax)
    for _ in range(70):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if f(np.exp(mid)) > target else (lo, mid)
    return int(np.ceil(np.exp(hi) - 1e-6))


# ------------------------------------------------------------------ spatial correlation
def layout(N, s):
    """Positions (mm) of copy j of design i: cells of pitch s on a near-square grid, designs on a 3 x 3 sub-grid."""
    nc = int(np.ceil(np.sqrt(N)))
    cells = np.array([(j % nc, j // nc) for j in range(N)], float) * s
    sub = np.array([(k % 3, k // 3) for k in range(len(UNITS))], float) * s / 3
    return cells, sub                        # position of (i, j) = cells[j] + sub[i]


def ring_equivalent(bl):
    """Merge each MZI phase row into the ring n_eff row (parallel rows; exact, Sec. 'What the MZI adds'):
    one n_eff row with the combined noise. Returns blocks with rows (n_eff, n_g)."""
    out = []
    for b in bl:
        k = b['kinds']
        if 'phase' not in k:
            out.append(b)
            continue
        i_n, i_g, i_p = k.index('neff'), k.index('ng'), k.index('phase')
        r = b['J'][i_p, 0] / b['J'][i_n, 0]
        s2 = 1.0 / (1.0 / b['S'][i_n, i_n] + r ** 2 / b['S'][i_p, i_p])
        out.append(dict(uid=b['uid'], kinds=['neff', 'ng'], J=b['J'][[i_n, i_g]], D=b['D'][[i_n, i_g]],
                        S=np.diag([s2, b['S'][i_g, i_g]])))
    return out


def bound_corr(params, bl, N, L, s, mode, sd=MOD):
    """CRLB with exp(-d/L) correlated local offsets (field shared by all structures, or same design only)."""
    bl = ring_equivalent(bl)
    cells, sub = layout(N, s)
    copies = [(i, j) for j in range(N) for i in range(len(bl))]
    P = np.array([cells[j] + sub[i] for i, j in copies])
    nr = bl[0]['J'].shape[0]
    cidx = np.repeat(np.arange(len(copies)), nr)
    d = np.sqrt(((P[cidx][:, None, :] - P[cidx][None, :, :]) ** 2).sum(-1))
    K = np.exp(-d / L) if L > 0 else (d == 0).astype(float) * (cidx[:, None] == cidx[None, :])
    if mode == 'same_design':
        ids = np.array([i for i, _ in copies])[cidx]
        K = K * (ids[:, None] == ids[None, :])
    J = np.vstack([bl[i]['J'] for i, _ in copies])
    Dr = np.vstack([bl[i]['D'] for i, _ in copies])
    C = np.diag(np.concatenate([np.diag(bl[i]['S']) for i, _ in copies]))
    for p in range(3):
        if sd[p]:
            C += sd[p] ** 2 * np.outer(Dr[:, p], Dr[:, p]) * K
    Jt = np.vstack([rf.whiten(J, C)] + rf.prior_rows(params))
    return ed.marginal_std(Jt, params.index('a'))


def correlation():
    params, bl = rf.blocks(UNITS)
    base = rf.bound(params, bl)
    res = dict(shared=base, target_2=2 * base, target_1p5=1.5 * base, L_mm=list(L_MM), s_mm=list(S_MM), sd=dict(width_nm=MOD[0], thickness_nm=MOD[1], angle_deg=MOD[2]),
               check_independent_N16=bound_corr(params, bl, 16, 0.0, 0.5, 'field'), check_independent_N16_ref=rf.bound(params, bl, 16, MOD), cases={})
    for mode in ('field', 'same_design'):
        for L in L_MM:
            for s in S_MM:
                f = lambda N: bound_corr(params, bl, int(N), L, s, mode)
                grid = {str(N): f(N) for N in (1, 4, 16, 64, 256)}
                out = dict(bound=grid)
                for lab, tgt in (('N_2', 2 * base), ('N_1p5', 1.5 * base)):
                    if grid['1'] <= tgt:
                        out[lab] = 1
                    elif grid['256'] > tgt:
                        out[lab] = None                     # more than NMAX copies
                    else:
                        lo, hi = 1, 256
                        while hi - lo > 1:
                            mid = (lo + hi) // 2
                            lo, hi = (mid, hi) if f(mid) > tgt else (lo, mid)
                        out[lab] = hi
                res['cases'][f'{mode}|L={L:g}|s={s:g}'] = dict(mode=mode, L_mm=L, s_mm=s, **out)
                print(mode, L, s, {k: round(v, 3) for k, v in grid.items()}, out['N_2'], out['N_1p5'], flush=True)
    return res


# ------------------------------------------------------------------ prediction
def prediction():
    params, bl = rf.blocks(UNITS)
    devices = ['strip_TE_500', 'rib_TE_500_slab50']
    res = {}
    for dev in devices:
        c = ed.columns(dev, ra.D[dev][1.0], ra.STAGE)
        pm = rf.LAM / ra.ng0(dev) * 1e12
        g = np.array([c[p][0] if p in c else 0.0 for p in params]) * pm
        akey = 'a' if dev.startswith('strip') else [q for q in c if q.startswith('ar')][0]
        own_cols = np.array([c['w'][0], c['t'][0], c[akey][0]]) * pm
        own = float(np.sqrt(np.sum((own_cols * np.array(MOD)) ** 2)))
        r = dict(own_offsets_pm=own, own_by_param_pm=dict(width=float(abs(own_cols[0]) * MOD[0]), thickness=float(abs(own_cols[1]) * MOD[1]), angle=float(abs(own_cols[2]) * MOD[2])))
        _, Jt = bound_sys(params, bl, None, None, 1, (0, 0, 0), return_cov=True)
        Cs = np.linalg.pinv(Jt.T @ Jt, rcond=1e-14)
        r['shared_ideal_pm'] = float(np.sqrt(g @ Cs @ g))
        r['by_N'] = {}
        for N in (1, 4, 16, 34, 61, 64, 256):
            _, Jt = bound_sys(params, bl, None, None, N, MOD, return_cov=True)
            C = np.linalg.pinv(Jt.T @ Jt, rcond=1e-14)
            mean = float(np.sqrt(g @ C @ g))
            r['by_N'][str(N)] = dict(process_mean_pm=mean, new_device_pm=float(np.hypot(mean, own)))
        res[dev] = r
        print(dev, json.dumps(r, default=float), flush=True)
    return res


def main():
    res = dict(units=UNITS, systematic=systematic(), prediction=prediction(), correlation=correlation())
    json.dump(res, open(os.path.join(OUT, 'systematic.json'), 'w'), indent=1, default=float)


if __name__ == '__main__':
    main()
