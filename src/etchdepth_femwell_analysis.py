"""
Etch-depth / nuisance-parameter study on femwell r12 derivatives: sidewall-angle Cramér-Rao bounds
as unknowns are added in stages.

Run from src/ (main environment):  python etchdepth_femwell_analysis.py
Reads  results_femwell/cache/solves_r12.jsonl          (w, t, a solves)
       results_etchdepth_femwell/cache/solves_r12.jsonl (slab, n_Si, n_ox solves; etchdepth_femwell_batch.py)
Writes results_etchdepth_femwell/data/{derivatives,step_check,stages}.json, step_check.csv

Parameters and units (eigenvectors are in these units):
  w   width [nm] (one lithography bias shared by every structure)
  t   SOI thickness [nm]
  a   sidewall angle of the full (strip) etch [deg]
  d70, d90   partial-etch depth of the 70 / 90 nm-slab rib [nm]; slab = t - d, so the rib's t column is
             dn/dt at fixed d = dn/dt|slab + dn/dslab, and dn/dd = -dn/dslab. The two ribs need different
             etch depths, so they come from different etch steps and get separate unknowns.
  ar70, ar90 sidewall angle of each partial etch [deg]
  nsi, nox   silicon and oxide index [units of 1e-3]; one oxide index for BOX and top cladding
  T   temperature [K], through dn/dT of Si and SiO2 (below)
Stages:
  A   w, t, a; the rib's partial etch shares the strip angle; slab thickness known (rib t column at fixed
      slab, as in the femwell summary)
  A*  as A but etch depth known (rib t column at fixed d); shown because B builds on it
  B   + d (per rib slab)
  C   + separate partial-etch angle (per rib slab); the strip angle a no longer enters the rib
  D   + nsi, nox
  E   + T. T's column is an exact combination of the nsi and nox columns, so it is not identifiable
      next to them and E equals D; E' = C + T (indices known at the reference temperature) shows what
      temperature alone costs.
Noise unchanged: sigma_neff 6.5e-6, sigma_ng 1e-3, sigma_phase 1e-3 rad with L = 100 um (all
placeholders, see fullvector_sensitivity.py); rib multiplier 1.0 (conservative) or 0.2 (optimistic).
"""
import csv
import json
import os

import numpy as np

import fullvector_sensitivity as fv

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, 'results_etchdepth_femwell')
CACHES = [os.path.join(REPO, 'results_femwell', 'cache', 'solves_r12.jsonl'), os.path.join(OUT, 'cache', 'solves_r12.jsonl')]
H = dict(w=5.0, t=2.0, a=1.0, s=2.0, nsi=0.002, nox=0.002)
FACTORS = (0.5, 1.0, 2.0)
WL0, DLAM = 1.55, 0.005
# thermo-optic coefficients at 1.55 um, 295 K (linear interpolation between tabulated wavelengths)
DNDT_SI = 1.86e-4   # 1.87e-4 (1.5 um), 1.78e-4 (2.0 um): Frey, Leviton, Madison, Proc. SPIE 6273, 62732J (2006), Table 4
DNDT_OX = 8.4e-6    # 8.45e-6 (1.5 um), 8.29e-6 (1.6 um): Leviton, Frey, arXiv:0805.0091, Table 6 (bulk Corning 7980)
STRIPS = [('strip', 'TE', w, None) for w in (400, 450, 500, 550, 600)] + [('strip', 'TM', 500, None)]
RIBS = [('rib', 'TE', w, s) for s in (70, 90) for w in (300, 400, 500, 600)]
GEO, MAT = ('w', 't', 'a', 's', 't_d'), ('nsi', 'nox')


def load():
    rows = {}
    for p in CACHES:
        with open(p, encoding='utf-8') as f:
            for line in f:
                if line.strip():
                    r = json.loads(line)
                    rows[r['key']] = r
    return rows


def nng(rows, u, p, d):
    n0, npl, nmi = rows[f'{u[0]}|{u[1]}|{u[2]}|{u[3]}|oxide|r12|{p}|{d:+.4f}']['n']
    return np.array([n0, n0 - WL0 * (npl - nmi) / (2 * DLAM)])


T_STRIP = 220.0  # nm


def impose_strip_symmetry(Dall):
    """Replace every strip's angle column by its exact first-order value, -t (pi/180) x width column.
    An oxide-clad strip is mirror-symmetric about mid-height, so a wall tilt about the bottom edge acts
    at first order exactly like a width change of -t tan(a); the finite-difference residual R - 1 is
    numerical noise and otherwise passes for angle information (results_montecarlo/SUMMARY.md, step 2).
    Default for every CRB that includes strips from 2026-10-27; results committed before then used the
    computed columns (strip_symmetry=False)."""
    out = {}
    for uid, Du in Dall.items():
        out[uid] = {}
        for f, D in Du.items():
            if uid.startswith('strip') and isinstance(D, dict) and 'a' in D:
                D = dict(D)
                D['a'] = -T_STRIP * np.radians(1.0) * np.asarray(D['w'])
            out[uid][f] = D
    return out


def derivatives(rows, strip_symmetry=True):
    """D[uid][f][p] = array([dn_eff, dn_g]) per unit of the study parameter."""
    out = {}
    for u in STRIPS + RIBS:
        uid = f'{u[0]}_{u[1]}_{u[2]}' + (f'_slab{u[3]}' if u[3] else '')
        ps = ('w', 't', 'a', 'nsi', 'nox') + (('s',) if u[0] == 'rib' else ())
        out[uid] = {'nominal': nng(rows, u, 'nom', 0.0).tolist(), 'unit': u}
        for f in FACTORS:
            D = {}
            for p in ps:
                h = f * H[p]
                D[p] = (nng(rows, u, p, h) - nng(rows, u, p, -h)) / (2 * h)
                if p in MAT:
                    D[p] = D[p] * 1e-3            # per 1e-3 of index
            if 's' in D:
                D['t_d'] = D['t'] + D['s']        # thickness at fixed etch depth (slab moves with t)
            out[uid][f] = D
    return impose_strip_symmetry(out) if strip_symmetry else out


def step_check(Dall):
    """2% between h/2, h, 2h; near-zero rule (|D(h)| < 10% of the largest in the same group, pass if the
    spread is within 2% of that largest) applied separately to the geometric and material groups,
    since their units differ."""
    rows = []
    for uid, Du in Dall.items():
        for k, obs in enumerate(('neff', 'ng')):
            for grp in (GEO, MAT):
                ps = [p for p in grp if p in Du[1.0]]
                big = max(abs(Du[1.0][p][k]) for p in ps)
                for p in ps:
                    v = [Du[f][p][k] for f in FACTORS]
                    dev = max(abs(x - v[1]) for x in v)
                    rel = dev / abs(v[1]) if v[1] else np.inf
                    near = abs(v[1]) < 0.1 * big
                    verdict = 'PASS' if rel <= 0.02 else ('PASS_NEAR_ZERO_ABS' if near and dev <= 0.02 * big else 'FAIL')
                    rows.append(dict(structure=uid, param=p, output=obs, D_half=v[0], D_h=v[1], D_2h=v[2],
                                     max_rel_dev=rel, frac_of_largest=abs(v[1]) / big, verdict=verdict))
    return rows


def columns(uid, D, stage):
    """Map one structure's derivatives (dict p -> [dneff, dng]) onto the stage's global parameters."""
    rib = uid.startswith('rib')
    slab = uid.split('slab')[-1] if rib else None
    c = {'w': D['w']}
    if stage == 'A':
        c['t'], c['a'] = D['t'], D['a']
        return c
    c['t'] = D['t_d'] if rib else D['t']
    if rib and stage != 'A*':
        c[f'd{slab}'] = -D['s']
    c[f'ar{slab}' if rib and stage in ('C', 'D', 'E', "E'") else 'a'] = D['a']
    if stage in ('D', 'E'):
        c['nsi'], c['nox'] = D['nsi'], D['nox']
    if stage in ('E', "E'"):
        c['T'] = (DNDT_SI * D['nsi'] + DNDT_OX * D['nox']) / 1e-3   # D per 1e-3 index -> per K
    return c


def rows_for(items, Dall, f, stage, rib_mult):
    """items: list of (uid, readouts) with readouts a subset of ('ring', 'mzi'); returns (params, J)."""
    k = 2 * np.pi / WL0 * fv.L_MZI
    rows = []
    for uid, readouts in items:
        m = rib_mult if uid.startswith('rib') else 1.0
        c = columns(uid, Dall[uid][f], stage)
        if 'ring' in readouts:
            rows.append({p: v[0] / (fv.SIGMA_NEFF * m) for p, v in c.items()})
            rows.append({p: v[1] / (fv.SIGMA_NG * m) for p, v in c.items()})
        if 'mzi' in readouts:
            rows.append({p: k * v[0] / (fv.SIGMA_PHASE * m) for p, v in c.items()})
    slabs = (50, 70, 90, 110, 130, 150)
    order = ['w', 't', 'a'] + [f'd{s}' for s in slabs] + [f'ar{s}' for s in slabs] + ['nsi', 'nox', 'T']
    params = [p for p in order if any(p in r for r in rows)]
    J = np.array([[r.get(p, 0.0) for p in params] for r in rows])
    return params, J


RANK_TOL = 1e-9  # relative singular-value cut on the column-normalized Jacobian: physical near-degeneracies
                 # (R = 1.000x) sit near 1e-4, exact constructions (T from nsi/nox) near 1e-16


def marginal_std(J, i):
    """sqrt of the CRB for parameter i with all others as nuisances. Columns are normalized first so the
    rank decision does not depend on units; the parameter is identifiable only if no null-space direction
    of J has a component along it (otherwise inf). Its variance is then the pseudo-inverse diagonal over
    the numerical rank, which is the CRB of an estimable parameter even when some nuisances are not."""
    norms = np.linalg.norm(J, axis=0)
    _, s, Vt = np.linalg.svd(J / norms, full_matrices=True)
    r = int(np.sum(s > RANK_TOL * s[0]))
    if r < Vt.shape[0] and np.max(np.abs(Vt[r:, i])) > 1e-6:
        return float('inf')
    return float(np.sqrt(np.sum(Vt[:r, i] ** 2 / s[:r] ** 2)) / norms[i])


def analyse(items, Dall, stage, rib_mult=1.0):
    out = {}
    for f in FACTORS:
        params, J = rows_for(items, Dall, f, stage, rib_mult)
        F = J.T @ J
        ev, evec = np.linalg.eigh(F)
        std = {p: marginal_std(J, i) for i, p in enumerate(params) if p == 'a' or p.startswith('ar')}
        out[f] = dict(params=params, angle_std=std, eig_min=float(ev[0]), eig_max=float(ev[-1]),
                      eigvec_min={p: float(x) for p, x in zip(params, evec[:, 0])},
                      cond=float(ev[-1] / ev[0]) if ev[0] > 0 else float('inf'))
    return out


def sets():
    te = [f'strip_TE_{w}' for w in (400, 450, 500, 550, 600)]
    strip = [(u, ('ring', 'mzi')) for u in te]
    tm = [('strip_TM_500', ('ring',))]
    rib = lambda w, s, r=('ring', 'mzi'): [(f'rib_TE_{w}_slab{s}', r)]
    two = [('strip_TE_550', ('ring',)), ('strip_TE_600', ('ring',))]
    S = {
        'strip TE (5 widths, ring+MZI)': strip,
        'strip TE + TM ring': strip + tm,
        'strip + rib70': strip + rib(500, 70),
        'strip + rib90': strip + rib(500, 90),
        'strip + both ribs': strip + rib(500, 70) + rib(500, 90),
        'two strip rings + rib70 ring': two + rib(500, 70, ('ring',)),
        'two strip rings + rib90 ring': two + rib(500, 90, ('ring',)),
    }
    # rescue options: a TM ring, or a second width of the same rib (same etch, so no new unknowns)
    for name in ('strip + rib70', 'strip + rib90', 'strip + both ribs', 'two strip rings + rib70 ring', 'two strip rings + rib90 ring'):
        S[name + ' + TM ring'] = S[name] + tm
    for s in (70, 90):
        for w2 in (300, 400, 600):
            S[f'strip + rib{s} + rib{s} {w2} nm'] = S[f'strip + rib{s}'] + rib(w2, s)
            S[f'two strip rings + rib{s} ring + rib{s} {w2} nm ring'] = S[f'two strip rings + rib{s} ring'] + rib(w2, s, ('ring',))
    return S


def main():
    rows = load()
    Dall = derivatives(rows)
    os.makedirs(os.path.join(OUT, 'data'), exist_ok=True)
    with open(os.path.join(OUT, 'data', 'derivatives.json'), 'w') as f:
        json.dump({u: {('h' if k == 1.0 else 'h/2' if k == 0.5 else '2h' if k == 2.0 else k):
                       ({p: list(v) for p, v in d.items()} if isinstance(d, dict) else d)
                       for k, d in Du.items()} for u, Du in Dall.items()}, f, indent=1)
    steps = step_check(Dall)
    with open(os.path.join(OUT, 'data', 'step_check.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(steps[0].keys()))
        w.writeheader()
        w.writerows(steps)
    fails = [r for r in steps if r['verdict'] == 'FAIL']
    print(f'step check: {len(steps)} derivatives, {sum(r["verdict"] == "PASS" for r in steps)} PASS, '
          f'{sum(r["verdict"] == "PASS_NEAR_ZERO_ABS" for r in steps)} near-zero PASS, {len(fails)} FAIL')
    for r in fails:
        print(f"  FAIL {r['structure']:22s} d{r['output']}/d{r['param']}  D(h/2,h,2h)=({r['D_half']:.4e}, {r['D_h']:.4e}, "
              f"{r['D_2h']:.4e})  rel dev {r['max_rel_dev']:.2%}  |D|/largest {r['frac_of_largest']:.3f}")
    res = {}
    for stage in ('A', 'A*', 'B', 'C', 'D', 'E', "E'"):
        res[stage] = {}
        for name, items in sets().items():
            res[stage][name] = {'cons': analyse(items, Dall, stage, 1.0), 'opt': analyse(items, Dall, stage, 0.2)}
    with open(os.path.join(OUT, 'data', 'stages.json'), 'w') as f:
        json.dump(dict(dndt_si=DNDT_SI, dndt_ox=DNDT_OX, stages={s: {n: {m: {str(k): v for k, v in x.items()}
                                                                          for m, x in d.items()} for n, d in st.items()}
                                                                  for s, st in res.items()}), f, indent=1)
    for stage, st in res.items():
        print(f'\n=== stage {stage} (rib conservative 1.0x; value at h [min..max over h/2, h, 2h]) ===')
        for name, d in st.items():
            c = d['cons']
            angs = c[1.0]['angle_std']
            txt = '  '.join(f"{p} {angs[p]:.3g}° [{min(c[f]['angle_std'][p] for f in FACTORS):.3g}..{max(c[f]['angle_std'][p] for f in FACTORS):.3g}]"
                            for p in angs)
            print(f'  {name:55s} {txt}')


if __name__ == '__main__':
    main()
