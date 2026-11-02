"""
Sensitivity of the headline strip-angle CRBs to each assumed input, varied one at a time (x0.5 and x2, or the
stated alternative), at the primary noise level (sourced 1 pm) and stage D, realistic priors.

Headline sets (read from results_final/sym_search_1pm/data):
  best   best one-etch set of 8 structures over all slabs (search.json, 'D|realistic', one_etch_overall['8'])
  two    strip TE 400-600 + rib 300 + 500 nm at the best slab for a two-width set (slab_pairs.json, 'D|realistic')

Varied (each alone, everything else at its base value):
  sigma_lambda 1 pm      -> 0.1 (best case), 0.5, 2 pm                    (noise level; sigma_neff, sigma_ng, sigma_phase scale with it)
  L = 100 um             -> 50, 200 um                   (ring/MZI length in the sigma_ng, sigma_phase conversions)
  n_g = 4.26             -> 3.37 (smallest n_g of the headline structures, from the solves)
  MZI fringe precision   -> sigma_phase x 2, x 0.5
  rib noise factor 1.0   -> 0.5, 2
  d prior 2 nm           -> 1, 4 nm
  n_Si prior 1e-4        -> 0.5e-4, 2e-4
  n_ox prior 1e-3        -> 0.5e-3, 2e-3
  derivative step h      -> h/2, 2h                      (existing stencils)
  strip symmetry imposed -> raw femwell strip derivatives (checks the R = 1 assumption)
Not varied (would need new solves; listed as NOT EVALUATED): wavelength step of the n_g difference (5 nm),
nominal n_Si / n_ox, trapezoid wall shape (no rounding), pivot location of the angle, material dispersion.

Run from src/ (main environment):  python assumption_sensitivity.py
Writes results_final/data/assumption_sensitivity.json
"""
import copy
import json
import os

import numpy as np

import etchdepth_femwell_analysis as ed
import fullvector_sensitivity as fv
import priors_search as ps

REPO = ed.REPO
SRC = os.path.join(REPO, 'results_final', 'sym_search_1pm', 'data')
OUT = os.path.join(REPO, 'results_final', 'data')
RM = ('ring', 'mzi')
BASE = dict(sl=1e-12, L=100e-6, ng=4.26, phase_mult=1.0)


def sigmas(sl, L, ng, phase_mult):
    lam = 1.55e-6
    return dict(SIGMA_NEFF=ng * sl / lam, SIGMA_NG=np.sqrt(2) * sl * ng ** 2 * L / lam ** 2,
                SIGMA_PHASE=phase_mult * 2 * np.pi * ng * L * sl / lam ** 2)


def crb(units, Dall, noise=None, prior=None, rib_mult=1.0, factor=1.0):
    for k, v in sigmas(**{**BASE, **(noise or {})}).items():
        setattr(fv, k, v)
    e = ps.evaluate([(u, RM) for u in units], Dall, 'D', prior or ps.PRIORS['realistic'], rib_mult=rib_mult,
                    factors=(factor,))
    return e[factor]['a']


def headline_sets():
    S = json.load(open(os.path.join(SRC, 'search.json')))
    best = S['D|realistic']['one_etch_overall']['8']
    P = json.load(open(os.path.join(SRC, 'slab_pairs.json')))['D|realistic']
    pairs = {k: v for k, v in P.items() if 'rib_TE_300' in k and 'rib_TE_500' in k}
    k2 = min(pairs, key=lambda k: pairs[k]['a_deg_h'])
    two = [f'strip_TE_{w}' for w in (400, 450, 500, 550, 600)] + k2.split(' + ')
    return {'best': best['set'], 'two': two}, k2


def main():
    os.makedirs(OUT, exist_ok=True)
    rows = ed.load()
    Dall = ps.derivatives(rows)
    Draw = ps.derivatives(rows, strip_symmetry=False)
    sets, two_key = headline_sets()
    real = ps.PRIORS['realistic']
    cases = [
        ('sigma_lambda', '0.1 pm (best case)', dict(noise=dict(sl=0.1e-12))),
        ('sigma_lambda', '0.5 pm', dict(noise=dict(sl=0.5e-12))), ('sigma_lambda', '2 pm', dict(noise=dict(sl=2e-12))),
        ('L (sigma_ng, sigma_phase conversion)', '50 um', dict(noise=dict(L=50e-6))),
        ('L (sigma_ng, sigma_phase conversion)', '200 um', dict(noise=dict(L=200e-6))),
        ('n_g in the conversion', '3.37', dict(noise=dict(ng=3.37))),
        ('MZI fringe precision', 'sigma_phase x0.5', dict(noise=dict(phase_mult=0.5))),
        ('MZI fringe precision', 'sigma_phase x2', dict(noise=dict(phase_mult=2.0))),
        ('rib noise factor', '0.5', dict(rib_mult=0.5)), ('rib noise factor', '2', dict(rib_mult=2.0)),
        ('d prior', '1 nm', dict(prior={**real, 'd': 1.0})), ('d prior', '4 nm', dict(prior={**real, 'd': 4.0})),
        ('n_Si prior', '0.5e-4', dict(prior={**real, 'nsi': 0.05})), ('n_Si prior', '2e-4', dict(prior={**real, 'nsi': 0.2})),
        ('n_ox prior', '0.5e-3', dict(prior={**real, 'nox': 0.5})), ('n_ox prior', '2e-3', dict(prior={**real, 'nox': 2.0})),
        ('derivative step', 'h/2', dict(factor=0.5)), ('derivative step', '2h', dict(factor=2.0)),
        ('strip symmetry', 'not imposed (raw femwell)', dict(D=Draw)),
    ]
    res = dict(base=dict(sigma_lambda_pm=1.0, L_um=100, n_g=4.26, phase_mult=1.0, rib_mult=1.0, priors=real, stage='D'),
               sets=sets, two_width_pair=two_key, results={},
               not_evaluated=['wavelength step of the n_g central difference (5 nm)', 'nominal n_Si 3.48 / n_ox 1.44',
                              'trapezoid wall shape (no corner rounding)', 'pivot location of the angle (bottom of wall)',
                              'material dispersion (none modelled)'])
    for name, units in sets.items():
        base = crb(units, Dall)
        r = dict(base=base, variations=[])
        for what, level, kw in cases:
            kw = dict(kw)
            v = crb(units, kw.pop('D', Dall), **kw)
            r['variations'].append(dict(assumption=what, level=level, a_deg=v, rel_change=v / base - 1))
        res['results'][name] = r
        print(f'{name}: base {base:.4g}')
        for x in r['variations']:
            print(f"   {x['assumption']:38s} {x['level']:26s} {x['a_deg']:.4g} ({x['rel_change']:+.1%})")
    with open(os.path.join(OUT, 'assumption_sensitivity.json'), 'w') as f:
        json.dump(res, f, indent=1)


if __name__ == '__main__':
    main()
