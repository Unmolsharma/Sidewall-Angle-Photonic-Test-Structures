"""
Headline CRBs recomputed with literature-sourced measurement noise instead of the placeholders.

Source (verified against the arXiv full text, v1):
  S. Janz, S. Dedyulin, D.-X. Xu, M. Vachon, S. Wang, R. Cheriton, J. Weber, "Measurement Accuracy in
  Silicon Photonic Ring Resonator Thermometers: Identifying and Mitigating Intrinsic Impairments",
  arXiv:2306.15643 [physics.optics] (2023), DOI 10.48550/arXiv.2306.15643. Journal version: UNVERIFIED.
  - optimistic: "the resonance wavelength uncertainty in Fig. 6 is dominated by the ±0.1 pm repeatability
    of the laser, and by the 0.1 pm wavelength step size used to acquire the spectra" -> sigma_lambda 0.1 pm
    (temperature-bath-stabilized ring, ~10 averaged spectra; instrument floor)
  - conservative: "Wavelength changes down to 1 pm or less can be routinely measured using calibrated
    tunable lasers" (Introduction) -> sigma_lambda 1 pm

Conversions (our assumptions, not from the source): resonance shift d(lambda)/lambda = d(n_eff)/n_g, so
  sigma_neff  = n_g sigma_lambda / lambda
  sigma_ng    = sqrt(2) sigma_lambda n_g^2 L / lambda^2   (n_g from one FSR, two resonances, ring length L)
  sigma_phase = 2 pi n_g L sigma_lambda / lambda^2        (MZI fringe position read to sigma_lambda)
with lambda = 1.55 um, L = 100 um (the length already used for the MZI; also the ring length in Janz et al.)
and n_g = 4.26, the largest group index among the headline structures (conservative). The MZI value
assumes the fringe position is read as precisely as a ring resonance; a broad MZI fringe is likely worse
(UNVERIFIED). Rib noise multiplier unchanged (1.0x, conservative).

Run from src/ (main environment):  python noise_sourced.py
Writes results_montecarlo/noise_sourced/data/noise_sourced.json (log logs/noise_sourced.log)
"""
import json
import os

import numpy as np

import etchdepth_femwell_analysis as ed
import fullvector_sensitivity as fv
import priors_search as ps

OUT = os.path.join(ed.REPO, 'results_montecarlo', 'noise_sourced', 'data')
LAM, L, NG = 1.55e-6, 100e-6, 4.26
SOURCE = dict(authors='S. Janz, S. Dedyulin, D.-X. Xu, M. Vachon, S. Wang, R. Cheriton, J. Weber',
              title='Measurement Accuracy in Silicon Photonic Ring Resonator Thermometers: Identifying and '
                    'Mitigating Intrinsic Impairments',
              venue='arXiv:2306.15643 [physics.optics]', year=2023, doi='10.48550/arXiv.2306.15643',
              journal_version='UNVERIFIED')
LEVELS = {'sourced optimistic (0.1 pm)': 0.1e-12, 'sourced conservative (1 pm)': 1e-12}
RM = ('ring', 'mzi')
SETS = {
    'best one-etch set (8), D realistic': (['strip_TE_400', 'strip_TE_600', 'strip_TM_400', 'strip_TM_600'] +
                                           [f'rib_TE_{w}_slab90' for w in (300, 400, 500, 600)], 'D', 'realistic'),
    'best one-etch set (8), C none': (['strip_TE_400', 'strip_TM_500', 'strip_TM_550', 'strip_TM_600'] +
                                      [f'rib_TE_{w}_slab90' for w in (300, 400, 500, 600)], 'C', 'none'),
    'two-width rib 300 + 500, slab 50, D realistic': ([f'strip_TE_{w}' for w in (400, 450, 500, 550, 600)] +
                                                      ['rib_TE_300_slab50', 'rib_TE_500_slab50'], 'D', 'realistic'),
    'strip + rib70 500 + 300, E realistic': ([f'strip_TE_{w}' for w in (400, 450, 500, 550, 600)] +
                                             ['rib_TE_500_slab70', 'rib_TE_300_slab70'], 'E', 'realistic'),
}


def sigmas(sl):
    return dict(SIGMA_NEFF=NG * sl / LAM, SIGMA_NG=np.sqrt(2) * sl * NG ** 2 * L / LAM ** 2,
                SIGMA_PHASE=2 * np.pi * NG * L * sl / LAM ** 2)


def main():
    os.makedirs(OUT, exist_ok=True)
    Dall = ps.derivatives(ed.load())           # strip symmetry imposed (default)
    placeholder = dict(SIGMA_NEFF=fv.SIGMA_NEFF, SIGMA_NG=fv.SIGMA_NG, SIGMA_PHASE=fv.SIGMA_PHASE)
    levels = {'placeholder': placeholder, **{k: sigmas(v) for k, v in LEVELS.items()}}
    res = dict(source=SOURCE, assumptions=dict(lambda_um=1.55, L_um=100, n_g=NG, rib_noise_mult=1.0),
               noise=levels, sets={})
    for k, s in levels.items():
        print(f"{k:30s} sigma_neff {s['SIGMA_NEFF']:.3g}  sigma_ng {s['SIGMA_NG']:.3g}  sigma_phase {s['SIGMA_PHASE']:.3g} rad")
    for name, (units, stage, prior) in SETS.items():
        res['sets'][name] = dict(units=units, stage=stage, prior=prior)
        for k, s in levels.items():
            for a, v in s.items():
                setattr(fv, a, v)
            v, lo, hi = ps.rng(ps.evaluate([(u, RM) for u in units], Dall, stage, ps.PRIORS[prior]))
            res['sets'][name][k] = dict(a_deg_h=v, a_deg_min=lo, a_deg_max=hi)
            print(f'{name:48s} {k:30s} {v:.4g} [{lo:.4g}-{hi:.4g}]')
        for a, v in placeholder.items():
            setattr(fv, a, v)
    with open(os.path.join(OUT, 'noise_sourced.json'), 'w') as f:
        json.dump(res, f, indent=1)


if __name__ == '__main__':
    main()
