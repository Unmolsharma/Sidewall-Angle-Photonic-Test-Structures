"""
Step 6b: what a sidewall-angle change does to a width extracted from strip
rings that ignore angle, and how wrong that makes a rib prediction.

Model (linear, femwell derivatives at step h, oxide clad, noise unchanged):
  - a true sidewall-angle change da (1 and 2 deg) shifts every strip TE ring's
    n_eff and n_g by J_a da (5 strip TE widths, 400-600 nm)
  - the extraction fits the strip ring data with angle held at nominal, by
    noise-weighted least squares (sigma_neff = 6.5e-6, sigma_ng = 1e-3):
      (A) width only (thickness assumed known)
      (B) width and thickness
  - the fitted "effective width" (and thickness) is used to predict the 500 nm
    rib TE (90 nm slab); the error is predicted minus true (true = J_rib,a da)
  - equivalent resonance-wavelength error of a ring near 1550 nm:
      d lambda = lambda * dn_eff / n_g   (fixed ring length and order)

Run from src/:  python femwell_impact.py
Writes results_femwell/data/impact_sidewall_angle.json
"""
import json
import os

import numpy as np

import fullvector_sensitivity as fv

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D = os.path.join(REPO, 'results_femwell', 'data')


def main():
    recs = {r['id']: r for r in json.load(open(os.path.join(D, 'sensitivity_db.json')))}
    strips = [recs[f'strip_TE_{w}'] for w in (400, 450, 500, 550, 600)]
    rib = recs['rib_TE_500']
    sig = []
    J = {p: [] for p in 'wta'}
    for r in strips:
        for obs, s in (('neff', fv.SIGMA_NEFF), ('ng', fv.SIGMA_NG)):
            sig.append(s)
            for p in 'wta':
                J[p].append(r[f'd{obs}_d{p}'])
    sig = np.array(sig)
    Jw, Jt, Ja = (np.array(J[p]) / sig for p in 'wta')
    out = dict(assumptions=dict(sigma_neff=fv.SIGMA_NEFF, sigma_ng=fv.SIGMA_NG, strip_widths_nm=[400, 450, 500, 550, 600],
                                rib='rib_TE_500 (90 nm slab)', wavelength_um=fv.WL0, n_g_rib=rib['ng']))
    for da in (1.0, 2.0):
        y = Ja * da                                   # whitened strip data caused by the angle change
        o = out[f'da={da:g}deg'] = {}
        for name, A in (('width_only', Jw[:, None]), ('width_and_thickness', np.stack([Jw, Jt], axis=1))):
            x, *_ = np.linalg.lstsq(A, y, rcond=None)
            resid = y - A @ x
            dw = float(x[0])
            dt = float(x[1]) if len(x) > 1 else 0.0
            res = {}
            for obs in ('neff', 'ng'):
                pred = rib[f'd{obs}_dw'] * dw + rib[f'd{obs}_dt'] * dt
                true = rib[f'd{obs}_da'] * da
                res[obs] = dict(predicted=pred, true=true, error=pred - true)
            dlam_nm = 1000 * fv.WL0 * res['neff']['error'] / rib['ng']
            o[name] = dict(effective_dw_nm=dw, effective_dt_nm=dt, strip_fit_chi2=float(resid @ resid),
                           strip_fit_chi2_without_fit=float(y @ y), rib=res, resonance_error_pm=1000 * dlam_nm)
            print(f"da = {da:g} deg, fit {name:20s}: effective dw = {dw:+.3f} nm, dt = {dt:+.4f} nm, "
                  f"strip chi2 {resid @ resid:.3g} (vs {y @ y:.3g} unfitted) | rib dn_eff error {res['neff']['error']:+.3e} "
                  f"(true change {res['neff']['true']:+.3e}), dn_g error {res['ng']['error']:+.3e} "
                  f"(true {res['ng']['true']:+.3e}) | resonance error {1000*dlam_nm:+.1f} pm")
    with open(os.path.join(D, 'impact_sidewall_angle.json'), 'w') as f:
        json.dump(out, f, indent=1)


if __name__ == '__main__':
    main()
