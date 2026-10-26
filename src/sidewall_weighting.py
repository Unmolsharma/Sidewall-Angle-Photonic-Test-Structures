"""
Step 3: symmetry argument from the boundary-perturbation formula.

Per unit length along an etched sidewall the perturbation weight is
    w(y) = [ (eps1 - eps2)|E_par|^2 - (1/eps1 - 1/eps2)|D_perp|^2 ] / (4P)
Width moves the wall uniformly (V = 1/2 per unit width), the tilt moves it by
-(y - y_pivot) tan(a). Hence, exactly to first order,
    (dn/da) / (dn/dw) = -2 tan(1 deg) (y_c - y_pivot),   y_c = Int y w dy / Int w dy
and with h the wall height
    R = (dn/da) / (-h tan(1 deg) dn/dw) = (y_c - y_pivot) / (h/2).
So the angle derivative is a fixed multiple of the width derivative (R = 1)
whenever the weighting along the wall is symmetric about mid-height (more
precisely: whenever its centroid is at mid-height).

This script evaluates w(y) from EMpy fields on boundary-aligned grids
(perturbation_derivatives.Structure) for strip TE oxide-clad, strip TE
air-clad, strip TM, rib TE 500, at several grid refinements, and reports the
centroid, R, and an asymmetry measure A = Int (w(y) - w(h - y)) dy / Int w dy
(0 = symmetric). The absolute boundary values from EMpy fields are not
reliable (see compare_derivatives.py), so only the SHAPE of w(y) is used here,
and R is checked against femwell's ratio of FD derivatives.

Run from src/:  python sidewall_weighting.py
Writes results_anisotropic/data/sidewall_weighting.json
"""
import json
import os

import numpy as np

import perturbation_derivatives as pd

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CASES = [('strip', 'TE', 'oxide'), ('strip', 'TE', 'air'), ('strip', 'TM', 'oxide'), ('rib', 'TE', 'oxide')]
REFINES = (1.0, 2.0, 3.0)


def analyse(geom, pol, clad, refine):
    st = pd.Structure(geom, 500, pol, refine, clad)
    m = st.solve()
    _, prof = pd.dn_dp(st, m, profile=True)
    p = [q for q in prof if q['x'] > 0][0]           # right wall (both walls are mirror images)
    y, dy, wgt = np.array(p['y']), np.array(p['dy']), np.array(p['weight'])
    h = st.t - st.y_pivot
    yc = np.sum(y * wgt * dy) / np.sum(wgt * dy)
    u = (y - st.y_pivot) / h                          # 0 at the pivot (bottom), 1 at the top
    wm = np.interp(1 - u, u, wgt)                     # weight mirrored about mid-height
    A = np.sum((wgt - wm) * dy) / np.sum(wgt * dy)
    lower = np.sum((wgt * dy)[u < 0.5]) / np.sum(wgt * dy)
    return dict(neff=float(m.neff.real), y_pivot=st.y_pivot, wall_height=h, y_centroid=float(yc),
                R=float((yc - st.y_pivot) / (h / 2)), ratio_nm_per_deg=float(-2 * np.tan(np.radians(1)) * (yc - st.y_pivot) * 1000),
                asymmetry=float(A), share_lower_half=float(lower),
                profile=dict(u=u.tolist(), weight_normalized=(wgt / np.max(wgt)).tolist()))


if __name__ == '__main__':
    out = {}
    for geom, pol, clad in CASES:
        for r in REFINES:
            k = f'{geom}_{pol}_500_{clad}|refine{r:g}'
            out[k] = a = analyse(geom, pol, clad, r)
            print(f"{k:30s} y_c - pivot = {1000*(a['y_centroid']-a['y_pivot']):6.2f} nm of h = {1000*a['wall_height']:.0f} nm   "
                  f"R = {a['R']:.4f}   (dn/da)/(dn/dw) = {a['ratio_nm_per_deg']:+.3f} nm/deg   "
                  f"weight in lower half {100*a['share_lower_half']:.1f}%", flush=True)
    with open(os.path.join(REPO, 'results_anisotropic', 'data', 'sidewall_weighting.json'), 'w') as f:
        json.dump(out, f, indent=1)
