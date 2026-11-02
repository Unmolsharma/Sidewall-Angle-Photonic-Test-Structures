"""
Mode-field centroid of the sidewall perturbation, for the closed-form angle/width ratio
    R_closed = (y_c - y_pivot) / (h/2),   h = height of the tilted wall,
compared in crosschecks.py with the femwell ratio R = (dn/da) / (dn/dw) / (-h pi/180).

A tilt about the pivot moves the wall at height y by tan(a) (y - y_pivot); to first order this acts like a
uniform wall shift of tan(a) (y_c - y_pivot), where y_c is the centroid of the first-order wall-perturbation
weight along the wall. The weight is the standard boundary-perturbation integrand (Johnson et al., PRE 65,
066611 (2002) -- citation UNVERIFIED in this repo), evaluated with the core-side fields:
    w(y) = (eps1 - eps2) |E_par|^2 + (eps1 - eps2) (eps1 / eps2) |E_perp|^2,
with E_perp = E_x (normal to a vertical wall), E_par = (E_y, E_z), eps1 = n_Si^2, eps2 = n_ox^2.
It is sampled at the quadrature points of the core elements within DELTA of the right wall
(x in [w/2 - DELTA, w/2], y on the wall), DELTA = 6 nm (half the r12 core mesh size). The
variant without E_z is also reported, since the scaling of femwell's stored E_z is not checked here.
For a strip, mirror symmetry about mid-height gives y_c = t/2 and R_closed = 1 exactly; the field value
is a numerical check of that.

Runs in the femwell virtual environment, from src/:
    C:/Users/unmol/.venvs/femwell/Scripts/python.exe crosscheck_fields.py
Writes results_final/data/crosscheck_fields.json
"""
import json
import os

import numpy as np

import femwell_batch as fb
from femwell.maxwell.waveguide import compute_modes

OUT = os.path.join(fb.REPO, 'results_final', 'data')
DELTA = 0.006
UNITS = ([('strip', 'TE', w, None) for w in (400, 500, 600)] + [('strip', 'TM', 500, None)] +
         [('rib', 'TE', w, s) for s in (50, 90) for w in (300, 400, 500, 600)])


def centroid(u):
    geom, pol, W, S = u
    w, t = W / 1000, 0.22
    slab = (S if S is not None else 90) / 1000
    basis0, eps = fb.build(geom, w, t, 0.0, slab, 'oxide', 'r12')
    guess = {'strip': 2.4 if pol == 'TE' else 1.75, 'rib': 2.5}[geom]
    m = fb.pick(compute_modes(basis0, eps, wavelength=fb.WL0, num_modes=2, order=2, n_guess=guess), pol)
    (et, et_b), (ez, ez_b) = m.basis.split(m.E)
    Et = et_b.interpolate(et).value              # (2, nel, nqp)
    Ez = ez_b.interpolate(ez).value              # (nel, nqp)
    X = et_b.global_coordinates().value          # (2, nel, nqp)
    dx = et_b.dx                                 # (nel, nqp)
    y0 = 0.0 if geom == 'strip' else slab
    core = np.zeros(m.basis.mesh.nelements, bool)
    core[basis0.get_dofs(elements='core0')] = True   # for a rib the slab is in the core too; y >= slab keeps the wall
    sel = core[:, None] & (X[0] >= w / 2 - DELTA) & (X[0] <= w / 2) & (X[1] >= y0) & (X[1] <= t)
    e1, e2 = fb.N_SI ** 2, fb.N_OX ** 2
    ex2, ey2, ez2 = np.abs(Et[0]) ** 2, np.abs(Et[1]) ** 2, np.abs(Ez) ** 2
    out = {}
    for label, wt in (('with_Ez', (e1 - e2) * (ey2 + ez2) + (e1 - e2) * (e1 / e2) * ex2),
                      ('without_Ez', (e1 - e2) * ey2 + (e1 - e2) * (e1 / e2) * ex2)):
        ww = (wt * dx)[sel]
        yc = float(np.sum(ww * X[1][sel]) / np.sum(ww))
        out[label] = dict(y_c_um=yc, R_closed=(yc - y0) / ((t - y0) / 2))
    return dict(n_eff=float(m.n_eff.real), te_fraction=float(m.te_fraction), y_pivot_um=y0, wall_height_um=t - y0,
                n_points=int(sel.sum()), delta_um=DELTA, **out)


if __name__ == '__main__':
    os.makedirs(OUT, exist_ok=True)
    res = {}
    for u in UNITS:
        uid = f'{u[0]}_{u[1]}_{u[2]}' + (f'_slab{u[3]}' if u[3] else '')
        res[uid] = centroid(u)
        print(uid, {k: (round(v['R_closed'], 4) if isinstance(v, dict) else v) for k, v in res[uid].items()})
    with open(os.path.join(OUT, 'crosscheck_fields.json'), 'w') as f:
        json.dump(res, f, indent=1)
