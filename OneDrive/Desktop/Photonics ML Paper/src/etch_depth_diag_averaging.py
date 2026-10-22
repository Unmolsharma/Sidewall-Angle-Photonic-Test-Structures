"""
DIAGNOSTIC ONLY (not used by any result): why the rib step-size check failed.

The existing epsfunc factories (fullvector_sensitivity.epsfunc_factory,
rib_waveguide.epsfunc_factory_rib) average each cell over a window of size
dx = np.diff(x_).mean(), dy = np.diff(y_).mean() -- the MEAN spacing of the
whole nonuniform grid (42.9nm x 36.7nm on the production grid) instead of
each cell's true size (26.7nm x 20nm in the core). Windows of neighbouring
cells therefore overlap unevenly, so the eps response to an interface moving
is 1.09x (y) / 1.24x (x) of the true rate at most positions and 0.54x /
0.62x at others, and every nominal slab value used (50..150nm) lands in a
0.54x dip. Derivatives then depend on step size and on where the nominal
interface sits relative to the grid.

This file builds the same rib geometry with per-cell bounds taken from the
actual grid (cell i spans [X[i], X[i+1]], which is what EMpy's centre
points xc = (X[:-1]+X[1:])/2 represent), and re-runs the h/2, h, 2h check at
one structure. If the diagnosis is right, the corrected averaging passes.

Run:  python etch_depth_diag_averaging.py
"""
import os
import json

import numpy as np

import fullvector_sensitivity as fv
from rib_waveguide import make_rib_grid, epsfunc_factory_rib
from etch_depth_sims import H0, nominal, shifted, OUT_DATA


def epsfunc_rib_cell_exact(X, Y, rib_width_um, thick_um, angle_deg, slab_h_um):
    tan_a = np.tan(np.radians(angle_deg))
    x_lo_e, x_hi_e = X[:-1], X[1:]
    y_lo_e, y_hi_e = Y[:-1], Y[1:]

    def epsfunc(x_, y_):
        assert len(x_) == len(X) - 1 and len(y_) == len(Y) - 1
        xl, yl = np.meshgrid(x_lo_e, y_lo_e, indexing='ij')
        xh, yh = np.meshgrid(x_hi_e, y_hi_e, indexing='ij')
        dx, dy = xh - xl, yh - yl
        yc = (yl + yh) / 2
        oy_slab = np.clip(np.minimum(yh, slab_h_um) - np.maximum(yl, 0.0), 0, None) / dy
        oy_rib = np.clip(np.minimum(yh, thick_um) - np.maximum(yl, slab_h_um), 0, None) / dy
        local_hw = np.clip(rib_width_um / 2 - tan_a * (yc - slab_h_um), 0, None)
        ox = np.clip(np.minimum(xh, local_hw) - np.maximum(xl, -local_hw), 0, None) / dx
        fill = np.clip(oy_slab + oy_rib * ox, 0, 1)
        return fill * fv.N_SI ** 2 + (1 - fill) * fv.N_OX ** 2
    return epsfunc


def volume_slope(make_eps, X, Y, s_nm):
    """d(effective slab thickness)/d(slab) in a far-slab column: ideal 1."""
    xc, yc, dyc = (X[:-1] + X[1:]) / 2, (Y[:-1] + Y[1:]) / 2, np.diff(Y)

    def eff(s):
        e = make_eps(s / 1000)(xc, yc)[5]
        return ((e - fv.N_OX ** 2) / (fv.N_SI ** 2 - fv.N_OX ** 2) * dyc).sum() * 1000
    return (eff(s_nm + 0.5) - eff(s_nm - 0.5))


def step_check(averaging, width_nm=500, slab_nm=90, pol='TE'):
    X, Y = make_rib_grid(1.4, 1.2, 30, 18)

    def solve(p, wl):
        slab = (p['t'] - p['d']) / 1000
        if averaging == 'existing':
            ep = epsfunc_factory_rib(p['w'] / 1000, p['t'] / 1000, p['a'], slab)
        else:
            ep = epsfunc_rib_cell_exact(X, Y, p['w'] / 1000, p['t'] / 1000, p['a'], slab)
        return fv.solve_mode_from_eps(ep, wl, pol, 4, X, Y).neff.real

    def neff_ng(p):
        n0, npl, nmi = (solve(p, fv.WL0 + s * 0.005) for s in (0, 1, -1))
        return n0, n0 - fv.WL0 * (npl - nmi) / 0.01

    p0 = nominal(width_nm, slab_nm)
    out = {}
    for k in ('w', 't', 'a', 'd'):
        for f in (0.5, 1.0, 2.0):
            h = f * H0[k]
            (n1, g1), (n2, g2) = neff_ng(shifted(p0, k, h)), neff_ng(shifted(p0, k, -h))
            out[(k, f)] = ((n1 - n2) / (2 * h), (g1 - g2) / (2 * h))
    return out


if __name__ == '__main__':
    X, Y = make_rib_grid(1.4, 1.2, 30, 18)
    print('volume-conservation slope at slab = 90nm (ideal 1.00):')
    print(f"  existing averaging : {volume_slope(lambda s: epsfunc_factory_rib(0.5, 0.22, 0.0, s), X, Y, 90):.3f}")
    print(f"  cell-exact         : {volume_slope(lambda s: epsfunc_rib_cell_exact(X, Y, 0.5, 0.22, 0.0, s), X, Y, 90):.3f}")

    results = {}
    for avg in ('existing', 'cell_exact'):
        res = step_check(avg)
        print(f'\n[{avg} averaging] rib 500nm, slab 90nm, production grid M0')
        print(f"{'param':>5} {'obs':>6} {'D(h/2)':>12} {'D(h)':>12} {'D(2h)':>12} {'max dev':>9}")
        for k in ('w', 't', 'a', 'd'):
            for j, obs in enumerate(('neff', 'ng')):
                d = [res[(k, f)][j] for f in (0.5, 1.0, 2.0)]
                dev = max(abs(d[0] - d[1]), abs(d[2] - d[1])) / abs(d[1])
                print(f"{k:>5} {obs:>6} {d[0]:>12.4e} {d[1]:>12.4e} {d[2]:>12.4e} {dev:>8.1%}")
                results[f'{avg}|{k}|{obs}'] = dict(D_half=d[0], D_h=d[1], D_2h=d[2], max_rel_dev=dev)
    with open(os.path.join(OUT_DATA, 'diag_averaging.json'), 'w') as f:
        json.dump(results, f, indent=1)
