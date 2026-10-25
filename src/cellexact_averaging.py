"""
Cell-exact sub-pixel dielectric averaging, built ALONGSIDE the original
epsfunc factories (fullvector_sensitivity.epsfunc_factory,
rib_waveguide.epsfunc_factory_rib, slot_waveguide.epsfunc_factory_slot),
which are left unchanged so old and new can be run side by side.

The bug being fixed: the original factories average every cell over a window
of dx = np.diff(x_).mean(), dy = np.diff(y_).mean() -- the MEAN spacing of
the whole nonuniform grid (~43 x 37 nm on the production grid) instead of
each cell's own size (~27 x 20 nm in the core). Neighbouring windows then
overlap unevenly, so moving an interface by 1 nm changes the modelled silicon
by anything from ~0.5 to ~1.25 nm worth, depending on where the interface
sits relative to the grid. That biases every finite-difference derivative,
by a different amount per parameter.

What this module does instead: EMpy evaluates eps at the centres of the cells
[X[i], X[i+1]] x [Y[j], Y[j+1]] (VFDModeSolver.build_matrix: xc = (x[:-1] +
x[1:])/2). Each factory here is given the actual grid X, Y, and fills every
cell with the EXACT area fraction of silicon inside that cell:

    fill = area(cell intersect silicon) / area(cell)
    eps  = fill*n_core^2 + (1 - fill)*n_clad^2      (same arithmetic mean as before)

For the trapezoidal sidewall the x-overlap of a cell with the core is a
piecewise-linear function of y, so it is integrated exactly with the
trapezoid rule between its breakpoints (not evaluated at the cell centre).
The total silicon area of the model therefore equals the true cross-section
area to rounding error, and moving any boundary by dh changes it by exactly
(boundary length) x dh.

Only the averaging window changes. Geometry conventions (trapezoid law, rib
slab, slot rail/gap coupling, SWG effective-medium core), grids, mode solve
and TEfrac classification are all reused unchanged from the original modules.
"""
from functools import partial

import numpy as np

import fullvector_sensitivity as fv
import rib_waveguide as rw
import slot_waveguide as sw
import swg_waveguide as swg

N_SI, N_OX, WL0 = fv.N_SI, fv.N_OX, fv.WL0


def _cell_edges(X, Y):
    xl, yl = np.meshgrid(X[:-1], Y[:-1], indexing='ij')
    xh, yh = np.meshgrid(X[1:], Y[1:], indexing='ij')
    return xl, xh, yl, yh


def band_area(X, Y, cx, hw0, tan_a, y0, y1, y_base):
    """Exact silicon area inside every cell for the band
        y0 <= y <= y1,  |x - cx| <= max(hw0 - tan_a*(y - y_base), 0).
    Returns an (len(X)-1, len(Y)-1) array of areas (um^2)."""
    xl, xh, yl, yh = _cell_edges(X, Y)
    ya = np.maximum(yl, y0)
    yz = np.minimum(yh, y1)
    valid = yz > ya
    yz = np.where(valid, yz, ya)

    # breakpoints of the overlap length L(y): where the local half-width hits
    # one of the cell's x-edges (from either side) or zero
    pts = [ya, yz]
    if tan_a != 0.0:
        for v in (xh - cx, xl - cx, cx - xl, cx - xh, np.zeros_like(xl)):
            pts.append(np.clip(y_base + (hw0 - v) / tan_a, ya, yz))
    pts = np.sort(np.stack(pts, axis=-1), axis=-1)

    hw = np.maximum(hw0 - tan_a * (pts - y_base), 0.0)
    L = np.clip(np.minimum(xh[..., None], cx + hw) - np.maximum(xl[..., None], cx - hw), 0.0, None)
    area = (0.5 * (L[..., 1:] + L[..., :-1]) * np.diff(pts, axis=-1)).sum(axis=-1)
    return np.where(valid, area, 0.0)


def slab_area(X, Y, y0, y1):
    """Exact area of the full-width layer y0 <= y <= y1 inside every cell."""
    xl, xh, yl, yh = _cell_edges(X, Y)
    return (xh - xl) * np.clip(np.minimum(yh, y1) - np.maximum(yl, y0), 0.0, None)


def cell_areas(X, Y):
    xl, xh, yl, yh = _cell_edges(X, Y)
    return (xh - xl) * (yh - yl)


def _make_epsfunc(X, Y, fill, n_core):
    xc_ref = (X[:-1] + X[1:]) / 2
    yc_ref = (Y[:-1] + Y[1:]) / 2
    eps = fill * n_core**2 + (1 - fill) * N_OX**2

    def epsfunc(x_, y_):
        # guard: these per-cell fills are only valid on the grid they were built for
        assert len(x_) == len(xc_ref) and len(y_) == len(yc_ref), "epsfunc called on a different grid"
        assert np.allclose(x_, xc_ref) and np.allclose(y_, yc_ref), "epsfunc called on a different grid"
        return eps
    return epsfunc


# ---------------------------------------------------------------------------
# Fill fractions per geometry (same conventions as the original factories)
# ---------------------------------------------------------------------------
def fill_strip(X, Y, width_um, thick_um, angle_deg):
    tan_a = np.tan(np.radians(angle_deg))
    return band_area(X, Y, 0.0, width_um / 2, tan_a, 0.0, thick_um, 0.0) / cell_areas(X, Y)


def fill_rib(X, Y, rib_width_um, thick_um, angle_deg, slab_h_um):
    tan_a = np.tan(np.radians(angle_deg))
    a = slab_area(X, Y, 0.0, slab_h_um)
    a = a + band_area(X, Y, 0.0, rib_width_um / 2, tan_a, slab_h_um, thick_um, slab_h_um)
    return np.clip(a / cell_areas(X, Y), 0, 1)


def fill_slot(X, Y, rail_width_um, gap_um, thick_um, angle_deg):
    tan_a = np.tan(np.radians(angle_deg))
    a = np.zeros((len(X) - 1, len(Y) - 1))
    for sign in (+1, -1):
        cx = sign * (gap_um / 2 + rail_width_um / 2)
        a = a + band_area(X, Y, cx, rail_width_um / 2, tan_a, 0.0, thick_um, 0.0)
    return np.clip(a / cell_areas(X, Y), 0, 1)  # same overlap guard as the original


# ---------------------------------------------------------------------------
# epsfunc factories
# ---------------------------------------------------------------------------
def epsfunc_strip_ce(X, Y, width_um, thick_um, angle_deg, n_core=N_SI):
    return _make_epsfunc(X, Y, fill_strip(X, Y, width_um, thick_um, angle_deg), n_core)


def epsfunc_rib_ce(X, Y, rib_width_um, thick_um, angle_deg, slab_h_um=rw.SLAB_H):
    return _make_epsfunc(X, Y, fill_rib(X, Y, rib_width_um, thick_um, angle_deg, slab_h_um), N_SI)


def epsfunc_slot_ce(X, Y, rail_width_um, gap_um, thick_um, angle_deg):
    return _make_epsfunc(X, Y, fill_slot(X, Y, rail_width_um, gap_um, thick_um, angle_deg), N_SI)


# ---------------------------------------------------------------------------
# Solvers with the same signatures as the originals, so they plug into
# fullvector_sensitivity.neff_ng / sensitivities via solve_fn=
# ---------------------------------------------------------------------------
def solve_neff_strip_ce(width_um, thick_um, angle_deg, wl=WL0, polarization='TE', n_modes=4,
                        X=fv.X, Y=fv.Y):
    ep = epsfunc_strip_ce(X, Y, width_um, thick_um, angle_deg)
    return fv.solve_mode_from_eps(ep, wl, polarization, n_modes, X, Y).neff.real


def solve_neff_rib_ce(width_um, thick_um, angle_deg, wl=WL0, polarization='TE', n_modes=4,
                      X=rw.X_RIB, Y=rw.Y_RIB, slab_h_um=rw.SLAB_H):
    ep = epsfunc_rib_ce(X, Y, width_um, thick_um, angle_deg, slab_h_um)
    return fv.solve_mode_from_eps(ep, wl, polarization, n_modes, X, Y).neff.real


def solve_neff_slot_ce(width_um, thick_um, angle_deg, wl=WL0, polarization='TE', n_modes=6,
                       X=sw.X_SLOT, Y=sw.Y_SLOT):
    rail_width_um, gap_um = sw.rail_gap_from_width(width_um)
    ep = epsfunc_slot_ce(X, Y, rail_width_um, gap_um, thick_um, angle_deg)
    return fv.solve_mode_from_eps(ep, wl, polarization, n_modes, X, Y).neff.real


def solve_neff_swg_ce(width_um, thick_um, angle_deg, wl=WL0, polarization='TE', n_modes=4,
                      duty_cycle=swg.DUTY_CYCLE):
    ep = epsfunc_strip_ce(fv.X, fv.Y, width_um, thick_um, angle_deg,
                          n_core=swg.n_eff_ordinary(duty_cycle))
    return fv.solve_mode_from_eps(ep, wl, polarization, n_modes, fv.X, fv.Y).neff.real


# Side-by-side switch: averaging -> geometry -> solve_fn(width_um, thick_um, angle_deg, wl, pol)
SOLVERS = {
    'original': dict(strip=fv.solve_neff, rib=rw.solve_neff_rib,
                     slot=sw.solve_neff_slot, swg=swg.solve_neff_swg),
    'cell_exact': dict(strip=solve_neff_strip_ce, rib=solve_neff_rib_ce,
                       slot=solve_neff_slot_ce, swg=solve_neff_swg_ce),
}


def solve_fn(averaging, geometry, slab_h_um=None):
    f = SOLVERS[averaging][geometry]
    if slab_h_um is not None:
        assert geometry == 'rib'
        f = partial(f, slab_h_um=slab_h_um)
    return f
