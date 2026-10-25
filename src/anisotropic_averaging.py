"""
Anisotropic (tensor) sub-pixel dielectric averaging, built ALONGSIDE the
original factories and cellexact_averaging.py, both left unchanged, so all
three averagings run side by side (see SOLVERS / solve_fn below).

Why: cell-exact averaging fills every partly filled cell with the exact
silicon fraction f but uses the arithmetic mean
    eps_par = f*eps_core + (1 - f)*eps_clad
for every field component. That is right for E parallel to the interface but
wrong for E normal to it, which wants the harmonic mean
    eps_perp = 1 / (f/eps_core + (1 - f)/eps_clad).
With the wrong mean, n_eff(boundary position) carries a sawtooth whose period
is one cell and whose SLOPE does not shrink with the cell size, so the
derivatives are not mesh-converged.

What this module does (Farjadpour et al., Opt. Lett. 31, 2972 (2006);
Kottke et al., PRE 77, 036611 (2008); the Meep scheme):
    eps_tensor = eps_par * (I - P) + eps_perp * P,      eps_zz = eps_par
where P = n n^T is the projector on the local interface normal n (in the
x-y cross-section; z is the propagation direction, always parallel to the
interfaces). Equivalently eps^-1 is averaged for the normal component and
eps for the parallel ones.

Choices made here:
  - Material fractions per cell are EXACT (same band/slab area integrals as
    cellexact_averaging), so silicon is still conserved exactly.
  - The normal is analytic, not sampled: each polygon edge of the geometry is
    clipped to each cell, and P is the length-weighted mean of n n^T over the
    edge pieces inside the cell. This is sign-free, exact for a cell crossed
    by one straight interface, a smooth blend at corners, and continuous in
    the geometry parameters (an edge leaving a cell has its in-cell length go
    to zero continuously), so it adds no new kinks to the derivatives.
  - More than two materials per cell (only used for the air-clad check:
    silicon, oxide below y = 0, air above) are handled with the same formula,
    eps_par = sum f_i eps_i and eps_perp = 1 / sum(f_i / eps_i), with the
    oxide/air plane included as an interface.

Geometry conventions (trapezoid law, rib slab, slot rail/gap coupling, SWG
effective-medium core), grids, mode solve and TEfrac classification are all
reused unchanged from the original modules. EMpy's VFDModeSolver takes the
tensor as an (Nx, Ny, 5) array (eps_xx, eps_xy, eps_yx, eps_yy, eps_zz).
"""
from functools import partial

import numpy as np

import fullvector_sensitivity as fv
import rib_waveguide as rw
import slot_waveguide as sw
import swg_waveguide as swg
import cellexact_averaging as ce

N_SI, N_OX, WL0 = fv.N_SI, fv.N_OX, fv.WL0
N_AIR = 1.0
BIG = 1e3  # um; edges that run to "infinity" (slab, oxide/air plane) end far outside every grid


# ---------------------------------------------------------------------------
# Interface edges and the per-cell normal projector
# ---------------------------------------------------------------------------
def _clip_lengths(X, Y, p, q):
    """Length of segment p->q inside every cell [X[i],X[i+1]] x [Y[j],Y[j+1]]
    (Liang-Barsky, vectorized over cells)."""
    xl, xh, yl, yh = ce._cell_edges(X, Y)
    d = np.subtract(q, p, dtype=float)
    seg = np.hypot(*d)
    lo = np.zeros_like(xl)
    hi = np.ones_like(xl)
    for k, (a, b) in enumerate(((xl, xh), (yl, yh))):
        if d[k] != 0.0:
            s1, s2 = (a - p[k]) / d[k], (b - p[k]) / d[k]
            lo = np.maximum(lo, np.minimum(s1, s2))
            hi = np.minimum(hi, np.maximum(s1, s2))
        else:
            outside = (p[k] < a) | (p[k] > b)
            hi = np.where(outside, -1.0, hi)
    return np.clip(hi - lo, 0.0, None) * seg


def normal_projector(X, Y, edges):
    """Length-weighted mean of n n^T over the edge pieces inside each cell.
    Returns (Pxx, Pxy, Pyy); zero in cells no edge crosses."""
    shape = (len(X) - 1, len(Y) - 1)
    L = np.zeros(shape)
    Pxx, Pxy, Pyy = np.zeros(shape), np.zeros(shape), np.zeros(shape)
    for p, q in edges:
        l = _clip_lengths(X, Y, p, q)
        dx, dy = np.subtract(q, p, dtype=float) / np.hypot(*np.subtract(q, p, dtype=float))
        nx, ny = dy, -dx
        L += l
        Pxx += l * nx * nx
        Pxy += l * nx * ny
        Pyy += l * ny * ny
    nz = L > 0
    for P in (Pxx, Pxy, Pyy):
        P[nz] /= L[nz]
    return Pxx, Pxy, Pyy


def trapezoid_edges(cx, hw0, tan_a, y0, y1):
    """Edges of the trapezoid |x - cx| <= hw0 - tan_a*(y - y0), y0 <= y <= y1."""
    hwt = max(hw0 - tan_a * (y1 - y0), 0.0)
    bl, br, tr, tl = (cx - hw0, y0), (cx + hw0, y0), (cx + hwt, y1), (cx - hwt, y1)
    return [(bl, br), (br, tr), (tr, tl), (tl, bl)]


def strip_edges(width_um, thick_um, angle_deg, cx=0.0):
    return trapezoid_edges(cx, width_um / 2, np.tan(np.radians(angle_deg)), 0.0, thick_um)


def rib_edges(rib_width_um, thick_um, angle_deg, slab_h_um):
    if slab_h_um <= 0.0:
        return strip_edges(rib_width_um, thick_um, angle_deg)
    tan_a = np.tan(np.radians(angle_deg))
    hw0 = rib_width_um / 2
    hwt = max(hw0 - tan_a * (thick_um - slab_h_um), 0.0)
    s, t = slab_h_um, thick_um
    return [((-BIG, 0.0), (BIG, 0.0)),            # slab bottom
            ((-BIG, s), (-hw0, s)), ((hw0, s), (BIG, s)),   # slab top, outside the rib
            ((hw0, s), (hwt, t)), ((-hw0, s), (-hwt, t)),   # rib sidewalls
            ((-hwt, t), (hwt, t))]                # rib top


def slot_edges(rail_width_um, gap_um, thick_um, angle_deg):
    e = []
    for sign in (+1, -1):
        e += strip_edges(rail_width_um, thick_um, angle_deg, cx=sign * (gap_um / 2 + rail_width_um / 2))
    return e


# ---------------------------------------------------------------------------
# Tensor assembly
# ---------------------------------------------------------------------------
def tensor_eps(fracs, eps_vals, P):
    """fracs: list of per-cell area fractions (sum to 1), eps_vals: matching
    permittivities, P: (Pxx, Pxy, Pyy). Returns the (Nx, Ny, 5) EMpy tensor."""
    e_par = sum(f * e for f, e in zip(fracs, eps_vals))
    e_perp = 1.0 / sum(f / e for f, e in zip(fracs, eps_vals))
    de = e_perp - e_par
    Pxx, Pxy, Pyy = P
    return np.stack([e_par + Pxx * de, Pxy * de, Pxy * de, e_par + Pyy * de, e_par], axis=-1)


def _make_epsfunc(X, Y, eps):
    xc_ref = (X[:-1] + X[1:]) / 2
    yc_ref = (Y[:-1] + Y[1:]) / 2

    def epsfunc(x_, y_):
        # guard: these per-cell values are only valid on the grid they were built for
        assert len(x_) == len(xc_ref) and len(y_) == len(yc_ref), "epsfunc called on a different grid"
        assert np.allclose(x_, xc_ref) and np.allclose(y_, yc_ref), "epsfunc called on a different grid"
        return eps
    return epsfunc


def _two_material(X, Y, fill, edges, n_core):
    return tensor_eps([fill, 1 - fill], [n_core**2, N_OX**2], normal_projector(X, Y, edges))


# ---------------------------------------------------------------------------
# epsfunc factories (same geometry conventions as the original factories)
# ---------------------------------------------------------------------------
def epsfunc_strip_an(X, Y, width_um, thick_um, angle_deg, n_core=N_SI):
    fill = ce.fill_strip(X, Y, width_um, thick_um, angle_deg)
    return _make_epsfunc(X, Y, _two_material(X, Y, fill, strip_edges(width_um, thick_um, angle_deg), n_core))


def epsfunc_rib_an(X, Y, rib_width_um, thick_um, angle_deg, slab_h_um=rw.SLAB_H):
    fill = ce.fill_rib(X, Y, rib_width_um, thick_um, angle_deg, slab_h_um)
    return _make_epsfunc(X, Y, _two_material(X, Y, fill, rib_edges(rib_width_um, thick_um, angle_deg, slab_h_um), N_SI))


def epsfunc_slot_an(X, Y, rail_width_um, gap_um, thick_um, angle_deg):
    fill = ce.fill_slot(X, Y, rail_width_um, gap_um, thick_um, angle_deg)
    return _make_epsfunc(X, Y, _two_material(X, Y, fill, slot_edges(rail_width_um, gap_um, thick_um, angle_deg), N_SI))


def strip_air_clad_fracs(X, Y, width_um, thick_um, angle_deg):
    """Silicon / oxide (y < 0) / air (y > 0, outside the core) area fractions."""
    area = ce.cell_areas(X, Y)
    f_si = ce.fill_strip(X, Y, width_um, thick_um, angle_deg)
    f_ox = ce.slab_area(X, Y, -BIG, 0.0) / area
    return f_si, f_ox, np.clip(1 - f_si - f_ox, 0.0, None)


def epsfunc_strip_air_clad(X, Y, width_um, thick_um, angle_deg, averaging='anisotropic'):
    """Strip on oxide with AIR top cladding. 'anisotropic' or 'cell_exact'
    (isotropic arithmetic mean with the same exact fractions)."""
    f_si, f_ox, f_air = strip_air_clad_fracs(X, Y, width_um, thick_um, angle_deg)
    eps_vals = [N_SI**2, N_OX**2, N_AIR**2]
    if averaging == 'cell_exact':
        return _make_epsfunc(X, Y, f_si * eps_vals[0] + f_ox * eps_vals[1] + f_air * eps_vals[2])
    hw0 = width_um / 2
    edges = strip_edges(width_um, thick_um, angle_deg) + [((-BIG, 0.0), (-hw0, 0.0)), ((hw0, 0.0), (BIG, 0.0))]
    return _make_epsfunc(X, Y, tensor_eps([f_si, f_ox, f_air], eps_vals, normal_projector(X, Y, edges)))


# ---------------------------------------------------------------------------
# Solvers with the same signatures as the originals (plug into
# fullvector_sensitivity.neff_ng / sensitivities via solve_fn=)
# ---------------------------------------------------------------------------
def solve_neff_strip_an(width_um, thick_um, angle_deg, wl=WL0, polarization='TE', n_modes=4,
                        X=fv.X, Y=fv.Y):
    ep = epsfunc_strip_an(X, Y, width_um, thick_um, angle_deg)
    return fv.solve_mode_from_eps(ep, wl, polarization, n_modes, X, Y).neff.real


def solve_neff_rib_an(width_um, thick_um, angle_deg, wl=WL0, polarization='TE', n_modes=4,
                      X=rw.X_RIB, Y=rw.Y_RIB, slab_h_um=rw.SLAB_H):
    ep = epsfunc_rib_an(X, Y, width_um, thick_um, angle_deg, slab_h_um)
    return fv.solve_mode_from_eps(ep, wl, polarization, n_modes, X, Y).neff.real


def solve_neff_slot_an(width_um, thick_um, angle_deg, wl=WL0, polarization='TE', n_modes=6,
                       X=sw.X_SLOT, Y=sw.Y_SLOT):
    rail_width_um, gap_um = sw.rail_gap_from_width(width_um)
    ep = epsfunc_slot_an(X, Y, rail_width_um, gap_um, thick_um, angle_deg)
    return fv.solve_mode_from_eps(ep, wl, polarization, n_modes, X, Y).neff.real


def solve_neff_swg_an(width_um, thick_um, angle_deg, wl=WL0, polarization='TE', n_modes=4,
                      duty_cycle=swg.DUTY_CYCLE, X=fv.X, Y=fv.Y):
    ep = epsfunc_strip_an(X, Y, width_um, thick_um, angle_deg, n_core=swg.n_eff_ordinary(duty_cycle))
    return fv.solve_mode_from_eps(ep, wl, polarization, n_modes, X, Y).neff.real


def solve_neff_strip_air(width_um, thick_um, angle_deg, wl=WL0, polarization='TE', n_modes=4,
                         X=fv.X, Y=fv.Y, averaging='anisotropic'):
    ep = epsfunc_strip_air_clad(X, Y, width_um, thick_um, angle_deg, averaging)
    return fv.solve_mode_from_eps(ep, wl, polarization, n_modes, X, Y).neff.real


# Side-by-side switch: averaging -> geometry -> solve_fn(width_um, thick_um, angle_deg, wl, pol)
SOLVERS = dict(ce.SOLVERS)
SOLVERS['anisotropic'] = dict(strip=solve_neff_strip_an, rib=solve_neff_rib_an,
                              slot=solve_neff_slot_an, swg=solve_neff_swg_an)


def solve_fn(averaging, geometry, slab_h_um=None):
    f = SOLVERS[averaging][geometry]
    if slab_h_um is not None:
        assert geometry == 'rib'
        f = partial(f, slab_h_um=slab_h_um)
    return f
