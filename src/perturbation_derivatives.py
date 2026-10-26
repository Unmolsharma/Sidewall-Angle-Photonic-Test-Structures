"""
Boundary-perturbation derivatives of n_eff and n_g (option 1), built
alongside the finite-difference derivatives (nothing existing is changed).

Formula (Johnson et al., PRE 65, 066611 (2002), moving dielectric boundaries in
high index contrast), at fixed frequency, for a boundary displaced by V_n per
unit parameter along the normal pointing from silicon (eps1) into the cladding
(eps2), i.e. V_n > 0 grows the silicon:

    dn_eff/dp = 1/(4P) * Int_boundary V_n [ (eps1 - eps2) |E_par|^2
                                            - (1/eps1 - 1/eps2) |D_perp|^2 ] dl
    P = 1/2 Re Int (E x H*) . z dA

in EMpy's units (D = eps E, D = curl H / (i k), so the vacuum impedance drops
out). E_par includes E_z. The normalization is checked numerically in
calibrate(): the volume form dn/d(eps_core) = Int_core |E|^2 dA / (4P) is
compared with a central difference of n_eff in eps_core (a smooth perturbation,
no moving boundary).

Fields: the nominal geometry (angle 0, so every cross-section is a union of
rectangles) is solved on a grid whose lines pass exactly through every
silicon boundary, so no cell is partly filled and no averaging is involved.
EMpy returns E (and D = eps E) at cell centres; E_par and D_perp on the
boundary are each extrapolated linearly from the two nearest cell centres on
each side and the two one-sided values averaged.

Boundary displacement per parameter (same conventions as the geometry code):
  width w   : every etched sidewall moves outward by dw/2 (strip: both walls;
              rib: the rib walls between slab top and rib top; slot: all four
              rail walls -- rail centres stay fixed because gap = GAP0 - (w - RAIL0))
  thickness : the top surface moves up by dt (rib: slab held fixed, so only
              the rib top moves)
  angle a   : top narrower, pivot at the BOTTOM of the tilted wall (y = 0 for
              strip/slot/SWG, y = slab top for rib); the wall moves inward by
              tan(a)(y - y_pivot), so per degree V_n = -(pi/180)(y - y_pivot)
n_g derivatives: dn/dp at lambda and lambda +- 5 nm on the same grid, then
dn_g/dp = D(lambda) - lambda (D(lambda+) - D(lambda-)) / (2 dlambda), the same
n_g definition as the FD code.
"""
import os
os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('MKL_NUM_THREADS', '1')

import numpy as np

import fullvector_sensitivity as fv
import rib_waveguide as rw
import slot_waveguide as sw
import swg_waveguide as swg
import cellexact_averaging as ce
import anisotropic_averaging as an

WL0, DLAM = fv.WL0, 0.005
DEG = np.pi / 180


# ---------------------------------------------------------------------------
# Grid with lines through every boundary
# ---------------------------------------------------------------------------
def _axis(breaks, core, half, step, n_clad):
    """Nonuniform axis like fv.make_grid (n_clad points per cladding segment,
    core region [core0, core1]) but with the core region split at `breaks`
    and each piece filled with the smallest uniform spacing <= step."""
    c0, c1 = core
    pts = sorted({c0, c1, *[b for b in breaks if c0 < b < c1]})
    ax = [np.linspace(-half, c0, n_clad, endpoint=False)]
    for a, b in zip(pts[:-1], pts[1:]):
        n = max(int(np.ceil((b - a) / step - 1e-9)), 1)
        ax.append(np.linspace(a, b, n, endpoint=False))
    ax.append(np.linspace(c1, half, n_clad))
    return np.concatenate(ax)


def aligned_grid(xbreaks, ybreaks, refine=1.0, slot=False):
    """Production-like spacing (strip/rib: 26.7 x 20 nm, slot: 20 x 20 nm in the
    core region) divided by `refine`; same domain and cladding points."""
    if slot:
        hx, hy, ncl, core_x, dx = 1.6, 1.2, 20, (-0.5, 0.5), 1.0 / 50
    else:
        hx, hy, ncl, core_x, dx = 1.4, 1.2, 18, (-0.4, 0.4), 0.8 / 30
    dy = 0.6 / 30
    ncl = int(round(ncl * min(refine, 2.0)))
    return (_axis(xbreaks, core_x, hx, dx / refine, ncl),
            _axis(ybreaks, (-0.1, 0.5), hy, dy / refine, ncl))


# ---------------------------------------------------------------------------
# Structures: nominal eps on the aligned grid + moving boundaries
# ---------------------------------------------------------------------------
class Structure:
    """geom in strip/rib/slot/swg; W in nm. Boundaries are tuples
    (orientation, coordinate, silicon side, span lo, span hi, eps_si, eps_clad, kind)
    with orientation 'v' (x = coord) or 'h' (y = coord), silicon side -1/+1 along
    the normal axis, kind 'side' (etched sidewall) or 'top'."""

    def __init__(self, geom, W, pol='TE', refine=1.0, clad='oxide', slab_nm=None, wl=WL0):
        self.geom, self.W, self.pol, self.clad = geom, W, pol, clad
        w, t = W / 1000, fv.H0
        self.t = t
        self.slab = (rw.SLAB_H if slab_nm is None else slab_nm / 1000) if geom == 'rib' else 0.0
        n_core = swg.n_eff_ordinary() if geom == 'swg' else fv.N_SI
        e_si = n_core ** 2
        e_cl = (an.N_AIR if clad == 'air' else fv.N_OX) ** 2
        self.e_si, self.e_cl = e_si, e_cl
        if geom == 'slot':
            rail, gap = sw.rail_gap_from_width(w)
            walls = [s * (gap / 2) for s in (-1, 1)] + [s * (gap / 2 + rail) for s in (-1, 1)]
            self.X, self.Y = aligned_grid(walls, [0.0, t], refine, slot=True)
            fill = ce.fill_slot(self.X, self.Y, rail, gap, t, 0.0)
            rails = [(-(gap / 2 + rail), -gap / 2), (gap / 2, gap / 2 + rail)]
        else:
            ybr = [0.0, t] + ([self.slab] if geom == 'rib' else [])
            self.X, self.Y = aligned_grid([-w / 2, w / 2], ybr, refine)
            fill = (ce.fill_rib(self.X, self.Y, w, t, 0.0, self.slab) if geom == 'rib'
                    else ce.fill_strip(self.X, self.Y, w, t, 0.0))
            rails = [(-w / 2, w / 2)]
        frac = fill[(fill > 1e-9) & (fill < 1 - 1e-9)]
        assert frac.size == 0, f'grid not aligned: {frac.size} partly filled cells'
        self.fill = fill
        y_pivot = self.slab
        self.y_pivot = y_pivot
        self.bounds = []
        for lo, hi in rails:
            self.bounds += [('v', lo, +1, y_pivot, t, 'side'), ('v', hi, -1, y_pivot, t, 'side'),
                            ('h', t, -1, lo, hi, 'top')]
        # eps per cell: silicon / cladding above y=0 / oxide below y=0
        yc = (self.Y[:-1] + self.Y[1:]) / 2
        e_bg = np.where(yc[None, :] < 0, fv.N_OX ** 2, e_cl) * np.ones((len(self.X) - 1, 1))
        self.eps = fill * e_si + (1 - fill) * e_bg

    def solve(self, wl=WL0):
        eps = self.eps
        X, Y = self.X, self.Y

        def epsfunc(x_, y_):
            assert len(x_) == len(X) - 1 and len(y_) == len(Y) - 1
            return eps
        n_modes = 6 if self.geom == 'slot' else 4
        return fv.solve_mode_from_eps(epsfunc, wl, self.pol, n_modes, X, Y)


def _centered(H):
    return (H[:-1, :-1] + H[1:, :-1] + H[:-1, 1:] + H[1:, 1:]) / 4


def power(mode, X, Y):
    Ex, Ey = mode.Ex, mode.Ey
    S = 0.5 * np.real(Ex * np.conj(_centered(mode.Hy)) - Ey * np.conj(_centered(mode.Hx)))
    A = np.outer(np.diff(X), np.diff(Y))
    return float(np.sum(S * A))


def _extrap(f1, f2, c1, c2, c0):
    return f1 + (f1 - f2) * (c0 - c1) / (c1 - c2)


def boundary_samples(st, mode, b):
    """For one boundary: positions s along it, segment lengths, |E_par|^2 and
    |D_perp|^2 on the boundary (each the mean of the two one-sided extrapolations),
    plus the one-sided disagreement as an accuracy indicator."""
    orient, c, side, lo, hi, kind = b
    X, Y = st.X, st.Y
    xc, yc = (X[:-1] + X[1:]) / 2, (Y[:-1] + Y[1:]) / 2
    E = dict(x=mode.Ex, y=mode.Ey, z=mode.Ez)
    eps = st.eps
    if orient == 'v':
        i = int(np.argmin(abs(X - c)))
        assert abs(X[i] - c) < 1e-9
        js = [j for j in range(len(yc)) if Y[j] >= lo - 1e-12 and Y[j + 1] <= hi + 1e-12]
        s, ds = yc[js], np.diff(Y)[js]
        par, nrm, cc, idx = ('y', 'z'), 'x', xc, lambda k: (k, js)
        cells = {-1: (i - 1, i - 2), +1: (i, i + 1)}
    else:
        j = int(np.argmin(abs(Y - c)))
        assert abs(Y[j] - c) < 1e-9
        js = [k for k in range(len(xc)) if X[k] >= lo - 1e-12 and X[k + 1] <= hi + 1e-12]
        s, ds = xc[js], np.diff(X)[js]
        par, nrm, cc, idx = ('x', 'z'), 'y', yc, lambda k: (js, k)
        cells = {-1: (j - 1, j - 2), +1: (j, j + 1)}
    vals = {}
    for sd, (k1, k2) in cells.items():
        def ext(F):
            return _extrap(F[idx(k1)], F[idx(k2)], cc[k1], cc[k2], c)
        Ep = [ext(E[p]) for p in par]
        Dn = ext(E[nrm] * eps)
        vals[sd] = (sum(abs(v) ** 2 for v in Ep), abs(Dn) ** 2)
    Epar2 = 0.5 * (vals[-1][0] + vals[1][0])
    Dn2 = 0.5 * (vals[-1][1] + vals[1][1])
    spread = (np.sum(abs(vals[-1][0] - vals[1][0]) * ds), np.sum(abs(vals[-1][1] - vals[1][1]) * ds))
    return s, ds, Epar2, Dn2, spread


def velocity(st, b, s, param):
    """Outward normal displacement of boundary b at positions s, per nm (w, t) or per degree (a)."""
    orient, c, side, lo, hi, kind = b
    if param == 'w':
        return np.full_like(s, 0.5e-3) if kind == 'side' else np.zeros_like(s)
    if param == 't':
        return np.full_like(s, 1e-3) if kind == 'top' else np.zeros_like(s)
    if param == 'a':
        return -DEG * (s - st.y_pivot) if kind == 'side' else np.zeros_like(s)
    raise ValueError(param)


def dn_dp(st, mode, params=('w', 't', 'a'), profile=False):
    P = power(mode, st.X, st.Y)
    out, prof = {}, []
    samples = [(b,) + boundary_samples(st, mode, b) for b in st.bounds]
    d_eps, d_inv = st.e_si - st.e_cl, 1 / st.e_si - 1 / st.e_cl
    for p in params:
        tot = 0.0
        for b, s, ds, Ep2, Dn2, _ in samples:
            tot += np.sum(velocity(st, b, s, p) * (d_eps * Ep2 - d_inv * Dn2) * ds)
        out[p] = tot / (4 * P)
    if profile:
        for b, s, ds, Ep2, Dn2, spread in samples:
            if b[5] == 'side':
                prof.append(dict(x=b[1], y=s.tolist(), dy=ds.tolist(),
                                 weight=((d_eps * Ep2 - d_inv * Dn2) / (4 * P)).tolist(),
                                 one_sided_spread=list(spread)))
        return out, prof
    return out


def derivatives(geom, W, pol='TE', refine=1.0, clad='oxide', slab_nm=None, params=('w', 't', 'a')):
    """n_eff, n_g and dn_eff/dp, dn_g/dp for p in params by boundary perturbation."""
    st = Structure(geom, W, pol, refine, clad, slab_nm)
    Ds, ns = [], []
    for wl in (WL0, WL0 + DLAM, WL0 - DLAM):
        m = st.solve(wl)
        ns.append(m.neff.real)
        Ds.append(dn_dp(st, m, params))
    n0 = ns[0]
    ng = n0 - WL0 * (ns[1] - ns[2]) / (2 * DLAM)
    res = dict(neff=n0, ng=ng)
    for p in params:
        res[f'dneff_d{p}'] = Ds[0][p]
        res[f'dng_d{p}'] = Ds[0][p] - WL0 * (Ds[1][p] - Ds[2][p]) / (2 * DLAM)
    return res


def calibrate(geom='strip', W=500, pol='TE', d=0.01):
    """Normalization check: volume form dn/d(eps_core) vs central difference in eps_core."""
    st = Structure(geom, W, pol)
    m = st.solve()
    P = power(m, st.X, st.Y)
    A = np.outer(np.diff(st.X), np.diff(st.Y))
    E2 = abs(m.Ex) ** 2 + abs(m.Ey) ** 2 + abs(m.Ez) ** 2
    pert = np.sum(st.fill * E2 * A) / (4 * P)
    base = st.eps.copy()
    vals = []
    for s in (1, -1):
        st.eps = base + s * d * st.fill
        vals.append(st.solve().neff.real)
    st.eps = base
    fd = (vals[0] - vals[1]) / (2 * d)
    return dict(perturbation=pert, finite_difference=fd, rel_diff=(pert - fd) / fd)
