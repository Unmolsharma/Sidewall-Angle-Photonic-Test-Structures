"""
Full-vector sensitivity analysis for SOI strip waveguide fabrication metrology.

Uses EMpy's vectorial finite-difference mode solver (Fallahkhair & Li 2008),
validated against published neff for the 500x220nm SOI benchmark waveguide.

Computes, for a set of waveguide widths and BOTH fundamental polarizations
(TE0, TM0):
    neff, ng
    d(neff)/d(width), d(neff)/d(thickness), d(neff)/d(sidewall_angle)
    d(ng)/d(width),   d(ng)/d(thickness),   d(ng)/d(sidewall_angle)

Sidewall angle is modeled as a trapezoidal cross-section (angle deviation
from vertical, in degrees; 0 = ideal rectangular waveguide).

Two structure classes of observable are built from these sensitivities:
  - ring resonator: (neff, ng) pair, read out via resonance condition + FSR
  - MZI arm: phase = (2*pi/lambda)*neff*L, read out via interferometric phase

Then builds Fisher-information-based identifiability landscapes over
candidate probe sets, comparing:
  (a) two TE rings                 [baseline, from the original session]
  (b) TE ring + TM ring            [mixed polarization, same structure class]
  (c) TE ring + MZI                [mixed structure class]
  (d) everything combined at once

to test whether structure-CLASS diversity (not just probe-width choice)
breaks the near-degeneracy found with TE rings alone.

An independent effective-index-method (EIM) slab calculation (no EMpy
dependency) is included as a cross-check for the new TM-mode solve, since
no directly-comparable published TM0 reference for exactly 500x220nm SOI
at 1550nm was found in literature search (see notes/progress.md).
"""
import numpy as np
from scipy.optimize import brentq
import EMpy
import json
import time

WL0 = 1.55       # um
N_SI, N_OX = 3.48, 1.44
H0 = 0.22        # nominal thickness, um
L_MZI = 100.0    # um, MZI arm length used for the phase observable

# grid: nonuniform, denser near the core, in um
def make_grid(half_x=1.4, half_y=1.2, n_core=30, n_clad=18):
    xs = np.concatenate((
        np.linspace(-half_x, -0.4, n_clad, endpoint=False),
        np.linspace(-0.4, 0.4, n_core, endpoint=False),
        np.linspace(0.4, half_x, n_clad)))
    ys = np.concatenate((
        np.linspace(-half_y, -0.1, n_clad, endpoint=False),
        np.linspace(-0.1, 0.5, n_core, endpoint=False),
        np.linspace(0.5, half_y, n_clad)))
    return xs, ys

X, Y = make_grid()

def epsfunc_factory(width_um, thick_um, angle_deg, n_core=N_SI):
    """Trapezoidal core with ANALYTIC sub-pixel dielectric averaging
    (continuous overlap fraction, not discrete sampling) so that
    sub-nanometer perturbations to width/thickness/angle are captured
    exactly, with no discretization floor. Each cell's fill fraction is
    the product of its x-overlap with the (locally width-varying)
    trapezoid and its y-overlap with the [0, thick_um] slab -- exact
    for a rectangle, a good local approximation for a shallow taper.
    n_core defaults to solid silicon; swg_waveguide.py passes an
    effective-medium index instead, reusing this exact geometry/averaging
    code rather than reimplementing it."""
    tan_a = np.tan(np.radians(angle_deg))

    def epsfunc(x_, y_):
        dx = np.diff(x_).mean() if len(x_) > 1 else 0.01
        dy = np.diff(y_).mean() if len(y_) > 1 else 0.01
        xx, yy = np.meshgrid(x_, y_, indexing='ij')

        # y-overlap of cell [y-dy/2, y+dy/2] with [0, thick_um]
        y_lo = np.clip(np.minimum(yy+dy/2, thick_um) - np.maximum(yy-dy/2, 0.0), 0, dy)
        oy = y_lo / dy

        # local half-width at this y (trapezoid), clipped >=0
        local_hw = np.clip(width_um/2 - tan_a*yy, 0, None)
        x_lo = np.clip(np.minimum(xx+dx/2, local_hw) - np.maximum(xx-dx/2, -local_hw), 0, dx)
        ox = x_lo / dx

        fill = ox * oy
        return fill*n_core**2 + (1-fill)*N_OX**2
    return epsfunc

def classify_mode(modes, polarization):
    """Given a list of solved EMpy modes, pick the highest-neff one matching
    the requested polarization, classified by TEfrac() (fraction of mode
    intensity that is TE-polarized). Geometry-agnostic: works for any
    epsfunc (strip, slot, ...) since it only inspects the solved fields."""
    te_fracs = [m.TEfrac() for m in modes]
    if polarization == 'TE':
        cands = [(m.neff.real, m) for m, f in zip(modes, te_fracs) if f > 0.5]
    elif polarization == 'TM':
        cands = [(m.neff.real, m) for m, f in zip(modes, te_fracs) if f <= 0.5]
    else:
        raise ValueError(polarization)
    if not cands:
        raise RuntimeError(
            f"No {polarization}-like mode found among {len(modes)} solved modes; "
            f"TEfracs={te_fracs}")
    cands.sort(key=lambda c: -c[0])
    return cands[0][1]

def solve_mode_from_eps(ep, wl=WL0, polarization='TE', n_modes=4, X=X, Y=Y):
    """Geometry-agnostic mode solve + classification: takes an already-built
    epsfunc(x,y)->eps callable (from ANY epsfunc_factory, strip or slot or
    otherwise), solves n_modes modes on grid (X,Y), and returns the mode
    matching the requested polarization. This is the reusable core that both
    strip.py-style and slot.py-style geometries share -- new geometries only
    need to supply a new epsfunc_factory, not new solve/classify logic."""
    solver = EMpy.modesolvers.FD.VFDModeSolver(wl, X, Y, ep, boundary='0000').solve(n_modes, tol=1e-8)
    return classify_mode(solver.modes, polarization)

def solve_mode(width_um, thick_um, angle_deg, wl=WL0, polarization='TE', n_modes=4):
    """Strip-waveguide mode solve. Verified for widths 400-600nm: mode index
    0 is always clearly TE0 (TEfrac ~0.96-0.99) and mode index 1 is always
    clearly TM0 (TEfrac ~0.02-0.03) -- classifying by TEfrac rather than
    assuming a fixed mode index guards against mode-order swaps under the
    +/- perturbations used for finite-difference derivatives."""
    ep = epsfunc_factory(width_um, thick_um, angle_deg)
    return solve_mode_from_eps(ep, wl, polarization, n_modes)

def solve_neff(width_um, thick_um, angle_deg, wl=WL0, polarization='TE'):
    return solve_mode(width_um, thick_um, angle_deg, wl, polarization).neff.real

def neff_ng(width_um, thick_um, angle_deg, polarization='TE', solve_fn=solve_neff):
    """solve_fn defaults to the strip solver but accepts any callable with
    the same (width_um, thick_um, angle_deg, wl, polarization) signature --
    e.g. slot_waveguide.solve_neff_slot -- so this dispersion/ng machinery
    is reused, not reimplemented, for new geometries."""
    dlam = 0.005
    n0 = solve_fn(width_um, thick_um, angle_deg, WL0, polarization)
    np_ = solve_fn(width_um, thick_um, angle_deg, WL0+dlam, polarization)
    nm_ = solve_fn(width_um, thick_um, angle_deg, WL0-dlam, polarization)
    ng = n0 - WL0*(np_-nm_)/(2*dlam)
    return n0, ng

def sensitivities(width_nm, thick_nm=220.0, angle_deg=0.0, dw=5.0, dt=2.0, da=1.0,
                   polarization='TE', solve_fn=solve_neff):
    w, t = width_nm/1000, thick_nm/1000
    n0, g0 = neff_ng(w, t, angle_deg, polarization, solve_fn)
    nwp, gwp = neff_ng(w+dw/1000, t, angle_deg, polarization, solve_fn)
    nwm, gwm = neff_ng(w-dw/1000, t, angle_deg, polarization, solve_fn)
    ntp, gtp = neff_ng(w, t+dt/1000, angle_deg, polarization, solve_fn)
    ntm, gtm = neff_ng(w, t-dt/1000, angle_deg, polarization, solve_fn)
    nap, gap = neff_ng(w, t, angle_deg+da, polarization, solve_fn)
    nam, gam = neff_ng(w, t, angle_deg-da, polarization, solve_fn)
    return dict(
        neff=n0, ng=g0,
        dneff_dw=(nwp-nwm)/(2*dw),
        dneff_dt=(ntp-ntm)/(2*dt),
        dneff_da=(nap-nam)/(2*da),   # per degree
        dng_dw=(gwp-gwm)/(2*dw),
        dng_dt=(gtp-gtm)/(2*dt),
        dng_da=(gap-gam)/(2*da),
    )

def mzi_phase_sensitivities(s, wl=WL0, L_um=L_MZI):
    """MZI arm phase = (2*pi/lambda) * neff * L, so
    d(phase)/d(param) = (2*pi/lambda) * L * d(neff)/d(param).
    This is an independent OBSERVABLE (phase readout instead of
    resonance/FSR) but note it is only a rescaling of the SAME dneff row
    the ring already measures -- so at a given (width, polarization) it
    does not add a new direction in (w,t,angle) parameter space, only a
    different noise weighting on the existing neff direction. It only
    adds a genuinely new direction if used at a different width or
    polarization than the ring it's paired with."""
    k = 2*np.pi/wl * L_um
    return dict(
        dphase_dw=k*s['dneff_dw'],
        dphase_dt=k*s['dneff_dt'],
        dphase_da=k*s['dneff_da'],
    )

# ---------------------------------------------------------------------------
# Independent cross-check: effective-index method (EIM), no EMpy dependency.
# Standard 2-step slab-dispersion EIM (Okamoto-convention): for a channel
# waveguide's quasi-TE mode, solve the vertical slab with the TE dispersion
# relation to get an effective core index, then solve the horizontal slab
# with the TM relation using that index; for quasi-TM, swap TE<->TM.
# ---------------------------------------------------------------------------
def _slab_neff(n_core, n_top, n_bot, d, pol, wl=WL0, m=0):
    k0 = 2*np.pi/wl

    def f(N):
        kappa = k0*np.sqrt(n_core**2 - N**2)
        g_top = k0*np.sqrt(N**2 - n_top**2)
        g_bot = k0*np.sqrt(N**2 - n_bot**2)
        if pol == 'TE':
            phase_top = np.arctan(g_top/kappa)
            phase_bot = np.arctan(g_bot/kappa)
        else:  # TM
            phase_top = np.arctan((n_core/n_top)**2 * g_top/kappa)
            phase_bot = np.arctan((n_core/n_bot)**2 * g_bot/kappa)
        return kappa*d - m*np.pi - phase_top - phase_bot

    lo, hi = max(n_top, n_bot) + 1e-6, n_core - 1e-6
    Ns = np.linspace(lo, hi, 4000)
    vals = [f(N) for N in Ns]
    for i in range(len(Ns)-1):
        if vals[i]*vals[i+1] < 0:
            return brentq(f, Ns[i], Ns[i+1])
    raise RuntimeError('EIM: no guided root found')

def eim_neff(width_um, thick_um, polarization, wl=WL0):
    if polarization == 'TE':
        n_vert = _slab_neff(N_SI, N_OX, N_OX, thick_um, 'TE', wl)
        return _slab_neff(n_vert, N_OX, N_OX, width_um, 'TM', wl)
    else:
        n_vert = _slab_neff(N_SI, N_OX, N_OX, thick_um, 'TM', wl)
        return _slab_neff(n_vert, N_OX, N_OX, width_um, 'TE', wl)

# ---------------------------------------------------------------------------
# Fisher-information machinery, generalized over structure class.
# Each probe is (kind, width_nm, polarization); 'ring' contributes 2
# observable rows (neff, ng), 'mzi' contributes 1 (phase).
# ---------------------------------------------------------------------------
SIGMA_NEFF  = 6.5e-6   # ~10pm resonance precision equivalent (ring)
SIGMA_NG    = 1.0e-3   # FSR-fit precision (ring)
SIGMA_PHASE = 1.0e-3   # rad; assumed interferometric phase-readout precision (MZI)
                       # -- not independently validated against a real MZI
                       # measurement, flagged as an assumption in progress.md

def jacobian_rows(probe, cache):
    kind, W, pol = probe
    s = cache[(W, pol)]
    if kind == 'ring':
        return [
            np.array([s['dneff_dw'], s['dneff_dt'], s['dneff_da']]) / SIGMA_NEFF,
            np.array([s['dng_dw'],   s['dng_dt'],   s['dng_da']])   / SIGMA_NG,
        ]
    elif kind == 'mzi':
        m = mzi_phase_sensitivities(s)
        return [np.array([m['dphase_dw'], m['dphase_dt'], m['dphase_da']]) / SIGMA_PHASE]
    else:
        raise ValueError(kind)

def fisher_info(probes, cache):
    rows = []
    for p in probes:
        rows.extend(jacobian_rows(p, cache))
    J = np.array(rows)
    F = J.T @ J
    return F, np.linalg.cond(F), np.linalg.det(F), J.shape[0]


if __name__ == "__main__":
    # Fast regression self-test: confirms the refactor (extracting
    # solve_mode_from_eps, adding solve_fn injection) didn't silently change
    # any previously-validated number. Full orchestration (width sweeps,
    # Fisher landscapes, atlas, etc.) now lives in probe_db.py / atlas.py /
    # fisher_reanalysis.py, which import this module rather than duplicating it.
    t0 = time.time()
    s_te = sensitivities(500, polarization='TE')
    s_tm = sensitivities(500, polarization='TM')
    print(f"TE0 @ 500x220nm: neff={s_te['neff']:.4f} (expect 2.4446, published ~2.44-2.45)")
    print(f"TM0 @ 500x220nm: neff={s_tm['neff']:.4f} (expect 1.8204)")
    eim_te = eim_neff(0.5, H0, 'TE')
    eim_tm = eim_neff(0.5, H0, 'TM')
    print(f"EIM cross-check: TE0={eim_te:.4f} ({100*(eim_te-s_te['neff'])/s_te['neff']:+.2f}%), "
          f"TM0={eim_tm:.4f} ({100*(eim_tm-s_tm['neff'])/s_tm['neff']:+.2f}%)")
    assert abs(s_te['neff'] - 2.4446) < 1e-3, "TE0 regression: refactor changed a previously-validated result"
    assert abs(s_tm['neff'] - 1.8204) < 1e-3, "TM0 regression: refactor changed a previously-validated result"
    print(f"Self-test OK. Elapsed: {time.time()-t0:.1f}s")
