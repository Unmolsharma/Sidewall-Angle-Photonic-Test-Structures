"""
Session 4, Task 4: subwavelength grating (SWG) waveguide -- the second
candidate structure tried against the PC3 gap.

IMPORTANT LIMITATION (stated up front, not buried): EMpy's solver here is a
2D cross-section eigenmode solver; it cannot represent a z-periodic Bloch
mode directly. A true SWG waveguide (periodic Si segments/gaps ALONG the
propagation direction) needs 3D/Bloch simulation to capture exactly. What's
implemented here is the standard zeroth-order effective-medium
approximation (EMA / Rytov approximation) used in the SWG literature as a
first-pass design tool before full 3D verification: at a sub-wavelength
period Lambda (below the Bragg condition), the periodic Si/cladding
segments are replaced by a single homogeneous anisotropic medium. For the
transverse (ordinary) field components that dominate a weakly-vectorial
guided mode, the Wiener/arithmetic-mean bound applies:
    n_eff_ordinary^2 = f*n_Si^2 + (1-f)*n_clad^2
where f is the silicon duty cycle. This scalar index is then used as the
CORE index of an otherwise-ordinary strip cross-section, reusing
fullvector_sensitivity.epsfunc_factory's n_core parameter (added this
session specifically for this reuse) and the same solve/classify core.

Consequence of the approximation: this captures the right ballpark physics
(a genuinely different, LOWER-index-contrast structure, since
n_eff_ordinary < n_Si) but does NOT capture higher-order Bloch effects
(dispersion near the band edge, segment-to-segment field discontinuities).
Results here should be read as "does a lower-index-contrast homogenized
structure reach new sensitivity directions," not as a quantitatively exact
SWG simulation -- flagged the same way the slot waveguide's lack of a
literature reference was flagged, so this isn't overclaimed.

Duty cycle f and period Lambda are treated as FIXED nominal design
constants (not fabrication-sensitivity DOF), consistent with how slot's
gap was NOT an independent free parameter either (Session 3 coupled it to
width) -- this keeps every structure comparable within the same
(width, thickness, angle) 3-parameter framework used throughout.
"""
import numpy as np
from fullvector_sensitivity import (
    WL0, N_SI, N_OX, epsfunc_factory, solve_mode_from_eps, X, Y,
    neff_ng as _neff_ng_generic, sensitivities as _sensitivities_generic,
)

DUTY_CYCLE = 0.5     # fraction of period that is silicon
PERIOD_UM = 0.25     # um; well below Bragg condition ~ wl/(2*n_eff) ~ 0.32um here

def n_eff_ordinary(f=DUTY_CYCLE):
    return np.sqrt(f*N_SI**2 + (1-f)*N_OX**2)

N_SWG = n_eff_ordinary()

def solve_neff_swg(width_um, thick_um, angle_deg, wl=WL0, polarization='TE', n_modes=4,
                    duty_cycle=DUTY_CYCLE):
    n_core = n_eff_ordinary(duty_cycle)
    ep = epsfunc_factory(width_um, thick_um, angle_deg, n_core=n_core)
    mode = solve_mode_from_eps(ep, wl, polarization, n_modes, X, Y)
    return mode.neff.real

def neff_ng_swg(width_um, thick_um, angle_deg, polarization='TE'):
    return _neff_ng_generic(width_um, thick_um, angle_deg, polarization, solve_fn=solve_neff_swg)

def sensitivities_swg(width_nm, thick_nm=220.0, angle_deg=0.0, dw=5.0, dt=2.0, da=1.0, polarization='TE'):
    return _sensitivities_generic(width_nm, thick_nm, angle_deg, dw, dt, da,
                                   polarization, solve_fn=solve_neff_swg)


if __name__ == "__main__":
    import time
    from fullvector_sensitivity import solve_neff as strip_solve_neff
    t0 = time.time()

    print("--- SWG (effective-medium) waveguide validation ---")
    print(f"n_eff_ordinary (f={DUTY_CYCLE}) = {N_SWG:.4f}  (vs solid Si {N_SI}, oxide {N_OX})")
    bragg_limit = WL0/(2*2.4)
    print(f"Period {PERIOD_UM*1000:.0f}nm vs approx Bragg limit {bragg_limit*1000:.0f}nm "
          f"({'OK, sub-wavelength regime' if PERIOD_UM < bragg_limit else 'WARNING: too close to Bragg condition'})")

    print("\nDuty cycle -> 1 consistency check (SWG should reduce EXACTLY to solid strip):")
    strip_ref = strip_solve_neff(0.5, 0.22, 0.0, WL0, 'TE')
    swg_full = solve_neff_swg(0.5, 0.22, 0.0, duty_cycle=1.0)
    print(f"  strip TE0(500nm)              = {strip_ref:.4f}")
    print(f"  SWG(duty_cycle=1) TE0(500nm)  = {swg_full:.4f}  "
          f"(diff = {swg_full-strip_ref:+.4f}, should be ~0)")

    print(f"\nNominal TE0 (duty_cycle={DUTY_CYCLE}, width=500nm) = {solve_neff_swg(0.5, 0.22, 0.0):.4f}")
    print(f"Elapsed: {time.time()-t0:.1f}s")
