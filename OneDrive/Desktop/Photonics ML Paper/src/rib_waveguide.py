"""
Session 4, Task 4: rib waveguide -- one of two candidate structures tried
against the PC3 gap identified in Session 3 (a specific width/sidewall-angle
combination that strip TE/TM and slot TE don't reach).

Geometry: a full-height rib of width `rib_width_um` sitting on a partially-
etched slab of fixed height SLAB_H, total silicon height `thick_um` at the
rib center (same "thickness" convention as strip/slot: total Si height).
Sidewall angle applies only to the etched rib sidewalls (the slab region is
untapered, since it isn't itself an etched sidewall in the same sense).
Reuses the exact per-cell analytic sub-pixel overlap-averaging technique
from the strip/slot epsfuncs (same style, written out explicitly here since
the two-height-region geometry doesn't reduce to a single call into the
existing strip epsfunc_factory) and the same solve_mode_from_eps/TEfrac
classification core.

Physical consistency check (in place of a literature comparison): as
SLAB_H -> 0, a rib waveguide IS a strip waveguide, so solving with SLAB_H=0
must reproduce the already-validated strip solver's neff. This is checked
before any sensitivity sweep is trusted.
"""
import numpy as np
from fullvector_sensitivity import (
    WL0, N_SI, N_OX, solve_mode_from_eps, neff_ng as _neff_ng_generic,
    sensitivities as _sensitivities_generic,
)

SLAB_H = 0.09  # nominal partial-etch slab height, um (~90nm, common SOI rib PDK value)

def make_rib_grid(half_x=1.4, half_y=1.2, n_core=30, n_clad=18):
    xs = np.concatenate((
        np.linspace(-half_x, -0.4, n_clad, endpoint=False),
        np.linspace(-0.4, 0.4, n_core, endpoint=False),
        np.linspace(0.4, half_x, n_clad)))
    ys = np.concatenate((
        np.linspace(-half_y, -0.1, n_clad, endpoint=False),
        np.linspace(-0.1, 0.5, n_core, endpoint=False),
        np.linspace(0.5, half_y, n_clad)))
    return xs, ys

X_RIB, Y_RIB = make_rib_grid()

def epsfunc_factory_rib(rib_width_um, thick_um, angle_deg, slab_h_um=SLAB_H):
    tan_a = np.tan(np.radians(angle_deg))

    def epsfunc(x_, y_):
        dx = np.diff(x_).mean() if len(x_) > 1 else 0.01
        dy = np.diff(y_).mean() if len(y_) > 1 else 0.01
        xx, yy = np.meshgrid(x_, y_, indexing='ij')

        # slab: silicon for all x, y in [0, slab_h_um]
        y_lo_slab = np.clip(np.minimum(yy+dy/2, slab_h_um) - np.maximum(yy-dy/2, 0.0), 0, dy)
        oy_slab = y_lo_slab / dy

        # rib: silicon within local_hw(y) for y in [slab_h_um, thick_um]
        y_lo_rib = np.clip(np.minimum(yy+dy/2, thick_um) - np.maximum(yy-dy/2, slab_h_um), 0, dy)
        oy_rib = y_lo_rib / dy
        local_hw = np.clip(rib_width_um/2 - tan_a*(yy - slab_h_um), 0, None)
        x_lo = np.clip(np.minimum(xx+dx/2, local_hw) - np.maximum(xx-dx/2, -local_hw), 0, dx)
        ox = x_lo / dx

        fill = np.clip(oy_slab + oy_rib*ox, 0, 1)
        return fill*N_SI**2 + (1-fill)*N_OX**2
    return epsfunc

def solve_neff_rib(width_um, thick_um, angle_deg, wl=WL0, polarization='TE', n_modes=4,
                    X=X_RIB, Y=Y_RIB, slab_h_um=SLAB_H):
    ep = epsfunc_factory_rib(width_um, thick_um, angle_deg, slab_h_um)
    mode = solve_mode_from_eps(ep, wl, polarization, n_modes, X, Y)
    return mode.neff.real

def neff_ng_rib(width_um, thick_um, angle_deg, polarization='TE'):
    return _neff_ng_generic(width_um, thick_um, angle_deg, polarization, solve_fn=solve_neff_rib)

def sensitivities_rib(width_nm, thick_nm=220.0, angle_deg=0.0, dw=5.0, dt=2.0, da=1.0, polarization='TE'):
    return _sensitivities_generic(width_nm, thick_nm, angle_deg, dw, dt, da,
                                   polarization, solve_fn=solve_neff_rib)


if __name__ == "__main__":
    import time
    from fullvector_sensitivity import solve_neff as strip_solve_neff
    t0 = time.time()

    print("--- Rib waveguide validation ---")
    print("\nSLAB_H -> 0 consistency check (rib should reduce EXACTLY to strip):")
    strip_ref = strip_solve_neff(0.5, 0.22, 0.0, WL0, 'TE')
    rib_zero_slab = solve_neff_rib(0.5, 0.22, 0.0, slab_h_um=0.0)
    print(f"  strip TE0(500nm)          = {strip_ref:.4f}")
    print(f"  rib(slab_h=0) TE0(500nm)  = {rib_zero_slab:.4f}  "
          f"(diff = {rib_zero_slab-strip_ref:+.4f}, should be ~0)")

    print(f"\nGrid convergence @ rib=500nm, slab_h={SLAB_H*1000:.0f}nm:")
    configs = [
        ('baseline 1.4x1.2 30/18', 1.4, 1.2, 30, 18),
        ('larger domain 2.0x1.8', 2.0, 1.8, 30, 18),
        ('finer core 1.4x1.2 45/18', 1.4, 1.2, 45, 18),
    ]
    vals = []
    for name, hx, hy, nc, ncl in configs:
        X, Y = make_rib_grid(hx, hy, nc, ncl)
        n = solve_neff_rib(0.5, 0.22, 0.0, X=X, Y=Y)
        vals.append(n)
        print(f"  {name:28s} TE0={n:.5f}")
    print(f"  spread: {max(vals)-min(vals):.4f} ({100*(max(vals)-min(vals))/np.mean(vals):.2f}% of mean)")

    print(f"\nNominal TE0 (production grid, rib=500nm) = {solve_neff_rib(0.5, 0.22, 0.0):.4f}")
    print("Literature: rib waveguides are extremely standard (partial-etch SOI PDK "
          "structures); the SLAB_H->0 exact-match-to-strip check above is treated as "
          "the primary validation since it directly confirms the geometry/solve is "
          "correct at a point where the answer is already independently known.")
    print(f"\nElapsed: {time.time()-t0:.1f}s")
