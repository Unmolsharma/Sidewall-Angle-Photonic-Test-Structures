# Derivative methods vs the femwell reference

Generated from `logs/derivative_method_comparison.log` (`src/compare_derivatives.py`), data in
`data/derivative_method_comparison.json`. Structures: strip TE 500 x 220 nm and rib TE 500 nm on a
90 nm slab, oxide clad. All derivatives at step h (5 nm width, 2 nm thickness, 1 deg angle).

## Methods

| Method | What it is | Source |
|---|---|---|
| femwell (reference) | Conforming FEM mesh (r12: 12 nm at the silicon, order-2 elements), every perturbed geometry remeshed, central differences | `src/femwell_crosscheck.py`, `data/femwell_derivatives_r12.json` |
| Perturbation | Boundary-perturbation integral (Johnson et al. 2002) from EMpy fields on a boundary-aligned grid; refine 1 and 2 | `src/perturbation_derivatives.py` |
| Aligned FD | EMpy FD on grids whose lines follow the moving boundary (fixed cell counts); refine 1-3; width and thickness only (a tilted wall cannot be grid-aligned) | `src/compare_derivatives.py` |
| Grid-shift N | EMpy FD averaged over N grid offsets (k/N of a core cell, diagonal), production mesh; N1 = plain FD; also N8 on the 45/18 mesh | `src/gridshift_derivatives.py`, `data/gridshift_validation.json` |

The femwell r8 mesh check was not available when this comparison ran (the "r8 +0.00%" entries in
the log are the r12 values reused); it is reported in `results_femwell/SUMMARY.md`.

## femwell reference

| | Strip TE 500 | Rib TE 500 |
|---|---|---|
| n_eff / n_g | 2.44871 / 4.0607 | 2.57070 / 3.7698 |
| dn_eff/dw (per nm) | +1.5861e-3 | +8.4426e-4 |
| dn_eff/dt (per nm) | +3.7238e-3 | +2.7956e-3 |
| dn_eff/da (per deg) | -6.0916e-3 | -1.5134e-3 |
| dn_g/dw | -1.7222e-3 | -4.7339e-5 (near zero) |
| dn_g/dt | +8.1151e-4 | +2.4670e-3 |
| dn_g/da | +6.6160e-3 | +1.8710e-4 (near zero) |
| Step check h/2, h, 2h | all 6 pass, max 0.67% | all 6 pass, max 1.17% (near-zero n_g entries) |
| (dn_eff/da)/(dn_eff/dw) | -3.8406 nm/deg, **R = 1.0001** | -1.7926 nm/deg, **R = 0.7900** |

R = (dn/da) / (-h tan 1 deg dn/dw), h = height of the tilted wall (220 nm strip, 130 nm rib).

## Deviation from femwell, % (n_eff derivatives)

| Derivative | Perturbation r1 / r2 | Aligned FD r1 / r2 / r3 | Cell-exact N1 / N4 / N8 (45/18) | Anisotropic N1 / N4 / N8 (45/18) |
|---|---|---|---|---|
| strip dn/dw | +39.3 / +70.5 | -4.3 / -2.5 / -1.8 | -5.5 / -3.7 / -3.4 (-2.6) | -13.1 / -0.6 / -0.6 (-0.8) |
| strip dn/dt | +24.8 / +51.3 | -1.1 / -0.9 / -0.7 | -1.6 / -0.1 / -0.0 (-0.2) | -4.4 / -0.2 / +0.2 (-0.5) |
| strip dn/da | -39.4 / -70.5 | n/a | +5.6 / +3.7 / +3.4 (+2.2) | +13.3 / +0.7 / +0.7 (+0.6) |
| rib dn/dw | +100.5 / +115.2 | -5.0 / -3.3 / -2.7 | -1.8 / -1.9 / -1.7 (-1.3) | -18.2 / +0.9 / +1.5 (+1.0) |
| rib dn/dt | +17.4 / +35.8 | -0.6 / -0.8 / -0.7 | +0.5 / +1.6 / +1.7 (+0.9) | -3.4 / +0.9 / +1.3 (+0.3) |
| rib dn/da | -64.6 / -95.0 | n/a | -1.7 / -1.0 / -1.6 (-1.2) | +11.4 / -3.6 / -3.8 (-1.3) |

## Deviation from femwell, % (n_g derivatives)

| Derivative | Aligned FD r1 / r2 / r3 | Cell-exact N1 / N4 / N8 (45/18) | Anisotropic N1 / N4 / N8 (45/18) |
|---|---|---|---|
| strip dn_g/dw | +14.7 / +8.4 / +6.3 | +10.4 / +11.4 / +11.0 (+7.8) | +0.7 / +8.3 / +9.3 (+5.9) |
| strip dn_g/dt | +34.8 / +17.9 / +13.5 | +72.6 / +35.5 / +34.7 (+21.7) | +108.6 / +39.6 / +35.0 (+25.8) |
| strip dn_g/da | n/a | -10.8 / -11.8 / -11.1 (-6.5) | -1.3 / -9.3 / -9.9 (-5.3) |
| rib dn_g/dw (near zero) | +174 / +98 / +74 | -19 / +41 / +41 (+19) | -198 / +61 / +81 (+33) |
| rib dn_g/dt | +5.5 / +2.4 / +1.6 | +9.5 / +4.1 / +4.0 (+1.8) | +18.1 / +6.9 / +6.1 (+3.9) |
| rib dn_g/da (near zero) | n/a | +48 / +11 / +16 (+22) | +149 / +6 / -4 (-12) |

The perturbation n_g derivatives are off by 59% to several thousand percent and are not listed.

## Angle/width ratio (dn_eff/da)/(dn_eff/dw), nm/deg

| Method | Strip TE 500 | Rib TE 500 |
|---|---|---|
| femwell | -3.8406 (R 1.0001) | -1.7926 (R 0.7900) |
| Perturbation r1 / r2 | -3.8422 / -3.8406 | -1.4714 / -1.6241 |
| Cell-exact N1 / N4 / N8 / N8 45/18 | -3.837 / -3.838 / -3.842 / -3.858 | -1.858 / -1.846 / -1.854 / -1.839 |
| Anisotropic N1 / N4 / N8 / N8 45/18 | -3.835 / -3.837 / -3.839 / -3.847 | -1.941 / -1.841 / -1.832 / -1.799 |

## What this shows

1. **femwell is a clean reference.** Every derivative passes the h/2, h, 2h check (within 0.67% for
   strip, 1.17% for rib's near-zero n_g entries), and n_eff is mesh-converged to 1e-5.
2. **Boundary perturbation from EMpy fields fails on absolute values**: +17% to +115% on n_eff
   derivatives, growing with grid refinement instead of converging. The near-boundary fields
   that EMpy reconstructs (E and D at cell centres, extrapolated to the wall) are not accurate enough.
   Its *ratios* are right only where the field weighting is symmetric (strip: -3.8422 / -3.8406 vs
   femwell -3.8406); for the rib it gives -1.47 / -1.62 vs -1.79. Usable only as qualitative support.
3. **Grid-shift averaging works for n_eff derivatives** with anisotropic averaging and N >= 4:
   within 0.2-3.8% of femwell on the production mesh, 0.3-1.3% on 45/18. Plain FD (N1) with
   anisotropic averaging is off by up to 18%, which matches the pilot. Cell-exact grid-shift is within
   0-3.7%.
4. **No EMpy method reproduces femwell's n_g derivatives**: grid-shift is off by 6-40% on the
   non-near-zero n_g derivatives, and aligned FD is converging toward femwell only slowly (strip
   dn_g/dt +34.8% -> +13.5% from r1 to r3).
5. **Aligned FD** (boundary-conforming EMpy grids) converges smoothly toward femwell for width and
   thickness (strip dn/dw -4.3 -> -2.5 -> -1.8%), which supports the diagnosis that the EMpy
   derivative errors come from boundaries crossing cells, not from the eigensolver.
6. **The key ratio.** femwell gives R = 1.0001 for strip and 0.790 for rib. Grid-shift EMpy gives
   strip -3.84 nm/deg (within 0.1%) and rib -1.80 to -1.85 nm/deg (0.4-3.4% from femwell).
