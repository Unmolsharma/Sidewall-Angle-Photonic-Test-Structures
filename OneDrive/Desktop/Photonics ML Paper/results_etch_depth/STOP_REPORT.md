# Etch-depth nuisance study — stopped at Step 2 (derivative checks failed)

Date: 2026-10-20. Steps 3 and 4 (Fisher analysis, figures, `report_etch_depth.md`) were **not run**.
The Step 2 step-size check failed. The cause is a bug in the existing dielectric-averaging
code, and it affects the existing strip and rib Jacobians as well as the new ones.
No existing file in `src/`, `data/`, `figures/` or `notes/` was modified.

## Step 0 — how the current code parameterizes the problem

| Item | Current value / method | Where | Source |
|---|---|---|---|
| Unknowns | width w [nm], thickness t [nm], sidewall angle a [deg] | `fullvector_sensitivity.sensitivities` | — |
| SWA convention | a = deviation from vertical, positive means the top is narrower. SWA = 90° − a | `epsfunc_factory` | — |
| **Nominal SWA** | **a = 0, i.e. SWA = 90° (vertical). Not 85° as in your spec.** ±1° central differences go up to 1° re-entrant | `sensitivities(angle_deg=0.0)` | — |
| Nominal t | 220 nm | `H0` | standard SOI |
| Widths | strip 400–600 nm, rib 400–600 nm (TE only), slot rails 150–210 nm | `probe_db.py`, `extend_db.py` | — |
| **Rib slab** | **fixed constant `SLAB_H = 90 nm` (so d = 130 nm)** | `rib_waveguide.py` | "common SOI rib PDK value", no citation |
| **Rib t-derivative** | **t is perturbed with the slab held fixed**, so the etch depth moves with t | `sensitivities_rib` | — |
| Is d / slab treated as known? | **Yes.** Slab thickness is a perfectly known constant. The slab sweep (`rib_slab_sweep.py`) changes the nominal value but never makes it an unknown | | |
| σ_neff | 6.5e-6 ("~10 pm resonance precision equivalent") | `fullvector_sensitivity.py` | **no source (placeholder)** |
| σ_ng | 1.0e-3 ("FSR-fit precision") | same | **no source (placeholder)** |
| σ_phase | 1.0e-3 rad, MZI L = 100 µm | same | **explicitly assumed in code** |
| Rib noise multiplier | 1.0× (conservative) / 0.2× (optimistic), applied to all three σ | `noise_sourcing.py` | "Record low-loss hybrid rib/wire waveguides…" (wide-rib caveat noted there) |
| Slot noise multiplier | 4× | `noise_sourcing.py` | Baehr-Jones et al., Opt. Express 19, 11529 (2011); Opt. Express 19, 26275 (2011) |
| FD steps | dw = 5 nm, dt = 2 nm, da = 1°, central differences | `sensitivities` | — |
| n_g | n_eff − λ·dn_eff/dλ, central difference, dλ = 5 nm. Constant N_Si = 3.48 and N_ox = 1.44, so **material dispersion is ignored** | `neff_ng` | — |
| Production mesh | half_x 1.4 µm, half_y 1.2 µm, 30 core / 18 cladding points per segment. Core cells 26.7 nm (x) × 20 nm (y). PEC boundary `'0000'`, 4 modes, tol 1e-8 | `make_grid`, `make_rib_grid` | — |

Observation on σ_neff (not changed, only flagged): 10 pm / 1550 nm = 6.45e-6 is Δλ/λ. For a ring,
Δλ_res/λ = Δn_eff/n_g, so a 10 pm resonance precision corresponds to σ_neff ≈ n_g × 6.45e-6 ≈ 2.6e-5.
Either the existing value means something else, or it is about 4× optimistic.

## What was built (new files only)

- `src/etch_depth_sims.py`: rib solve in the (w, t, a, d) parameterization, with slab = t − d.
  It reuses `epsfunc_factory_rib` and `solve_mode_from_eps` unchanged.
- `src/etch_depth_checks.py`: step, mesh and dλ check tables.
- `src/etch_depth_fisher.py`: Steps 3–4 (Case A/B Fisher matrices, priors, structure sets,
  greedy E/A design, figures). **Written but deliberately not run.**
- `src/etch_depth_diag_averaging.py`: the diagnostic below.
- Outputs are in `results_etch_depth/data/`.

Checks that passed:
- n_g in the new module matches `fv.neff_ng` to 1.4e-11. That residual is ARPACK random-start
  repeatability, measured at ~1e-14 on n_eff (`ng_method_check.json`).
- The new rib derivatives reproduce the stored database. dn/dw and dn/da match exactly, and the
  stored slab-fixed dn_eff/dt = dt|_d + dd = 3.128e-3 (n_g: 2.854e-3 vs stored 2.855e-3).
  So the reparameterization is correct.

## Step 2 failure

Step check on all 30 ribs (5 widths × slabs 50–150 nm, h/2 / h / 2h): **34 of 240 entries pass
the 2% criterion** (`data/step_check_M0.csv`). dn_eff/dd and dn_g/dd deviate by 36–41% at
**every** structure. The mesh-convergence run was stopped unfinished because it is contaminated
by the same cause.

### Cause: the sub-pixel averaging does not conserve silicon volume

`epsfunc_factory` and `epsfunc_factory_rib` average each cell over a window of
`dx = np.diff(x_).mean()` and `dy = np.diff(y_).mean()`. That is the **mean spacing of the whole
nonuniform grid** (42.9 nm × 36.7 nm), not each cell's own size (26.7 nm × 20 nm in the core).
The windows of neighbouring cells overlap unevenly, so when an interface moves, the effective
silicon changes at the wrong rate, and that rate depends on where the interface sits:

| Interface | Rate of effective Si change per nm of interface motion (ideal 1.00) |
|---|---|
| horizontal (t, slab) | 1.09 at most positions, **0.54** within ±2 nm of a cell centre |
| vertical (w, and SWA via the taper) | 1.24 or **0.62**, alternating every ~13 nm |

Cell centres in y fall at 10, 30, 50, 70, 90, … nm, so **every nominal slab (50, 70, 90, 110,
130, 150 nm) sits exactly in a 0.54× dip**. A step of h = 1–2 nm stays in the dip, while 2h = 4 nm
partly leaves it. That produces the ~40% step dependence. The strip top at 220 nm sits on a
1.09× plateau, so the existing t-derivatives are ~9% high. Strip and rib widths of 400/500/600 nm
sit in 0.62× dips; 450 and 550 straddle a transition.

The Session-1 fix ("analytic overlap averaging fixed the zero thickness derivative") removed the
zero, but left this position-dependent error of up to ~2×.

### Diagnostic: cell-exact averaging, rib 500 nm / slab 90 nm, production grid

| param | obs | existing D(h) | existing max dev | cell-exact D(h) | cell-exact max dev |
|---|---|---|---|---|---|
| w | n_eff | 5.59e-4 | 16.2% | 8.29e-4 | 0.3% |
| w | n_g | 3.3e-6 | 1183% | −5.62e-5 | 5.0% |
| t (fixed d) | n_eff | 4.20e-3 | 10.2% | 4.45e-3 | 0.7% |
| t (fixed d) | n_g | 1.20e-3 | 55.0% | 1.64e-4 | **49.9%** |
| a | n_eff | −9.34e-4 | 11.9% | −1.524e-3 | 0.0% |
| a | n_g | 1.44e-4 | 13.2% | 2.77e-4 | 2.7% |
| d | n_eff | −1.07e-3 | 39.7% | −1.640e-3 | 0.7% |
| d | n_g | 1.66e-3 | 39.8% | 2.534e-3 | 0.8% |

What this shows:
1. The diagnosis is right. With true per-cell bounds the volume slope is exactly 1.000, and every
   n_eff derivative and both d-derivatives are stable to <1%.
2. The existing derivatives are **biased, not just noisy**. dn_eff/da is 63% larger once fixed,
   dn_eff/dd 53%, dn_eff/dw 48%, dn_eff/dt 6%. The per-parameter biases differ, so the
   **direction** of each sensitivity vector changes. The rib result was explicitly a direction
   effect (the "PC3 gap" finding, and the 571× / 11–51× improvements), so those numbers cannot be
   trusted until recomputed. Strip derivatives carry the same bug, with different biases at each width.
3. Even after the fix, two small n_g derivatives still fail: dn_g/dt|_d at 50% and dn_g/dw at 5%.
   dn_g/dt|_d is ~15× smaller than dn_g/dd (near-cancellation between the rib-top and slab
   contributions), so a relative 2% test on it is ill-conditioned. It still needs either a
   Richardson / smaller-step study or an agreed criterion for near-zero entries. I have not
   decided this on your behalf.

## Decisions needed before Steps 3–4 can run

1. **Fix the averaging in the existing geometry code** (strip, rib, slot, SWG all share the
   pattern). This changes previously validated numbers: the TE0 = 2.4446 regression assert, the
   whole `sensitivity_db`, and every Session 2–5 headline. Options: (a) fix it in place and
   regenerate everything; (b) add corrected factories alongside the old ones and rebuild into a
   new data folder, leaving the old results as a record.
2. Agree how to judge near-zero derivatives (for example, 2% relative or an absolute floor tied to
   noise-weighted Fisher contribution), or require a Richardson check for them.
3. Nominal SWA: keep the existing 90° (vertical), or move to your 85°? 85° would also avoid
   stepping into re-entrant geometry.

Context for the prior sweep: a secondary source reports demonstrated final etch-depth variation
of 2.5–3.3 nm (ScienceDirect topic page "Silicon Waveguide"; not verified against the primary
paper). A Boning (MIT) DFM talk attributes radial wafer variation mainly to SOI thickness and
dry-etch depth. These are context only and were not used in any calculation.
