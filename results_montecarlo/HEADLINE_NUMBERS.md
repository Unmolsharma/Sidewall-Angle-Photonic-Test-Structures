# Headline numbers (state as of 2026-10-29)

These are the numbers we would currently put in the paper, each with its source file (paths relative
to the repository root). The strip angle is the full-etch sidewall angle, as a 1-sigma Cramér-Rao bound
in degrees at derivative step h, with [min–max over h/2, h, 2h].

All current CRBs that include strips impose the exact strip symmetry (R = 1); see
`results_montecarlo/SUMMARY.md`. **Noise values are placeholders** (σ_neff 6.5e-6, σ_ng 1e-3,
σ_phase 1e-3 rad, L = 100 µm; none sourced), so every precision below scales with them.

Stages:
- B: + etch depth;
- C: + separate rib angle;
- D: + n_Si, n_ox;
- E: + T.

"Realistic" priors: d 2 nm, n_Si 1e-4, n_ox 1e-3, T 1 K.

## Current numbers

### Physics of the degeneracy

| Number | Value | Source |
|---|---|---|
| Strip R exact by symmetry, (π/180)/tan 1° | analytic, = 0.9999 to 4 digits (mirror symmetry; see `results_montecarlo/SUMMARY.md`, step 2) | `results_femwell/data/crosscheck_and_ratios.json` (femwell strip R 0.99987–1.0084) |
| Strip R, femwell residual about exact | TE −2.7e-4 to +1.9e-3; TM up to 8.5e-3 | `results_montecarlo/data/strips_info.json` |
| Rib R (n_eff), TE 400–600 nm, 90 nm slab | 0.776–0.799 | `results_femwell/data/crosscheck_and_ratios.json` |
| Rib R (n_eff), all slabs 50–150 nm, 300–600 nm | 0.712–0.873 | `results_femwell/data/crosscheck_and_ratios.json` |
| Strips-only sets, any size | ∞ (no angle information) | `results_montecarlo/sym_search/data/search.json` |
| 1° tilt looks like a strip width change of (220 nm·π/180) | 3.84 nm; fitted Monte Carlo slope 3.83–3.84 nm/deg | `results_montecarlo/data/mc_valley.json` |
| Ring resonance error from a 1° strip tilt treated as width | −711 pm | `results_femwell/data/impact_sidewall_angle.json` |
| Mesh check r12 vs r8, every derivative | within 0.30% | `results_femwell/data/meshcheck_r12_vs_r8.json` |

### Structure search (symmetry imposed)

All from `results_montecarlo/sym_search/data/summary_tables.md`, built from `search.json` and
`slab_pairs.json`.

| Number | Value |
|---|---|
| Best one-etch set, 8 structures, D realistic | **0.162 [0.118–0.171]**, 90 nm slab: strip TE 400, TE 600, TM 400, TM 600 + rib 300/400/500/600 nm |
| Same, C no priors | 0.144 [0.105–0.153] |
| Best set by size, D realistic, k = 4/5/6/7/8 | 0.413 / 0.329 / 0.174 / 0.167 / 0.162 |
| Best 8 per slab, D realistic, 50/70/90/110/130/150 nm | 0.359 / 0.187 / 0.162 / 0.795 / 1.21 / 3.81 |
| Two-width rib set (strip TE 5 + rib 300 + 500), D realistic, 50/70/90/110/130/150 nm | **0.367** / 0.425 / 0.561 / 0.802 / 1.22 / 3.83 → **best slab 50 nm** |
| Same, C no priors, 50 nm slab | 0.359 [0.359–0.369] |
| Single rib (strip TE 5 + rib 500), C or D, no priors | ∞ at every slab |
| Single rib, realistic priors | 24.3–43.1 (set by the priors) |

### Priors (strip + rib70)

All from `results_montecarlo/sym_search/data/priors_table.csv`.

| Number | + rib 500 nm | + rib 500 + 300 nm |
|---|---|---|
| B, no priors (shared angle) | 0.334 | 0.00433 |
| C, no priors | ∞ | 0.419 |
| D, no priors | ∞ | 0.451 |
| D, realistic | 25.1 | 0.425 |
| E, realistic | 25.2 | 0.428 |
| C, prior ar − a ≤ 1° | 1.24 | 0.404 |
| D, realistic + ar − a ≤ 2° | 2.41 | 0.421 |

### Monte Carlo (stage D, realistic priors, 2,000 sets per case)

All from `results_montecarlo/data/mc_results.json` and `results_montecarlo/data/mc_fits.csv`.

| Number | (a) best multi-width, slab 90 | (b) strip + rib 300 + 500, slab 50 | (c) strip + rib 500, slab 50 |
|---|---|---|---|
| MC std of fitted strip angle [95% CI] | **0.163 [0.155–0.172]** | **0.374 [0.363–0.385]** | 2.47 (bounded by the ±3° box) |
| CRB (femwell derivatives) | 0.162 | 0.367 | 29.4 |
| MC std / CRB | 1.01 (1.12 vs surrogate CRB) | 1.02 | – |
| Bias ± SE | +0.0097 ± 0.0037 | −0.015 ± 0.008 | −0.17 ± 0.06 |
| at bound / ambiguous / failed (of 2000) | 11 / 0 / 0 | 0 / 0 / 0 | 1163 / 11 / 6 |

### Surrogate accuracy and femwell noise

All from `results_montecarlo/data/surrogate_accuracy.json`.

| Number | Value |
|---|---|
| Held-out n_eff RMS, strips | 0.04–0.37 σ_phase |
| Held-out n_eff RMS, rib 400–600 nm (slab 90), rib 500 nm (slab 50) | 1.26–2.50 σ_phase |
| Held-out n_eff RMS, rib 300 nm, slab 90 / slab 50 | 8.1 / 30.1 σ_phase (femwell remeshing noise, up to ~1e-4) |

### Step check

| Number | Value | Source |
|---|---|---|
| All 372 derivatives | 366 PASS, 3 near-zero PASS, 3 marginal FAIL (2.04–2.25%, all n_g) | `results_montecarlo/sym_search/data/step_check_all.csv` |

### Sourced measurement noise (old placeholders vs sourced)

Source, verified against the arXiv full text (v1): S. Janz, S. Dedyulin, D.-X. Xu, M. Vachon, S. Wang,
R. Cheriton, J. Weber, "Measurement Accuracy in Silicon Photonic Ring Resonator Thermometers:
Identifying and Mitigating Intrinsic Impairments", arXiv:2306.15643 [physics.optics] (2023),
DOI 10.48550/arXiv.2306.15643. Journal version: UNVERIFIED.
- Optimistic, 0.1 pm: "the resonance wavelength uncertainty in Fig. 6 is dominated by the ±0.1 pm
  repeatability of the laser, and by the 0.1 pm wavelength step size used to acquire the spectra"
  (temperature-bath-stabilized ring, ~10 averaged spectra).
- Conservative, 1 pm: "Wavelength changes down to 1 pm or less can be routinely measured using
  calibrated tunable lasers" (Introduction).

The conversions are ours, not the source's:
- σ_neff = n_g·σ_λ/λ;
- σ_ng = √2·σ_λ·n_g²·L/λ², from one FSR of a ring of length L;
- σ_phase = 2π·n_g·L·σ_λ/λ², an MZI fringe read as precisely as a ring resonance (UNVERIFIED; a broad
  fringe is likely worse).

Inputs: λ = 1.55 µm, L = 100 µm, and n_g = 4.26, the largest among the headline structures. No
published values were found and verified for an FSR/group-index precision or an MZI phase-readout
precision as such, so those two entries are derived from the same resonance-wavelength precision.
All numbers below are from `results_montecarlo/noise_sourced/data/noise_sourced.json`.

| Noise | σ_neff | σ_ng | σ_phase (rad) |
|---|---|---|---|
| Placeholder (old) | 6.5e-6 | 1.0e-3 | 1.0e-3 |
| Sourced, 0.1 pm | 2.75e-7 | 1.07e-4 | 1.11e-4 |
| Sourced, 1 pm | 2.75e-6 | 1.07e-3 | 1.11e-3 |

| Strip-angle CRB (°) [h/2–2h] | Placeholder | Sourced 0.1 pm | Sourced 1 pm |
|---|---|---|---|
| Best one-etch set (8), D realistic | 0.162 [0.118–0.171] | **0.0139** [0.0101–0.0147] | **0.139** [0.101–0.147] |
| Best one-etch set (8), C no priors | 0.144 [0.105–0.153] | 0.0124 [0.0089–0.0131] | 0.124 [0.089–0.131] |
| Two-width rib 300 + 500, slab 50, D realistic | 0.367 [0.367–0.377] | **0.0408** [0.0408–0.0418] | **0.389** [0.389–0.400] |
| Strip + rib70 500 + 300, E realistic | 0.428 [0.428–0.444] | 0.0479 [0.0479–0.0496] | 0.455 [0.455–0.472] |

The placeholders sit close to the 1 pm level. At 0.1 pm the bounds shrink about tenfold, but they
then demand about 1 mK temperature stability (Janz et al. quote ~80 pm/K) and derivatives accurate far
beyond the remeshing noise of the rib 300 nm solves (see the next table). The Monte Carlo was run at
the placeholder noise only.

### Robustness to the rib 300 nm derivatives

From `results_montecarlo/robustness/data/robustness.json`, written by `src/robustness_rib300.py`.
Each rib 300 nm derivative entry is multiplied by (1 + u), u uniform in ±p, over 2,000 draws. You
asked for ±3% and ±5%; ±0.3% and ±1% are added because they match the size of the mesh errors.

| Strip-angle CRB (°), median [5–95%] | Best one-etch set (8), D realistic | Two-width 300 + 500, slab 50, D realistic |
|---|---|---|
| Nominal | 0.162 | 0.367 |
| ±0.3% | 0.160 [0.128–0.215] | 0.367 [0.343–0.396] |
| ±1% | 0.159 [0.0828–0.357] | 0.363 [0.293–0.484] |
| ±3% | 0.111 [0.0415–0.347] | 0.366 [0.215–1.34] |
| ±5% | 0.0706 [0.0267–0.322] | 0.353 [0.167–2.56] |
| **Fixed mesh** (no remeshing across the stencil) | not finished (18 of 25 solves; paused at the battery floor) | **0.372 [0.370–0.380] (+1.3%)** |

**Fixed mesh:** for the rib 300 nm on the 50 nm slab, the geometric derivatives without remeshing
agree with the remeshed ones within 0.13% (width, thickness, angle and slab; n_eff and n_g; h/2, h,
2h). The two-width number moves by only +1.3%. So the remeshing noise that spoils the surrogate fits
(up to ~1e-4 in n_eff between unrelated geometries) mostly cancels in the finite differences. Within
one stencil, the remeshed solves differ from the fixed-mesh ones by up to 2e-5
(`results_montecarlo/fixedmesh/cache/solves_r12.jsonl` vs `results_femwell/cache/solves_r12.jsonl`).

**Fragility:** the angle information comes from small differences between nearly parallel rib
columns.
- The two-width number tolerates the measured ~0.1% derivative error (±0.3% gives −7% to +8%).
- It falls apart at a few percent.
- The best one-etch set is more sensitive: random ±0.3% errors already give −21% to +33%.

The fixed-mesh and mesh-refinement checks (0.13%, 0.30%) are the relevant derivative accuracies. Real
mesh errors are systematic and correlated across structures, not independent random errors, which
should matter less than this test. Still, the 0.162° figure should carry a stated derivative-accuracy
caveat.

### Scope notes

- **Slab range:** 50 nm is the edge of the tested range (50–150 nm). The two-width optimum is at that
  edge (0.367° at 50 nm, rising steadily with slab thickness), so a thinner slab might be better still.
  Adding 30 and 40 nm slabs for the 300/500 nm ribs needs the full stencil (w, t, a, slab, n_Si, n_ox
  at h/2, h, 2h: 37 solves per structure, 148 in all), estimated at ~70 min with 2 workers at the
  measured 58 s per rib solve and worker (`results_montecarlo/cache/solves_r12.jsonl`), over the 1 h
  limit. It was not run.
- **Slots and SWG are dropped from the main results.** Oxide-clad slots and SWG waveguides are
  mirror-symmetric about mid-height like strips, so they are blind to the sidewall angle in the same
  way (R = 1 by symmetry; femwell slot R 0.9991–1.0015, SWG 0.9988–1.0005,
  `results_femwell/data/crosscheck_and_ratios.json`). Their earlier numbers are listed as RETIRED
  below and were not rechecked with the symmetry imposed.

## Conditional numbers (only under stated assumptions; not recomputed with the symmetry)

These hold only if the rib is etched in the same step as the strip, so the angle is shared, and the
depth is known (Stage A). That is not the usual process, so none of them is a headline number.

| Number | Value | Source |
|---|---|---|
| Strip + 5 ribs, conservative, Stage A | 0.000606° | `results_femwell/data/robustness_checks.json` |
| Two strip rings + one rib, Stage A | 0.00824° | `results_femwell/data/robustness_checks.json` |
| Strip + rib70 500 nm, Stage A | 0.00146° | `results_etchdepth_femwell/data/stages.json`, `results_etchdepth_femwell/logs/analysis.log` |

## RETIRED

| Number | Where it appeared | Reason |
|---|---|---|
| Strip-alone angle CRB 0.612° (0.34–0.612) | `results_femwell/SUMMARY.md`, `COMPARISON.md` | Strip symmetry: true value is ∞; the finite number is derivative noise in R |
| Strip-alone 1.33°, 1.39° (stages A–E) | `results_etchdepth_femwell/SUMMARY.md`, `results_priors_search/SUMMARY.md` | Same: ∞ under the symmetry |
| Best strips-only set 0.61–0.68° (0.614°) | `results_priors_search/SUMMARY.md`, `results_montecarlo/SUMMARY.md` step 2 | Same: ∞ for every size |
| Strip + TM ring 1.17°, 1.36° | `results_etchdepth_femwell/SUMMARY.md` | Strips only: ∞ |
| Two strip rings 118° (90–427°), 788° | `results_femwell/SUMMARY.md`, `results_etchdepth_femwell/SUMMARY.md` | Strips only: ∞ |
| Stage C/D strip + one rib = 1.33° / 1.39° ("same as strip alone") | `results_etchdepth_femwell/SUMMARY.md` | Now ∞ without priors (25.1° with realistic priors) |
| "900x" (1.33° → 0.0015°) and "~1000x", 1,009x (560–1,009x), 2,270x | `results_etchdepth_femwell/SUMMARY.md`, `results_femwell/SUMMARY.md`, `results_priors_search/SUMMARY.md` | Divide by the strip-alone baseline (∞), and assume Stage A |
| 14,352x (two rings + rib) | `results_femwell/SUMMARY.md` | Divides by the 118° strip baseline |
| 23x / 51x / 11–14x rib improvement factors, and 320x / 703x | `results_femwell/COMPARISON.md`, `results/CORRECTED_VS_ORIGINAL.md` | Ratios to the strip baseline |
| "about 4x", "3x", 3.3x, 4.11x, 307.5x, 1.5x, 1.16x | `results_etchdepth_femwell/SUMMARY.md`, `results_priors_search/SUMMARY.md` (T1, T2) | Ratios to strip alone; the denominator is ∞. Compare rib sets with each other instead |
| Slab gains 1.42e6x, 1,301,591x, 789x and the ordering 70 > 90 > 50 > 110 > 130 > 150 | `results_femwell/SUMMARY.md`, `COMPARISON.md` | Ratios with a singular denominator; replaced by the symmetry-imposed slab tables (two widths: 50 best; best 8: 90 best) |
| Slot set 0.308° (0.108–0.308), slot gains 33x / 4.6x / 1.7x, 2.6x, 6.8x, 2.0x | `results_femwell/SUMMARY.md`, `COMPARISON.md`, `results/CORRECTED_VS_ORIGINAL.md` | Oxide-clad slots (and SWG) are also mirror-symmetric about mid-height, so the same artifact applies; not rechecked with the symmetry imposed |
| Best one-etch sets 0.145–0.151° (C) and 0.167–0.176° (D realistic) | `results_priors_search/SUMMARY.md` | Superseded by the symmetry-imposed 0.144° and 0.162–0.174° |
| Step-2 table 0.400° / 0.414° (strip + rib70 500 + 300) and 0.145 / 0.167 | `results_montecarlo/SUMMARY.md` step 2 | Computed without the symmetry; use 0.419° (C) / 0.428° (E realistic) and 0.144 / 0.162 |
| B-stage strip + one rib 0.324° (rib70), 0.32–0.35° | `results_etchdepth_femwell/SUMMARY.md` | Superseded by 0.334° with the symmetry |
