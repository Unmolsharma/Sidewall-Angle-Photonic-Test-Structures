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
| Strip R = (dn/da)/(dn/dw)/(−t·π/180), exact by symmetry | 0.99990 (= (π/180)/tan 1°) | `results_montecarlo/SUMMARY.md` (step 2) |
| Strip R, femwell residual about exact | TE −2.7e-4 to +1.9e-3; TM up to 8.5e-3 | `results_montecarlo/data/strips_info.json` |
| Rib R | 0.776–0.799 | `results_femwell/SUMMARY.md` |
| Strips-only sets, any size | ∞ (no angle information) | `results_montecarlo/sym_search/data/search.json` |
| 1° tilt looks like a strip width change of | 3.84 nm | `results_femwell/SUMMARY.md`; `results_montecarlo/data/mc_fits.csv` (slope 3.84 nm/deg) |
| Ring resonance error from a 1° strip tilt treated as width | −711 pm | `results_femwell/SUMMARY.md` |
| Mesh check r12 vs r8 | agree to 0.3% | `results_femwell/SUMMARY.md` |

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

## Conditional numbers (only under stated assumptions; not recomputed with the symmetry)

These hold only if the rib is etched in the same step as the strip, so the angle is shared, and the
depth is known (Stage A). That is not the usual process, so none of them is a headline number.

| Number | Value | Source |
|---|---|---|
| Strip + 5 ribs, conservative, Stage A | 0.000606° | `results_femwell/SUMMARY.md` |
| Two strip rings + one rib, Stage A | 0.00824° | `results_femwell/SUMMARY.md` |
| Strip + rib70 500 nm, Stage A | 0.00146° | `results_etchdepth_femwell/SUMMARY.md` |

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
