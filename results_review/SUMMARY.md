# Review analyses (reviewer critique of v1, 2026-11-03/04)

Every number below is read from the data file named in its section. Stage D, realistic priors (etch depth 2 nm,
n_Si 1e-4, n_ox 1e-3, T 1 K), 1 pm resonance-wavelength precision, ring and MZI per structure unless stated.
"Two-width set" = strip TE rings 400-600 nm + ribs 300 and 500 nm on a 50 nm slab; "best set" = the best
ring-and-MZI set of eight on a 90 nm slab. Strip symmetry imposed throughout. Bounds are 1-sigma CRLB of the strip
angle, in degrees.

Scripts (main environment unless stated): `src/review_analyses.py`, `src/review_practical_context.py`,
`src/literature_verify.py`, `src/make_literature.py`, `src/review_surrogate_check.py` (design / analyze),
`src/review_surrogate_batch.py` (femwell environment), `src/review_feasibility.py`, `src/review_systematic.py`. All rows are also in `results_final/HEADLINE_NUMBERS.md`
sections L-Z.

## 1. Correlated noise, calibration offset, temperature

Data: `results_review/data/review_analyses.json` (key `noise`).

- Baseline: two-width 0.3891, best 0.139.
- Correlation rho between the n_eff and n_g errors of one ring. Derived value -1/sqrt(2) (n_eff from one of the two
  resonances that define the FSR): two-width 0.3365, best 0.1142. rho = 0.9: 0.2996 and 0.0795; rho = -0.9: 0.2985
  and 0.0802. The independent-error baseline is conservative.
- Common wavelength-calibration offset (one extra unknown): unknown, two-width 1.07 and best 0.1393; prior 10 pm,
  0.3892 and 0.139; prior 1 pm, 0.3891 and 0.139. The two-width set needs an absolute reference of about 10 pm.
- Temperature as a shared unknown: prior 0.1 K, two-width 0.3891; 0.5 K, 0.3898; 1 K, 0.3917; unknown, 0.4122. The
  best set stays at 0.139.
- Independent temperature drift per structure: 1 mK, two-width 0.3892 and best 0.1396; 10 mK, 0.3935 and 0.1877;
  100 mK, 0.4599 and 0.4993.

## 2. MZI readout and ring-only sets

Data: `results_review/data/review_analyses.json` (key `mzi`).

- Two-width set: rings only 0.3932; MZI at the derived noise 0.3891; MZI 3x noisier 0.3925; 10x noisier 0.3932.
- Best set: rings only 0.1901; MZI 0.139; 3x 0.1815; 10x 0.1892.
- Exhaustive ring-only search (at least one strip and one rib): best sets of 4-8 give 0.4356, 0.3606, 0.2029, 0.1941
  and 0.1895. The best of eight is TE strips 400, 550, 600 nm, TM strip 400 nm and ribs 300-600 nm on a 90 nm slab.
- Ring-only two-width set by slab: 30 nm 0.4202, 40 nm 0.4079, 50 nm 0.3932, 70 nm 0.4554, 90 nm 0.6001,
  110 nm 0.8576, 130 nm 1.298, 150 nm 4.094.
- Decision: the MZI adds about 1% for the two-width set, so the paper proposes ring-only structures.

## 3. Per-structure departures from the shared parameters

Data: `results_review/data/review_analyses.json` (key `variation`).

- Bias gain (RMS strip-angle bias per unit of independent per-structure offset): two-width 219.7 deg/deg and
  57.62 deg/nm; best 194.1 deg/deg and 74.41 deg/nm.
- Random angle offsets of std 0.25 deg: bias std 54.92 (two-width), 48.52 (best). Random width offsets of 0.5 nm:
  28.81 and 37.21.
- Linear angle trend 0.25 deg per 100 nm of width: bias 6.869 (two-width), -36.93 (best).
- Random-effects CRLB (per-structure deviations as unknowns with a Gaussian prior):
  - angle spread 0.01 deg: two-width 0.5502, best 0.5204; 0.1 deg: 2.295 and 1.024;
  - width spread 0.1 nm: 1.139 and 0.7154; 0.5 nm: 3.38 and 2.107;
  - thickness spread 0.1 nm: 2.544 and 2.047.
- Spread that doubles the CRLB: two-width 0.01749 deg (angle), 0.03938 nm (width), 0.01705 nm (thickness); best
  0.001739 deg, 0.00365 nm, and below the 0.001 nm grid for thickness.
- A common offset of all structures is absorbed by the shared parameters; only differences between structures bias
  the fit.

## 4. Practical value: resonance-prediction error

Data: `results_review/data/review_analyses.json` (key `practical`), `results_review/data/practical_context.json`.

- 1-sigma predicted-resonance error, all parameter uncertainties propagated, 1 pm: two-width set, strip TE 500 nm
  ring 0.4569 pm and rib 500 nm (50 nm slab) 0.7733 pm; best set, strip 0.5595 pm and rib 500 nm (90 nm slab)
  0.5214 pm.
- With rho = -0.707: two-width 0.3599 and 0.6412 pm. With an unknown calibration offset: two-width 982.2 pm, best
  39.68 pm. At 0.1 pm: two-width 0.04631 and 0.07733 pm.
- Published scale (verified abstracts): Zortman et al. 105 GHz within a die = 841.5 pm at 1550 nm, 1 THz across a
  wafer = 8.014 nm; Jayatilleka et al. trimming to +-32 pm; Selvaraja et al. linewidth control 0.6 nm within chip,
  2 nm chip-to-chip.

## 5. Literature

Data: `results_review/data/literature_records.json`; output `paper/LITERATURE.md`, `paper/refs.bib`.

- 50 verified entries (Crossref or arXiv record, OpenAlex abstract, patent page, or full text read). Not verified
  and not cited: the arXiv version of Klein et al. and an electrical sidewall-angle inference preprint (arXiv lookups
  rate-limited).

## 6. Femwell solves at perturbed geometries and Monte Carlo on femwell data

Data: `results_review/data/surrogate_check_design.json`, `results_review/cache/solves_r12.jsonl`,
`results_review/data/surrogate_check.json`.

- Design: 20 truths per set drawn from the CRLB covariance; 160 solves for the best set, 140 for the two-width set.
  All 300 solved (r12 mesh, 0 failed; battery-paused runs resumed from the cache on 2026-11-06/07). 500 fits per set:
  25 noise draws at 1 pm per truth, data = femwell values, model = surrogate, ring and MZI.
- Surrogate vs femwell, n_eff in MZI-phase noise sigmas: two-width set, rib 300 nm (50 nm slab) RMS 5.166, largest
  12.42; strips largest 0.8627; n_g largest 0.0524 n_g sigmas. Best set, rib 300 nm (90 nm slab) RMS 1.519, largest
  3.366.
- Two-width set: std 0.4034 deg (95% bootstrap 0.3781-0.4273), 1.037 x CRB 0.3891; mean 0.02463 +- 0.01804; 500 ok.
- Best set: std 0.2067 deg (0.1859-0.2253), 1.487 x CRB 0.139; 14 fits at the box edge, 0 failed.
- Reading: the two-width bound holds with an independent forward model; the best set does not reach its bound. Section
  11 shows that the 300 nm rib data alone do not explain that gap.

## 7. Realistic local variation (one copy of each structure)

Data: `results_review/data/feasibility.json` (key `variation`). Script: `src/review_feasibility.py`. Literature inputs: `results_final/data/citations.json`, `results_review/data/practical_context.json`.

Sources: width within-chip control better than 0.6 nm (Selvaraja et al., abstract); per-die std of width and thickness 1.26 and 0.30 nm and local maxima 3 and 0.5 nm between devices less than a few hundred microns apart (Xing et al., full text). These describe copies of one design. No verified measurement of a local sidewall-angle spread was found: **all angle levels are ASSUMED**; width levels below 0.6 nm and thickness 0.1 nm are ASSUMED.

Two-width set, ring and MZI, shared-parameter bound 0.3891 deg. CRLB (deg):

| width (nm) | t 0.1: angle 0.1 | 0.3 | 0.5 | 1 | t 0.3: angle 0.1 | 0.3 | 0.5 | 1 |
|---|---|---|---|---|---|---|---|---|
| 0.1 | 3.067 | 4.333 | 4.958 | 5.421 | 5.16 | 5.624 | 6.079 | 6.636 |
| 0.3 | 3.419 | 4.419 | 4.98 | 5.423 | 5.249 | 5.673 | 6.1 | 6.639 |
| 0.6 | 4.081 | 4.634 | 5.042 | 5.429 | 5.491 | 5.812 | 6.163 | 6.648 |
| 1 | 4.687 | 4.915 | 5.146 | 5.441 | 5.851 | 6.042 | 6.279 | 6.668 |

Best set of 8, ring and MZI, shared-parameter bound 0.139 deg. CRLB (deg):

| width (nm) | t 0.1: angle 0.1 | 0.3 | 0.5 | 1 | t 0.3: angle 0.1 | 0.3 | 0.5 | 1 |
|---|---|---|---|---|---|---|---|---|
| 0.1 | 2.23 | 2.832 | 3.407 | 4.318 | 4.282 | 4.524 | 4.843 | 5.593 |
| 0.3 | 2.468 | 2.953 | 3.467 | 4.331 | 4.371 | 4.586 | 4.882 | 5.607 |
| 0.6 | 2.967 | 3.261 | 3.637 | 4.37 | 4.598 | 4.755 | 4.997 | 5.653 |
| 1 | 3.544 | 3.69 | 3.914 | 4.448 | 4.932 | 5.034 | 5.209 | 5.748 |

## 8. Replication: N copies of each structure, fitted jointly

Data: `results_review/data/feasibility.json` (key `replication`).

Shared process parameters, Gaussian offsets independent between copies. N for 1.5x and 2x refers to the shared-parameter bound of one copy.

Two-width set, ring and MZI (rings only in parentheses):

| width nm, angle deg, thickness nm | N=1 | 4 | 16 | 64 | 256 | N for 1.5x | N for 2x |
|---|---|---|---|---|---|---|---|
| 0.1, 0.1, 0.1 | 3.067 | 1.556 | 0.7911 | 0.4086 | 0.2177 | 31 (30) | 17 (17) |
| 0.3, 0.3, 0.1 | 4.419 | 2.215 | 1.113 | 0.5659 | 0.2978 | 61 (59) | 34 (33) |
| 0.6, 0.5, 0.3 | 6.163 | 3.119 | 1.581 | 0.8001 | 0.409 | 123 (120) | 68 (67) |
| 1, 1, 0.3 | 6.668 | 3.404 | 1.745 | 0.8864 | 0.4499 | 151 (147) | 84 (82) |

Best set of 8, ring and MZI (rings only in parentheses):

| width nm, angle deg, thickness nm | N=1 | 4 | 16 | 64 | 256 | N for 1.5x | N for 2x |
|---|---|---|---|---|---|---|---|
| 0.1, 0.1, 0.1 | 2.23 | 1.116 | 0.559 | 0.2799 | 0.1401 | 116 (62) | 65 (35) |
| 0.3, 0.3, 0.1 | 2.953 | 1.485 | 0.752 | 0.3816 | 0.1924 | 218 (116) | 122 (65) |
| 0.6, 0.5, 0.3 | 4.997 | 2.528 | 1.285 | 0.6652 | 0.346 | 730 (384) | 404 (211) |
| 1, 1, 0.3 | 5.748 | 2.918 | 1.482 | 0.7663 | 0.4004 | 995 (522) | 550 (286) |

A part of each offset common to all copies of a design does not average. Two-width set, ring and MZI, offsets 0.3 nm, 0.3 deg, 0.1 nm: per-design fraction, bound at N = 256, limit as N grows:

- 0.01: 0.299, 0.02536.
- 0.03: 0.3079, 0.06685.
- 0.1: 0.3662, 0.1758.
- 0.3: 0.5963, 0.4926.

## 9. Layout (linear gradient, ring pitch 50 um ASSUMED)

Data: `results_review/data/feasibility.json` (key `layout`).

Two-width set, offsets 0.3 nm, 0.3 deg, 0.1 nm, N = 16: angle bias per nm/mm of thickness gradient and bound with the gradients fitted:

- blocked: 2.337 deg; 1.205 deg.
- interleaved: 0.5842 deg; 1.114 deg.
- mirrored: 0 deg; 1.113 deg.

## 10. Common wavelength offset

Data: `results_review/data/feasibility.json` (key `calibration`).

- Two-width set, ring and MZI: no offset 0.3891; unknown offset 1.0697; differential formulation 1.0697. Rings only: 0.39325; 1.0846; 1.0846.
- The differential formulation equals the unknown-offset bound: differencing removes the same information.
- Absolute reference, ring and MZI: 0.1 pm 0.3891; 1 pm 0.38911; 5 pm 0.38914; 10 pm 0.38924; 20 pm 0.38963; 50 pm 0.39239; 100 pm 0.40198; 200 pm 0.43695; 500 pm 0.59644.
- Reference std at which the bound is 1.1x, 1.5x, 2x: 178.6 pm, 476.7 pm, 902 pm.
- Independent reference error per structure: 0.1 pm 0.3892; 1 pm 0.3963; 5 pm 0.4411.

## 11. Which structure causes the femwell gap of the best set

Data: `results_review/data/feasibility_ablation.json`.

500 fits per variant, same noise draws as section 6; intervals are 95% bootstrap.

| data | best set std (deg) | std/CRB | two-width std (deg) | std/CRB |
|---|---|---|---|---|
| femwell, all structures | 0.2067 [0.1859-0.2253] | 1.487 | 0.4034 | 1.037 |
| surrogate for the 300 nm rib | 0.1988 [0.1785-0.218] | 1.43 | 0.4023 | 1.034 |
| femwell for the 300 nm rib only | 0.1664 [0.1453-0.1857] | 1.197 | 0.4034 | 1.037 |
| surrogate, all structures | 0.1689 [0.1469-0.1886] | 1.215 | 0.4022 | 1.034 |

Reading: the 300 nm rib data do not explain the gap. Part of it is present with surrogate data everywhere (perturbed truths), the rest needs femwell data for the other structures. The test does not cover the effect of the 300 nm rib derivatives on the bound itself.

## 12. Systematic width offsets (third review)

Data: `results_review/data/systematic.json` (key `systematic`). Script: `src/review_systematic.py`.

Two-width set, ring and MZI. Width offset as a function of drawn width, x = (W - mean W)/100 nm, coefficients fitted with no prior or the prior std shown (nm). CRLB (deg):

| model | none | 0.01 | 0.03 | 0.1 | 0.3 | 1 |
|---|---|---|---|---|---|---|
| lin_all | 0.4518 | 0.4029 | 0.434 | 0.4496 | 0.4516 | 0.4518 |
| quad_all | 0.6747 | 0.4741 | 0.4971 | 0.5602 | 0.6463 | 0.6717 |
| lin_strip | 0.3899 | 0.3893 | 0.3896 | 0.3898 | 0.3899 | 0.3899 |
| quad_strip | 1.5 | 0.4436 | 0.5509 | 1.001 | 1.397 | 1.49 |
| rib_offset | ∞ | 0.3891 | 0.3892 | 0.39 | 0.3969 | 0.4682 |
| rib_offset_lin | ∞ | 0.4101 | 0.5494 | 1.315 | 3.105 | 4.724 |
| separate_lin | ∞ | 0.4103 | 0.5497 | 1.315 | 3.109 | 4.741 |
| per_design | ∞ | 0.4445 | 0.5481 | 1.139 | 2.512 | 4.474 |

Reading: smooth offsets shared by strips and ribs, or on the strips alone, are absorbed; a rib-only offset (relative to the strips, or between the two rib widths) or a free offset per design makes the angle unidentifiable without a prior.

Fitted systematic term plus random offsets 0.3 nm / 0.3 deg / 0.1 nm, copies needed for 2x (1.5x) of the shared bound:

- none: N = 1 4.419 deg; N = 34 0.7686 deg; 34 (61) copies.
- lin_all: N = 1 5.154 deg; N = 34 0.9163 deg; 48 (85) copies.
- quad_all: N = 1 5.544 deg; N = 34 0.9934 deg; 56 (100) copies.
- quad_strip: N = 1 4.426 deg; N = 34 0.7694 deg; 34 (61) copies.

## 13. Replication with spatially correlated offsets

Data: `results_review/data/systematic.json` (key `correlation`).

rho(d) = exp(-d/L), L ASSUMED (no published correlation length found); cells of pitch s, one copy of each design per cell. Copies needed for 2x (1.5x); bound at N = 256 (deg):

| model | L mm | s mm | copies | N = 256 |
|---|---|---|---|---|
| field | 0.05 | 0.15 | 29 (51) | 0.2801 |
| field | 0.05 | 0.5 | 33 (60) | 0.2973 |
| field | 0.05 | 1.5 | 34 (61) | 0.2978 |
| field | 0.2 | 0.15 | 17 (31) | 0.2287 |
| field | 0.2 | 0.5 | 27 (49) | 0.2745 |
| field | 0.2 | 1.5 | 33 (59) | 0.2961 |
| field | 1 | 0.15 | 8 (15) | 0.2027 |
| field | 1 | 0.5 | 15 (26) | 0.2139 |
| field | 1 | 1.5 | 23 (41) | 0.2563 |
| field | 5 | 0.15 | 3 (6) | 0.2579 |
| field | 5 | 0.5 | 7 (12) | 0.2135 |
| field | 5 | 1.5 | 12 (21) | 0.2013 |
| same_design | 0.05 | 0.15 | 36 (66) | 0.3089 |
| same_design | 0.05 | 0.5 | 34 (61) | 0.2978 |
| same_design | 0.05 | 1.5 | 34 (61) | 0.2978 |
| same_design | 0.2 | 0.15 | 71 (136) | 0.4387 |
| same_design | 0.2 | 0.5 | 38 (69) | 0.3165 |
| same_design | 0.2 | 1.5 | 34 (61) | 0.2979 |
| same_design | 1 | 0.15 | >256 (>256) | 0.914 |
| same_design | 1 | 0.5 | 99 (197) | 0.5225 |
| same_design | 1 | 1.5 | 47 (87) | 0.3515 |
| same_design | 5 | 0.15 | >256 (>256) | 1.392 |
| same_design | 5 | 0.5 | >256 (>256) | 1.061 |
| same_design | 5 | 1.5 | 169 (>256) | 0.6681 |

Check: L = 0, N = 16 gives 1.113 deg, the independent-offset value.

## 14. Resonance prediction under local variation

Data: `results_review/data/systematic.json` (key `prediction`).

- strip_TE_500: ideal case 0.4569 pm; own local offsets 734.5 pm; process mean at N = 1, 34, 256: 329.3, 56.5, 20.61 pm; new device at N = 34: 736.7 pm.
- rib_TE_500_slab50: ideal case 0.7733 pm; own local offsets 440.3 pm; process mean at N = 1, 34, 256: 356.6, 63.15, 23.17 pm; new device at N = 34: 444.8 pm.
