# Review analyses (reviewer critique of v1, 2026-11-03/04)

Every number below is read from the data file named in its section. Stage D, realistic priors (etch depth 2 nm,
n_Si 1e-4, n_ox 1e-3, T 1 K), 1 pm resonance-wavelength precision, ring and MZI per structure unless stated.
"Two-width set" = strip TE rings 400-600 nm + ribs 300 and 500 nm on a 50 nm slab; "best set" = the best
ring-and-MZI set of eight on a 90 nm slab. Strip symmetry imposed throughout. Bounds are 1-sigma CRLB of the strip
angle, in degrees.

Scripts (main environment unless stated): `src/review_analyses.py`, `src/review_practical_context.py`,
`src/literature_verify.py`, `src/make_literature.py`, `src/review_surrogate_check.py` (design / analyze),
`src/review_surrogate_batch.py` (femwell environment). All rows are also in `results_final/HEADLINE_NUMBERS.md`
sections L-P.

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
- Status: see the paper section "Monte Carlo on femwell data" and `paper/REVIEW.md`.
