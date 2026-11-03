# Review of the draft (step 4)

Reread as a skeptical reviewer on 2026-11-03. Numbers checked by `python src/provenance_audit.py --paper`
(report: `results_final/audit/PAPER_AUDIT.md`).

## Audit result

- 164 numbers with a decimal point or exponent in `main.tex` and `sections/*.tex` (comments, `\todo` notes, labels,
  references and file names excluded): **0 untraced, 0 mismatches**, traced against `results_final/HEADLINE_NUMBERS.md`,
  `CLAIMS.md`, `ASSUMPTIONS.md`, `CITATIONS.md` and `results_final/data/validation.json`.
- Caveat: the chance-match rate of this pool is 28% (a random 3-digit number would "trace" that often), so the
  audit proves that every number exists in a source file at the shown precision, not that it was taken from the
  intended row. Each number was placed by hand from the row noted in `paper/NOTES.md`.
- Integers (counts such as 1190 of 2000, 7 of 32, widths, slabs) are not audited by the script; they were checked
  by hand against HEADLINE_NUMBERS.md and CLAIMS.md.
- Figures read the same data files as HEADLINE_NUMBERS.md. Fig. 1 asserts that its recomputed centroids equal
  `crosscheck_fields.json` to 1e-9 um; Fig. 4 prints the Monte Carlo and CRLB values from `mc_results_1pm.json`.
- Compilation: `latexmk -pdf main.tex` runs without errors, undefined references or undefined citations. Five overfull
  boxes of 0.94-1.45 pt remain; they come from empty boxes (not text) and stay with `\finaltrue`.

## Claims vs evidence

- No claim is stronger than its CLAIMS.md status. "Air cladding barely matters" (UNSUPPORTED) is not claimed; air
  cladding is a stated limitation. "Thinner slab is better" (PARTIAL) is written as "the optimum among the slabs
  simulated is 50 nm".
- Wording follows the decisions of 2026-11-03: the EMpy lesson is split (EMpy derivatives disagree and fail step
  checks; derivative noise in any solver fakes strip information), and the single-rib priors sentence uses four numbers
  at one stage and noise level (stage C, 1 pm, 70 nm slab), without "rescue".
- Simulation-only language throughout ("simulations indicate", no "we demonstrate").

## Weakest points (what a reviewer will push on)

1. **No measurement.** Everything is simulation; the noise model is cited only for the resonance wavelength
   (Janz et al.), and the n_g and MZI phase noise are our derivation. The bounds scale linearly with
   sigma_lambda and move by tens of percent with the rib noise factor, the ring length and the n_g used in the
   conversion (ASSUMPTIONS.md).
2. **Numerical uncertainty of the headline numbers.** The two-width bound moves +5.2% when the domain is enlarged
   (derivatives change by <= 0.35%); the best set moves -15.3% with fixed-mesh rib derivatives and -27.5% at half
   the step. A reviewer may ask for a second, independent solver for the derivatives, or a mesh/domain refinement
   study of the bound itself rather than of individual derivatives.
3. **Observation model.** We treat each ring as an absolute n_eff measurement (known resonance order and ring
   length) and every structure in a set as sharing the same width offset, thickness and angles. Real dies have
   local variation between structures and unknown ring lengths; both would add unknowns. Now stated as
   limitations, but not quantified.
4. **Monte Carlo surrogate.** The Monte Carlo fits a smooth surrogate, not femwell itself. For the 300 nm ribs the
   surrogate misses held-out femwell solves by up to 27 sigma_phase (femwell remeshing scatter). The estimator is
   therefore tested on a smoothed model; the best-set spread is 0.952 of the femwell CRLB because the surrogate's
   own bound is 0.1249 deg. A check with a handful of direct femwell fits would close this.
5. **Untested modeling choices:** ideal trapezoid walls, pivot at the bottom of the wall, nominal indices without
   dispersion, the 5 nm wavelength step for n_g, TE ribs only, 1.55 um only (NOT EVALUATED in ASSUMPTIONS.md).
6. **Validation.** Only against the earlier EMpy solutions (n_eff within 0.79%/0.15%) and one published EIM value
   (2.602; femwell 1.2% lower). The "2.44-2.45 published" strip value was dropped because no source could be found.
7. **Literature coverage is thin.** Three `\todo{cite}` remain in the introduction; the sidewall-angle variability
   and AFM sidewall-metrology statements need verified references, and the claim that optical extraction fixes the
   angle needs the full text of Xing et al. (2018) (our fetch failed on a TLS error).

## Open TODOs in the source

- `introduction.tex`: cite measured sidewall-angle variation in Si photonic etch processes; cite 3D/CD-AFM sidewall
  metrology on Si waveguides; verify that Xing et al. 2018 fix the angle (search snippet: 85 deg from a foundry
  handbook).
- `main.tex`: e-mail, affiliation details / corresponding author / dates, co-author last name, acknowledgments.

## What would make it stronger before arXiv

1. Close the three citation TODOs with verified sources.
2. A direct femwell Monte Carlo check (e.g. 50-100 sets fitted with femwell in the loop) for the two-width set.
3. A domain + mesh refinement of the two-width bound itself (not only derivatives), to state a numerical error bar.
4. An observation model with unknown ring length (or relative resonance shifts) and a local-variation term, to show
   the two-width conclusion survives.
5. Air-clad sets: compute the bound instead of only R (would decide whether the oxide-clad restriction matters).
6. Wall-shape sensitivity (corner rounding, pivot at mid-height or top) for the two-width set.
7. An independent derivative check with a second full-vector solver on two or three structures.

## A 2-page CLEO 2027 summary (due 23 Nov 2026)

Keep: the symmetry result with Fig. 1 (cut to panels a and c, or a two-panel single-column version), the counting
rule in one paragraph, the single-rib sentence (stage C numbers), the two-width set (0.39 deg at 1 pm, 0.0408 deg
at 0.1 pm) with the Monte Carlo agreement (Fig. 4b), and the 711 pm cost of ignoring the angle. Mention the best set
(0.12-0.14 deg) in one sentence with its caveat.
Cut: the EMpy comparison, the full stage-by-stage discussion (keep one sentence and maybe Fig. 3 as the second
figure instead of Fig. 4), the slab study (one sentence: optimum at 50 nm), the priors sensitivity table, the fixed-mesh
and domain checks (one sentence), the surrogate details. Two figures maximum; about 8-10 references.
