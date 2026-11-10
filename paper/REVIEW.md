# Review of the draft

v1 (reviewed 2026-11-03) is archived in `paper/archive/v1/`. This file covers v2 and the changes after the external
reviews of the 2026-11-08 and 2026-11-10 drafts. Numbers are checked by `python src/provenance_audit.py --paper`
(report: `results_final/audit/PAPER_AUDIT.md`); the review analyses are summarized in `results_review/SUMMARY.md` and
are rows L-Z of `results_final/HEADLINE_NUMBERS.md`.

## Audit result (2026-11-10)

- Paper: 594 numbers with a decimal point or exponent, **0 untraced**. The pool's chance-match rate is 62%, so a TRACED
  verdict alone is weak evidence; the numbers written in the last two rounds (67 + 30) were also checked one by one
  against their own JSON paths: **0 mismatches**.
- HEADLINE_NUMBERS.md: 3964 numbers, 0 untraced; CLAIMS.md: 39, 0 untraced; results_review/SUMMARY.md: 474, 0 untraced.
- Compilation: `latexmk -pdf main.tex` runs without errors, undefined references or citations (14 pages); eight overfull
  boxes of 0.94-1.45 pt remain (empty boxes, as before). Placeholders are yellow `	odo` boxes (set `inaltrue` in
  main.tex to hide them).
- References: 46 cited, from the 51 verified entries of refs.bib.

## Points of the 2026-11-10 review (third round)

| # | Request | Status | Where |
|---|---|---|---|
| 1 | Systematic vs random variation | done | Sec. "Random and systematic offsets", Table V. Smooth width-dependent offsets shared by strips and ribs are absorbed (linear 0.452 deg, quadratic 0.675 deg, no prior); a width offset of the ribs relative to the strips, or between the two rib widths, makes the angle unidentifiable without a prior (1 nm prior on the rib-strip offset: 0.468 deg; 0.03 nm on the rib-rib difference: 0.549 deg). Why replication averages one and not the other is stated with the copy-mean equation. |
| 2 | Spatial correlation | done | Table IV: exp(-d/L), L = 0.05-5 mm ASSUMED (Chrostowski 2014 and Lu 2017 describe correlated variation but give no length), copy spacing 0.15-1.5 mm, two limiting models. Shared field: 3-34 copies for 2x; same-design correlation: 34 to more than 256. |
| 3 | Prediction under realistic variation | done | Sec. "Resonance prediction": 0.457 / 0.773 pm labelled ideal case; process mean 329 pm (N = 1), 56.5 pm (N = 34); a new device's own offsets add 735 pm (strip) and 440 pm (rib) |
| 4 | Pair 0.39 deg with the realistic result | done | abstract and conclusion: ideal 0.39 deg, 4.42 deg for one copy, 34 / 61 copies, random-vs-systematic caveat; conclusion keeps the best set secondary with 0.2067 deg |
| 5 | Placeholders as 	odo | done | author names, dates, affiliations and e-mails, corresponding author, acknowledgments |
| r4 | Independent expert check of the perturbation argument and Fisher implementation | open | not something this repo can do; see open issues |

### Decisions left to you

- The abstract is now 231 words: within the IEEE limit of 250 but above the earlier 200-word target, because of the pairing
  and caveat you asked for. Cut the first sentence's second clause or the etch-bias sentence if you want it shorter.
- The rib-to-strip width offset is the new weak point. The model has always assumed that strips and ribs share one width
  offset; the analysis shows the angle depends on that assumption. Decide how prominently to state it (now: Results K,
  Discussion B and C, abstract).
- Which correlation model to lead with. The paper presents both as bounds; a test chip would decide.

## Points of the 2026-11-08 review

| # | Point | Status | Where |
|---|---|---|---|
| A | Eight-structure result needs caution | addressed, with a changed attribution (see below) | Table I, Sec. "Multi-width rib sets", "Monte Carlo on femwell data". CRLB 0.139 deg shown beside 0.2067 deg on femwell data; two-width set is the headline in abstract, introduction and conclusion |
| B | Interpret variability; connect to layout | addressed (model), open (data) | Sec. "Realistic local variation and replication", "Placement on the chip", Discussion "What it would take" |
| C | Common calibration offset | addressed | Sec. "Measurement-model variants": differential formulation gives no gain (1.070 deg, same as an unknown offset); reference of 0.1-5 pm returns the baseline; NIST gas cell cited |
| D | Measurement model needs grounding | partly addressed | No measured noise for the complete structures exists here; stated as a limitation |
| 1 | Make the two-width design primary | done | abstract, contribution list, conclusion |
| 2 | Stronger closest-work comparison | open | no new work read; see "Sections you should rewrite yourself" |

### Changes from the requested text

1. **The femwell gap is not attributed to the 300 nm rib.** The request was to attribute the gap between 0.139 and 0.207 deg to
   the 300 nm rib derivative accuracy. A new test (500 fits per variant, same noise draws; rows U) does not support it: surrogate data
   for that rib only gives 0.1988 deg, femwell data for that rib only gives 0.1664 deg, and surrogate data everywhere gives 0.1689 deg
   (1.215 x CRLB). Part of the gap is present with no model error at these perturbed truths, and the rest needs femwell data
   for the other structures. The paper states this. Rib derivative accuracy remains the stated uncertainty of the bound itself
   (-15.3% with fixed-mesh derivatives), which this test does not address.
2. **"Needs a reference of about 10 pm" was wrong as written and is corrected.** The angle bound rises by 10% only at about 179 pm,
   and a predicted true wavelength inherits the reference error (10.01 pm for a 10 pm reference).
3. **"The two-width set is better because it tolerates ten times more variation" is removed.** At realistic spreads the best set has
   the lower bound (2.95 against 4.42 deg at 0.3 nm, 0.3 deg, 0.1 nm). The two-width set stays the proposal because its bound is the one reached
   with an independent forward model.

### What the replication rule does not cover

The rule assumes offsets independent between copies. A part common to all copies of a design does not average: at a tenth of the local
spread the two-width bound tends to 0.176 deg, and at three tenths of the 0.6 nm / 0.5 deg / 0.3 nm spreads no number of copies
reaches twice the shared bound. Neighbouring copies that share an offset count as fewer copies, and per-structure etch-depth offsets were not modeled.

### Source status of the variation levels

- Width: 0.6 nm (Selvaraja et al., abstract) and 1.26 nm per-die std (Xing et al., full text) bracket the upper levels; 0.1 and 0.3 nm are ASSUMED.
- Thickness: 0.30 nm per-die std (Xing et al.) matches the 0.3 nm level; 0.1 nm is ASSUMED.
- Sidewall angle: **all levels ASSUMED**. Two searches found no measured local spread of the angle.
- These sources describe copies of one design, and Xing et al.'s local values are maxima. Xing et al. give no die size, so gradients are quoted per unit (nm/mm).
- Xing et al. 2022 (ACS Photonics), Lu et al. 2017 and the ResearchGate pages could not be read (HTTP 403), so no numbers were taken from them.

## Target venue: what would need to change

The manuscript is 13 pages in IEEEtran two-column format.

**Journal of Lightwave Technology** (IEEE Photonics Society page and IEEE Author Center, read 2026-11-09)
- Format already matches: "standard two-column, single-spaced". The page asks for "a MS Word document and a PDF"; confirm whether LaTeX source is accepted.
- Length: "8 pages in length ... before incurring mandatory over length page charges". At 13 pages the paper is 5 over. A figure on the page reported $260 per extra page, but that came from an automated summary, so confirm it. Cutting to 8 would mean moving the Monte Carlo detail, the slab study and Fig. 6 or Fig. 7 to a supplement.
- Abstract: a single paragraph of up to 250 words (now about 198), self-contained, "without abbreviations, footnotes, or references"; complies.
- Keywords: 3-5 (now 8), cut three. Author biographies are optional.

**Optics Express**
- Read: "accepts manuscripts submitted in MS Word or LaTeX"; a universal template exists, and its visual styling is optional.
- **Not verified**: the Optica style guide and publication-fee pages blocked automated access, so abstract length, page charges and reference format must be checked by hand. One older template suggests an abstract of about 100 words; if that holds, this abstract needs halving.
- Move to the Optica template (single column), so regenerate the figures for its text width in `paper/scripts/style.py`; `table*` and `figure*` become ordinary floats; replace IEEEtran-only elements (`\IEEEkeywords`, `\thanks`) and the bibliography style.
- Length is less constrained in a single-column format, so the full analysis could stay in the main text; confirm when the fee page is checked.

**Both**: fill the co-author surname, e-mail, affiliations and acknowledgments, and check data- and code-availability requirements.

## Earlier reviewer points (2026-11-03 critique)

| # | Point | Status | Where |
|---|---|---|---|
| 1 | Monte Carlo not independent of the surrogate | addressed | Sec. "Monte Carlo on femwell data", Table I column "fw.": two-width set 0.4034 deg = 1.037 x CRB; best set 0.2067 deg = 1.487 x CRB |
| 2 | Correlated noise | addressed | Fig. 7: correlation lowers the bound |
| 2 | MZI precision separate from ring | addressed | Sec. "Ring and MZI readout"; proposal is ring-only |
| 2 | Temperature | addressed | shared unknown and independent drift per structure |
| 2 | Lab vs production noise | partly addressed | Xing et al. fitting errors compared in Theory; no production-tester noise source found |
| 3 | Assumed uniformity; width-dependent etch | addressed (analysis), open (data) | Sec. "Departures from shared process parameters" and replication section |
| 4 | Identifiability, not practical metrology | addressed | title, abstract, four-finding framing |
| 5 | Practical value of 0.39 deg | addressed | Sec. "Resonance prediction", Discussion A |
| s | Literature positioning; search; material assumptions; ring geometry | addressed or partly | Related work; Sec. "Test Structures and Design Search"; dispersion, doping and unknown ring length still not modeled |
| s | Generalizability | open | air cladding, TM ribs, other platforms |
| w | Wording and placeholders | addressed | except e-mail, affiliation details, acknowledgments and co-author surname |

## Sections you should rewrite yourself

1. **Introduction, paragraph 1** (the concrete problem). Make it sound like your lab's motivation.
2. **Related work.** Xing 2018, the NRC patent and the SRM 2519a certificate are read in full; the rest is from verified abstracts. Read at least Lu 2017, Butt 2024, Zhang 2024 and Dwivedi 2015, and strengthen the comparison with the closest extraction work.
3. **Discussion: "What it would take".** New and written to the model's limits. Your judgement of what a foundry can place, and of which per-design offsets to expect, is better than mine.
4. **Monte Carlo on femwell data.** Decide whether to keep the best set of eight in Table I and Fig. 8, given the unresolved gap.
5. **Next steps.** Replace with the experiment you can actually do.

## Final open-issues list

**Physics and model**
1. The angle's identifiability rests on strips and ribs sharing one width offset. A rib-to-strip offset must be known to about 1 nm, and a difference between the two rib widths to a few hundredths of a nanometre. No source gives either.
2. No measured local spread of the sidewall angle and no correlation length were found, so both are ASSUMED. The difference between adjacent waveguides of different widths has not been reported.
3. Per-structure etch-depth offsets, non-Gaussian offsets, curvature and pattern-density effects in the layout model, and the unknown length of each ring are not modeled.
4. The best set of eight: the gap between 0.139 and 0.2067 deg is only partly located.
5. Air cladding, TM ribs, wall shape, dispersion, doping and other platforms remain open.

**Numerics**
6. There is no second full-vector solver for the derivatives, and the bound's domain dependence is 5.2% for the two-width set.
7. An independent expert should check the perturbation argument and the Fisher implementation, especially the observables (n_eff from the resonance with known order and length; n_g from one FSR) and the noise covariance (a 1 pm resonance precision converted with n_g = 4.26 and L = 100 um).

**Experiment**
8. No device was fabricated or measured. The decisive test is a chip with the two-width set replicated and interleaved, adjacent and separated copies, cross-sections of the same dies, and a gas-cell reference.

**Manuscript**
9. Fill the `	odo` placeholders: co-author surname, dates, affiliations, e-mails, corresponding author and acknowledgments.
10. Venue formatting: JLT allows 8 pages before overlength charges, and the paper is now 14. The Optics Express limits are still unverified (see above).
11. The sections listed under "Sections you should rewrite yourself" are yours to rewrite. "What it would take" and "Random and systematic offsets" are new.
12. Commits are on branch `practical-feasibility`; `main` does not have the last two review rounds until you merge.
