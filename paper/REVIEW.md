# Review of the draft

v1 (reviewed 2026-11-03) is archived in `paper/archive/v1/`. This file covers v2 and the changes after the external
review of the 2026-11-08 draft. Numbers are checked by `python src/provenance_audit.py --paper` (report:
`results_final/audit/PAPER_AUDIT.md`); the review analyses are summarized in `results_review/SUMMARY.md` and are rows
L-V of `results_final/HEADLINE_NUMBERS.md`.

## Audit result (2026-11-09)

- Paper: 509 numbers with a decimal point or exponent, **0 untraced**. The pool now holds so many values that its
  chance-match rate is 55%, so a TRACED verdict alone is weak evidence. The 67 numbers newly written into the prose were
  therefore also checked one by one against their own JSON paths: **0 mismatches**.
- HEADLINE_NUMBERS.md: 3123 numbers, 0 untraced; CLAIMS.md: 39, 0 untraced; results_review/SUMMARY.md: 342, 0 untraced.
- Compilation: `latexmk -pdf main.tex` runs without errors, undefined references or citations (13 pages); eight overfull
  boxes of 0.94-1.45 pt remain (empty boxes, as before).
- References: 46 cited, from the 51 verified entries of refs.bib (new: NIST SRM 2519a, read from the certificate).

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

## Still open

- No measurement; no measured per-design offset between adjacent waveguides of different widths; no measured local angle spread.
- No air-clad bound; TM ribs; wall-shape sensitivity; unknown length per ring; a second full-vector solver for derivatives; the bound's domain dependence (5.2% for the two-width set).
- The cause of the best set's femwell gap is only partly located. A smoother forward model and more truths would be needed to resolve it.
- Layout model: linear gradients only; no curvature, short-range correlation or pattern density.
