# Review of the draft

v1 (reviewed 2026-11-03) is archived in `paper/archive/v1/`. This file covers v2, the rewrite after the external
critique of 2026-11-03. Numbers are checked by `python src/provenance_audit.py --paper` (report:
`results_final/audit/PAPER_AUDIT.md`); the review analyses are summarized in `results_review/SUMMARY.md` and are rows
L-P of `results_final/HEADLINE_NUMBERS.md`.

## Audit result (v2)

- Paper: 262 numbers with a decimal point or exponent, **0 untraced, 0 mismatches**. Chance-match rate 34%, so
  every new number was also placed by hand from its row (sections L-P of HEADLINE_NUMBERS.md; Xing et al. values
  from results_final/CITATIONS.md). The audit now also reads numbers printed as "1.1 × 10−5" in quoted papers.
- HEADLINE_NUMBERS.md: 562 numbers, 0 untraced; results_review/SUMMARY.md: 119 numbers, 0 untraced.
- Compilation: `latexmk -pdf main.tex` runs without errors, undefined references or citations (11 pages); seven
  overfull boxes of 0.94-1.45 pt remain (empty boxes, as in v1).
- References: 45 cited, all from the 50 verified entries of refs.bib. Xing et al. 2018 is now verified from the full
  text (quotes with pages in results_final/CITATIONS.md); the PDF is kept in literature_pdfs/, which git ignores.
- Framing decision (2026-11-03, "Reframe"): lead with identifiability. 0.39 deg is the bound under a shared-parameter
  model, and the sensitivity to local variation is a central result that sets a co-location requirement.

## Reviewer points

| # | Point | Status | Where |
|---|---|---|---|
| 1 | Monte Carlo not independent of the surrogate | addressed | Sec. "Monte Carlo on femwell data", Table I column "fw.": 300 femwell solves at 20 perturbed geometries per set, 500 fits per set on femwell data. Two-width set 0.4034 deg = 1.037 x CRB, 0 fits at the box edge. Best set 0.2067 deg = 1.487 x CRB, 14 at the edge: it does not reach its bound, and the paper now says so (Results, Limitations). |
| 2 | Noise model optimistic: correlated errors | addressed | Sec. "Measurement-model variants", Fig. 7: rho sweep, derived rho = -1/sqrt(2); correlation lowers the bound |
| 2 | MZI precision separate from ring | addressed | Sec. "Ring and MZI readout": MZI 3x and 10x noisier, ring-only search; proposal is ring-only |
| 2 | Calibration | addressed | common offset: unknown 1.07 deg (two-width), 10 pm prior recovers the baseline |
| 2 | Temperature | addressed | shared unknown (0.1-1 K, unknown) and independent drift per structure (1-100 mK) |
| 2 | Lab vs production noise; 0.1 pm idealized | partly addressed | 1 pm is primary; linear scaling stated. New: measured MZI fitting errors of Xing et al. (n_eff 1.1e-5, n_g 1.0e-3) compared in Theory; our n_g noise matches, our n_eff noise is smaller, stated as a limitation. No production-tester noise source found. |
| 3 | Assumed uniformity; width-dependent etch / loading | addressed (analysis), open (data) | Sec. "Departures from shared process parameters", Fig. 6; bias gains, linear trend, random-effects bounds, doubling spreads. Discussion B now compares with Xing et al.'s measured local width variation (per-die std 1.26 nm, local up to 3 nm between copies of one design), far above the 0.039 nm doubling spread; differences between adjacent waveguides of different widths remain unreported. |
| 4 | Identifiability, not practical metrology | addressed | title, abstract, intro paragraph 3 ("computational study of test-structure design and identifiability") |
| 5 | Practical value of 0.39 deg | addressed | Sec. "Resonance prediction" and Discussion A: 0.457 / 0.773 pm predicted resonance vs 710.8 pm/deg, Zortman 841.5 pm, trimming +-32 pm |
| s | Literature positioning, novelty | addressed | Related work (4 paragraphs, 45 verified references cited), Contribution subsection; `paper/LITERATURE.md` lists the 5 closest papers |
| s | Explain the search | addressed | Sec. "Test Structures and Design Search": candidates, one-etch constraint, size <= 8, objective, exhaustive + greedy check |
| s | Material assumptions | partly addressed | n_Si, n_ox unknown with priors (stage D); dispersion and doping still not modeled |
| s | Ring geometry (known lengths and orders) | partly addressed | stated; common calibration offset tested; unknown length per ring not modeled |
| s | Width dependence | addressed | linear angle trend vs width (6.869 deg bias per 0.25 deg/100 nm, two-width) |
| s | Temperature excluded from primary results | addressed | temperature variants in Sec. "Measurement-model variants"; stage E in Fig. 3 |
| s | Generalizability | open | air cladding (R = 0.9782, no bound), TM ribs, other platforms listed as limitations/next steps |
| s | MZI redundancy | addressed | ring-only proposal |
| w | Weak literature review | addressed | see above |
| w | Bring the motivation and 711 pm forward | addressed | first paragraph of the introduction |
| w | Overly definitive language | addressed | "to first order", "for the oxide-clad cross-section", "simulations indicate", shared-parameter qualifier on every headline bound |
| w | Abstract too numeric | addressed | under 200 words, three numbers (0.39 deg, 1 pm, 0.039 nm), order problem, discovery, one rib, fix, precision and meaning, method |
| w | Placeholders | addressed | all `\todo` removed except e-mail, affiliation details and acknowledgments; co-author surname still `[LAST NAME]` |


## Sections you should rewrite yourself

1. **Introduction, paragraph 1** (the concrete problem). It is the reader's first contact; make it sound like your
   lab's motivation. Keep the three numbers (710.8 pm, 841.5 pm, +-32 pm) and the sources.
2. **Related work.** Xing 2018 and the NRC patent are now read in full; the rest is written from verified abstracts.
   Read at least Lu 2017, Butt 2024, Zhang 2024 and Dwivedi 2015 in full and correct anything the abstract-level
   summary gets wrong; add anything your advisor considers standard.
3. **Contribution list.** Decide how strongly you want to claim the co-location result; it is the part most likely to
   draw a reviewer's objection, because it argues against the practicality of the proposal.
4. **Discussion A and B** (foundry meaning, shared-parameter breakdown). These carry the paper's argument for a real
   test chip; your judgement of what a foundry would actually place is better than mine. The new comparison with
   Xing et al.'s local width variation (1.26 nm std, up to 3 nm local, against a 0.039 nm doubling spread) is the
   strongest argument against practicality in the paper. Decide how you want to frame it: as the requirement a
   co-located design must beat, or as the reason the measurement model needs per-structure terms.
5. **Monte Carlo on femwell data** (Results). Short and factual now; decide whether to keep the best set of eight
   in Table I at all, given that it does not reach its bound on femwell data.
6. **Next steps.** Replace with the experiment you can actually do (chip, process, metrology access).

## Still open

- No measurement; no air-clad bound; TM ribs; wall-shape sensitivity; unknown length per ring; a second full-vector
  solver for derivatives; a mesh/domain refinement of the bound itself (the domain check moves the two-width bound
  by 5.2%).
- One candidate reference (the arXiv version of Klein et al.) could not be verified and is not cited.
- The best set of eight does not reach its bound with femwell data (1.487 x CRB); a smoother forward model (fixed-mesh
  solves for the 300 nm ribs) would be needed to say whether that is model error or a property of the set.
