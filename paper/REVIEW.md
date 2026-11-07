# Review of the draft

v1 (reviewed 2026-11-03) is archived in `paper/archive/v1/`. This file covers v2, the rewrite after the external
critique of 2026-11-03. Numbers are checked by `python src/provenance_audit.py --paper` (report:
`results_final/audit/PAPER_AUDIT.md`); the review analyses are summarized in `results_review/SUMMARY.md` and are rows
L-P of `results_final/HEADLINE_NUMBERS.md`.

## Audit result (v2)

- Paper: 230 numbers with a decimal point or exponent, **0 untraced, 0 mismatches**. Chance-match rate 34%, so
  every new number was also placed by hand from its row (sections L-P of HEADLINE_NUMBERS.md).
- HEADLINE_NUMBERS.md: 508 numbers, 0 untraced; results_review/SUMMARY.md: 101 numbers, 0 untraced.
- Compilation: `latexmk -pdf main.tex` runs without errors, undefined references or citations; seven overfull boxes of
  0.94-1.45 pt remain (empty boxes, as in v1).
- Framing decision (2026-11-03, "Reframe"): lead with identifiability. 0.39 deg is the bound under a shared-parameter
  model, and the sensitivity to local variation is a central result that sets a co-location requirement.

## Reviewer points

| # | Point | Status | Where |
|---|---|---|---|
| 1 | Monte Carlo not independent of the surrogate | partly addressed | Sec. "Monte Carlo on femwell data"; 20 femwell-perturbed geometries per set, 500 fits on femwell data. See the status line below. |
| 2 | Noise model optimistic: correlated errors | addressed | Sec. "Measurement-model variants", Fig. 7: rho sweep, derived rho = -1/sqrt(2); correlation lowers the bound |
| 2 | MZI precision separate from ring | addressed | Sec. "Ring and MZI readout": MZI 3x and 10x noisier, ring-only search; proposal is ring-only |
| 2 | Calibration | addressed | common offset: unknown 1.07 deg (two-width), 10 pm prior recovers the baseline |
| 2 | Temperature | addressed | shared unknown (0.1-1 K, unknown) and independent drift per structure (1-100 mK) |
| 2 | Lab vs production noise; 0.1 pm idealized | partly addressed | 1 pm is primary; linear scaling stated; no production-tester noise source found |
| 3 | Assumed uniformity; width-dependent etch / loading | addressed (analysis), open (data) | Sec. "Departures from shared process parameters", Fig. 6; bias gains, linear trend, random-effects bounds, doubling spreads. No published value for adjacent-structure differences was found. |
| 4 | Identifiability, not practical metrology | addressed | title, abstract, intro paragraph 3 ("computational study of test-structure design and identifiability") |
| 5 | Practical value of 0.39 deg | addressed | Sec. "Resonance prediction" and Discussion A: 0.457 / 0.773 pm predicted resonance vs 710.8 pm/deg, Zortman 841.5 pm, trimming +-32 pm |
| s | Literature positioning, novelty | addressed | Related work (4 paragraphs, 50 verified references), Contribution subsection; `paper/LITERATURE.md` lists the 5 closest papers |
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

Status of point 1: see `results_review/SUMMARY.md` section 6.

## Sections you should rewrite yourself

1. **Introduction, paragraph 1** (the concrete problem). It is the reader's first contact; make it sound like your
   lab's motivation. Keep the three numbers (710.8 pm, 841.5 pm, +-32 pm) and the sources.
2. **Related work.** I wrote it from the verified abstracts only. Read at least Xing 2018, Lu 2017, Butt 2024, Zhang
   2024 and the NRC patent in full and correct anything the abstract-level summary gets wrong; add anything your
   advisor considers standard.
3. **Contribution list.** Decide how strongly you want to claim the co-location result; it is the part most likely to
   draw a reviewer's objection, because it argues against the practicality of the proposal.
4. **Discussion A and B** (foundry meaning, shared-parameter breakdown). These carry the paper's argument for a real
   test chip; your judgement of what a foundry would actually place is better than mine. Check the sentence that no
   published value for adjacent-structure differences exists.
5. **Next steps.** Replace with the experiment you can actually do (chip, process, metrology access).

## Still open

- No measurement; no air-clad bound; TM ribs; wall-shape sensitivity; unknown length per ring; a second full-vector
  solver for derivatives; a mesh/domain refinement of the bound itself (the domain check moves the two-width bound
  by 5.2%).
- Two candidate references could not be verified (arXiv rate limits) and are not cited.
