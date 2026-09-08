# Progress Log

## Session 1 (baseline, reconstructed from project_brief.md — see note below)

**Note on this entry:** this repo had no `notes/progress.md` when this session
started, even though `project_brief.md` references it twice. The numbers below
are reconstructed from project_brief.md's "Progress so far" section, which is
the only surviving record of the prior session. If a more detailed original
log exists elsewhere, it was not found in this repo and should replace this
reconstruction.

Built and validated a full-vector finite-difference mode solver (EMpy,
Fallahkhair & Li 2008 algorithm) against published neff~2.44-2.45 for the
500x220nm SOI benchmark waveguide (TE0). Computed d(neff)/d(w,t,angle) and
d(ng)/d(w,t,angle) for TE0 across widths 400-600nm. Built a 3-parameter
(width, thickness, sidewall angle) Fisher-information landscape over
two-TE-ring probe pairs.

**Key finding:** two same-type TE ring probes are badly conditioned
(cond(F) ~ 1e5-1e9) for jointly resolving width/thickness/angle — sidewall
angle sensitivity barely changes across ring widths, so it's nearly
degenerate with width/thickness combinations when using TE rings alone.
Condition number still swings ~30,000x depending on which two ring widths
are chosen — probe design matters even within an under-constrained setup.

**Known prior bug (fixed before, per project_brief.md):** grid quantization
initially killed the thickness derivative until sub-pixel/analytic overlap
averaging was added to `epsfunc_factory`.

---

## Session 2 — 2026-07-09: does structure-class diversity break the degeneracy?

### Housekeeping
- Repo had no `src/` or `notes/` directories; `fullvector_sensitivity.py` and
  `project_brief.md` were sitting at repo root. Restructured into
  `src/fullvector_sensitivity.py` and `notes/{project_brief,progress}.md`
  per this session's instructions.
- Environment: `EMpy` (the photonics package) was not installed, and `pip show EMpy`
  resolves to an unrelated PyPI package (`empy`, a templating engine) due to a
  case-insensitive name collision. The correct install is
  `pip install ElectromagneticPython`, which provides the `EMpy` import. Flagging
  this because a naive `pip install EMpy` silently installs the wrong package.
- The original script had a hardcoded output path (`/home/claude/...`), a
  leftover from a different (Linux/cloud) environment. Changed to a relative
  path so the script runs standalone on this Windows machine.
- Re-ran the untouched TE-only baseline first, before changing anything, to
  confirm the refactor doesn't silently change prior results: reproduced
  neff(500nm)=2.4446 and the exact same TE-TE condition-number landscape
  (best 240744, worst 7.32e9, 30388x spread) reported in Session 1. Refactor
  is confirmed behavior-preserving for the TE-only path.

### 1. TM-mode sensitivity (new)
Added polarization-aware mode selection: solve N=4 modes per geometry, classify
each by EMpy's `TEfrac()` (fraction of mode intensity that is TE-polarized),
and pick the highest-neff mode matching the requested polarization. Checked
across widths 400-600nm: mode index 0 is always clearly TE0 (TEfrac 0.96-0.99)
and mode index 1 is always clearly TM0 (TEfrac 0.02-0.03) — stable mode
ordering, but classifying by TEfrac rather than assuming a fixed index guards
against order swaps under the +/- finite-difference perturbations.

TM0 sensitivities (same finite-difference machinery as TE, reused unchanged):

| W(nm) | neff(TM) | ng(TM) | dneff/dw | dneff/dt | dneff/da | dng/dw | dng/dt | dng/da |
|---|---|---|---|---|---|---|---|---|
| 400 | 1.7374 | 3.483 | 6.22e-4 | 6.33e-3 | -2.39e-3 | 2.2e-3 | 3.1e-2 | -8.4e-3 |
| 450 | 1.7821 | 3.634 | 5.18e-4 | 6.90e-3 | -2.02e-3 | 1.4e-3 | 3.0e-2 | -5.6e-3 |
| 500 | 1.8204 | 3.739 | 4.97e-4 | 7.31e-3 | -2.00e-3 | 1.2e-3 | 2.8e-2 | -4.6e-3 |
| 550 | 1.8531 | 3.810 | 5.75e-4 | 7.60e-3 | -2.20e-3 | 1.2e-3 | 2.7e-2 | -4.5e-3 |
| 600 | 1.8798 | 3.857 | 6.17e-4 | 7.81e-3 | -2.36e-3 | 1.1e-3 | 2.6e-2 | -4.2e-3 |

**Verified numerically (task requirement, not assumed):** the
d(neff)/d(angle) ÷ d(neff)/d(thickness) ratio genuinely differs between
polarizations, by roughly a constant factor across widths:

| W(nm) | TE ratio | TM ratio | TM/TE |
|---|---|---|---|
| 400 | -1.563 | -0.378 | 0.242 |
| 450 | -1.150 | -0.293 | 0.254 |
| 500 | -0.964 | -0.274 | 0.284 |
| 550 | -0.968 | -0.290 | 0.299 |
| 600 | -0.966 | -0.303 | 0.314 |

TM's angle sensitivity is consistently ~3-4x smaller relative to its
thickness sensitivity than TE's is. This is a real, non-degenerate
difference in sensitivity-direction geometry, not noise (it's monotonic and
consistent across all 5 widths) — but see section 3 for why this alone
doesn't fix the Fisher conditioning problem.

### 2. MZI phase observable (new)
`d(phase)/d(param) = (2*pi/lambda) * d(neff)/d(param) * L`, L=100um, computed
directly from the already-validated dneff sensitivities (no new mode solves
needed). Example at 500nm:

| W(nm) | pol | dphase/dw (rad/nm) | dphase/dt (rad/nm) | dphase/da (rad/deg) |
|---|---|---|---|---|
| 500 | TE | 0.4177 | 1.6561 | -1.5962 |
| 500 | TM | 0.2014 | 2.9615 | -0.8107 |

**Important caveat found during implementation, not assumed:** at a *fixed*
(width, polarization), the MZI phase row is *exactly* a scalar multiple of
the ring's own dneff row (both come from the same dneff/d(w,t,a) vector,
just rescaled by `2*pi*L/lambda` and a different noise weight). It is **not**
an independent direction in parameter space unless paired with a different
width or polarization than the ring it's compared against. Confirmed this
numerically: every same-width ring+MZI Fisher matrix in section 3(c) below
is exactly rank-deficient (det ≈ 0 up to floating-point noise, including
small *negative* "determinants" like -2.4e-3, which is floating-point noise
around an exact zero for what must be a PSD matrix — flagging this
explicitly per the instruction to flag suspicious-looking numbers, but it is
explained, not a bug: J^T J with linearly dependent rows is exactly singular).

### 3. Mixed-structure-class Fisher-information comparison

All condition numbers below use the same noise model as Session 1's ring
observables (sigma_neff=6.5e-6, sigma_ng=1.0e-3, both applied identically to
TE and TM rings since both are read out the same way — resonance + FSR).
The MZI phase observable uses an **assumed, not independently validated**
sigma_phase=1e-3 rad (typical order for interferometric photodetection, but
not calibrated against a real measurement — flagged as an open assumption).
Swept sigma_phase from 1e-6 to 1e-2 rad to check robustness (see below) —
the qualitative conclusion below does not depend on this choice.

| Config | Best cond(F) | Worst cond(F) |
|---|---|---|
| (a) TE ring + TE ring [baseline] | 2.41e5 | 7.32e9 |
| (b) TE ring + TM ring | 3.46e8 | 1.54e12 |
| (c) TE ring + MZI | 7.41e5 | ~1e16 (exactly singular at matching width/pol, see above) |
| (d) 3-probe: TE ring + TM ring + MZI-TE, all @ 500nm | 1.25e9 | — |
| (d) everything: TE ring + TM ring + MZI-TE, all 5 widths (25 observables) | 5.69e5 | — |

**Follow-up robustness check:** starting from the single best (a) pair found
(TE ring @550nm + TE ring @600nm, cond=2.41e5) and adding one more probe on
top (rather than building fresh 2-probe combos), every single addition
tested — TM ring at any width, MZI at any width/polarization, or even a
*third TE ring* — makes the condition number *worse*, not better. Best
addition found was a third TE ring (cond 3.03e5), not a different structure
class.

### KEY FINDING: structure-class diversity, as tested, does NOT fix the
### condition-number degeneracy — and the reason is informative

Answering the direct question: **no**, adding a TM ring or an MZI phase
observable does not improve on a well-chosen pair of TE rings, by the
condition-number metric. Every mixed-class configuration tested has a worse
best-case condition number than the TE-only baseline; some (same-width
ring+MZI) are exactly singular.

**But this is not simply "diversity doesn't help" — a follow-up eigenvalue
decomposition shows the condition-number metric itself is likely the wrong
diagnostic.** Checking lambda_min(F) (the quantity that actually bounds the
achievable precision on the worst-determined parameter, via Cramér-Rao)
instead of cond(F) = lambda_max/lambda_min:

| Addition to TE(550)+TE(600) base | lambda_min before | lambda_min after | cond(F) after |
|---|---|---|---|
| + TE ring @400nm (same class) | 6.33 | 7.64 | 3.03e5 (worse than base's 2.41e5) |
| + TM ring @400nm | 6.33 | 6.96 | 3.61e5 (worse) |
| + MZI(600nm, TE) | 6.33 | **10.21** | 6.61e5 (worse) |

Every addition *increases* lambda_min (as it must — Fisher information is
additive and PSD, so lambda_min can only grow), meaning every added probe
*does* genuinely improve the worst-case parameter precision. But cond(F)
gets worse anyway, because lambda_max grows even faster (added probes are
directionally correlated with the already-strong direction, since all
observables tested are neff/ng-based measurements of the same mode family).
**By lambda_min, MZI is actually the best of the three additions tested** —
the opposite ranking from what cond(F) suggests.

This means: (1) the original condition-number diagnostic conflates "well
distributed sensitivity directions" with "large amount of information along
the strong direction," and probably should be replaced or supplemented by
an E-optimality-style criterion (maximize lambda_min) or A-optimality
(minimize trace(F^-1)) — standard alternatives from optimal-experiment-design
theory that this project hasn't used yet; (2) structure-class diversity
*may* still be the right lever, but the paper needs the right objective to
see it — cond(F) alone would have led to discarding MZI/TM entirely, even
though MZI gave the largest lambda_min improvement of anything tested.

### 4. TM-mode solve validation
- **Grid/domain convergence:** TM0 @ 500x220nm stable at 1.80-1.82 across
  domain half-widths (1.4-2.0um) and mesh densities (30-45 core points) —
  same order of stability as TE0 (2.43-2.45 across the same sweep).
- **Independent cross-check (effective-index method, no EMpy/FD dependency):**
  implemented the standard 2-step slab-dispersion EIM from scratch
  (raw numpy/scipy, `_slab_neff` + `eim_neff` in the script) as a fully
  independent computation.

  | W(nm) | FD TE0 | EIM TE0 | %diff | FD TM0 | EIM TM0 | %diff |
  |---|---|---|---|---|---|---|
  | 400 | 2.2308 | 2.3210 | +4.04% | 1.7374 | 1.7903 | +3.05% |
  | 500 | 2.4446 | 2.4946 | +2.05% | 1.8204 | 1.8459 | +1.40% |
  | 600 | 2.5620 | 2.5966 | +1.35% | 1.8798 | 1.8864 | +0.35% |

  TE0's EIM-vs-FD gap (1.3-4%, decreasing with width) matches the
  well-documented systematic EIM overestimate for high-index-contrast strip
  waveguides. TM0 shows the *same sign and comparable magnitude* of gap,
  which is the expected signature if the TM solve is correct — if the TM
  mode-picking logic were buggy (e.g. silently returning a higher-order
  mode), there would be no reason for its EIM discrepancy to track TE's so
  closely.
- **Literature search — flagged, not fully resolved:** could not find a
  directly comparable published TM0 value for exactly 500x220nm SOI at
  1550nm; both full-text sources found via search (an MDPI paper, an IEEE
  paper) returned 403/no-access. One indirect data point surfaced via search
  snippet: TM0 ≈ 1.74 for a *460nm-wide* (not 500nm) 220nm-thick guide at
  1550nm — in the right ballpark and consistent with our width-trend
  (TM0 increasing with width), but not independently full-text-verified, so
  **treated as suggestive, not confirmatory**. The EIM cross-check above is
  the primary validation evidence for the TM solve; a cleaner published
  500x220nm TM0 reference would still be worth locating if time allows.

### Next steps
1. **Re-run the Fisher landscape with lambda_min (or trace(F^-1), or D-optimality
   det(F)) as the comparison metric instead of / alongside cond(F).** This is
   the most important next step — the condition-number framing used so far
   may have been hiding a real, positive effect from structure diversity.
2. If lambda_min-based reranking shows MZI/TM genuinely helping, redo the
   full-width-pair landscape scan (like section 3 above) but reporting
   lambda_min instead of cond(F), to see whether the width/class choice that
   maximizes lambda_min is stable and interpretable.
3. Get an actual sigma_phase estimate from realistic MZI photodetection
   noise (shot noise / TIA noise budget) rather than the assumed 1e-3 rad —
   swept 1e-6 to 1e-2 and the cond(F) conclusion was robust to this choice,
   but a lambda_min-based re-analysis should redo this sweep too.
4. Still open: try a structure with a genuinely different field distribution
   (slot waveguide, sub-wavelength grating) rather than just a different
   polarization/readout of the same strip-waveguide mode family — since TE
   ring, TM ring, and MZI here are all still "neff of a solid strip
   waveguide," just measured three ways, which may explain why they're more
   directionally correlated than hoped.
5. Locate or compute a cleaner independent TM0 @ 500x220nm@1550nm reference
   to fully close out the validation flagged in section 4.

---

## Session 3 — 2026-07-12: sensitivity database, atlas, multi-metric Fisher, slot waveguide, information map

**Scope discipline for this session:** no optimization, Bayesian optimization,
inverse design, or ML was written or run. Everything below is diagnostic:
does a rich, structured information landscape actually exist before any of
that is justified.

### 0. Refactor (Task 1 infrastructure)
`src/fullvector_sensitivity.py` was refactored to separate the
geometry-agnostic mode-solving core from the strip-specific geometry:
- `solve_mode_from_eps(ep, wl, polarization, n_modes, X, Y)` -- new, takes
  any already-built `epsfunc` and does the solve+TEfrac-classify step. This
  is the reused core; `solve_mode` (strip) is now a thin wrapper around it.
- `neff_ng` and `sensitivities` now accept a `solve_fn=` parameter
  (default: strip's `solve_neff`), so the exact same finite-difference
  derivative machinery is reused by any new geometry rather than
  reimplemented. `slot_waveguide.py` uses this.
- Confirmed behavior-preserving: re-ran a regression self-test
  (`fullvector_sensitivity.py`'s new trimmed `__main__`) reproducing
  TE0=2.4446, TM0=1.8204 @ 500nm exactly, before building anything on top.
- Repo output layout: `data/` (CSV+JSON tables) and `figures/` (PNGs) added
  at repo root, alongside `src/`.

### 1. Sensitivity database (Task 1)
`src/probe_db.py` builds two tables:
- `data/sensitivity_db.{csv,json}` -- **15 raw records**, one per
  (geometry, polarization, width): strip TE x5 widths, strip TM x5,
  slot TE x5 (400-600nm strip / 150-210nm slot rail, see section 5). Each
  record carries neff, ng, all 6 first derivatives, a normalized
  sensitivity vector (raw / fab-tolerance-scaled / unit, see below),
  simulation settings, convergence status, and validation status as
  strings (not just numbers) so provenance travels with the data.
- `data/probe_db.{csv,json}` -- **30 derived records**: each raw record
  wrapped by a `ring` readout (2 Fisher rows, reusing Session 2's
  SIGMA_NEFF/SIGMA_NG) and an `mzi` readout (1 Fisher row, SIGMA_PHASE),
  so `fisher_reanalysis.py` and future work can assemble any probe
  combination without re-deriving Jacobians by hand.

**Normalization convention (a documented modeling choice):** each raw
vector (dneff/dw, dneff/dt, dneff/da) is scaled by assumed typical
fabrication tolerances -- sigma_w=5nm, sigma_t=2nm, sigma_a=1deg (the same
magnitudes already used as finite-difference step sizes, for consistency)
-- before unit-normalizing, so the three axes are comparable despite
different physical units. This choice was not swept/stress-tested this
session; it directly affects the atlas geometry in section 2 and is the
most important assumption to revisit if atlas conclusions look sensitive
to it later.

### 2. Sensitivity atlas (Task 2) -- candidate Figure 1
`src/atlas.py`. For all 15 raw records: pairwise cosine similarity, angle,
and Euclidean distance (the last computed on the *scaled*, not unit,
vectors, so it reflects magnitude differences too -- deliberately
different information than angle/cosine). Saved as full matrices
(`data/atlas_{cosine_similarity,angle_deg,euclidean_distance}.csv`) plus
`figures/fig1_sensitivity_atlas.png` (3 pairwise-axis projections + a
hierarchically-clustered cosine-similarity heatmap).

**Quantified, not eyeballed:**
- within-class pairs (same geometry+polarization, different width only):
  n=30, angle 0.1-13.9 deg, mean 4.0 deg -- essentially redundant directions,
  confirming Session 1/2's qualitative finding numerically.
- cross-class pairs (different geometry or polarization): n=75, angle
  20.9-73.7 deg, mean 46.4 deg -- genuinely different directions.
- most redundant pair: strip-TM-450 <-> strip-TM-550 (0.11 deg).
- most orthogonal pair: strip-TM-500 <-> slot-TE-210 (73.7 deg).
- hierarchical clustering leaf order groups perfectly into 3 blocks
  (all slot-TE together, all strip-TM together, all strip-TE together) --
  no mixing, visible directly in the heatmap and confirmed by the leaf order.

### 3. Multi-metric Fisher re-analysis (Task 3)
`src/fisher_reanalysis.py`. Recomputed every Session-2-style 2-probe
landscape using cond(F), lambda_min (E-optimality), det(F) (D-optimality),
and trace(F^-1) (A-optimality), plus a new (d) config (TE ring + slot TE
ring) and an isolated slot-contribution test.

**Numerical hygiene note:** the exactly-rank-deficient same-width
ring+MZI points identified analytically in Session 2 reappeared here as
near-singular matrices whose numerical trace(F^-1) can go *negative* from
floating-point roundoff around an exact zero eigenvalue (impossible for a
real PSD matrix). These points (cond(F) > 1e12) were excluded from ranking
correlations and reported separately, rather than silently corrupting the
correlation numbers -- 5 such points in config (c), 1 in config (b), 0 in
(a)/(d).

**Ranking agreement across metrics (Spearman correlation, well-posed
subset only):**

| Config | cond~lambda_min | cond~det | cond~trace_inv | det~trace_inv |
|---|---|---|---|---|
| (a) TE ring + TE ring | +1.000 | +0.465 | +0.995 | +0.469 |
| (b) TE ring + TM ring | +0.997 | +0.888 | +0.997 | +0.903 |
| (c) TE ring + MZI | +0.995 | +0.547 | +0.989 | +0.562 |
| (d) TE ring + slot TE ring | **+0.859** | +0.688 | +0.859 | +0.927 |

**Key finding: metric agreement is high within same-or-similar-class
combos and degrades specifically when a genuinely new direction (slot) is
introduced** -- config (d) is the only one where cond(F) and lambda_min
disagree on the *best* combo (cond picks stTE550+slTE165; lambda_min picks
stTE400+slTE165). This is a self-consistent, non-arbitrary pattern: metrics
mostly agree when there's nothing genuinely new to disagree about, and
diverge specifically at the point where something informationally new
enters -- which is itself evidence that slot is doing something the
others weren't.

**Isolated slot contribution** (apples-to-apples: same strip-only probe
set, with/without slot's 15 additional Fisher rows added on top):

| | n observables | cond(F) | lambda_min | trace(F^-1) |
|---|---|---|---|---|
| strip only (TE+TM, ring+mzi, all widths) | 30 | 5.79e5 | 1.45e2 | 6.88e-3 |
| + slot TE (ring+mzi, all widths) | 45 | 4.09e4 | 1.25e4 | 8.00e-5 |
| **improvement** | | **14.2x** | **86.0x** | **86.0x** |

This is a much larger, more consistent improvement (across 3 of 4 metrics)
than anything TM or MZI produced in Session 2. Note this specific
comparison uses a much larger probe budget (30/45 observables) than the
2-probe landscapes above -- it answers "does slot help on top of an
already-large strip set," a different (complementary) question from "what's
the best 2-structure choice," not a contradiction of it.

**Which metric is actually right for this problem?** Cramér-Rao says
lambda_min(F) bounds the worst-case achievable parameter precision, and
trace(F^-1) bounds the average -- both map directly onto "can we identify
width, thickness, AND angle," which is the actual scientific question. det(F)
(confidence-ellipsoid volume) rewards raw information accumulation somewhat
independent of *shape*, which is why it diverges from the other three most
sharply (Spearman as low as +0.465 vs cond). cond(F) alone -- what Session 2
used throughout -- conflates "how informative is the strongest direction"
with "how well-distributed is the information," and Session 2's own
follow-up eigenvalue check already caught it giving a misleading answer.
**Recommendation: use lambda_min (E-optimality) as the primary metric going
forward**, with trace(F^-1) as a secondary check (they agree well, Spearman
+0.859 to +1.000 across all 4 configs) and det(F) reported for context but
not optimized for; cond(F) should be treated as a numerical-conditioning
diagnostic, not a design objective.

### 4. Does the atlas explain the Fisher results? (Task 4)
`src/atlas_vs_fisher.py`. Computed the actual Spearman correlation between
atlas pairwise angle and each pair's Fisher metrics (ring readout, n=105
pairs), rather than assuming "more orthogonal = more informative."

**Counterintuitive result, confirmed not a bug:**
```
Spearman(angle, cond)       = +0.478   (wrong sign: larger angle -> worse conditioning)
Spearman(angle, lambda_min) = -0.369   (wrong sign: larger angle -> worse worst-direction precision)
Spearman(angle, det)        = +0.558
Spearman(angle, trace_inv)  = +0.369   (wrong sign)
```
**Mechanism (verified numerically, not hand-waved):** printed each raw
record's unit vector's angle-axis component directly:

| family | angle-axis component (unit vector) |
|---|---|
| strip-TE | -0.38 to -0.48 |
| strip-TM | -0.13 to -0.18 (mostly thickness: 0.97-0.98) |
| slot-TE | **-0.60 to -0.62** (largest in the dataset) |

strip-TM's large mutual angle from strip-TE is mostly explained by TM
rotating toward the *thickness* axis (which TE already partially covers),
not the angle axis -- so despite "looking different" by mutual angle, TM
doesn't bring much new angle-axis information, exactly matching Session 2's
finding that TM gave the smallest lambda_min improvement of the three
additions tested. slot-TE has the single largest angle-axis component of
anything simulated, exceeding even the best strip-TE width -- this is the
mechanistic reason slot's Fisher contribution (section 3) was so much larger.

Confirmed quantitatively: summed *individual* angle-axis content
(|unit_a_i|+|unit_a_j|) correlates with lambda_min at **+0.117** (correct
sign) vs. mutual pairwise angle's **-0.369** (wrong sign) -- individual
axis-alignment is a better predictor than mutual dissimilarity, though
still weak/incomplete, because (a) the atlas as built uses only the dneff
direction (not dng, which also feeds the ring's Fisher rows) and (b) it
doesn't encode absolute noise-weighted magnitude. **Answering Task 4's
questions directly:**
- *Are clustered probes redundant?* Yes -- within-class pairs are
  near-parallel (mean 4 deg) and Session 2's cond(F) landscape shows only
  modest (~30,000x) spread from width choice alone.
- *Are orthogonal probes consistently more informative?* **No** -- mutual
  angle is a weak and even wrong-signed predictor. What matters is each
  probe's own alignment with the specific poorly-constrained axis (sidewall
  angle here), not how different it looks from its pairing partner.
- *Does the atlas predict Fisher performance?* Only partially, and not via
  the naive "bigger angle is better" reading -- the atlas's per-axis
  component decomposition is informative (it correctly flagged slot as the
  angle-rich structure before Fisher analysis confirmed it), but pairwise
  angle/cosine alone is not a reliable proxy, and this was the most
  important negative sub-result of the session precisely because it would
  have been easy to over-read the "candidate Figure 1" plot naively.
- *Unoccupied regions?* See section 6 (this was answered more precisely
  via PCA than via the 2D-projection wedge scan below, but for completeness:
  the empty-wedge scan shows every 2D projection is occupied by only a
  narrow angular sector -- e.g. width-vs-angle spans essentially one line
  (357.8 deg empty) since dneff/dw > 0 and dneff/da < 0 for every structure
  tested -- consistent with section 6's finding that all tested structures
  collapse onto a thin sub-manifold of the full 3D space.

### 5. Slot waveguide (Task 5) -- the one new structure class
`src/slot_waveguide.py`. Two silicon rails (nominal 180nm each), 220nm
thick, separated by a nominal 120nm gap, reusing the strip's exact
per-cell analytic sub-pixel dielectric averaging (applied twice, once per
rail) and the exact same `solve_mode_from_eps`/TEfrac classification core
-- no new discretization or mode-solving logic was written.

**Fabrication-parameter modeling choice (documented, not the only possible
one):** the single scalar "width" parameter is modeled as a coupled
lithographic CD bias -- rail_width(w) = w, gap(w) = GAP0 - (w - RAIL0) --
so a positive bias simultaneously widens both rails AND narrows the gap,
matching how a real dose/etch bias erodes or grows every silicon boundary
the same way. A simpler rail-only model (gap held fixed) was deliberately
rejected because it would have missed the specific reason slot waveguides
are interesting: their well-known strong gap-sensitivity.

**Validation** (same discipline as Session 2's TM validation):
- *Grid/domain convergence:* 1.40% neff spread across domain half-width
  1.6-2.2um, core mesh 50-70pts -- same order as strip/TM's convergence checks.
- *Physical consistency check (used in place of an analytic EIM, which
  doesn't have a simple closed form for coupled-rail slot modes):* as gap
  grows (120nm to 1000nm), slot TE0 converges monotonically toward the
  independently-validated strip TE0(rail_width) solver's value, agreeing to
  within 0.13% at 1000nm gap -- exactly the physically-required decoupled
  limit, using already-validated code as ground truth rather than a new
  external method.
- *Literature comparison:* attempted (Almeida et al. 2004, an MDPI paper, a
  Nature Scientific Reports paper, Wikipedia) -- **no full-text-accessible
  published neff table found for this exact geometry** (paywalled/blocked
  in every case tried). Computed neff (~1.6-1.9) is in the qualitatively
  expected range commonly cited for slot waveguides at telecom dimensions,
  but this is flagged, per the working-style instruction, as not
  independently literature-confirmed -- exactly like TM0's status in
  Session 2.
- *Sanity-checked a nonlinearity, not assumed it was a bug:* dneff_dw and
  dneff_da jump noticeably at rail>=195nm (gap<=105nm) rather than
  following the smoother trend at smaller rail widths. Re-verified on a
  finer grid (70x24 vs 50x20 mesh) -- the jump persists essentially
  unchanged (dneff_dw: 0.01034 vs 0.01034; dneff_da: -0.03852 vs -0.03852),
  ruling out a discretization artifact. Attributed to the well-known
  nonlinear field-enhancement scaling of slot waveguides as the gap narrows
  (E-field enhancement scales roughly as 1/gap) -- genuine physics, not noise.

Slot's sensitivity vectors were added to the database and atlas (sections
1-2) alongside strip TE/TM. **Quantified, not eyeballed, that it occupies a
genuinely different direction:** mean angle to strip-TE = ~48 deg, to
strip-TM = ~57 deg (both comfortably inside the "cross-class" distribution,
section 2); its angle-axis component (-0.60 to -0.62) is the largest of any
structure simulated. Slot TM mode was confirmed to exist (mode1, TEfrac~0.02
at nominal dimensions) but was not simulated this session -- cheap follow-up.

### 6. Information map (Task 6)
`src/information_map.py`. PCA (from-scratch via SVD, no sklearn dependency)
on the 15 raw records' fab-tolerance-scaled sensitivity vectors.
**t-SNE/UMAP deliberately not used** -- both are designed to unfold local
neighborhood structure in high-dimensional, large-n data, and are known to
produce unstable, over-interpretable embeddings for n=15 in a space that's
already only 3D; using them here would manufacture apparent structure
rather than reveal it. This is a documented methodological choice, not a
skipped step.

```
PC1: 95.3% variance, loadings (width, thickness, angle) = (+0.785, -0.191, -0.589)
PC2:  4.7% variance, loadings = (-0.170, -0.981, +0.091)
PC3:  0.009% variance, loadings = (-0.596, +0.029, -0.803)
```

**PC3 is the headline result.** Every one of the 15 structures simulated
this session -- across 3 geometry/polarization families and 5 widths each
-- has a PC3 score under 0.0005 in magnitude (vs. PC1 scores up to 0.05).
That means all currently-tested probes collapse onto a thin, essentially
2D sub-manifold of the full 3D fabrication-sensitivity space, and there is
a specific, quantified, entirely unpopulated direction -- a
width/sidewall-angle combination distinct from the one slot itself
provides (PC3's loading is (-0.596, +0.029, -0.803), vs. slot's own
dominant direction which sits along PC1) -- that no structure tried so far
reaches, regardless of geometry class or width.

Family centroids are well-separated relative to their own intra-family
spread (slot<->TE: 6.3x, slot<->TM: 7.9x, TE<->TM: 12.7x), confirming
geometry families occupy small, distinct, non-overlapping regions --
answering Task 6's central question: **yes, current process monitors
occupy only a small region of Fisher-information space, and the gap
(PC3) is now identified and quantified, not just asserted.**

Figure: `figures/fig2_information_map.png`.

### 7. Scientific interpretation (Task 7)

**What has been learned.** Three independent lines of evidence this
session -- the atlas's cross-class vs. within-class angle statistics
(section 2), the multi-metric Fisher re-analysis and isolated slot
contribution (section 3), and the PCA information map (section 6) -- all
converge on the same conclusion from different directions: structure-class
diversity is real and quantifiable in this problem, but not all diversity
is equally useful, and the earlier metric (condition number) used to judge
it in Session 2 was specifically bad at detecting the useful kind.

**Session 2's negative conclusion is revised, not contradicted.** Session 2
found that condition number got *worse* whenever a probe was added, no
matter how different its structure class, and initially read this as
"structure-class diversity doesn't help." Session 2's own final eigenvalue
check already put a crack in that reading (lambda_min always improved).
This session shows why more completely: cond(F) and lambda_min/trace(F^-1)
agree well on similar-class combos and specifically diverge when a
genuinely new direction (slot) is introduced (section 3's table). The
"diversity doesn't help" conclusion was an artifact of the metric, not a
property of the physics -- once lambda_min is used, slot delivers an 86x
improvement over the strip-only probe set. This is exactly the kind of
revision the working style asks to document plainly rather than paper over.

**Is there a meaningful optimization landscape?** The strongest evidence
for "yes" is section 6's PC3 result: a specific, unoccupied,
sub-1%-of-variance direction exists that none of 15 structures across 3
families and 10 widths reach. That is concrete, falsifiable evidence of
unexploited design space, not a hand-wavy "there's probably more out
there." The strongest evidence for caution is that this is still a small
sample (n=15, 3 families) and the normalization/noise-model assumptions
(section 1, section 3's ring-noise-model caveat) haven't been stress-tested
-- the landscape's *shape* is credible, its exact *contours* are not yet
pinned down.

**Are existing structures redundant?** Within a class, largely yes (mean
4 deg pairwise angle, modest condition-number spread from width alone).
Across classes, no -- and this is the more important finding, since a
"process monitor suite" built purely from same-class-different-width
probes (which is close to current foundry PCM practice, per the project
brief's framing) is leaving real, quantifiable information on the table.

**Do new structure classes provide complementary information?** Yes,
demonstrably for slot -- and mechanistically explained, not just
empirically observed: slot's sensitivity direction is unusually
angle-rich (section 4), which is specifically the axis that was hardest to
constrain with rings/MZIs alone in Session 2. This is a much stronger claim
than "slot has a different cond(F)" -- it's "slot is different in the
specific way that matters for this problem's known weak point."

**Is inverse design still justified?** Cautiously yes, with two important
refinements to how it should now be framed relative to the original project
brief: (1) the objective function should be E-optimality (lambda_min) or
A-optimality (trace(F^-1)), not condition number -- optimizing cond(F)
directly, as originally implied, would likely have led inverse design
toward misleading targets, given section 3's evidence that cond(F)
disagrees with the theoretically-motivated metrics specifically in the
mixed-structure-class regime where inverse design would need to operate;
(2) the PC3 gap (section 6) gives a concrete, physically-interpretable
target -- an inverse-design search (or even further hand-designed
candidates) should be judged partly by whether it populates that direction,
not just by whether it improves an aggregate scalar. Absent this session's
diagnostic work, an inverse-design phase would have been optimizing a
metric already shown to be misleading, against a design space whose
unoccupied regions were unknown.

**Honest limitations carried forward:**
- Atlas normalization (assumed fab tolerances sigma_w/t/a) is an untested
  modeling choice; atlas geometry (and therefore some Task 4/6 conclusions)
  could shift under a different choice.
- Ring/MZI/slot-ring noise model reuses Session 2's strip-ring
  SIGMA_NEFF/SIGMA_NG for every geometry and polarization, including slot --
  an unverified simplification; a slot ring's real achievable Q and
  measurement precision could differ materially from a strip ring's.
- The atlas direction vector uses only dneff, not dng -- shown in section 4
  to make it an incomplete predictor of full ring-Fisher behavior. A fuller
  atlas including dng direction and per-observable noise-weighting is a
  natural, well-motivated next step now that plain dneff-direction has been
  shown insufficient.
- No literature-confirmed neff reference was found for the slot geometry
  simulated (search attempted, consistently paywalled) -- flagged exactly
  like TM0's equivalent gap in Session 2, not glossed over.
- Slot TM mode exists (confirmed) but wasn't simulated this session.
- n=15 raw structures is still a small sample for the PCA/clustering claims
  in section 6; the qualitative shape (thin sub-manifold, clear family
  separation) is well-supported by the numbers, but exact PC-variance
  percentages should be expected to shift as more structures are added.

### Next steps
1. Add dng-direction and noise-weighting into the atlas vector definition
   (section 4's identified gap) and re-run Task 4's correlation analysis --
   check whether the wrong-signed angle/cond correlation resolves once the
   atlas better reflects what actually drives Fisher rows.
2. Simulate slot TM (confirmed to exist, cheap addition) and re-run the
   full atlas/Fisher/PCA pipeline with 20 (not 15) raw structures.
3. Design or search (still not "optimize" -- more candidate hand structures)
   for something that populates the PC3 gap identified in section 6, as a
   concrete falsification test of whether that direction is reachable by
   any realistic photonic structure at all.
4. Stress-test the fab-tolerance normalization (sigma_w=5, sigma_t=2,
   sigma_a=1) and the slot ring noise-model assumption -- both are
   currently single-point, undocumented-sensitivity choices that feed
   directly into this session's headline numbers.
5. Only after (1)-(4): begin scoping what an E-optimality/A-optimality-based
   inverse design objective would actually look like, now that condition
   number has been ruled out as the right target.

---

## Session 4 — 2026-07-13: external validation, dng-aware atlas, stress tests, slot TM, rib/SWG search for the PC3 gap

**Scope discipline maintained:** no optimization, Bayesian optimization,
inverse design, or ML written or run this session either. Everything below
is still diagnostic -- specifically, closing out the gaps Session 3 flagged
before any of that would be justified.

### External validation (added to notes/project_brief.md, not repeated here)
A literature-grounding section was appended to `project_brief.md` this
session: sidewall angle is a real, spatially-varying fabrication error
(etch loading effects), it's independently reported as performance-limiting
(a photonic-crystal-cavity Q-factor study), and -- the direct hit -- a 2023
process-variation-aware compact-modeling paper found their standard optical
extraction fitting procedure's error-function minimum stayed roughly flat
across sidewall angles 85-90deg, i.e. published, independent confirmation
that standard optical extraction can't resolve sidewall angle, matching this
project's own Fisher/PCA finding. Current industry practice is destructive
SEM cross-sectioning, not a smarter optical technique. This reframes the
PC3 gap from Session 3 as specifically "find a structure enabling
non-destructive, all-optical sidewall-angle characterization," not just
"fill an abstract empty direction" -- see `project_brief.md` for the full
section (this doesn't change any of this session's numbers, only their
framing/motivation).

### Task 1: dng-aware atlas vector, re-run Task 4 correlation
Session 3's atlas used only the dneff direction; `src/atlas_v2.py` rebuilds
the per-structure vector as the concatenated, ALREADY noise-weighted ring
Jacobian rows -- `[dneff_dw,dneff_dt,dneff_da]/SIGMA_NEFF` stacked with
`[dng_dw,dng_dt,dng_da]/SIGMA_NG`, unit-normalized in 6D -- i.e. exactly
what a ring's Fisher contribution is actually built from, not a new
modeling choice layered on top.

**Result: v1 and v2 angles are nearly identical** (e.g. strip-TE-400 vs
strip-TM-500: 39.8deg v1 vs 42.5deg v2; full comparison table in the script
output). This is expected once you note SIGMA_NEFF (6.5e-6) vs SIGMA_NG
(1.0e-3) differ by ~150x, so the neff-row dominates the 6D vector's norm
and dng barely rotates it.

**Re-running Session 3's correlation with v2 angle (`src/atlas_v2_vs_fisher.py`):**
```
Spearman(v2 angle, cond)       = +0.468   [v1 was +0.478]
Spearman(v2 angle, lambda_min) = -0.381   [v1 was -0.369]
Spearman(v2 angle, det)        = +0.518   [v1 was +0.558]
Spearman(v2 angle, trace_inv)  = +0.381   [v1 was +0.369]
```
**The wrong-signed correlation does NOT resolve -- it's essentially
unchanged.** This is itself a useful negative result: it confirms Session
3's mechanistic explanation was the right one. The problem was never "the
atlas is missing dng information" -- it's the more fundamental geometric
fact that mutual pairwise angle between two vectors doesn't tell you
whether their shared subspace happens to align with the specific
poorly-constrained axis (sidewall angle). Adding dng didn't fix it because
dng wasn't the missing ingredient; individual axis-alignment (established
in Session 3) remains the better (if still incomplete) predictor.

### Task 2: stress-testing two flagged assumptions
`src/stress_tests.py`.

**(a) Fab-tolerance normalization sweep** (sigma_w, sigma_t, sigma_a
scaling used to build atlas vectors) -- swept 7 configurations from
(2,2,1) to (5,5,1) to (1,1,1, i.e. no rescaling at all, an extreme control
case):

| (sigma_w,sigma_t,sigma_a) | within-mean | cross-mean | PC1% | PC2% | PC3% |
|---|---|---|---|---|---|
| (5,2,1) [Session 3 baseline] | 4.0 | 46.4 | 97.2 | 2.8 | 0.01 |
| (10,2,1) | 3.7 | 42.2 | 97.4 | 2.6 | 0.01 |
| (2,2,1) | 4.0 | 47.3 | 97.7 | 2.3 | 0.01 |
| (5,5,1) | 4.3 | 45.9 | 98.6 | 1.4 | 0.01 |
| (5,2,0.5) | 4.0 | 47.0 | 97.4 | 2.6 | 0.01 |
| (5,2,2) | 3.8 | 43.7 | 97.2 | 2.8 | 0.01 |
| (1,1,1) [no rescaling, control] | 3.8 | 44.9 | 97.1 | 2.8 | 0.00 |

**Robust:** within/cross-class separation and the thin-2D-manifold finding
(PC3 < 0.01% everywhere) hold across every configuration tested, including
the no-rescaling control. Slot's |angle-axis component| was checked
against strip-TE and strip-TM in every configuration and **was the largest
of the three in all 7** -- the qualitative "slot is the angle-rich
structure" finding is not an artifact of the specific sigma choice.

**(b) Slot-ring noise-model stress test** -- Session 3's headline "slot
gives an 86x lambda_min improvement" assumed slot rings have the SAME
measurement precision (SIGMA_NEFF, SIGMA_NG) as strip rings. Slot
waveguides are well documented to have higher propagation loss than strip
(sidewall-roughness scattering is proportionally worse given field
enhancement at the narrow gap's Si-oxide interfaces), which would mean
lower ring Q and worse (larger) sigma in practice. Re-ran the isolated
contribution test under pessimistic slot-noise multipliers:

| slot noise vs strip | n_obs | cond(F) | lambda_min | lambda_min improvement |
|---|---|---|---|---|
| same as strip [Session 3 assumption] | 45 | 4.09e4 | 1.25e4 | 86.0x |
| 3x worse | 45 | 6.88e4 | 1.60e3 | 11.0x |
| 5x worse | 45 | 1.29e5 | 7.07e2 | 4.9x |
| 10x worse | 45 | 2.83e5 | 3.03e2 | 2.1x |
| 30x worse (extreme) | 45 | 5.13e5 | 1.64e2 | 1.1x |

**This is NOT robust.** The 86x headline number requires an optimistic (and
unverified) noise assumption. Under a plausible pessimistic assumption
(5-10x worse, consistent with typical reported slot-vs-strip loss ratios),
the improvement is a much more modest 2-5x. Important nuance: it never
reverses sign (lambda_min never gets WORSE than the strip-only baseline,
consistent with Fisher information being strictly additive/PSD) -- so
"slot helps" survives as a qualitative claim, but "slot helps by ~86x"
does not, and should not be repeated without an actual measured or
literature-sourced slot-ring Q estimate.

### Task 3: slot TM added (20 raw structures)
`src/extend_db.py` (idempotent -- checks existing IDs before appending, a
bug caught and fixed this session after an initial careless re-run
silently duplicated 5 records; flagging the fix, not hiding it). Slot TM
sensitivities computed across the same 5 rail widths.

**Finding, plainly negative:** slot-TM-210 sits only **2.05 degrees** from
strip-TE-400 -- essentially redundant, the smallest cross-class angle found
in the whole dataset (smaller than some within-class angles). Slot TM's
PCA centroid also sits unusually close to strip-TE's (centroid distance
0.0044, ratio only 5.2x vs. 6-15x for other cross-class pairs) --
**adding a new geometry/polarization combination does not automatically
mean new information; here it happens to land almost exactly where an
existing, already-tested structure already sits.** This is worth
remembering as a caution against assuming "different structure = different
direction" as a blanket rule -- Session 3 already showed strip-TM was
similarly less useful than its "looks different" mutual angle suggested;
slot-TM extends that caution to a second case.

With slot TM included (20 structures), PC3 variance = 0.0109% (up
marginally from Session 3's 0.009% with 15 structures, still essentially
zero) -- the thin-manifold finding holds with more data, not just the
original small sample.

### Task 4: searching for the PC3 gap -- rib waveguide and SWG
Two candidates, `src/rib_waveguide.py` and `src/swg_waveguide.py`, both
reusing `solve_mode_from_eps`/TEfrac classification unchanged.

**Rib waveguide:** full-height rib (width = fab parameter) on a fixed
90nm partial-etch slab, same trapezoid-taper sidewall-angle convention as
strip/slot. *Validated* via an exact limiting-case check: as slab
height -> 0, a rib IS a strip, and the solver reproduces strip's neff
exactly (diff = 0.0000 @ 500nm). Grid convergence: 0.49% spread (better
than strip's own).

**SWG (subwavelength grating) waveguide, EFFECTIVE-MEDIUM APPROXIMATION
ONLY -- stated up front, not buried:** EMpy solves 2D cross-sections, not
z-periodic Bloch modes, so a true SWG needs 3D simulation this toolchain
cannot do. What's implemented is the standard 0th-order Rytov/effective-
medium approximation: a homogeneous core of index
`n_eff = sqrt(f*n_Si^2 + (1-f)*n_clad^2)` (f=duty cycle=0.5, period=250nm,
comfortably sub-Bragg) replaces the periodic segments, reusing
`epsfunc_factory`'s `n_core` parameter (generalized this session,
backward-compatible default=solid silicon, regression-checked against
the exact TE0/TM0 self-test before anything was built on top). *Validated*
the EMA implementation (not the physical accuracy of the approximation
itself, which is a stated limitation) via an exact duty_cycle->1 check
(reduces to solid strip, diff=0.0000 @ 500nm).

**Quantitative test of "does either reach the PC3 gap"** (projected each
candidate's unit sensitivity vector onto the PC3 direction identified from
the 20-structure set, i.e. the gap BEFORE either candidate was added):

| structure | \|projection onto PC3-gap direction\| |
|---|---|
| strip-TE-500 / strip-TM-500 (existing) | 0.027 / 0.031 |
| slot-TE-180 / slot-TM-180 (existing) | 0.018 / 0.025 |
| **rib-TE-400 -> rib-TE-600** | **0.154 -> 0.087** (3-6x larger than anything existing) |
| swg-TE-400 -> swg-TE-600 | 0.026 -> 0.031 (indistinguishable from strip) |

**Confirmed with an isolated Fisher-contribution test** (same method as
Session 3's slot test, strip-only 30-observable baseline):

| addition | n_obs | lambda_min improvement | cond improvement |
|---|---|---|---|
| + slot (TE+TM, this session's full set) | 60 | 90.8x | 13.4x |
| **+ rib (TE only, 5 widths)** | 45 | **571.7x** | **510.5x** |
| + swg (TE only, 5 widths) | 45 | 1.1x | 1.0x |

**Ruled out a magnitude artifact before trusting this:** rib's raw
dneff_dw/dt/da at 500nm (5.6e-4, 3.1e-3, -9.3e-4) are actually SMALLER in
magnitude than strip's (1.0e-3, 4.1e-3, -3.9e-3) at the same width -- the
571x improvement is a genuine DIRECTION effect (a different dw:dt:da
ratio, specifically closer to the previously-unreached combination), not
inflated sensitivity magnitude.

**Answering Task 4 directly: rib waveguide is a real, substantial hit --
by far the best candidate tried against the PC3 gap, outperforming even
slot on a per-observable basis (571x from just 15 new rows vs slot's 91x
from 60 mixed TE+TM rows). SWG (at least in this effective-medium
approximation) is not -- it contributes essentially nothing beyond what
strip/slot already provide**, consistent with it being, in this
approximation, just a uniformly-rescaled-index version of the same solid
rectangular-core topology rather than a structurally new field
distribution. Even with rib added, PC1+PC2 still capture 99.9% of variance
(PC3 up from ~0.01% to ~0.1%, a genuine ~10x increase but still far from
"gap filled") -- rib moves meaningfully toward the gap without fully
closing it.

Full 30-structure atlas and PCA figures regenerated
(`figures/fig1_sensitivity_atlas.png`, `figures/fig2_information_map.png`;
color palette extended to 6 families, a rendering bug from the first pass
this session where 3 new families all rendered as indistinguishable gray).

### Updated overall picture after Session 4
The core Session 3 conclusions survive scrutiny, with important
refinements:
- The dneff-only-vs-dng-aware atlas distinction (Task 1) doesn't matter
  much in practice -- reassuring for Session 3's use of the simpler v1
  atlas, and the wrong-signed mutual-angle correlation is now understood
  more solidly (not an artifact of missing dng).
- The atlas's qualitative shape is robust to the fab-tolerance
  normalization choice (Task 2a) -- another assumption that turned out not
  to be load-bearing.
- The slot noise-model assumption (Task 2b) IS load-bearing, and Session
  3's 86x number should be cited as "86x under an optimistic, unverified
  noise assumption; 2-11x under more realistic pessimistic assumptions" --
  a real, honest downgrade of a headline number, not a reversal.
- Not every new geometry/polarization adds new information -- slot TM
  (Task 3) turned out to sit almost exactly on top of strip TE, a second
  concrete instance (after Session 3's strip TM) of "looks different,
  isn't."
- The PC3 gap (Session 3's central positive finding) is not yet closed,
  but Task 4 found a real, substantial, mechanistically-sound step toward
  it (rib), and a real negative control that clarifies what DOESN'T help
  (SWG's index-only change) -- exactly the kind of contrast that makes the
  positive result more credible, not less.

### Next steps
1. Rib clearly deserves the SAME validation depth slot got in Session 3 --
   push toward literature comparison (rib is an extremely standard PDK
   structure, a published neff reference should be findable) and a finer
   convergence study specifically in the regime where rib's PC3-projection
   is largest (narrower rib widths, <400nm, where the effect was strongest).
2. Sweep rib's slab height (fixed at 90nm this session) -- if slab height
   itself is what's driving the new direction, that's an important
   mechanistic detail worth nailing down before any future inverse-design
   phase treats "rib" as a single fixed recipe.
3. Get a real slot-ring Q/noise estimate (measured or literature-sourced)
   to resolve Task 2b's now-open range (2-86x) into a single trustworthy
   number.
4. Try a genuinely 3D-simulated SWG (or accept it's out of scope for this
   toolchain and drop it as a candidate) -- the effective-medium result
   here is a real negative finding for the EMA-level approximation
   specifically, not necessarily for true SWG physics, and that distinction
   should not be blurred in any future writeup.
5. Given rib's outsized per-observable impact, consider whether a rib-based
   probe (rather than slot) should become the paper's primary
   "complementary structure" example -- slot was the one recommended in
   the original project brief, but this session's evidence now points more
   strongly at rib.

---

## Session 5 — 2026-07-25: rib validation, sourced noise model, angle-specific CRB, slab-height mechanism

**Scope discipline maintained:** still no optimization/inverse-design/ML.
This session closes out the five gaps Session 4 left open before rib can
anchor the paper, per the working style's "validate before headlining"
discipline already applied to slot in Session 3.

### Task 1: rib validation (literature + narrow-width convergence)
`src/rib_validation.py`.

**(a) Literature/EIM cross-check.** A web search found an EXACT-geometry
hit: IntechOpen, "Silicon Waveguides and Temperature-Dependent
Performance" (online-first 1224995), reports for W=500nm/H=220nm/
slab=90nm/1.55um/TE-like -- the project's own nominal point --
neff_strip=2.494, neff_rib=2.602, computed via EIM in MATLAB (no citation
given for these specific numbers; read as the chapter authors' own
calculation, so treated as literature-ADJACENT, not gold-standard, exactly
Session 2's caveat level for its own EIM comparisons).

Rather than compare our FD number to a black-box external number, we
independently re-derived the SAME EIM value from scratch: a 3-layer rib
EIM (two vertical 2-layer slab solves for rib-region/slab-region effective
indices, one lateral 3-layer slab solve), built entirely from the
already-validated `_slab_neff` primitive (no new numerical method). Result:

| | published (MATLAB EIM) | ours (numpy/scipy EIM) | diff |
|---|---|---|---|
| strip @500nm | 2.494 | 2.4946 | +0.02% |
| rib @500nm | 2.602 | 2.6033 | +0.05% |

Two independent EIM implementations agree to ~0.1% on both geometries --
this both validates our 3-layer EIM generalization and lends credibility
to the published number. Then, using the SAME EIM-vs-FD bias direction
Session 1/2 already validated for strip (EIM systematically overestimates
neff for high-index-contrast Si waveguides, decreasing with width):

| W(nm) | rib FD | rib EIM(ours) | %diff |
|---|---|---|---|
| 400 | 2.4596 | 2.5176 | +2.35% |
| 450 | 2.5157 | 2.5644 | +1.94% |
| 500 | 2.5619 | 2.6033 | +1.61% |
| 550 | 2.5998 | 2.6355 | +1.37% |
| 600 | 2.6300 | 2.6622 | +1.22% |

(strip's own Session 1/2 baseline @500nm: +2.05%.) Rib's FD-vs-EIM gap is
the same sign and comparable order of magnitude as strip's already-trusted
gap -- the exact signature Session 2 used to argue the TM solve was
correct, now giving the same confidence to rib's FD solve. Combined with
Session 4's SLAB_H->0 exact-match-to-strip check (unchanged, still holds),
rib now has two independent pieces of validation evidence, matching the
depth slot got in Session 3.

**(b) Narrow-width (<400nm) convergence.** Session 4's PC3-projection was
largest at rib's narrowest tested width (400nm), but Session 4's grid
convergence check only tested the 500nm nominal point. Repeated the same
4-config domain/mesh convergence check at 200-390nm:

| W(nm) | spread | %ofmean |
|---|---|---|
| 200 | 0.0104 | 0.48% |
| 250 | 0.0145 | 0.64% |
| 300 | 0.0133 | 0.57% |
| 350 | 0.0117 | 0.49% |
| 390 | 0.0153 | 0.62% |

Worst case 0.64% vs. the nominal point's already-accepted 0.49% -- not
meaningfully worse. Mode-existence/TEfrac sanity check (guards against a
silent wrong-mode pick, the same failure class as Session 1/2's TE/TM
order-swap concern) confirms a clean, unambiguous TE0 (TEfrac 0.98-0.99)
down to 200nm. **The 571x claim's most-leveraged regime is not resting on
an under-converged grid.**

### Task 2: the load-bearing noise question, sourced explicitly
`src/noise_sourcing.py`. Per the session instructions, the sourced values
are reported here BEFORE their downstream use, not just cited in passing.

**Slot vs strip** -- two independent SAME-PROCESS, SAME-PAPER comparisons
(both avoid comparing loss numbers across unrelated fabrication runs):
- Baehr-Jones et al., "Reduced propagation loss in silicon strip and slot
  waveguides coated by atomic layer deposition," Opt. Express 19, 11529
  (2011): strip 2+/-1 dB/cm, slot 7+/-2 dB/cm (same TiO2-ALD-coated run)
  -> ratio 3.5x.
- A companion ALD-coating slot paper (Opt. Express 19, 26275, 2011): strip
  ~1 dB/cm, slot 3.7-5 dB/cm (100/150/220nm slot widths) -> ratio 3.7-5x.

Both sources agree on direction and land in a tight 3.5-5x band despite
different absolute loss numbers (different process maturity between
papers) -- the RATIO, not the absolute dB/cm, is the transferable
quantity. **Sourced slot noise multiplier: 4x** (central estimate),
applied identically to SIGMA_NEFF and SIGMA_NG (same simplifying
assumption Session 4's stress test used).

**Rib vs strip** -- one same-study comparison found: "Record low-loss
hybrid rib/wire waveguides for silicon photonic circuits" (Ghent/imec-style
hybrid rib-straight + wire-bend design): straight rib sections
0.27+/-0.012 dB/cm vs a 460nm-wide photonic wire (strip) at
1.36+/-0.06 dB/cm -> ratio 0.20x (rib LOWER-loss, opposite direction from
slot). **Caveat stated up front:** this "low-loss rib" literature class is
built to keep the mode away from the fully-etched sidewall, most effective
for WIDE ribs; the source paper's use case (straight low-loss sections
feeding compact wire bends) suggests a wide rib, not this project's narrow
400-600nm single-mode rib. No narrow-rib-specific loss measurement was
found. **Sourced rib noise multiplier reported as a bracket, not a single
point: 0.2x (optimistic, literature ratio) to 1.0x (conservative,
no-benefit-assumed floor; no literature found suggesting narrow rib is
WORSE than strip).**

Conversion from loss ratio to noise multiplier: Q_intrinsic =
2*pi*ng/(lambda*alpha), so sigma_target/sigma_strip = loss_target/loss_strip
(radius-independent, avoids assuming an absolute ring geometry).

**Re-running the isolated Fisher-contribution test (same method as Session
4's slot stress test) under these sourced multipliers:**

| rib noise multiplier | lambda_min improvement |
|---|---|
| same as strip [Session 4 default, unjustified] | 571.7x |
| optimistic (0.2x, sourced) | 4006.4x |
| conservative (1.0x, no benefit assumed) | 571.7x |

**Rib's headline number does NOT collapse under a realistic noise
assumption the way slot's did.** The conservative floor equals Session 4's
original (uncorrected) number exactly, and the sourced-optimistic case is
even larger. This is the opposite direction from slot's 86x -> 2-11x
collapse -- an honest, non-obvious asymmetry between the two structures'
noise sensitivity, reported plainly.

| slot noise multiplier | lambda_min improvement |
|---|---|
| same as strip [Session 4, optimistic/unsourced] | 90.8x |
| 3x worse [Session 4, unsourced round number] | 11.5x |
| **4x worse [THIS SESSION, sourced, 3.5-5x band]** | **7.1x** |
| 5x worse [Session 4, unsourced round number] | 5.1x |
| 10x worse [Session 4, unsourced round number] | 2.1x |

**Slot's Session-4-flagged open 2-86x range is resolved: cite 7.1x going
forward**, not the bracket. (The "same as strip" row recomputes to 90.8x
here, not Session 4's originally-reported 86.0x -- verified this is NOT a
bug: Session 4's stress_tests.py part_b ran before that same session's
Task 3 added slot TM to sensitivity_db.json, so it saw 5 slot records
(TE only); this session's `load_records()` sees all 10 (TE+TM, since
extend_db.py's addition is now permanent in the saved database). More
slot information genuinely available now -> a somewhat larger "same as
strip" baseline number. Confirmed by directly counting slot records
(10, both polarizations) vs. the 5 TE-only that existed when Session 4's
number was produced. Doesn't affect the resolved 7.1x conclusion, which
uses this session's (correct, current) data throughout.)

### Task 3: the angle-specific metric
`src/angle_crb.py`. Every headline number through Session 4 was
lambda_min(F) -- the worst-case direction blended across width, thickness,
AND angle. The paper's actual claim is narrower: sidewall angle
specifically goes from unresolvable to resolvable. Computed the marginal
Cramer-Rao bound on angle alone: (F^-1)_angle,angle, i.e. the variance
after marginalizing out width/thickness (not the much-too-optimistic
1/F_angle,angle).

**Three framings, all reported (not just the most favorable one):**

**(A) Each geometry alone** (5 widths, ring+mzi, TE, marginal angle std in
degrees): strip 0.0315, rib-conservative 0.0157, rib-optimistic 0.0031,
slot-sourced 0.0105.

**(B) Strip baseline + addition** (isolated-contribution framing, matching
Task 2's method): strip-only baseline angle std = 0.0210 deg; + rib
(conservative) -> 0.0009 deg (23.2x improvement); + rib (optimistic) ->
0.0004 deg (51.4x); + slot (sourced) -> 0.0081 deg (2.6x). **Rib gives the
largest angle-specific improvement of anything tested, by a wide margin
over slot** -- consistent with, and now more precisely targeted than, the
aggregate lambda_min story.

**FLAGGED (not glossed over):** both (A) and (B) show strip ALONE already
achieving sub-0.1deg angle precision -- comfortably below the 1deg
fab-tolerance bar used throughout this project's atlas normalization. This
looks inconsistent with Session 1's "badly conditioned... nearly
degenerate" framing for angle. Investigated with:

**(C) Minimal, Session-1/2-realistic probe budget** (exactly 2 TE rings --
Session 2's own best-conditioned pair, TE@550nm+TE@600nm, cond=2.41e5 --
NOT the 15-row all-5-widths-plus-MZI set used in (A)/(B)): angle std =
0.188 deg. Still technically under 1deg, but the mechanism is now visible
by decomposing the full marginal covariance:

```
full marginal std (width[nm], thickness[nm], angle[deg]) = [0.373, 0.141, 0.188]
worst-eigenvalue eigenvector (w, t, a) = [-0.923, -0.121, -0.365]  (eigval=6.33, matches Session 2's own reported 6.33 exactly)
```

**Mechanism, confirmed not asserted:** the worst (lambda_min) direction is
dominated by WIDTH (|-0.92|), not angle (|-0.36|). This is exactly why
aggregate lambda_min looked catastrophic (cond=2.4e5) while the
angle-MARGINAL CRB alone is comparatively well-behaved -- lambda_min was
never purely an angle problem, it was substantially a width-degeneracy
problem, and blending them into one scalar (as every headline number
through Session 4 did) obscures which parameter is actually hard to pin
down. This is the concrete mechanistic justification for why Task 3 asked
for the angle-specific metric rather than aggregate lambda_min.

Adding one rib probe (TE@400nm) to the minimal 2-ring strip base still
gives the largest angle-specific improvement (11.0-14.3x depending on the
noise bracket) vs. slot's 1.9x and a same-class strip control's 2.0x --
**rib's angle-specific advantage survives at the realistic minimal-budget
scale, not just the large-budget scale.**

**SCOPE CAVEAT, stated explicitly:** every CRB number here is an idealized
bound under uncorrelated Gaussian noise at the project's ASSUMED (not
independently measured) SIGMA_NEFF/SIGMA_NG -- a best-case
information-theoretic bound, not a real nonlinear-extraction result. The
literature-cited compact-modeling paper (project_brief.md) found a REAL
extraction procedure's fit residual stayed flat across 85-90deg sidewall
angle under real (correlated, non-Gaussian, model-uncertain) noise -- a
different, harder problem than this linearized CRB captures. **The
absolute sub-degree CRB numbers should NOT be read as "the angle problem
is already solved without rib"; the RELATIVE ranking (rib >> slot for
angle-specific information, at both minimal and large probe budgets, under
every noise assumption tested) is the robust, actionable finding.**

### Task 4: does slab height drive the effect, or is rib robust to it?
`src/rib_slab_sweep.py`. Session 4 fixed SLAB_H=90nm by convention and
never swept it. Swept slab height 50-150nm (realistic partial-etch range,
etch depth 70-170nm out of the fixed 220nm total) at widths 300-600nm,
re-deriving the PC3-gap direction from the strip+slot-only (pre-rib)
20-structure set exactly as Session 4 did, then checking both the PC3
projection and the isolated Fisher lambda_min improvement.

**Consistency check:** our re-derived projection at Session 4's exact
point (slab_h=90nm, W=400nm) = 0.1540, matching Session 4's reported 0.154
almost exactly -- confirms this session's re-implementation is computing
the same quantity Session 4 did.

| slab_h(nm) | proj @300nm | proj @400nm | proj @500nm | proj @600nm | lambda_min improvement |
|---|---|---|---|---|---|
| 50 | 0.184 | 0.102 | 0.077 | 0.051 | 788.9x |
| 70 | 0.243 | 0.135 | 0.100 | 0.073 | **961.4x (best tested)** |
| 90 [Session 4's value] | 0.284 | 0.154 | 0.111 | 0.087 | 834.3x |
| 110 | **0.300 (best proj)** | 0.157 | 0.111 | 0.093 | 585.3x |
| 130 | 0.288 | 0.143 | 0.100 | 0.089 | 351.6x |
| 150 | 0.243 | 0.115 | 0.080 | 0.077 | 181.5x (worst tested) |

**PLAINLY: slab height is NOT a free parameter the effect is robust to.**
PC3 projection @400nm spans 0.102-0.157 (1.54x) and is non-monotonic,
loosely peaking near 90-110nm slab height. The isolated lambda_min
improvement spans 181.5x-961.4x (5.30x) and is monotonically DECREASING
with slab height across the tested range (best at the shallowest slab
tested, 50-70nm; worst at the deepest, 150nm) -- a genuinely different
optimum than the PC3-projection metric's, another instance (like Session
3/4's cond-vs-lambda_min divergences) of two reasonable metrics disagreeing
on the exact optimum while agreeing on the qualitative story.

**Important reassurance:** the benefit never disappears across the tested
range -- even the worst slab height (181.5x) is dramatically larger than
slot's realistic sourced number (7.1x, Task 2), so Session 4's decision to
lead with rib is not fragile to this specific choice. But Session 4's
90nm value, chosen by convention (a common SOI PDK slab thickness) rather
than optimized, is neither the best-PC3-projection point (~90-110nm, close
by luck) nor the best-lambda_min point (~70nm) -- **rib should be treated
as a slab-height FAMILY with its own internal design freedom, not a single
fixed recipe, and any future inverse-design phase must expose slab height
as a free/optimized parameter rather than fixing it at 90nm by
convention.**

### Task 5: bookkeeping
`notes/project_brief.md` updated with a new section recording rib (not
slot) as the primary complementary-structure example, with the Session 4
evidence for the switch and this session's validation/noise/angle-specific/
slab-height follow-up summarized inline.

### Updated overall picture after Session 5
- Rib now has validation depth matching slot's (Session 3): an
  independently-reimplemented literature-adjacent EIM cross-check
  (agreeing to <0.1% with the published number) plus the pre-existing
  SLAB_H->0 exact-match-to-strip check, and confirmed-adequate convergence
  in the narrow-width regime the 571x claim depends on most.
- The noise-model risk that gutted slot's headline number (86x -> 2-11x,
  Session 4) does NOT similarly gut rib's -- sourced literature points the
  other way (rib may be lower-loss than strip), though with a real
  width-mismatch caveat on how directly that literature applies to this
  project's narrow rib. Slot's own open range is now resolved to a single
  citable number: 7.1x.
- The paper's real claim -- angle-specific resolvability, not blended
  3-parameter lambda_min -- now has its own metric, and rib wins on it
  decisively (11-51x angle-specific improvement vs slot's ~2-2.6x),
  including at the realistic minimal-probe-budget scale Session 1
  originally used to motivate the whole project. The mechanism (lambda_min
  being width-dominated, not angle-dominated, in the minimal-budget case)
  is now understood, not just observed.
- Rib is a slab-height family, not a fixed recipe: the effect's magnitude
  swings up to ~5x across a realistic partial-etch range, with the
  lambda_min-optimal and PC3-projection-optimal slab heights disagreeing
  with each other and with Session 4's arbitrary 90nm convention.
- Every number in this session carries the SAME idealized-CRB scope
  caveat as the project's existing Fisher machinery (linearized,
  Gaussian, assumed-not-measured noise floor) -- this is a pre-existing
  project-wide limitation, not new to Session 5, but Task 3's absolute
  sub-degree numbers make it more visible and worth stating explicitly
  going forward whenever an absolute (not relative) CRB number is quoted.

**Diagnostic phase is now essentially complete.** Rib's validation, its
noise-realistic Fisher number, its angle-specific resolvability claim, and
its slab-height mechanism are all now pinned down with the same rigor slot
received in Session 3. The next session should be ready to scope the
E-optimality-based (or angle-marginal-CRB-based, per Task 3's finding that
this may be the more on-point objective than blended lambda_min)
inverse-design objective -- the first actual optimization/ML work in this
project.

### Next steps
1. Scope the inverse-design objective: E-optimality (lambda_min) vs. the
   angle-marginal-CRB metric introduced this session -- Task 3's finding
   that lambda_min's worst direction can be width-dominated rather than
   angle-dominated suggests the angle-marginal metric may be the more
   correct objective to optimize directly, not just report alongside.
2. If pursuing rib further, slab height should be a free design variable
   in any optimizer, not fixed at 90nm -- Task 4 found the lambda_min- and
   PC3-projection-optimal slab heights differ from each other and from the
   Session 4 convention.
3. The rib-vs-strip noise multiplier is still a bracket (0.2x-1.0x), not a
   single sourced number, because no narrow-rib-specific (<600nm,
   single-mode) propagation-loss measurement was found in this session's
   search -- worth a more targeted literature search, or a first-principles
   scattering-loss estimate (e.g. Payne-Lacey roughness model) scaled by
   the narrow rib's actual modal overlap with its exposed sidewall, if a
   tighter number is needed before publication.
4. The angle-marginal CRB's idealized-noise scope caveat (Task 3) applies
   to every Fisher number in this project, not just this session's new
   ones -- worth a dedicated paragraph in the eventual paper's limitations
   section distinguishing "information-theoretic best case" from "what a
   real nonlinear extraction achieves," directly citing the compact-model
   paper's contrasting real-world finding.
5. slot TM's mode and the SWG EMA-approximation negative result (Session
   4) remain open items from that session's own next-steps list, not
   revisited this session (out of scope for the five tasks assigned).
