# Project Brief: Inverse-Designed Photonic Process Monitors

## The research question
"Can different classes of photonic process monitors (rings, MZIs, TM-mode
devices, slot waveguides, etc.) provide complementary information that
makes fabrication parameter estimation well-conditioned -- and can we
then INVERSE-DESIGN structures purpose-built to maximize that information,
rather than relying on hand-picked traditional structures?"

## Background / motivation
Silicon photonic chips are extremely sensitive to fabrication variation:
tiny (nm-scale) errors in waveguide width, thickness, and sidewall angle
during etching/lithography shift optical performance significantly.
Foundries and designers characterize this variation using "process control
monitor" (PCM) test structures -- rings, MZIs, waveguides -- placed on the
wafer and measured optically after fabrication (this is well-established:
e.g. Ghent/imec groups extract sub-nm-accurate width/thickness wafer maps
this way). These structures are currently chosen by hand/expert intuition,
NOT by any formal information-theoretic optimization.

## The core idea
Reframe photonic device design: instead of designing structures for
performance or for ROBUSTNESS (insensitivity to fabrication error, which
is the dominant existing research direction), design structures that are
DELIBERATELY MAXIMALLY SENSITIVE to specific fabrication parameters --
i.e. optimal sensors of the fabrication process itself. Use Fisher
information / optimal experiment design theory (well-established in
quantum metrology and sensing, but never applied to a device's OWN
manufacturing geometry as the target parameter) to:
1. Quantify how identifiable fabrication parameters (width, thickness,
   sidewall angle, ...) are from a given set of candidate probe structures.
2. Eventually: inverse-design or select structures/placements that
   maximize that identifiability under a chip-area budget.

## Why this is believed to be novel (after extensive literature review)
Multiple rounds of adversarial lit-review search across robust inverse
design, fabrication-aware optimization, photonic circuit yield prediction,
GNN surrogates for PIC yield, and active-learning calibration ideas found
that ALL of those specific framings are already substantially published
or productized (e.g. Luceda Circuit Analyzer for layout-aware yield MC,
BOSON^-1 for robust inverse design, PRISM for fabrication digital twins).
The one gap that survived repeated adversarial searching: nobody has
inverse-designed or Fisher-information-optimized a structure whose
PURPOSE is sensing its own fabrication process. Fisher-information
optimization of photonic structures is common in quantum sensing
(gyroscopes, gas sensors, interferometers) but always targets an EXTERNAL
physical signal, never the device's own manufacturing geometry.

## Progress so far (see progress.md for detailed log)
Built and validated (against published neff~2.44-2.45 @ 500x220nm SOI
waveguide) a full-vector finite-difference mode solver (via EMpy, pure
Python, Fallahkhair-Li algorithm -- NOT Meep, avoid Meep/MPB, it requires
heavy compilation that wastes time/disk in this environment). Computed
sensitivities d(neff)/d(width,thickness,sidewall_angle) and d(ng)/d(...)
for TE-mode ring resonators across widths 400-600nm.

KEY FINDING SO FAR: two same-type TE ring probes are badly conditioned
(Fisher-info condition number ~1e5-1e9) for jointly resolving THREE
fabrication parameters (width, thickness, sidewall angle) -- sidewall
angle sensitivity barely changes across ring widths, so it's nearly
degenerate with width/thickness combinations when using rings alone.
BUT: condition number still swings ~30,000x depending on which two
ring widths are chosen as the probe pair -- so probe design clearly
matters even within a fundamentally under-constrained setup.

## Immediate open question (what we're testing next)
Does adding a DIFFERENT structure class (TM-mode devices, MZIs, slot
waveguides) provide a sensitivity direction that's not degenerate with
what rings alone give you -- i.e. does structure-class DIVERSITY fix the
conditioning problem? This is a pure diagnostic step before any
optimization/inverse-design code gets written. If yes, that becomes the
central empirical motivation for the whole paper. If no, we need to
reconsider what actually breaks the degeneracy.

## Constraints / working style
- I'm an undergrad (Vanderbilt, ECE) working with access to a lab GPU
  server (not needed yet -- everything so far is CPU/scipy-scale) and
  possibly a Lumerical license later if simulation needs exceed what
  EMpy/analytic methods can do.
- Prefer lightweight, pip-installable tools (EMpy, numpy, scipy) over
  anything requiring compilation (Meep/MPB caused a prior session to
  waste an hour "thinking" doing nothing but a failed/slow build).
- I want honest, validated numbers -- if a result looks numerically
  suspicious (e.g. a derivative that's suspiciously exactly zero, or a
  quantity identical across parameter values that should differ), STOP
  and flag it rather than presenting it as a real finding. This has
  happened once already (grid quantization killed the thickness
  derivative until sub-pixel/analytic overlap averaging fixed it) --
  stay alert for the same class of bug elsewhere.
- Target venues: Optics Express / ACS Photonics / JOSA B for the full
  paper, ML4PS or a photonics-EDA workshop for a faster first placement.
- Update notes/progress.md after each session with what was done, what
  was found (including null/negative results), and what's next.

## External validation of the core motivation (2026-07-13)

Literature search confirms the PC3 gap (Session 3) isn't just a simulation
artifact -- it matches a directly published failure mode:

- Sidewall angle is a real, spatially-varying fabrication error caused by
  etch loading effects (etch rate depends on local pattern density), so it
  varies across a single die/wafer, not just chip-to-chip.
- It's performance-limiting: at least one published photonic-crystal-cavity
  study attributed underperforming measured Q-factor directly to nonvertical
  sidewall angle, confirmed via SEM.
- DIRECT HIT: a 2023 process-variation-aware compact-modeling paper found
  that their standard optical parameter-extraction fitting procedure's
  error-function minimum stayed roughly constant across sidewall angles
  85-90 degrees -- i.e. published, independent confirmation that a standard
  optical extraction method cannot resolve sidewall angle from its
  measurements, matching what our own Fisher/PCA analysis found.
- Current industry workaround is NOT a smarter optical technique -- it's
  dropping optical extraction and using destructive SEM cross-sectioning
  instead, which is slow/expensive and can't be done routinely on every
  wafer.

REFRAMED TARGET: the PC3 gap search is no longer just "fill an empty
direction in an abstract 3D space" -- it's specifically "find a structure
enabling non-destructive, all-optical sidewall-angle characterization,
as a real alternative to destructive SEM." This is the paper's
introduction-level motivating claim going forward; cite the compact-model
degeneracy finding directly when writing it up.

## Primary complementary-structure example updated: rib, not slot (2026-07-13/14)

As of Session 4/5, the paper's primary "complementary structure" example
is **rib waveguide**, not slot waveguide (slot was the structure named in
the original project framing above and in Session 3's headline result).
This supersedes any earlier text in this document that treats slot as the
flagship complementary-structure example -- slot remains a genuinely
useful, validated secondary example, but rib is now the stronger one to
lead with.

**Why the switch (Session 4 evidence):** searching for a structure that
reaches the PC3 gap (Session 3's identified unoccupied fabrication-
sensitivity direction), rib was tested alongside a subwavelength-grating
(SWG, effective-medium approximation) candidate. Rib's projection onto the
PC3-gap direction (0.087-0.154 across widths 400-600nm) was 3-6x larger
than anything previously simulated, including slot (0.018-0.031); SWG
under the EMA approximation showed essentially no gap-filling behavior
(0.026-0.031, indistinguishable from strip) and was ruled out as a
negative control. In an isolated Fisher-contribution test (strip-only
30-observable baseline, same method used to establish slot's original 86x
number), adding rib (TE only, 5 widths, 15 new rows) gave a 571.7x
lambda_min improvement, vs slot's 90.8x from 60 mixed TE+TM rows -- i.e.
rib's per-observable impact is far larger than slot's. A magnitude-artifact
check ruled out "just bigger sensitivities": rib's raw dneff/d(w,t,a) at
500nm are actually SMALLER in magnitude than strip's at the same width --
the improvement is a genuine DIRECTION effect (an unreached dw:dt:da
ratio), not inflated sensitivity.

**Session 5 follow-up (validation, noise-realism, and the actual
angle-specific claim):**
- Literature/EIM cross-check: an independently re-derived 3-layer rib EIM
  (built from the same validated slab-dispersion primitive already used
  for strip/TM) lands within 0.05% of a same-geometry literature EIM value
  (IntechOpen, W=500/H=220/slab=90nm, neff_rib=2.602), and the FD-vs-EIM
  gap for rib (1.2-2.4% across 400-600nm) is the same sign/magnitude as
  strip's already-validated FD-vs-EIM gap -- same validation signature
  Session 2 used to trust the TM solve. Narrow-width (<400nm) grid
  convergence, the regime the 571x claim leans on hardest, is comparable
  to the nominal-width convergence already checked in Session 4 (<=0.64%
  vs 0.49%) -- not a source of concern.
- The 571x headline uses strip's noise model unchanged for rib rings, the
  same unvalidated-carryover issue Session 4 flagged for slot. Sourcing a
  same-process rib-vs-strip loss comparison from literature (a hybrid
  rib/wire low-loss waveguide paper: rib 0.27dB/cm vs strip 1.36dB/cm,
  though for a WIDER rib than this project's 400-600nm geometry -- flagged
  as an optimistic upper bound, not a narrow-rib-specific measurement)
  gives a noise-multiplier bracket of 0.2x (optimistic) to 1.0x
  (conservative, no-benefit-assumed floor). Under that bracket, rib's
  isolated lambda_min improvement is 571.7x (conservative) to ~4000x
  (optimistic) -- i.e. rib's headline number does not collapse the way
  slot's did (86x -> 2-11x) under a realistic noise assumption; if
  anything it may be understated. Slot's own open 2-86x range (Session 4)
  is resolved to ~7x using a literature-sourced 4x noise multiplier (two
  same-process strip-vs-slot loss papers, ratio 3.5-5x, central estimate
  4x).
- The paper's actual claim is about SIDEWALL ANGLE specifically, not the
  blended 3-parameter lambda_min. The angle-marginal Cramer-Rao bound
  (mechanistically confirmed to be a materially different, more targeted
  metric than lambda_min -- lambda_min's worst-eigenvalue direction for a
  minimal 2-strip-ring probe set is dominated by WIDTH, not angle) shows
  rib gives the largest angle-specific improvement of any structure tested
  (11-51x depending on probe budget/noise assumption), vs slot's more
  modest ~2-2.6x. IMPORTANT CAVEAT carried into any future writeup: these
  are idealized Fisher/CRB numbers under the project's assumed (not
  independently measured) noise floor -- a best-case information bound,
  not a real nonlinear-extraction result like the compact-model paper's
  85-90deg-flat finding cited above. The RELATIVE ranking (rib > slot for
  angle-specific information) is the robust claim; the absolute sub-degree
  CRB values should not be read as "the problem is already solved."
- Rib's slab height (fixed at 90nm by convention in Session 4) genuinely
  matters: swept 50-150nm (realistic partial-etch range), both the PC3-gap
  projection (1.5x spread, non-monotonic, loosely peaking near 90-110nm)
  and the isolated lambda_min improvement (5.3x spread, peaking near
  70nm, monotonically decreasing toward deeper slabs) vary substantially
  across the range. Rib should be treated as a SLAB-HEIGHT FAMILY with its
  own internal design freedom, not a single fixed recipe -- relevant to
  how any future inverse-design phase parameterizes it. The benefit never
  disappears across the tested range (worst case still 181.5x, comfortably
  above slot's realistic ~7x), so the switch to rib as the flagship example
  is not fragile to this choice, but the specific 90nm value should not be
  assumed optimal without further sweeping/optimization.

See notes/progress.md Session 4 and Session 5 entries for full numeric
detail, methodology, and honestly-reported negative/uncertain results.
