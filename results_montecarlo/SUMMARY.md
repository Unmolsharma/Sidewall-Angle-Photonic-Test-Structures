# Strip symmetry, symmetry-imposed search, and the Monte Carlo check

Status: steps 1–3 done. Step 2 (strips-only angle information) is kept unchanged further down.
Headline numbers with sources, and retired numbers with reasons, are in `HEADLINE_NUMBERS.md`.
Noise values are the usual **placeholders** throughout (σ_neff 6.5e-6, σ_ng 1e-3, σ_phase 1e-3 rad,
L = 100 µm; none sourced).

## Strip symmetry is now the default

**From 2026-10-27 every Cramér-Rao bound that includes strips imposes the exact strip symmetry:**
each oxide-clad strip's angle column is replaced by −220 nm·(π/180)·(its width column), i.e.
R = 1 exactly. This is the default of `derivatives()` in `src/etchdepth_femwell_analysis.py` and
`src/priors_search.py` (`strip_symmetry=True`). In the nonlinear Monte Carlo model a strip sees the
angle only through its mid-height width (strips are solved with vertical walls). Results committed
before 2026-10-27 were computed without it. Only `src/strips_angle_info.py` (step 2, which studies
the raw residuals) and `src/step_check_all.py` (which checks the raw solves) switch it off.

## Step 1. Search and priors table rerun with the symmetry imposed

Sources: `sym_search/data/summary_tables.md` (tables S1–S3), `search.json`, `search_table.csv`,
`slab_pairs.json`, `priors_table.csv`, all written by `src/sym_search.py` (log
`logs/sym_search.log`). Strip-angle CRB, deg, at step h [min–max over h/2, h, 2h].

New solves:
- strip TM index solves: 48 solved, 0 failed (`../results_etchdepth_femwell/logs/striptm_r12.log`);
- rib 300 and 500 nm at 50/110/130/150 nm slabs (slab, n_Si, n_ox steps): 144 solved, 0 failed
  (`../results_etchdepth_femwell/logs/ribslab_r12.log`).

Step check over all 372 derivatives (`sym_search/data/step_check_all.csv`, from
`src/step_check_all.py`): 366 PASS, 3 PASS near zero, **3 marginal FAIL**, all n_g:
- strip TM 600, angle: 2.11% (replaced by the symmetry anyway);
- rib 300 nm, slab 110, slab: 2.25%;
- rib 600 nm, slab 110, rib angle: 2.04% (limit 2%).

### Best slab for a two-width rib set

Strip TE at 5 widths + rib 300 + 500 nm (table S1):

| Slab (nm) | C, none | D, realistic priors |
|---|---|---|
| **50** | **0.359 [0.359–0.369]** | **0.367 [0.367–0.377]** |
| 70 | 0.419 [0.419–0.435] | 0.425 [0.425–0.441] |
| 90 | 0.555 [0.555–0.597] | 0.561 [0.561–0.603] |
| 110 | 0.796 [0.743–0.865] | 0.802 [0.749–0.871] |
| 130 | 1.21 [1.06–1.3] | 1.22 [1.06–1.31] |
| 150 | 3.81 [2.51–3.81] | 3.83 [2.52–3.83] |

The shallowest slab is best for two widths, and the bound gets worse steadily as the slab thickens.

### Best one-etch set overall

Table S3, stage D, realistic priors, 8 structures: **0.162° [0.118–0.171]** at the 90 nm slab, from
strip TE 400 + TE 600 + TM 400 + TM 600 + rib 300/400/500/600 nm. Note the −27% step spread at h/2.
Stage C without priors gives 0.144 [0.105–0.153].

Best 8-structure set per slab (D, realistic):

| Slab (nm) | 50 | 70 | 90 | 110 | 130 | 150 |
|---|---|---|---|---|---|---|
| CRB (°) | 0.359 | 0.187 | 0.162 | 0.795 | 1.21 | 3.81 |

Slabs 70 and 90 offer four rib widths, the others only two. With 300 + 500 nm only, slabs 70 and 90
give 0.419 and 0.555, so for the same two widths the 50 nm slab is best.

Best set by size (D, realistic):

| Structures | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|
| CRB (°) | 0.413 (slab 50) | 0.329 | 0.174 | 0.167 | 0.162 (slab 90) |

**Strips only: ∞ for every set size** (no angle information under the symmetry).

### Single rib

Strip TE 5 + rib 500 nm (table S2):
- ∞ at stages C and D without priors.
- With realistic priors, 24–43°: 29.4 (slab 50), 25.1 (70), 24.3 (90), 26.3 (110), 32.1 (130),
  43.1 (150). The priors set the angle, not the data.

### Priors table

From `priors_table.csv`, strip + rib70. Strip alone is ∞ everywhere.

| Stage, prior | + rib 500 nm | + rib 500 + 300 nm |
|---|---|---|
| B, none (shared angle) | 0.334 | 0.00433 |
| C, none | ∞ | 0.419 |
| C, realistic | 25.1 | 0.419 |
| C, ar − a ≤ 1° | 1.24 | 0.404 |
| D, none | ∞ | 0.451 |
| D, realistic | 25.1 | 0.425 |
| D, realistic + ar − a ≤ 2° | 2.41 | 0.421 |
| E, realistic | 25.2 | 0.428 |

## Step 3. Monte Carlo (symmetry imposed)

Code and files:
- `src/montecarlo.py` (subcommands plan, fit, mc) and `src/montecarlo_batch.py` (femwell solves);
- design: `data/mc_design.json`;
- solves: `cache/solves_r12.jsonl`, 448 solved, 0 failed (`logs/mc_solves_r12.log`; 2 h 43 min on
  battery, 2 workers);
- surrogates: `data/surrogates.json`, accuracy in `data/surrogate_accuracy.{json,csv}` (log
  `logs/mc_fit.log`);
- Monte Carlo: `data/mc_results.json`, per-set fits in `data/mc_fits.csv` (log `logs/mc_run.log`).

### Cases

Stage D, realistic priors on d, n_Si, n_ox; truth = nominal.
- (a) best multi-width one-etch set: strip TE 400, TE 600, TM 400, TM 600 + rib 300/400/500/600 nm,
  90 nm slab (7 unknowns);
- (b) strip TE 400–600 (5) + rib 300 + 500 nm, 50 nm slab;
- (c) strip TE 400–600 (5) + rib 500 nm, 50 nm slab.

### Surrogate

Each structure has its own polynomial model of n_eff and n_g in scaled inputs:
- rib inputs: w, t, rib angle, slab, n_Si, n_ox;
- strip inputs: mid-height width, t, n_Si, n_ox.

Each model is fitted to a Latin-hypercube design (20 strip or 30 rib points), plus the existing
derivative solves that fall inside the box.
- **Box:** ±6 marginal CRB σ (floors 1 nm, 0.25°, 1e-4). Case c uses a wide box: w ±12 nm, strip
  angle ±3°, rib angle ±5°.
- **Basis:** full cubic, quadratic or reduced, chosen per structure by leave-one-out error on the
  training points only.

Accuracy on 6 (strip) or 8 (rib) held-out solves, as RMS in units of the noise (σ_phase as an n_eff
noise is 2.47e-6):

| Surrogate | n_eff RMS (σ_phase) | n_eff max | n_g RMS (σ_ng) |
|---|---|---|---|
| strips (7 structures) | 0.04–0.37 | ≤ 1.8e-6 | ≤ 0.004 |
| rib 400/500/600, slab 90 | 1.26–2.50 | ≤ 1.2e-5 | ≤ 0.004 |
| rib 500, slab 50 | 1.26 | 4.7e-6 | 0.013 |
| rib 500, slab 50, wide box | 5.32 | 2.8e-5 | 0.044 |
| **rib 300, slab 90** | **8.06** | 4.5e-5 | 0.035 |
| **rib 300, slab 50** | **30.1** | 1.1e-4 | 0.26 |

**The rib 300 nm misses are femwell noise, not surrogate error.** Three observations point to that:
- the errors are as large near the box centre as at its edges;
- the mode is the same in every solve (TE fraction 0.961–0.962);
- the errors do not follow the element count.

So they are remeshing noise of the weakly confined 300 nm rib, up to ~1e-4 in n_eff, far above
σ_phase. No smooth surrogate can do better.

This does not affect the spread test, because the Monte Carlo generates and fits its data with the
same surrogate. It does mean **femwell's own n_eff for the rib 300 nm on a 50 nm slab is only good to
~30 σ_phase.** That limits its finite-difference derivatives to a few percent, which is consistent
with the step check passing at 2%. How far the smooth model sits from femwell at nominal is the
"fit to femwell nominal" row in the table below.

### Monte Carlo

- 2,000 noisy sets per case: whitened noise N(0,1) per observable, and prior centres drawn from the
  priors.
- Nonlinear least squares (`scipy least_squares`, trust-region, bounds = the box), with 5 starts per
  set (nominal + 4 random).
- A set is *ambiguous* if another start reaches χ² within 1 of the best with a strip angle further
  away than min(CRB, 1°).
- A set is *at bound* if any parameter ends within 0.1% of the box edge.

| | (a) best multi-width | (b) strip + 2 rib widths | (c) strip + 1 rib |
|---|---|---|---|
| CRB, femwell derivatives | 0.162 | 0.367 | 29.4 |
| CRB, surrogate at nominal | 0.146 | 0.372 | 32.0 |
| **MC std of fitted strip angle** [95% bootstrap] | **0.163 [0.155–0.172]** | **0.374 [0.363–0.385]** | 2.47 (box-limited) |
| MC std / CRB (femwell; surrogate) | 1.01; 1.12 | 1.02; 1.01 | – |
| robust std (1.4826·MAD) | 0.143 | 0.364 | 3.75 |
| 5–95% of fitted angle | −0.232 to +0.271 | −0.622 to +0.599 | −3.00 to +2.49 |
| bias ± SE | +0.0097 ± 0.0037 | −0.015 ± 0.008 | −0.17 ± 0.06 |
| ok / at bound / ambiguous / failed | 1989 / 11 / 0 / 0 | 2000 / 0 / 0 / 0 | 820 / 1163 / 11 / 6 |
| mean χ² (dof) | 20.1 (20) | 16.6 (17) | 15.1 (14) |
| fit to noise-free femwell nominal | +0.021 | +0.005 | +0.001 |

**(b) behaves exactly as the CRB says.**
- The spread is 1.01–1.02 × CRB.
- There are no failed, ambiguous or at-bound fits, and the bias is under 0.05 σ.
- All other parameters also match their CRBs, e.g. width 1.43 vs 1.42 nm, rib angle 0.528 vs 0.524°.

**(a) matches the femwell CRB (1.01×) but is 12% above the surrogate's own CRB.**
- Its robust std (0.143) sits at that CRB. The core is Gaussian at the CRB, and a small tail inflates
  the std: 11 at-bound sets, and 1.35% of sets beyond 3 CRB.
- The surrogate CRB is 10% below the femwell CRB because its slopes come from the smooth fit through
  noisy rib 300 nm solves.
- The bias, +0.010°, is 2.7 SE but only 0.06 σ.
- The model error at nominal (+0.021°) is 0.13 σ.

**(c) shows the width/angle valley.**
- The data do not determine the angle (CRB ~30°).
- Fits wander to the box edge: 1163 of 2000 end at a bound, 11 are ambiguous, 6 failed.
- The fitted width and angle lie on the line w − 3.84·a = const (correlation 1.000, slope 3.84
  nm/deg, from `data/mc_fits.csv`), i.e. constant mid-height width.
- Its 2.47° spread is set by the ±3° box, not by information.

**Conclusion:** for the sets that are identifiable, the linearized CRB predicts the real
nonlinear-fit spread well (within about 1–12%). The 0.16° and 0.37° figures stand, at the
placeholder noise.

## Run notes

- **Rib-slab batch, on battery:**
  - It ran with a 25% floor, which is checked only when a new solve is submitted.
  - Two solves took ~3 h each (probably sleep or throttling), so the battery reached 5% before the
    floor check paused new solves (17:19, `../results_etchdepth_femwell/logs/ribslab_r12.log`).
  - A second stall of ~2 h 19 min followed, between 18:37 and 20:56.
  - The results are complete (144/144), but the floor check does not protect against a long
    in-flight solve.
- **Monte Carlo:** the first run was stopped by the system for low memory, with no results written.
  It was rerun with 2 worker processes (`NPROC = 2`) and took ~68 min.
- **Rib TM:** remains future work, as agreed.

## Step 2. Strips-only angle information (unchanged from the previous stop)

All numbers below come from `data/strips_info.json` (produced by `src/strips_angle_info.py`; log
`logs/strips_info.log`); the per-draw values are in `data/strips_info_perturbation.csv`. Angles are
1-sigma Cramér-Rao bounds on the strip (full-etch) sidewall angle, in degrees, at step h, with
conservative rib noise. Noise values are the usual placeholders (σ_neff 6.5e-6, σ_ng 1e-3, σ_phase
1e-3 rad; none sourced).

### Verdict (step 2)

**Strips alone carry no real sidewall-angle information.** For an oxide-clad strip the top and
bottom claddings are the same, so the mode is mirror-symmetric about mid-height. A wall tilt about the
bottom edge splits into two parts:
- a uniform wall shift, which is exactly a width change of −t·tan a;
- a rotation about mid-height, which is odd under that mirror and so changes n_eff and n_g only at
  second order.

At first order, therefore, dn/da = −t·(π/180)·dn/dw for every strip: TE and TM, n_eff and n_g alike.
Equivalently, R = (π/180)/tan 1° = 0.99990 exactly. The finite strips-only numbers reported so far
come entirely from numerical residuals in R:
- strip TE alone: 1.33°;
- best strips-only set: 0.61°;
- the 0.61-1.39° strip-alone baselines in `results_priors_search`;
- the strip-baseline numbers in `results_femwell`.

So every "rib vs strip alone" ratio has the wrong denominator: the true strip-alone precision is ∞,
and the right comparison is rib sets against each other. The rib-based results themselves are real:
the multi-width rib sets are almost unchanged when the exact symmetry is imposed.

### 1. Residuals are numerical noise (`R_minus_exact` in `data/strips_info.json`)

R − 0.99990 for the strips at h/2, h, 2h:
- TE: −2.7e-4 to +1.9e-3.
- TM: up to 8.5e-3 (TM 600 n_g at h), and 2.9e-2 for TM 600 n_g at 2h.

The sign and size change erratically with step size (for example strip TE 500 n_eff: −1.1e-4,
+2.4e-4, +6.2e-5), which is the signature of solve-to-solve mesh and solver noise. Only TM 600 n_g
grows steadily with the step. For comparison, the r12-vs-r8 mesh check agreed only to 0.3%, so these
residuals sit well below the derivatives' own accuracy.

### 2. Imposing the exact symmetry (`exact_symmetry` vs `nominal`)

Here each strip's angle column is replaced by −t·(π/180) × its width column.

| Set | Stage, priors | As computed | Exact symmetry |
|---|---|---|---|
| Strip TE, 5 widths | C, none | 1.33 | **∞** |
| Best strips-only set (8 structures) | C, none | 0.614 | **∞** |
| Strip TE 5 + rib70 500 nm | C, none | 1.33 | **∞** |
| Strip TE 5 + rib70 500 + 300 nm | C, none | 0.400 | 0.419 |
| Strip TE 5 + rib70 500 + 300 nm | E, realistic | 0.414 | 0.428 |
| Best one-etch set (8) | C, none | 0.145 | 0.147 |
| Best one-etch set (8) | D, realistic | 0.167 | 0.169 |

- **With one rib width:** if the rib has its own angle (Stage C), the strip angle cannot be
  determined at all.
- **With two or more widths of the same rib:** the rib widths pin the width and thickness, and the
  strips then turn the "effective width" they see into the angle. Those results change by at most 5%
  under exact symmetry.

### 3. Where the strips-only number comes from (`decomposition`, `strips_only_subsets`)

**Rows contributing to the angle's marginal information,** best strips-only set (Stage C):
- TM MZI phase 61%;
- TE MZI phase 27%;
- TM n_eff 9%;
- TE n_eff 4%;
- n_g: 0%.

**Sub-sets of that set:**

| Sub-set | Strip angle |
|---|---|
| All | 0.614 |
| TE structures only | 2.09 |
| TM structures only | 0.742 |
| n_eff + MZI phase only | 0.614 |
| n_g only | 45.4 |

So it comes almost entirely from n_eff-type rows (the MZI phase is n_eff with a smaller sigma),
mostly TM. That is where the residuals are largest, and the rows with the smallest noise amplify
them most.

**For the real rib-based sets, the information comes from the ribs:**
- Strip TE 5 + rib70 500 + 300 nm (Stage C): rib n_g 91%.
- Best one-etch set (Stage C): rib MZI phase 75%, rib n_eff 11%, rib n_g 9%; strips about 4%.

### 4. Random perturbation of the strip derivatives (`perturbation`; 2,000 draws per level)

Each strip derivative entry is multiplied by (1 + u), with u uniform in ±p. Rib derivatives are
unchanged. Each cell gives the median with the 5-95% range below it.

| Set | Nominal | ±0.1% | ±0.3% | ±1% | ±3% |
|---|---|---|---|---|---|
| Strip TE, 5 widths | 1.33 | 0.279<br>0.149-0.744 | 0.093<br>0.049-0.238 | 0.0275<br>0.015-0.073 | 0.0091<br>0.0049-0.0245 |
| Best strips-only (8) | 0.614 | 0.265<br>0.141-0.512 | 0.094<br>0.048-0.183 | 0.0285<br>0.014-0.057 | 0.0094<br>0.0048-0.019 |
| Strip TE 5 + rib70 500 + 300 (C) | 0.400 | 0.233<br>0.140-0.364 | 0.091<br>0.049-0.209 | 0.0274<br>0.015-0.072 | 0.0091<br>0.0049-0.0244 |
| Strip TE 5 + rib70 500 + 300 (E, realistic) | 0.414 | 0.280<br>0.173-0.411 | 0.120<br>0.063-0.270 | 0.042<br>0.021-0.112 | 0.017<br>0.007-0.052 |
| Best one-etch set (C, 8) | 0.145 | 0.135<br>0.106-0.151 | 0.093<br>0.048-0.145 | 0.035<br>0.015-0.096 | 0.012<br>0.005-0.039 |
| Best one-etch set (D, realistic, 8) | 0.167 | 0.162<br>0.136-0.199 | 0.135<br>0.087-0.269 | 0.070<br>0.033-0.218 | 0.028<br>0.012-0.106 |

**Reading:**
- The perturbations make things look *better*, because random errors break the exact
  proportionality and create apparent angle information in the strips.
- **Strips-only sets:** a 0.1% error already moves the result 2-5x, and at ±1% it is 20-50x. The
  number is set by the size of the derivative error, not by physics.
- **Best 8-structure rib sets:** move 3-7% (median) at ±0.1%. The mesh check suggests the
  derivatives are accurate to about 0.1-0.3%. At ±1% and above they inherit the same fake strip
  information.
- **Strip TE 5 + two rib widths:** moves 30-40% already at ±0.1% (0.400 → 0.233 in C, 0.414 →
  0.280 in E). Its real angle information is small, so even a little fake strip information competes
  with it. Its trustworthy value is the symmetrized one, 0.419° / 0.428°.
- **Practical consequence:** any CRB that includes strips should be computed with the strip angle
  columns symmetrized (Section 2). Otherwise derivative noise can pass for angle information.

### What this means for earlier results

These earlier numbers measure numerical residuals, not physics:
- strip-alone angle CRBs: 1.33°, 1.39°, 0.61°, 0.612° (strip baseline in `results_femwell`);
- two strip rings: 118°;
- every "x better than strip alone" ratio.

These stand: the Stage A–B rib results, the multi-width rib results, the best one-etch sets, and the
conclusion that a single rib with its own angle is useless. The last one in fact gets stronger: that
case is ∞, not "equal to strip alone". I have not changed any earlier files.
