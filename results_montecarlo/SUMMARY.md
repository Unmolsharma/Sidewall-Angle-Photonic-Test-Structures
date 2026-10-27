# Strips-only angle information, and the Monte Carlo check (stopped after step 2)

Status: **stopped after step 2, as instructed, because the strips-only baseline is not just fragile
but an artifact.** The Monte Carlo check (step 3) has not been started. Step 1 is partly done (see
the end).

All numbers below come from `data/strips_info.json` (produced by `src/strips_angle_info.py`; log
`logs/strips_info.log`); the per-draw values are in `data/strips_info_perturbation.csv`. Angles are
1-sigma Cramér-Rao bounds on the strip (full-etch) sidewall angle, in degrees, at step h, with
conservative rib noise. Noise values are the usual placeholders (σ_neff 6.5e-6, σ_ng 1e-3, σ_phase
1e-3 rad; none sourced).

## Verdict

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

## 1. Residuals are numerical noise (`R_minus_exact` in `data/strips_info.json`)

R − 0.99990 for the strips at h/2, h, 2h:
- TE: −2.7e-4 to +1.9e-3.
- TM: up to 8.5e-3 (TM 600 n_g at h), and 2.9e-2 for TM 600 n_g at 2h.

The sign and size change erratically with step size (for example strip TE 500 n_eff: −1.1e-4,
+2.4e-4, +6.2e-5), which is the signature of solve-to-solve mesh and solver noise. Only TM 600 n_g
grows steadily with the step. For comparison, the r12-vs-r8 mesh check agreed only to 0.3%, so these
residuals sit well below the derivatives' own accuracy.

## 2. Imposing the exact symmetry (`exact_symmetry` vs `nominal`)

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

## 3. Where the strips-only number comes from (`decomposition`, `strips_only_subsets`)

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

## 4. Random perturbation of the strip derivatives (`perturbation`; 2,000 draws per level)

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

## What this means for earlier results

These earlier numbers measure numerical residuals, not physics:
- strip-alone angle CRBs: 1.33°, 1.39°, 0.61°, 0.612° (strip baseline in `results_femwell`);
- two strip rings: 118°;
- every "x better than strip alone" ratio.

These stand: the Stage A–B rib results, the multi-width rib results, the best one-etch sets, and the
conclusion that a single rib with its own angle is useless. The last one in fact gets stronger: that
case is ∞, not "equal to strip alone". I have not changed any earlier files.

## Step 1 status

- **Strip TM index derivatives (48 solves, ~6 min on AC):** queued. Started at 11:41 but paused,
  because the laptop is on battery (88%). It runs as soon as the charger is connected. The Stage D
  search and the priors table have not been rerun yet.
- **Rib slab/n_Si/n_ox at 50/110/130/150 nm:**
  - Not run: 288 solves × ~65-75 s / 3 workers ≈ **1.7-2 h**, over your ~1 h limit.
  - Cheaper option: two widths (300 and 500 nm) at the 4 slabs, 144 solves, about 55 min. That is
    enough for the "best slab for a two-width rib set" question.
- **Rib TM:** future work, as requested.

## Step 3 (Monte Carlo): not started

Given the above, I recommend running it with the strip angle columns symmetrized, on:
1. the best one-etch multi-width set (Stage D, realistic priors);
2. strip + two rib widths;
3. strip + one rib, which is expected to be unidentifiable. It would show the nonlinear fit
   wandering along the width/angle valley.

A strips-only Monte Carlo would only test the symmetry, which is exact.
