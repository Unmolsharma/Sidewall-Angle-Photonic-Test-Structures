# femwell derivative run: summary and decision log

Unattended overnight run, started 2026-10-25 evening. Each decision is logged with its reason.
Numbers are filled in only from files in this repo; nothing is estimated.

## Decision log

1. **Engine.** femwell (conforming finite-element mesh, order-2 elements) is the main derivative
   engine, per instructions. EMpy grid-shift-averaged derivatives are used only as a cross-check of
   n_eff derivatives; boundary-perturbation results only as support for the symmetry argument.
2. **Separate Python environment.** femwell runs in `C:/Users/unmol/.venvs/femwell` (femwell 0.1.12,
   numpy 2.5.4, scipy 1.18.1, scikit-fem 12.0.2, gmsh 4.15.2). Installing it into the main Python
   could have upgraded numpy/scipy under EMpy and changed every FD result.
3. **Mode search.** Two modes near a guess (the EMpy cell-exact n_eff of the same structure); the
   highest-n_eff mode of the requested polarization is kept (te_fraction > 0.5 = TE). If none is
   found, or it is more than 3% from the guess, the solve is repeated with 6 modes and no guess. Any
   fallback is logged. Reason: one shift-invert mode is 3-4x faster than 4 modes, and the guess only
   steers the eigensolver, it is not used in any result.
4. **Resumability.** Every geometry solve (n_eff at 1.55 um and 1.55 +- 0.005 um) is appended to
   `cache/solves_<mesh>.jsonl` the moment it finishes; a restart skips cached solves.
5. **Battery rule.** Before every new solve the runner checks AC power (Windows
   GetSystemPowerStatus). On battery it submits nothing new, lets running solves finish, and
   re-checks every 5 minutes. All pauses are in `logs/*.log`. At 2026-10-25 ~23:30 the laptop was on
   battery (66%).
6. **Memory rule.** Worker count = min(4, (free RAM - 4 GB) / per-worker estimate), with ~2 GB per
   worker for r12 and ~3.5 GB for r8 (measured earlier: r8 rib workers used 2-3.5 GB). No new solve
   is submitted while free RAM is below one worker's estimate + 2 GB. Reason: earlier, 4 femwell r8
   workers plus 4 EMpy workers ran ~10x slower than normal, consistent with memory pressure.

## Geometry conventions (same as the EMpy code)

- lambda = 1.55 um; n_Si = 3.48, n_ox = 1.44, constant (no material dispersion); air = 1.0.
- Nominal thickness 220 nm. Width = bottom width of the etched wall.
- Sidewall angle a = deviation from vertical in degrees, top narrower. The tilt pivots about the
  **bottom of the tilted wall**: y = 0 for strip, slot and SWG; y = slab top for rib. The wall moves
  inward by tan(a)(y - y_pivot).
- Rib: slab (default 90 nm) across the whole domain, rib from slab top to t; the thickness derivative
  moves the rib top with the slab held fixed. Slab sweep: 50-150 nm, widths 300-600 nm.
- Slot: rails of width w, gap = max(120 nm - (w - 180 nm), 20 nm), rail centres fixed at
  +-(gap/2 + w/2) = +-150 nm; each rail tilts about its own bottom edge.
- SWG: zeroth-order effective medium: strip with core index sqrt(0.5 n_Si^2 + 0.5 n_ox^2) = 2.6631 (computed).
- Oxide below y = 0; top cladding oxide (air only for the cladding check).
- Domain x in [-1.4, 1.4] um (slot [-1.6, 1.6]), y in [-1.2, 1.2] um, as in EMpy.
- Steps h = 5 nm (w), 2 nm (t), 1 deg (a); central differences at h/2, h, 2h; every perturbed
  geometry remeshed. n_g = n_eff - lambda dn_eff/dlambda with dlambda = 5 nm on the same mesh.

## Steps 1-4: method comparison, mesh check, cladding cases, production mesh

7. **Step 1 done.** The derivative-method comparison finished (laptop on battery for most of it; it
   was a short job, not throttled on purpose) and is written up in
   `results_anisotropic/METHOD_COMPARISON.md`.
8. **Step 2 done (r8 mesh check, alone, 2 workers).** Started 23:09 on AC power; free RAM 12.2 GB, so
   the memory rule gave 2 workers at ~3.5 GB each. 51 solves in 10.5 min, none failed, no fallbacks.
   Result (`data/meshcheck_r12_vs_r8.json`, `logs/meshcheck_r12_vs_r8.log`): every derivative agrees
   within 0.30% between r12 and r8 (largest: strip dn_g/dt, -0.30%; rib dn_g/da -0.22%, near zero).
   n_eff differs by at most 2.2e-5. r8 passes every step check. Ratios r12 / r8:
   strip -3.8406 / -3.8397 nm/deg (R 1.0001 / 0.9999), rib -1.7926 / -1.7928 (R 0.7900 / 0.7901),
   air-clad strip -3.7570 / -3.7565 (R 0.9783 / 0.9782).
9. **Step 3 done.** The air-clad strip TE and rib TE cases were run in femwell at r12 and r8 (above).
   The EMpy weighting profiles had given R 0.86-0.91 (air) and 0.63-0.72 (rib), drifting with
   refinement; femwell gives 0.978 and 0.790, stable to 1e-4. The EMpy profile values are not used.
10. **Step 4 decision: r12 is production.** r12 and r8 agree within 0.30% (< 2%) on all n_eff and
   n_g derivatives, so the rule says r12. (For near-zero entries the same convention as the step check
   was used: difference within 2% of the largest derivative of that output. Only rib dn_g/dw and
   dn_g/da are near zero, and they agree within 0.22% relative anyway.)

## Step 5: full femwell batch (r12)

11. **Runtime estimate (before launch).** Timing sample (`logs/timing_r12.log`, 4 parallel workers,
   AC power, 23:21): seconds per perturbed geometry, three wavelengths included: strip TE 6, strip TM
   6, slot TE 5, slot TM 5, SWG 5, rib 19 (12,800 elements vs ~4,000). Batch = 54 units (30
   structures + 24 slab-sweep ribs) x 19 geometries = 1,026 geometry solves: 475 non-rib x ~5.5 s +
   551 rib x ~19 s = ~13,100 worker-seconds, **about 55 min with 4 workers**. Well within the 8 h
   budget.
12. **Workers.** Free RAM 12.7 GB at launch, ~2 GB per r12 worker -> 4 workers (memory rule above).
13. **Sanity of the sampled modes.** femwell n_eff vs EMpy cell-exact: strip TE 2.44871 / 2.45402,
   strip TM 1.76824 / 1.77788, slot TE 1.57381 / 1.58483, slot TM 1.57092 / 1.60261, rib 2.57070 /
   2.57161, SWG 1.69958 / 1.70697; te_fraction 0.98 / 0.05 / 0.97 / 0.04 / 0.99 / 0.97. All the
   requested polarization; slot TM differs most (2.0%), as expected for a mode whose dominant field is
   normal to the silicon walls.
14. **Batch stopped for low system memory (00:08, 26 Oct).** The femwell batch was
   stopped by the session that was supervising the run ("system is running low on memory while the session was
   idle"), after 857 of 1,026 geometry solves (its own log showed 7.8-8.3 GB free just before).
   The instruction that comes with that event is not to restart the command automatically, so **I did
   not restart it.** Every finished solve was cached; no work was lost. Complete: all 30 database
   structures; slab sweep 50, 70, 90 nm (all 4 widths). Incomplete: slab 110 nm (600 nm width
   missing), 130 and 150 nm (not started). To finish, run
   `C:/Users/unmol/.venvs/femwell/Scripts/python.exe src/femwell_batch.py full r12` (resumes from the
   cache, ~170 rib solves, about 15 min), then rerun the post-processing in item 16.
   **Resolved 26 Oct: see item 19.**
15. **Partial slab dropped.** Slab 110 nm had 3 of 4 widths. Its gain would come from fewer rib
   widths than the other slabs, so the whole slab is left out rather than compared with an unlike set
   (`data/skipped.json`). **Superseded by item 19: all six slabs are now complete and nothing is
   skipped.**
16. **Post-processing run on the complete units** (light, < 1 GB): `src/femwell_postprocess.py
   database r12`, `src/downstream_analysis.py --data ../results_femwell/data --slab
   rib_slab_sweep_db_femwell.json --fig ../figures/femwell --label femwell`, `src/femwell_impact.py`,
   `src/femwell_comparison_table.py`. Since 26 Oct also `src/femwell_crosscheck_ratios.py` (writes
   `data/crosscheck_and_ratios.json`; until then this was an ad hoc command, so its slab entries
   stopped at 90 nm) and `src/femwell_script_consistency.py` (item 20).
17. **Mode-search fallbacks.** 23 solves used the 6-mode fallback (22 slot TE at 165/195/210 nm, 1
   slot TM 180): the EMpy guess was more than 3% from the femwell value. All returned TE-like modes
   where TE was requested (te_fraction 0.97-0.98) and every slot derivative passes the step check, so
   the mode choice is consistent within each structure.
18. **One original script failed on the femwell data (fixed 26 Oct).** `fisher_reanalysis.py`
   crashed (`min() iterable argument is empty`) in section (c), TE ring + MZI: with femwell
   derivatives **all 25** ring + MZI combinations have condition number >= 1e12, the script's own
   "exactly rank-deficient" threshold, so its ranked list was empty. (The overnight report said this
   was the two-TE-ring case; that was wrong. Two TE rings: 8 of 20 combinations are above the
   threshold, the other 12 rank normally, cond 3.0e11-9.7e11.) This is a result, not a data error: a
   ring and an MZI of a symmetric strip see width and angle through the same exact degeneracy
   (R = 1.000). Fix: when every combination is rank-deficient the script now says so and skips the
   ranking; output on the cell-exact database is byte-identical to before. All 9 original scripts now
   run (`logs/original_scripts_on_femwell_db/failed_scripts.json` is empty).

19. **Slab sweep finished (26 Oct, on your instruction).** Resumed from the cache with
   `femwell_batch.py full r12 2 25`: at most 2 workers, and, overriding the battery rule for this
   run only, solving on battery down to 25% charge (it would have paused and logged below that; it
   never got there). The runner also now drops to 1 worker, and logs it, whenever free RAM falls
   below the per-worker estimate + 2 GB; free RAM stayed at 6.1-9.6 GB, so that never triggered.
   9.4 GB free at start. 169 solves (slab 110 nm at 600 nm, all of 130 and 150 nm), 0 failed,
   11:25-15:46 including a laptop sleep; on battery a rib solve took ~105 s per worker, against
   ~19 s on AC. The cache had no truncated lines after the sleep. Windows reported AC = 0 (battery)
   throughout, including after the charger was said to be connected. All 1,026 solves are cached;
   `data/skipped.json` is empty. Database structures are unaffected: every number outside the slab
   sweep is identical to the earlier run.
20. **Cross-script consistency.** `src/femwell_script_consistency.py` reads the same quantity from
   every output that reports it, plus an independent Fisher computation from `probe_db.json`:
   strip-baseline and two-ring angle CRBs, rib/slot CRBs and improvement factors (writeup_numbers vs
   angle_crb.log, printed digits), ring-pair condition numbers (writeup_numbers vs
   fisher_reanalysis.json), R for all 30 structures (robustness_checks vs crosscheck_and_ratios vs
   sensitivity_db), the 6b "true" rib changes, and the slab list. **84 checks, 0 disagree**
   (`data/script_consistency.json`). The same pipeline run on the committed cell-exact database
   reproduces all 330 numbers in `results_corrected/data/writeup_numbers_corrected.json` exactly.

## femwell mesh validation (r12 vs r8)

| Structure | n_eff r12 / r8 | Largest derivative difference | Ratio (dn/da)/(dn/dw) r12 / r8, nm/deg | R r12 / r8 |
|---|---|---|---|---|
| Strip TE 500, oxide | 2.448711 / 2.448706 | dn_g/dt -0.30% | -3.8406 / -3.8397 | 1.0001 / 0.9999 |
| Rib TE 500, oxide | 2.570704 / 2.570682 | dn_g/da -0.22% (near zero) | -1.7926 / -1.7928 | 0.7900 / 0.7901 |
| Strip TE 500, air top | 2.388615 / 2.388605 | dn_g/dw +0.02% | -3.7570 / -3.7565 | 0.9783 / 0.9782 |

Every n_eff and n_g derivative agrees within 0.30%; all r8 step checks pass.
Source: `data/meshcheck_r12_vs_r8.json`, `logs/meshcheck_r12_vs_r8.log`.

## Step-size check (h/2, h, 2h; 2% plus the near-zero rule)

| Structure type | femwell fails / 30 | Cell-exact fails / 30 |
|---|---|---|
| Strip TE | 0 | 6 |
| Strip TM | **1** (600 nm dn_g/da, 2.11%) | 16 |
| Slot TE | 0 | 13 |
| Slot TM | 0 | 12 |
| Rib TE | 0 | 5 |
| SWG TE | 0 | 0 |
| **Total, 30 structures** | **1 / 180** | 52 / 180 (original averaging: 103 / 180) |
| Slab sweep | **1 / 144** (all six slabs; slab 110 nm, 600 nm, dn_g/da, 2.04%, a near-zero entry at 10.5% of the largest) | 14 / 144 (all six slabs) |

All 179 database passes are plain 2% passes. In the slab sweep, one pass (slab 150 nm, 600 nm,
dn_g/da, 0.8% of the largest derivative) needed the near-zero rule.
Source: `data/step_check.csv`, `logs/postprocess_r12.log`, `figures/femwell/fig3_step_check.png`.

## Angle/width ratios (femwell r12)

R = (dn/da) / (-h tan 1 deg dn/dw), h = height of the tilted wall (220 nm; rib 220 nm - slab).
R = 1 means a tilt is exactly a width change. Source: `data/crosscheck_and_ratios.json`.

| Structure | (dn_eff/da)/(dn_eff/dw), nm/deg | R (n_eff), 5 widths | R (n_g), 5 widths |
|---|---|---|---|
| Strip TE, oxide | -3.8397 .. -3.8408 | 0.9999 .. 1.0002 | 0.9999 .. 1.0014 |
| Strip TM, oxide | -3.8397 .. -3.8465 | 0.9999 .. 1.0017 | 1.0012 .. 1.0084 |
| Slot TE, oxide | -3.8412 .. -3.8459 | 1.0003 .. 1.0015 | 0.9998 .. 1.0010 |
| Slot TM, oxide | -3.8366 .. -3.8419 | 0.9991 .. 1.0005 | 0.9999 .. 1.0009 |
| SWG TE, oxide | -3.8399 .. -3.8422 | 1.0000 .. 1.0005 | 0.9988 .. 0.9996 |
| Rib TE, 90 nm slab (400-600 nm) | -1.7607 .. -1.8124 | **0.7759 .. 0.7987** | 0.59 .. 1.74 (both n_g derivatives near zero; not meaningful) |
| Rib TE, 90 nm slab, 300 nm | -1.7020 | 0.7501 | |
| Rib TE, 70 nm slab (300-600 nm) | -2.0529 .. -2.1797 | 0.7841 .. 0.8325 | |
| Rib TE, 50 nm slab (300-600 nm) | -2.4665 .. -2.5911 | 0.8312 .. 0.8732 | |
| Rib TE, 110 nm slab (300-600 nm) | -1.4013 .. -1.4822 | 0.7298 .. 0.7719 | |
| Rib TE, 130 nm slab (300-600 nm) | -1.1214 .. -1.1759 | 0.7138 .. 0.7485 | |
| Rib TE, 150 nm slab (300-600 nm) | -0.8698 .. -0.8969 | 0.7119 .. 0.7340 | |
| **Strip TE 500, air top cladding** | -3.7570 | **0.9783** | 0.9534 |

The femwell ratios confirm the symmetry mechanism to 0.1-0.2% for every top-bottom-symmetric
structure, in both polarizations. With air on top the strip is no longer symmetric, but the effect
is small: R = 0.978, a 2.2% departure, against 20-25% for the rib. So the "tilt = width change"
conclusion depends on the cladding only weakly; the slab is what breaks it. The boundary-perturbation
formula gives the same reading: (dn/da)/(dn/dw) = -2 tan(1 deg)(y_c - y_pivot), with y_c the centroid
of the field weighting along the wall, so R = 1 exactly when that weighting is centred on the wall
(`src/sidewall_weighting.py`; its EMpy-based numbers are only qualitative, see
`results_anisotropic/METHOD_COMPARISON.md`).

## Cross-check of n_eff derivatives against grid-shift EMpy (step 6c, no new EMpy runs)

| Derivative | femwell | Anisotropic N4 / N8 | Cell-exact N4 / N8 |
|---|---|---|---|
| Strip TE 500 dn/dw | +1.5861e-3 | -0.6 / -0.6% | -3.7 / -3.4% |
| Strip TE 500 dn/dt | +3.7238e-3 | -0.2 / +0.2% | -0.1 / -0.0% |
| Strip TE 500 dn/da | -6.0916e-3 | +0.7 / +0.7% | +3.7 / +3.4% |
| Rib TE 500 dn/dw | +8.4426e-4 | +0.9 / +1.5% | -1.9 / -1.7% |
| Rib TE 500 dn/dt | +2.7956e-3 | +0.9 / +1.3% | +1.6 / +1.7% |
| Rib TE 500 dn/da | -1.5134e-3 | -3.6 / -3.8% | -1.0 / -1.6% |

The batch values equal the earlier femwell r12 reference run to < 0.001%.
Source: `logs/crosscheck_and_ratios.log`.

## Headline sidewall-angle precisions (1-sigma Cramér-Rao bound, femwell)

Value at step h; in brackets the min..max over steps h/2, h, 2h. Noise assumptions unchanged
(rib conservative = 1.0x strip noise, optimistic = 0.2x; slot 4x).

| Structure set | Angle precision | Note |
|---|---|---|
| Strip, 5 widths, TE+TM, ring+MZI | 0.612° (0.34 .. 0.612°) | near-singular; strip cannot determine angle |
| + 5 ribs, conservative | **0.000606°** (0.000606 .. 0.000607°) | stable to the step |
| + 5 ribs, optimistic | **0.000269°** (0.000269 .. 0.000270°) | stable |
| + 5 slots, sourced 4x noise | 0.308° (0.108 .. 0.308°) | essentially no help |
| Two strip TE rings (550 + 600 nm) | 118° (90.5 .. 427°) | cannot determine angle at all |
| Two strip rings + one rib (400 nm), conservative | **0.00824°** (0.00824 .. 0.00825°) | stable |
| Two strip rings + one rib (400 nm), optimistic | 0.00782° (0.00782 .. 0.00783°) | stable |
| Rib alone (5 widths, TE ring+MZI), conservative | 0.24° (0.231 .. 0.24°) | |

Improvement ratios over a strip baseline are not headline numbers: the baselines are near-singular
(0.34-0.612° and 90-427° depending on the step), so ratios such as 1,009x (560-1,009x) or 14,352x
(10,988-51,738x) mostly measure how singular the baseline is. All values: `COMPARISON.md`.

## Step 6b: what an ignored sidewall angle does to a rib prediction

Strip TE rings at 400-600 nm fitted with angle held at nominal (noise-weighted least squares, same
sigma as the Fisher analysis), applied to the 500 nm rib (90 nm slab); linear in the angle, femwell
derivatives at step h. Source: `data/impact_sidewall_angle.json`, `logs/impact_sidewall_angle.log`.

| True angle change | Effective width error from strip rings | Rib n_eff change: predicted / true | Rib n_eff error | Rib n_g error (true change) | Ring resonance error near 1550 nm |
|---|---|---|---|---|---|
| 1° | -3.840 nm | -3.242e-3 / -1.513e-3 | **-1.729e-3** | -5.3e-6 (+1.87e-4) | **-711 pm** |
| 2° | -7.681 nm | -6.484e-3 / -3.027e-3 | **-3.457e-3** | -1.06e-5 (+3.74e-4) | **-1,422 pm** |

- The strip rings absorb the tilt completely. Fitting a width change leaves chi^2 = 0.07 (1°) against
  6.4e6 with no fit, so the angle is invisible in strip data. Fitting thickness as well changes
  nothing (fitted dt = 0.0000 nm).
- The rib responds to the tilt only 79% as strongly as to the equivalent width change, so the
  effective width overpredicts the rib's n_eff shift by more than 2x per degree (-3.24e-3 predicted
  vs -1.51e-3 true). For a ring built from that rib this is a 0.71 nm resonance error per degree of
  sidewall angle. The n_g error is small (-5.3e-6 per degree) because the rib's n_g barely responds
  to width or angle.
- Resonance error = lambda x (n_eff error) / n_g(rib), n_g = 3.7698 (femwell).

## Which conclusions changed

Full numbers: `COMPARISON.md` (original | cell-exact | femwell, each with its h/2..2h range).

**Still true, and firmer:**
- Symmetric structures cannot separate sidewall angle from width: R = 0.999-1.002 for strip TE/TM,
  slot TE/TM and SWG (cell-exact 0.995-1.013); rib 0.776-0.799.
- The "empty direction": PC3 of the first 15 structures is 9.5e-6% (cell-exact 0.0011%, original
  0.009%); with the rib, 0.35%. Seven scalings: largest 1.0e-5% (the original data actually reached
  0.014%, not "below 0.01%"). PC3 loadings 0.609 / 0.000 / 0.793 (width / thickness / angle).
- Rib projection onto that direction 0.136-0.255, against <= 0.0008 for every other structure.
- Angle precision with the rib included is stable: 0.000606° (5 ribs, conservative), 0.00824° (two
  strip rings + one rib). Cell-exact gave 0.000651° and 0.00905°.
- Directions cluster cleanly by type; within-type mean angle 7.3° vs cross-type 44.2° (15 structures).
- Same-width ring + MZI is exactly singular (condition number 1.4e16).
- Slab-thickness ranking, all six slabs (smallest-eigenvalue gain of two strip TE rings + the 4-width
  rib set): **70 > 90 > 50 > 110 > 130 > 150 nm** at h, h/2 and 2h, and the same order as cell-exact.

  | Slab (nm) | femwell gain at h | h/2 .. 2h | cell-exact gain at h | R (n_eff), 300-600 nm |
  |---|---|---|---|---|
  | 50 | 1.30e6x | 4.01e5 .. 4.09e5x | 1.18e5x | 0.831 .. 0.873 |
  | 70 | **1.59e6x** | 4.91e5 .. 4.99e5x | 1.47e5x | 0.784 .. 0.832 |
  | 90 | 1.42e6x | 4.36e5 .. 4.45e5x | 1.34e5x | 0.750 .. 0.799 |
  | 110 | 1.02e6x | 3.12e5 .. 3.22e5x | 9.89e4x | 0.730 .. 0.772 |
  | 130 | 6.27e5x | 1.93e5 .. 1.98e5x | 6.23e4x | 0.714 .. 0.749 |
  | 150 | 3.25e5x | 1.00e5 .. 1.02e5x | 3.34e4x | 0.712 .. 0.734 |

  The absolute gains are not meaningful: the denominator is the near-singular two-ring baseline,
  which is why h and h/2..2h differ by ~3x. The order is what holds. Thicker slabs break the
  width/angle symmetry more (R falls from ~0.85 at 50 nm to ~0.72 at 150 nm), yet they help less,
  because the rib's derivatives shrink as the etch gets shallower (500 nm rib: dn_eff/da -3.12e-3
  per degree at 50 nm vs -3.00e-4 at 150 nm; dn_eff/dw 1.21e-3 vs 3.42e-4 per nm). 70 nm balances
  the two.

**Changed vs cell-exact:**
- Two TE rings: condition numbers 3.0e11-4.6e13 (cell-exact 9.7e9-3.1e10); best/worst spread 153x
  (cell-exact 3.1x, original 30,000x). "The choice of ring widths hardly matters" does not hold as
  stated, but every pair is near-singular, so the spread is not meaningful either.
- Additions to the 550+600 pair: the TM ring is still the best (smallest eigenvalue 2.76e-3 vs
  baseline 4.53e-6), but a third TE ring now also helps substantially (500 nm: 5.55e-4; 400 nm:
  9.85e-4; cell-exact ~1.7x). The MZI still adds nothing (4.54e-6). All ranges are wide.
- Spearman correlation: 0.090 (-0.009 .. 0.176) for 15 structures, -0.023 for 30. Effectively zero
  (cell-exact +0.21, original -0.37). "Looking different doesn't predict being useful" holds; neither
  the "wrong sign" nor a positive correlation does.
- Slot: gains 33x / 4.6x / 1.7x (same / 4x / 10x noise), far smaller than cell-exact (474x / 46.5x /
  9.5x) and below the original (90.8x / 7.1x / 2.1x). In the angle CRB the slot set gives 0.308° vs
  strip 0.612°. Slot TM vs nearest strip TE direction: 0.24° (cell-exact 1.6°).
- Two strip rings + a third strip ring: 14.7x (11-16x) vs cell-exact 1.7x. Not meaningful: the two-ring
  baseline is 118° (90-427°).
- Strip baseline angle CRB 0.612° (cell-exact 0.209°, original 0.021°).
- SWG gain 1.5x (1.2-1.6x) vs 1.0x: still a negative control.

**Changed vs the original writeup** (already flagged in the cell-exact comparison and still true with
femwell): the MZI is not the best addition; the 30,000x ring-pair spread; the "wrong sign"
correlation; the slot gains; sub-0.1° strip precision; the 23x / 51x / 11-14x rib improvement factors
(now baseline-limited and not meaningful); the "400 nm" third ring was the 500 nm ring.

## Strip TM

femwell strip TM derivatives pass 29 of 30 step checks (cell-exact 14 of 30); the one failure is
dn_g/da at 600 nm (2.11%). femwell resolves the top and bottom interfaces exactly, so this supports
the earlier explanation that cell-exact's TM problem came from arithmetic averaging of a field normal
to those interfaces. Strip TM R = 0.9999-1.0017: TM is as degenerate as TE.

## Surprising or worth checking

1. The supervising session stopped the batch for low memory although the runner logged 7.8-8.3 GB free; Chrome and
   other apps were also running. Nothing was lost, and the slab sweep was finished later (item 19).
2. With femwell derivatives, every strip TE ring + MZI pair and 8 of 20 strip TE ring pairs are
   rank-deficient by the analysis code's own threshold (>= 1e12). This is what an exact width/angle
   degeneracy looks like numerically (item 18).
3. Rib n_g derivatives with respect to width and angle are both near zero (-4.7e-5 and +1.9e-4 at
   500 nm), so the rib's n_g ratios (0.59-1.74) are meaningless; only n_eff ratios are quoted.
4. The third-strip-ring gain jumped from 1.7x (cell-exact) to 14.7x. This comes from arithmetic on
   near-singular matrices, not a physical improvement.
5. The EMpy guess was more than 3% off femwell for slot TE at 195/210 nm (item 17).
6. The batch ran in ~47 min for 845 solves, far faster than the earlier femwell timings; those were
   taken on battery with EMpy competing for memory.

## Files

- `data/`: sensitivity_db (30 structures, femwell r12), probe_db, rib_slab_sweep_db_femwell (50-150
  nm, all six slabs), step_check.csv, script_consistency.json, fisher_reanalysis.json, writeup_numbers_femwell.json, robustness_checks.json (h/2 and 2h variants, R),
  meshcheck_r12_vs_r8.json, crosscheck_and_ratios.json, impact_sidewall_angle.json, skipped.json, and
  the original scripts' outputs (atlas_*, ...).
- `COMPARISON.md`: every writeup number, original | cell-exact | femwell, with step ranges.
- `logs/`: batch logs (timing, mesh check, full run), post-processing, original scripts.
- `cache/solves_r8.jsonl`, `cache/solves_r12.jsonl`: every femwell solve.
- `figures/femwell/`: Figures 1-3 regenerated from the femwell database.
