# Etch-depth / nuisance-parameter study (femwell r12): summary

Question: how much of the sidewall-angle precision, and of the rib's advantage, survives when more
of the process is unknown? Unknowns are added in stages; at each stage the sidewall-angle
Cramér-Rao bound (1-sigma, all other unknowns marginalized, no priors) is reported at derivative
step h, with the min..max over h/2, h, 2h in brackets.

## Verdict (plain language)

- **The rib's advantage depends on two assumptions: that the partial etch has the same sidewall
  angle as the full etch, and that its depth is known.** With both (Stage A), adding one rib takes
  the angle from 1.33° (strip only) to 0.0015°, about 900x better, and the result is step-stable.
- **An unknown etch depth (Stage B) costs most of that:** strip + one rib gives 0.32-0.35°, only
  about 4x better than strip alone. With the optimistic rib noise it is 0.07°. The result is still
  step-stable (within ±6%).
- **A separate partial-etch angle (Stage C) removes the advantage completely.** Strip + one rib gives
  1.33°, identical to strip alone; the rib's own angle is 2.5-3.0°. The reason is simple
  counting. One rib width gives two independent observables, n_eff and n_g (the MZI phase is
  proportional to n_eff). It also brings two new unknowns, depth and angle, so it adds no net
  information on the shared parameters.
- **Rescue: a second width of the same rib helps; a TM ring does not.** With a second rib width
  (same etch, so no new unknowns), strip + rib70 at 500 and 300 nm gives 0.40° [0.38-0.43] for the
  strip angle and 0.66° for the rib angle, about 3x better than strip alone. A TM ring changes the
  strip angle only from 1.33° to 1.17°.
- **Index and temperature uncertainty cost little on top of that (Stages D, E).** The strip angle
  goes from 1.33° to 1.39° at h, but the step range widens to 1.16-5.94°, because these sets are
  already close to singular. Temperature alone, with the indices known at a reference temperature
  (E'), changes almost nothing (1.33°). Temperature with free indices (E) is exactly the same as D,
  because the temperature column is a fixed combination of the two index columns.
- **The minimal sets fail.** Two strip rings + one rib ring give 6 observables. From Stage D on
  they have 7 or more unknowns, so the angle is not identifiable at all (∞). Already at Stage C they
  are back to the two-ring baseline of 118°.
- **The worst-determined direction does not change:** about 0.97·w + 0.25·a, i.e. a 1° tilt looks
  like a 3.84 nm width change, the same R = 1 degeneracy as before. From Stage C on the rib angle
  joins that direction (w 0.84, ar90 0.49, a 0.22). The new unknowns (d, n_Si, n_ox) are not what
  is worst determined. They cost precision by removing the information that broke that direction.

So the earlier headline, "a rib breaks the width/angle degeneracy by about 1000x", holds only if
the rib is etched in the same step as the strip and its depth is known. If the partial etch is a
separate step (the usual case), a single rib is worthless for the strip angle. Two rib widths from
that etch give about 3x, not 1000x.

## Stage-by-stage table (rib noise conservative 1.0x; optimistic 0.2x in parentheses, at h)

Sidewall angle of the full (strip) etch, a, in degrees. From Stage C on, the rib's own angle (ar70 /
ar90) is listed on the next line.

| Set | A | A* (depth known) | B (+ depth) | C (+ rib angle) | D (+ n_Si, n_ox) | E (+ T) = D | E' (C + T) |
|---|---|---|---|---|---|---|---|
| Strip TE, 5 widths, ring+MZI | 1.33 [0.82-1.86] | 1.33 [0.82-1.86] | 1.33 [0.82-1.86] | 1.33 [0.82-1.86] | 1.39 [1.16-5.94] | 1.39 [1.16-5.94] | 1.33 [0.82-1.89] |
| Strip TE + TM 500 ring | 1.17 [0.82-1.74] | 1.17 | 1.17 | 1.17 | 1.36 [0.93-2.48] | 1.36 | 1.33 [0.82-1.86] |
| Strip + rib70 (500 nm) | 0.00146 [0.00146-0.00146] (0.00066) | 0.00184 (0.0013) | 0.324 [0.31-0.33] (0.067) | 1.33 [0.82-1.86]; ar70 2.46 [1.53-3.43] | 1.39 [1.16-5.94]; ar70 2.58 [2.17-10.9] | same as D | 1.33; ar70 2.46 |
| Strip + rib90 (500 nm) | 0.00147 [0.00147-0.00147] (0.00067) | 0.00192 (0.0014) | 0.352 [0.33-0.36] (0.073) | 1.33 [0.82-1.86]; ar90 3.02 [1.89-4.20] | 1.39 [1.16-5.94]; ar90 3.17 [2.68-13.3] | same as D | 1.33; ar90 3.02 |
| Strip + both ribs | 0.00112 [0.00112-0.00112] (0.00064) | 0.00162 (0.0013) | 0.242 [0.24-0.24] (0.049) | 1.33; ar70 2.46; ar90 3.02 | 1.39; ar70 2.58; ar90 3.17 | same as D | 1.33; ar70 2.46; ar90 3.02 |
| Two strip rings (550+600) + rib70 ring | 0.0053 [0.0053-0.0053] (0.0038) | 0.0084 (0.0076) | 0.334 [0.33-0.34] (0.067) | 118 [90-427]; ar70 217 | ∞ | ∞ | 788 [477-1710]; ar70 1450 |
| Two strip rings + rib90 ring | 0.0043 [0.0043-0.0043] (0.0023) | 0.0114 (0.0108) | 0.366 [0.37-0.37] (0.074) | 118 [90-427]; ar90 265 | ∞ | ∞ | 788 [477-1710]; ar90 1770 |

"Strip + both ribs" assumes two partial etches (70 and 90 nm slabs need different depths), each
with its own depth and, from Stage C, its own angle.

### Rescues

| Set | B | C | D |
|---|---|---|---|
| Strip + rib70 + TM ring | 0.321 | 1.17; ar70 2.17 | 1.36; ar70 2.51 |
| Strip + rib70 (500) + rib70 300 nm | 0.0043 | **0.40 [0.38-0.43]; ar70 0.66** | **0.43 [0.43-0.46]; ar70 0.71** |
| Strip + rib70 (500) + rib70 400 nm | 0.0116 | 0.60; ar70 0.99 | 0.67; ar70 1.09 |
| Strip + rib70 (500) + rib70 600 nm | 0.0327 | 0.84; ar70 1.36 | 0.92; ar70 1.49 |
| Strip + rib90 (500) + rib90 300 nm | 0.0045 | 0.51 [0.48-0.57]; ar90 1.02 | 0.55; ar90 1.08 |
| Two strip rings + rib70 ring + TM ring | 0.333 | 5.3 [5.3-102] | 60 [53-496] |
| Two strip rings + rib70 ring + rib70 300 nm ring | 0.0205 | 0.43 [0.43-0.45]; ar70 0.71 | 1.36; ar70 2.12 |

The widest width contrast (300 vs 500 nm) rescues the most, and those sets are step-stable to
within about ±10%. With a 600 nm second width the spread is up to ±25% (0.65-1.05° for rib90 in
Stage C). The TM ring does almost nothing once the rib needs
its own angle.

### Worst-determined direction (smallest eigenvalue of F, conservative, at h)

Units: w, t, d in nm; a, ar in degrees; n_Si, n_ox in 1e-3; T in K. The eigenvector is listed by its
largest components.

| Stage | Strip + rib90 | Strip + both ribs |
|---|---|---|
| A | 2.4e4: w −0.97, a −0.22, t +0.07 | 3.8e4: w −0.97, a −0.21, t +0.10 |
| A* | 1.4e4: w 0.97, a 0.22, t −0.07 | 1.8e4: w 0.97, a 0.21, t −0.09 |
| B | 0.48: w 0.94, a 0.25, d90 0.24 | 0.96: w −0.91, d70 −0.25, a −0.24, d90 −0.23 |
| C | 0.027: w −0.84, ar90 −0.49, a −0.22 | 0.023: w 0.78, ar90 0.46, ar70 0.37, a 0.20 |
| D | 0.025: w 0.84, ar90 0.50, a 0.22 | 0.021: same as C |
| E | ~0 (exact): T 0.98, n_Si −0.18 | ~0 (exact): T 0.98, n_Si −0.18 |
| E' | 0.027: same as C | 0.023: same as C |

The strip-only set has 0.036 (0.033 in D): w −0.97, a −0.25 in every stage except E, where the exact T direction is smallest.

## Decisions

1. **Engine:** femwell r12, the same mesh and conventions as `results_femwell` (oxide clad,
   λ = 1.55 µm, n_Si 3.48, n_ox 1.44, no material dispersion; width = bottom width; tilt pivots at
   the bottom of each wall, the slab top for the rib).
2. **The STOP_REPORT questions were taken as answered by the femwell work:**
   - Averaging bug: femwell uses a conforming mesh, so there is no averaging.
   - Near-zero rule: the existing rule (2%, or within 2% of the group's largest when |D| < 10% of it).
   - Nominal angle: kept at 90° (a = 0), as in every femwell result.
3. **Etch depth:** slab = t − d. The rib's t column is dn/dt at fixed d = dn/dt|slab + dn/dslab,
   and dn/dd = −dn/dslab. The only new geometric solve is the slab perturbation with the rib top
   fixed, step h = 2 nm. Stage A keeps the earlier convention (slab fixed and known), so it reproduces
   the femwell summary. A* shows depth-known as a reference for Stage B.
4. **Separate unknowns per partial etch:** the 70 nm and 90 nm slabs need different etch depths, so
   d70/d90 (Stage B) and ar70/ar90 (Stage C) are separate. A second width of the same rib shares its
   etch's unknowns.
5. **Index derivatives** come directly from femwell with the material index perturbed, h = 0.002 (h/2
   0.001, 2h 0.004), reported per 1e-3 of index. One oxide index covers BOX and top cladding. Real
   BOX (thermal) and cladding (deposited) oxides differ, which would add another unknown.
6. **Temperature,** through the thermo-optic coefficients only (295 K, linear interpolation to
   1.55 µm):
   - Si: dn/dT = 1.86e-4 /K (1.87e-4 at 1.5 µm, 1.78e-4 at 2.0 µm). Source: B. J. Frey,
     D. B. Leviton, T. J. Madison, "Temperature-dependent refractive index of silicon and germanium,"
     Proc. SPIE 6273, 62732J (2006), arXiv:physics/0606168, Table 4.
   - SiO2: dn/dT = 8.4e-6 /K (8.45e-6 at 1.5 µm, 8.29e-6 at 1.6 µm). Source: D. B. Leviton,
     B. J. Frey, "Temperature-dependent absolute refractive index measurements of synthetic fused
     silica," arXiv:0805.0091, Table 6.

   Both are bulk values, Corning 7980 for the silica; film oxides can differ. Thermal expansion and
   thermal stress are ignored. The T column is exactly 0.186·∂/∂n_Si + 0.0084·∂/∂n_ox (index in
   1e-3), so Stage E equals Stage D. E' (C + T) is the informative temperature case.
7. **Identifiability:** a parameter's CRB is ∞ when the noise-weighted Jacobian has a null
   direction with a component along it. Columns are unit-normalized first, and the rank cut is at
   1e-9 of the largest singular value. Physical near-degeneracies sit near 1e-4 and exact ones near
   1e-16. Otherwise the CRB is the pseudo-inverse diagonal over the numerical rank. A first version
   used a plain pseudo-inverse of the nuisance block. It gave finite but round-off-dependent values
   for unidentifiable sets, and E ≠ D. It was replaced; E now equals D to 1e-10.
8. **Sets:** strip TE 400-600 nm (ring + MZI), TM 500 nm (ring), rib TE 500 nm at slabs 70 and 90
   (ring + MZI). Second rib widths 300/400/600 nm. The minimal set is the 550 + 600 nm strip TE rings
   plus one rib ring. Rib widths 300-600 nm were cheap enough to include.
9. **Workers:** BLAS limited to 1 thread per worker. On this machine, today, the default 16 threads
   made one strip eigen-solve 1.4-5 s erratic and about 3x slower than 1 thread. The run used up to
   3 workers (11 GB free at start).

## Checks

- **Reproduction:**
  - Stage A, strip TE alone: 1.331°, the same as `crb_alone['strip x1.0']`.
  - Two strip rings alone: 118.23°.
  - Two rings + rib 400 nm (90 nm slab): 0.0082377°, the same as the femwell summary (0.008238°).
- **Step check (h/2, h, 2h) on all 172 raw and derived derivatives** (t_d = t|slab + slab):
  - 171 PASS (largest deviation 1.17%), 1 near-zero PASS, **0 FAIL**.
  - The near-zero pass is dn_g/dt at fixed d for the 600 nm rib on the 90 nm slab: a 4.5% spread, at
    2.3% of the largest derivative in its group.
  - The STOP_REPORT's unresolved dn_g/dt|_d (50% spread with EMpy) is therefore step-stable in
    femwell.
  - The near-zero rule is applied within the geometric group (w, t, a, slab, t_d) and the material
    group (n_Si, n_ox) separately, because their units differ.
  - The earlier w, t, a derivatives were step-checked in `results_femwell` (1/180 fail).
- **Mode picks:** all 216 solves found the mode near the guess, with no fallbacks. TE fraction ≥ 0.966
  for TE and ≤ 0.045 for TM.

## Runtimes

- **Timing sample** (5 solves, 2 workers, default BLAS threads): 29 s for a strip geometry, 103-131 s
  for a rib (3 wavelengths each). Estimate before launch: 2.6 h at those times.
- **Batch:** 211 solves, 3 workers, 1 BLAS thread each. Median 21 s per strip solve and 65 s per rib
  solve per worker; about 3 solves per minute overall.
  - Ran 20:54-21:22 (120/216 cached), then stopped on your request.
  - Resumed 21:34 but paused on battery until 22:24.
  - Finished 22:59 with 0 failures.
  - Memory dropped to 4.0 GB free once (22:40), so the runner briefly went to 1 worker.

## Noise assumptions (unchanged; placeholders flagged)

| Quantity | Value | Status |
|---|---|---|
| σ_neff (ring) | 6.5e-6 | **placeholder**, no source. STOP_REPORT notes that 10 pm resonance precision gives σ_neff ≈ n_g × 6.45e-6 ≈ 2.6e-5, about 4x larger |
| σ_ng (ring) | 1e-3 | **placeholder**, no source |
| σ_phase (MZI), L = 100 µm | 1e-3 rad | **placeholder**, explicitly assumed in code |
| Rib multiplier | 1.0 (conservative) / 0.2 (optimistic) | 0.2 sourced from wide-rib literature; caveat that it may not transfer to narrow ribs |
| Priors | none | all parameters free. Realistic priors (for example on n_Si) would soften Stages D-E |

## Files

- `data/derivatives.json`: n_eff, n_g derivatives per structure at h/2, h, 2h (index per 1e-3).
- `data/step_check.csv`
- `data/stages.json`: every stage × set × noise case × step, with angle CRBs, smallest
  eigenvalue/eigenvector and condition number.
- `logs/`: batch, timing and analysis logs.
- `cache/solves_r12.jsonl`: every new femwell solve.
- Code:
  - `src/etchdepth_femwell_batch.py`: solves, femwell venv.
  - `src/etchdepth_femwell_analysis.py`: stages, main env.
  - `src/femwell_batch.py`: now also perturbs slab, n_Si and n_ox. Old cache keys are unchanged.
