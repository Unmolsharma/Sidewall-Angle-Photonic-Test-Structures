# Paper notes (step 1)

Working notes for the arXiv draft. Sources: `results_final/HEADLINE_NUMBERS.md` (cited below as **H:Lnn**, the
line number in that file at commit b721f15), `results_final/CLAIMS.md` (**C#n**), `results_final/ASSUMPTIONS.md`,
`CITATIONS.md`, `REPRODUCTION.md`, `ENVIRONMENT.md`, and the data files these point to.

## 1. What was simulated

- **Solver:** femwell 0.1.12 (scikit-fem 12.0.2, gmsh 4.15.2), finite-element full-vector mode solver on a
  conforming triangular mesh (ENVIRONMENT.md). Mesh "r12": core resolution 0.012 um (**H:L213**). Mode search:
  2 modes near a guess, keep the highest-n_eff mode of the requested polarization (TE fraction > 0.5), fallback 6
  modes without guess (source in `results_final/data/environment.json`, `femwell_env.solver_source`).
- **Geometry:** silicon thickness t = 220 nm (**H:L212**), oxide below and oxide cladding above (n_ox = 1.44
  for both, **H:L209**), n_Si = 3.48 (**H:L208**), no dispersion, wavelength 1.55 um (**H:L210**).
  Width = bottom width of the wall. Sidewall angle = deviation from vertical, top narrower, pivot at the bottom
  of the tilted wall (y = 0 for strips, slab top for ribs). Rib: slab of thickness s across the domain, rib from the
  slab top to t. Domain x in [-1.4, 1.4] um, y in [-1.2, 1.2] um (femwell_batch.py `build`).
- **Group index:** n_g = n_eff - lambda dn_eff/dlambda by a central difference with dlambda = 5 nm (**H:L211**)
  on the same mesh.
- **Derivatives:** central differences at steps h (w 5 nm, t 2 nm, a 1 deg, slab 2 nm, n_Si and n_ox 0.002),
  each perturbed geometry remeshed; step check at h/2, h, 2h with a 2% rule.
- **Convergence checks available:**
  - step check over 428 derivatives: 422 PASS, 3 PASS near zero, 3 FAIL (marginal, largest spread 2.25% vs
    the 2% limit) (**H:L163-L167**);
  - r12 vs r8 mesh: largest derivative difference 0.305% (**H:L152**; strip and rib TE 500 nm only);
  - fixed-mesh stencil (nominal mesh moved onto each perturbed geometry, no remeshing) for the rib 300 nm
    structures: derivatives agree with the remeshed ones within 0.78% (90 nm slab, **H:L90**) and 0.13% (50 nm
    slab, **H:L103**);
  - **no domain-size convergence check exists** (gap; see section 6).

## 2. Structures and what they are sensitive to

| Structure | n_eff (TE unless noted) | R (femwell) | Source |
|---|---|---|---|
| Strip TE 400/500/600 nm | 2.2297 / 2.4487 / 2.5705 | 1.0001 / 1.0001 / 0.9999 | `results_final/data/crosscheck_fields.json`, `crosschecks.json` |
| Strip TM 500 nm | 1.7682 | 1.0012 | same |
| Rib 300-600 nm, 50 nm slab | 2.1044-2.6021 | 0.8312-0.8732 | same |
| Rib 300-600 nm, 90 nm slab | 2.3274-2.6405 | 0.7501-0.7987 | same |

R = (dn/da)/(dn/dw)/(-h_wall tan 1 deg). R = 1 means a tilt acts exactly like a width change. Closed form
R = (y_c - y_pivot)/(h/2), with y_c the centroid of the wall-perturbation weight from the mode field: strips
agree with femwell within 0.785% (**H:L149**), ribs within 1.25% median, 5.38% max (**H:L147-L148**). For
strips, R = 1 is exact by mirror symmetry; 220 nm tan 1 deg = 3.8401 nm (**H:L143**), femwell 3.8396-3.8725
nm/deg over all strips, TE and TM, n_eff and n_g (**H:L144-L146**).

Each ring gives two observables (n_eff via resonance position, n_g via FSR); an MZI gives one phase, which is
parallel to the ring's n_eff row (C#6: ring + MZI equals a ring with the combined n_eff noise).

Unknowns by stage: B = w, t, a + etch depth d; C = + separate rib angle ar; D = + n_Si, n_ox; E = + T.

## 3. Numbers planned for the paper

**Headline (two-width set, 50 nm slab, stage D, realistic priors, 1 pm):** 0.3891 deg at h, [0.3891-0.3997]
over steps (**H:L31, L37-L38**); fixed mesh 0.3941 (+1.293%, **H:L100, L102**); 0.1 pm best case 0.04077
(**H:L122**). Monte Carlo spread 0.3976 [0.385-0.4096], spread/CRB 1.022, bias -0.015 ± 0.0089, 0 failed,
0 ambiguous, 0 at bound (**H:L67-L76**).

**Best one-etch set (8 structures, 90 nm slab, stage D, realistic, 1 pm):** 0.139 deg at h [0.1007-0.1473]
(**H:L18-L20**); fixed mesh 0.1177, -15.28% (**H:L87, L89**), fixed-mesh max over steps 0.1556 (**H:L88**).
Planned wording: "about 0.12-0.15 deg, with a ~15% uncertainty from derivative accuracy". Stage C, no priors:
0.1237 (**H:L21**). By set size 4-8: 0.4353, 0.308, 0.1533, 0.1447, 0.139 (**H:L24-L28**). 0.1 pm best case
0.0139 (**H:L118**). Monte Carlo: spread 0.1323 [0.1253-0.1396], spread/CRB 0.952, bias 0.0075 ± 0.003, 0
failed, 0 ambiguous, 4 at bound (**H:L56-L65**).

**Slab dependence (two-width set):** 30: 0.414, 40: 0.4027, 50: 0.3891, 70: 0.4521, 90: 0.5972, 110: 0.8549,
130: 1.296, 150: 4.089 deg (**H:L29-L36**, order as in the file). Optimum at 50 nm; C#12 PARTIAL.

**Single rib:** infinite (unidentifiable) at every slab without priors, stages C and D (**H:L39-L40**); with
realistic priors 24.4-45.44 deg (**H:L41-L42**), set by the priors. Etch depth known to 0.01 nm: 5.032 deg
(**H:L49**). Rib-strip angle link of 1 deg: 1.247 deg (**H:L47**). Strip + rib70 500 nm, stage D realistic:
25.19 (**H:L48**). Strip + rib70 500 + 300 nm: stage C none 0.4477, stage D realistic 0.4521, stage E realistic
0.4552 (**H:L44-L46**).

**Strips only:** infinite for every set size 1-8 (**H:L43**).

**Valley (strip + one rib 500 nm, 50 nm slab, Monte Carlo at placeholder noise):** CRB 29.44 deg, 1163 of 2000
fits at the ±3 deg box edge, slope 3.842 nm/deg, correlation 1 (**H:L77-L80**).

**Cost of ignoring the angle:** a 1 deg strip tilt fitted as width (-3.8403 nm, **H:L151**) mispredicts the rib
500 nm (90 nm slab) resonance by -710.8 pm (**H:L150**), from d lambda = lambda dn_eff/n_g.

**Search:** greedy worse than exhaustive in 7 of 32 stage/prior/slab cases, never better, largest ratio 1.12
(C#7, `results_final/data/claims_checks.json`).

**Robustness (random derivative errors on the rib 300 nm structure):** two-width set ±0.3%: 0.3638-0.4198
(5-95%), ±1%: 0.3107-0.5131 (**H:L104-L107**); best set ±0.3%: 0.1088-0.1885, ±1%: 0.06996-0.354
(**H:L91-L94**). (±3%, ±5% rows **H:L95-L98, L108-L111**.)

**Methods (EMpy):** EMpy cell-exact vs femwell over 30 structures, median |relative difference| 0.0474-0.14 by
derivative; EMpy step check 373 PASS, 222 FAIL, 53 near-zero (C#10, `claims_checks.json`). Derivative noise can
create false angle information: the best strips-only set gives 0.6135 deg from raw femwell derivatives but is
infinite with the exact symmetry imposed (C#11, `results_montecarlo/data/strips_info.json`).

**Other mirror-symmetric guides (limitation/aside only):** slot R 0.99908-1.0015, SWG 0.9988-1.0005 (**H:L153-L156**);
slot results themselves are retired and will not be reported.

## 4. Noise model and assumptions

- **Cited:** resonance-wavelength precision 1 pm (primary) and 0.1 pm (best case), Janz et al., arXiv:2306.15643
  (**H:L130-L131**); ~80 pm/K (**H:L136**), so 0.1 pm needs about 1.25 mK (**H:L137**).
- **Derived (our conversion; say so):** sigma_neff = n_g sigma_lambda/lambda = 2.75e-6; sigma_ng = sqrt(2)
  sigma_lambda n_g^2 L/lambda^2 = 0.00107; sigma_phase = 2 pi n_g L sigma_lambda/lambda^2 = 0.00111 rad
  (**H:L132-L134**), with L = 100 um (**H:L206**) and n_g = 4.26 (**H:L207**). MZI fringe read as precisely as a
  ring resonance: ASSUMED.
- **Priors ("realistic", ASSUMED):** d 2 nm, n_Si 1e-4, n_ox 1e-3, T 1 K (**H:L201-L204**); n_Si prior based
  on the cited 8.67e-5 bulk-Si uncertainty (**H:L186**). dn/dT from Frey (Si) and Leviton (SiO2), interpolated
  (**H:L187-L192**).
- **Rib noise factor 1.0** (ASSUMED, **H:L205**).
- **Sensitivity (ASSUMPTIONS.md):** noise level scales the bound linearly; rib noise factor, L and n_g of the
  conversion change it by tens of percent; priors change it by < 1%; MZI precision matters for the best set
  (-35% / +24%) but not the two-width set (< 1%).
- **Not evaluated:** n_g wavelength step, nominal indices, trapezoid wall shape, pivot location, dispersion.

## 5. Main claim vs CLAIMS.md

| Supporting result | Status | Note |
|---|---|---|
| 1. Symmetry obstruction, R = 1 for strips, 3.84 nm/deg, rib R < 1 | C#1 SUPPORTED; closed form in **H:L147-L149** | ok |
| 2a. One rib with separate etch angle unidentifiable | C#3 SUPPORTED | ok |
| 2b. Etch-depth knowledge does not rescue it | C#4 SUPPORTED | ok |
| 2c. "Only a tight angle link between the etches does" | **not in CLAIMS.md** | evidence: 1 deg link gives 1.247 deg (**H:L47**), which is still well above the two-width set; no claim row says "only". See question 2. |
| 2c'. Single-rib sentence (replaces "only a tight angle link rescues"; decided 2026-11-03) | not a CLAIMS.md row | All at stage C, 1 pm, strip TE 400-600 + rib 500 nm on the 70 nm slab, from `results_final/sym_search_1pm/data/priors_table.csv` and `d_prior_sweep.json`: no prior: unidentifiable; realistic priors 25.18 deg (H row "Priors table, strip + rib70 (500 nm), stage C, realistic priors"); depth known to 0.01 nm 5.032 deg (**H:L49**); angle link 1 deg 1.247 deg (**H:L47**); two widths of the same rib (500 + 300 nm, same slab), no priors 0.4477 deg (**H:L44**). The user's draft quoted 0.39 deg for two widths, which is the stage D / 50 nm slab value; the stage C / 70 nm value keeps all four numbers at one stage. |
| 2d. Counting rule | reasoning, not a claims row | consistent with C#3 and the stage table; state as an argument, not a result |
| 3a. Multi-width sets work | C#5 SUPPORTED | ok |
| 3b. Monte Carlo matches the CRB | C#8 SUPPORTED | ok |
| 3c. Exhaustive beats greedy | C#7 SUPPORTED as "greedy can be worse" | 7 of 32 worse, 25 equal: write "never better, worse in 7 of 32" |
| 4. ~711 pm per degree | **H:L150** COMPUTED | ok (it is -710.8 pm for one specific rib) |
| 5a. EMpy derivatives unreliable | C#10 SUPPORTED | ok |
| 5b. "EMpy can fake information" | **evidence is from femwell residuals, not EMpy** (C#11) | see question 1 |
| 5c. femwell passed the step-size and mesh checks | partly: 3 marginal FAILs of 428 (**H:L164, L167**) | state the counts |
| Slab dependence "thinner is better" | C#12 **PARTIAL** | optimum at 50 nm; write it that way |
| Air cladding | C#2 **UNSUPPORTED** | oxide-clad only; limitation |

## 6. Gaps and TODOs

1. **Domain-size convergence:** no data in results_final. Section 4 of the requested outline asks for it.
2. **Strip n_eff validation:** the "2.44-2.45 published" value has no source (UNVERIFIED); the femwell-EMpy n_eff
   comparison is not in results_final (only derivative differences are). The only verified external value is the
   EIM rib value 2.602 (Jiao & Pikus, **H:L195**) vs femwell 2.5707 (rib 500 nm, 90 nm slab,
   `crosscheck_fields.json`); EIM is an approximation, so this is a sanity check only.
3. **Fig. 4 noise mismatch:** the valley Monte Carlo (case c) exists only at the placeholder noise; cases a and b
   at 1 pm.
4. **Fig. 1 field profile:** `crosscheck_fields.json` stores the centroid y_c, not the profile. The figure script
   would re-solve the nominal modes with femwell (same settings) to plot the field; the plotted y_c must equal
   the stored one.
5. **Citations still to verify:** literature on sidewall-angle effects, SEM/AFM/scatterometry metrology, optical
   extraction of width/thickness (Xing 2018 and Dwivedi 2015 are UNVERIFIED), Fisher/CRB references, femwell.
6. **Best-set range:** "0.12-0.15 deg" spans the fixed-mesh value 0.1177 and the remeshed 2h value 0.1473; the
   h/2 value 0.1007 and the fixed-mesh 2h value 0.1556 lie outside it. Proposed wording: "0.12-0.14 deg at the
   nominal step (fixed vs remeshed mesh), 0.10-0.16 deg over step sizes".

## 7. Reference verification (2026-11-01/03)

| Key | What it is cited for | Verified how |
|---|---|---|
| janz2023 | 1 pm / 0.1 pm wavelength precision, 80 pm/K | arXiv full text read; quotes and pages in `results_final/CITATIONS.md` |
| frey2006 | n_Si uncertainty 8.67e-5, n_Si 3.48314, dn/dT of Si | arXiv full text (Tables 1, 2, 4); CITATIONS.md |
| leviton2008 | n of fused silica 1.44502, dn/dT of SiO2 | arXiv full text (Tables 4, 6); CITATIONS.md |
| johnson2002 | boundary-perturbation formula (method) | APS abstract page + PDF title page; Crossref record |
| jiao2025 | EIM n_eff 2.602 (validation only) | IntechOpen page, Section 2.3 quoted; Crossref record (book title as in Crossref) |
| lu2017 | layout-dependent variability models | Crossref record (title-level claim) |
| orji2018 | cross-sectioning destructive; scatterometry model-dependent, non-unique | PMC author manuscript, Table 1 quoted by fetch agent; DOI confirmed by NIST/imec/PMC listings |
| raymond2005 | scatterometry in manufacturing (model-based) | Crossref record (title-level) |
| xing2018 | optical extraction of width/thickness | Crossref record; full text NOT read (TLS error), so the "fixed angle" statement is a \todo |
| dwivedi2015 | optical extraction of n_eff from interferometers | Crossref record (title-level) |
| kay1993 | CRLB | library catalogue / publisher listing (ISBN 0133457117) |
| gustafsson2020, geuzaine2009 | scikit-fem, Gmsh | Crossref records |
| femwell, empy | software | installed package metadata (version, repository URL) |

Dropped: Selvaraja et al. JSTQE (DOI did not resolve), Chrostowski & Hochberg book and Bogaerts et al. 2012 (verified
but not needed in the current text; still in refs.bib, not cited, so not printed).
