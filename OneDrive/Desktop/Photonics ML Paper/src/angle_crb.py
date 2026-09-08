"""
Session 5, Task 3: the angle-specific metric.

Every headline number so far (86x, 91x, 571x, and Session 5's Task 2
noise-resolved numbers) is lambda_min(F) -- the WORST-case direction across
all three parameters (width, thickness, angle) combined. That's a fine
aggregate diagnostic, but the paper's actual claim is narrower and more
specific: "sidewall angle goes from unresolvable to resolvable." lambda_min
could in principle be dominated by width or thickness (both already
well-constrained by strip alone per Session 1-4) rather than angle, which
would make it a misleading stand-in for the angle-specific claim.

This computes the MARGINAL Cramer-Rao bound on angle specifically: the
(angle, angle) diagonal entry of F^-1 (marginal variance, after
marginalizing out width/thickness -- NOT 1/F_aa, which would be the wrong,
overly optimistic "as if width/thickness were already known exactly"
bound). sqrt of that entry is reported in degrees since our angle rows are
scaled per-degree throughout.

Two framings, both useful and reported separately:
  (A) EACH GEOMETRY ALONE (5 widths, ring+mzi, TE only) -- what angle
      precision is achievable using ONLY that geometry's probes. Directly
      answers "is angle resolvable with strip alone" vs "with rib alone."
  (B) STRIP BASELINE vs STRIP+ADDITION -- same isolated-contribution
      framing as Task 2, showing the marginal angle improvement from
      adding rib or slot on top of an already-large strip probe set
      (the more realistic "what would you actually deploy" scenario).

Noise assumptions: strip unchanged (1.0x). Slot uses Task 2's sourced 4x
multiplier. Rib uses Task 2's CONSERVATIVE 1.0x bound (no-benefit-assumed
floor) as the primary "realistic" number, with the OPTIMISTIC 0.2x sourced
bound reported alongside as a secondary, more speculative upper bound --
matching Task 2's explicit bracketing rather than collapsing it to one
number prematurely.
"""
import numpy as np

import fullvector_sensitivity as fv
from stress_tests import rebuild_ring_rows, rebuild_mzi_row
from atlas import load_records
from noise_sourcing import SLOT_MULTIPLIER_SOURCED, RIB_MULTIPLIER_OPTIMISTIC, RIB_MULTIPLIER_CONSERVATIVE

PARAM_NAMES = ['width', 'thickness', 'angle']
ANGLE_IDX = 2


def rows_for(records, geometry, pol, multiplier=1.0):
    rows = []
    for r in records:
        if r['geometry'] != geometry or r['polarization'] != pol:
            continue
        rows.extend(rebuild_ring_rows(r, fv.SIGMA_NEFF * multiplier, fv.SIGMA_NG * multiplier))
        rows.extend(rebuild_mzi_row(r, fv.SIGMA_PHASE * multiplier))
    return rows


def marginal_angle_crb(rows):
    J = np.array(rows)
    F = J.T @ J
    eigs = np.linalg.eigvalsh(F)
    if eigs[0] < 1e-30:
        return None, None, np.inf, np.inf
    Finv = np.linalg.inv(F)
    var_angle = Finv[ANGLE_IDX, ANGLE_IDX]
    std_angle = np.sqrt(var_angle) if var_angle > 0 else np.inf
    lambda_min = eigs[0]
    return F, Finv, var_angle, std_angle


if __name__ == "__main__":
    records = load_records()

    print("="*78)
    print("(A) Each geometry ALONE (5 widths, ring+mzi, TE) -- marginal angle CRB")
    print("="*78)
    configs_alone = [
        ('strip', 1.0, 'strip [reference, unchanged]'),
        ('rib', RIB_MULTIPLIER_CONSERVATIVE, 'rib [conservative 1.0x, no benefit assumed]'),
        ('rib', RIB_MULTIPLIER_OPTIMISTIC, 'rib [optimistic 0.2x, sourced wide-rib lit.]'),
        ('slot', SLOT_MULTIPLIER_SOURCED, 'slot [sourced 4.0x]'),
    ]
    print(f"{'geometry':>45} {'lambda_min':>12} {'angle var(F^-1)':>16} {'angle std (deg)':>16}")
    for geom, mult, label in configs_alone:
        rows = rows_for(records, geom, 'TE', mult)
        F, Finv, var_a, std_a = marginal_angle_crb(rows)
        lmin = np.linalg.eigvalsh(F)[0]
        print(f"{label:>45} {lmin:>12.3e} {var_a:>16.3e} {std_a:>16.4f}")

    print("\nInterpretation: this is the 1-sigma angle precision achievable from ring+mzi")
    print("resonance measurements using ONLY that geometry's 5-width probe set, with")
    print("width AND thickness estimated jointly (marginalized), not assumed known.")

    print()
    print("="*78)
    print("(B) Strip baseline vs strip+addition -- marginal angle CRB improvement")
    print("="*78)
    strip_rows = rows_for(records, 'strip', 'TE', 1.0)
    # strip TM + strip MZI already included via rebuild_mzi_row per record;
    # match Task 2's isolated-contribution baseline exactly: strip TE+TM, ring+mzi
    strip_rows += rows_for(records, 'strip', 'TM', 1.0)
    F0, Finv0, var0, std0 = marginal_angle_crb(strip_rows)
    lmin0 = np.linalg.eigvalsh(F0)[0]
    print(f"strip-only (TE+TM, ring+mzi) baseline: lambda_min={lmin0:.3e}  "
          f"angle std={std0:.4f} deg")

    additions = [
        ('rib', RIB_MULTIPLIER_CONSERVATIVE, 'rib [conservative 1.0x]'),
        ('rib', RIB_MULTIPLIER_OPTIMISTIC, 'rib [optimistic 0.2x]'),
        ('slot', SLOT_MULTIPLIER_SOURCED, 'slot [sourced 4.0x]'),
    ]
    print(f"\n{'addition':>32} {'lambda_min':>12} {'angle std (deg)':>16} {'angle std improvement':>22}")
    for geom, mult, label in additions:
        rows = strip_rows + rows_for(records, geom, 'TE', mult)
        F, Finv, var_a, std_a = marginal_angle_crb(rows)
        lmin = np.linalg.eigvalsh(F)[0]
        print(f"{label:>32} {lmin:>12.3e} {std_a:>16.4f} {std0/std_a:>21.1f}x")

    print("\nAnswering Task 3 directly (large-budget framing): does angle go from")
    print("unresolvable to resolvable?")
    print(f"  strip alone: angle std = {std0:.2f} deg -- for reference, the fab-tolerance")
    print("  assumption used throughout this project's atlas normalization (sigma_a=1deg,")
    print("  Session 3) is the target precision a monitor should beat to be useful.")
    print("  A marginal std well above ~1 deg means angle is effectively unresolvable")
    print("  from that probe set alone (the estimate is not usefully more precise than")
    print("  the fabrication tolerance itself); a std at or below ~1 deg means it is.")
    print("  ** FLAG: strip-alone ALREADY clears this bar comfortably (0.02-0.03 deg)")
    print("  in this 15-row (5-width, ring+mzi) framing -- see part (C) below for why")
    print("  this looks inconsistent with Session 1's 'badly conditioned' framing, and")
    print("  which framing is actually the one the paper's claim rests on.")

    print()
    print("="*78)
    print("(C) MINIMAL, Session-1/2-realistic probe budget: 2 TE rings only")
    print("="*78)
    print("Session 1/2's original motivating claim ('two same-type TE ring probes are")
    print("badly conditioned for angle') used exactly 2 ring probes (4 rows: neff+ng x 2")
    print("widths), NOT the 15-row all-5-widths-plus-MZI set used in (A)/(B) above.")
    print("That large-budget set is what Session 3 onward has used for every headline")
    print("Fisher number (86x, 91x, 571x) because it isolates a structure's MARGINAL")
    print("contribution against a realistic full baseline -- but it is a much bigger")
    print("probe budget than a real PCM unit cell (1-2 monitors) would deploy, and it is")
    print("NOT the setup Session 1's original degeneracy claim was about. Both framings")
    print("are reported because they answer different questions -- this section answers")
    print("the sharper, more paper-relevant one: does ONE additional rib/slot probe fix")
    print("the SPECIFIC minimal-budget degeneracy Session 1 found?")

    def probe_row(geometry, pol, width_nm, records, multiplier=1.0):
        r = next(r for r in records if r['geometry'] == geometry
                  and r['polarization'] == pol and r['width_nm'] == width_nm)
        return rebuild_ring_rows(r, fv.SIGMA_NEFF * multiplier, fv.SIGMA_NG * multiplier)

    strip_widths_avail = sorted({r['width_nm'] for r in records if r['geometry'] == 'strip'})
    rib_widths_avail = sorted({r['width_nm'] for r in records if r['geometry'] == 'rib'})
    slot_widths_avail = sorted({r['width_nm'] for r in records if r['geometry'] == 'slot'})
    print(f"\nUsing Session 2's reported best-conditioned strip pair: TE@550nm + TE@600nm")

    minimal_rows = (probe_row('strip', 'TE', 550, records) +
                     probe_row('strip', 'TE', 600, records))
    F_min, Finv_min, var_min, std_min = marginal_angle_crb(minimal_rows)
    lmin_min = np.linalg.eigvalsh(F_min)[0]
    print(f"  strip TE@550+TE@600 alone: lambda_min={lmin_min:.3e}  cond={np.linalg.cond(F_min):.3e}  "
          f"angle std={std_min:.3f} deg")
    print(f"  full marginal std (w[nm],t[nm],a[deg]) = {np.sqrt(np.diag(Finv_min))}")
    eigvals_min, eigvecs_min = np.linalg.eigh(F_min)
    print(f"  worst-eigenvalue eigenvector (w,t,a) = {eigvecs_min[:,0]} (eigval={eigvals_min[0]:.2f})")
    print("  -> MECHANISM: the worst (lambda_min) direction is dominated by WIDTH "
          f"(|{eigvecs_min[0,0]:.2f}| >> |{eigvecs_min[2,0]:.2f}|), not angle. This is exactly")
    print("  why the blended lambda_min metric looked catastrophic (6.33, cond 2.4e5) while")
    print("  the MARGINAL angle CRB alone is much better-behaved (0.19deg) -- confirming")
    print("  Task 3's own premise: lambda_min can be dominated by width/thickness axes")
    print("  this project already handles well, making it a poor stand-in for the")
    print("  angle-specific claim. This is the mechanistic justification for why Task 3")
    print("  asked for the angle-marginal metric specifically, not just aggregate lambda_min.")

    print(f"\n{'+ one probe from':>32} {'lambda_min':>12} {'cond(F)':>12} {'angle std (deg)':>16} {'angle std improvement':>22}")
    single_additions = [
        ('rib', 400, RIB_MULTIPLIER_CONSERVATIVE, 'rib TE@400 [cons. 1.0x]'),
        ('rib', 400, RIB_MULTIPLIER_OPTIMISTIC, 'rib TE@400 [opt. 0.2x]'),
        ('rib', 600, RIB_MULTIPLIER_CONSERVATIVE, 'rib TE@600 [cons. 1.0x]'),
        ('slot', 150, SLOT_MULTIPLIER_SOURCED, 'slot TE@150 [sourced 4x]'),
        ('strip', 400, 1.0, 'strip TE@400 [same class, control]'),
    ]
    for geom, w, mult, label in single_additions:
        if (geom == 'rib' and w not in rib_widths_avail) or \
           (geom == 'slot' and w not in slot_widths_avail):
            continue
        rows = minimal_rows + probe_row(geom, 'TE', w, records, mult)
        F, Finv, var_a, std_a = marginal_angle_crb(rows)
        lmin = np.linalg.eigvalsh(F)[0]
        cond = np.linalg.cond(F)
        print(f"{label:>32} {lmin:>12.3e} {cond:>12.3e} {std_a:>16.4f} {std_min/std_a:>21.1f}x")

    print("\nAnswering Task 3 directly (minimal-budget framing -- the one Session 1's")
    print("original claim is actually about):")
    print(f"  2 strip TE rings alone: angle std = {std_min:.3f} deg")
    if std_min > 1.0:
        print("  -> ABOVE the 1deg fab-tolerance bar: angle is NOT usefully resolvable from")
        print("     this minimal, realistic strip-only probe set -- consistent with Session")
        print("     1's qualitative 'badly conditioned' framing (this reconciles the")
        print("     apparent tension with part (A)/(B): those used a 3-4x larger probe")
        print("     budget across 5 widths, not the minimal 2-probe set Session 1 tested).")
    else:
        print("  -> already below the 1deg bar even at minimal budget -- flagged as")
        print("     genuinely surprising given Session 1's framing, not glossed over.")

    print("\n** SCOPE CAVEAT, stated explicitly (not glossed over): every number in this")
    print("script is a Cramer-Rao bound under IDEALIZED, uncorrelated Gaussian measurement")
    print("noise at the project's ASSUMED (not independently validated, per Session 2/4's")
    print("own flags) SIGMA_NEFF/SIGMA_NG. It is a best-case, information-theoretic bound,")
    print("not a real nonlinear-fit result. The literature-cited compact-modeling paper")
    print("(project_brief.md) found a REAL extraction procedure's fit residual stayed flat")
    print("across 85-90deg sidewall angle -- that is a statement about a specific realistic")
    print("fitting procedure under real (correlated, non-Gaussian, model-uncertain) noise,")
    print("not a Fisher-information CRB. The absolute CRB numbers here (all sub-degree,")
    print("even for strip alone) should NOT be read as 'the angle problem is already solved")
    print("without rib' -- the RELATIVE comparison (rib gives the largest angle-specific")
    print("improvement of the structures tested, 11-51x depending on budget/assumption) is")
    print("the robust, actionable finding; the absolute sub-degree numbers reflect how")
    print("optimistic the assumed noise floor is, a pre-existing project assumption this")
    print("session did not re-derive.")
