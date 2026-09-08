"""
Session 5, Task 2: the load-bearing noise question.

Session 4's stress_tests.py showed slot's headline lambda_min improvement
(86x) swings to 2-11x under pessimistic-but-UNSOURCED noise multipliers
(3x/5x/10x/30x picked as round numbers, not derived from any specific
measurement). Rib's 571x headline (Session 4, Task 4) carries the exact
same unvalidated assumption -- it used strip's SIGMA_NEFF/SIGMA_NG
unchanged for rib rings, with no justification given either way.

This script:
  (a) sources REAL propagation-loss numbers (dB/cm) for strip, rib, and
      slot ring resonators from literature -- reported explicitly below,
      not just used silently -- and converts them to noise multipliers via
      the standard ring relation Q_intrinsic = 2*pi*ng/(lambda*alpha), so
      Q_target/Q_strip = alpha_strip/alpha_target = loss_strip/loss_target.
      Assuming resonance-fit precision (sigma_neff) scales ~1/Q at fixed
      SNR/fitting-fraction-of-linewidth (the same assumption Session 4's
      stress test used implicitly when it applied a flat noise multiplier),
      sigma_target/sigma_strip = loss_target/loss_strip. This lets a
      dB/cm ratio be applied directly as a SIGMA_NEFF/SIGMA_NG multiplier,
      exactly like Session 4's stress_tests.py machinery, without needing
      to assume a specific ring radius or absolute Q.
  (b) re-runs the isolated Fisher-contribution test for rib under the
      sourced multipliers (same method as Session 4's slot stress test).
  (c) resolves slot's Session-4-flagged 2-86x open range into one number
      using the sourced slot multiplier.

=====================================================================
SOURCED VALUES (reported explicitly, per instruction, before use)
=====================================================================

STRIP (reference baseline -- unchanged from Session 1-4, SIGMA_NEFF=6.5e-6,
SIGMA_NG=1.0e-3, multiplier defined as 1.0x by construction).

SLOT vs STRIP -- two independent SAME-PROCESS, SAME-PAPER comparisons
found (both avoid the apples-to-oranges problem of comparing loss numbers
from unrelated fabrication runs):
  1. Baehr-Jones et al., "Reduced propagation loss in silicon strip and
     slot waveguides coated by atomic layer deposition," Opt. Express 19,
     11529 (2011): strip 2+/-1 dB/cm, slot 7+/-2 dB/cm (TiO2-ALD-coated
     devices, same fabrication run) -> ratio 3.5x.
  2. A companion/related ALD-coating slot-waveguide paper (same search
     cluster, "Low-loss silicon slot waveguides and couplers fabricated
     with optical lithography and atomic layer deposition," Opt. Express
     19, 26275, 2011): strip ~1 dB/cm, slot 3.7-5 dB/cm (100/150/220nm
     slot widths) -> ratio 3.7-5x.
  Both sources agree on ORDER (slot worse) and land in a tight 3.5-5x
  band despite different absolute loss numbers (different process
  maturity between the two papers) -- the RATIO is the transferable
  quantity, the absolute dB/cm is not. SOURCED SLOT MULTIPLIER: 4x
  (central estimate of the 3.5-5x band), applied to both SIGMA_NEFF and
  SIGMA_NG identically (same simplifying assumption Session 4 used).

RIB vs STRIP -- one same-study comparison found:
  "Record low-loss hybrid rib/wire waveguides for silicon photonic
  circuits" (Ghent/imec group, hybrid rib-straight + wire-bend design):
  straight RIB sections 0.27+/-0.012 dB/cm vs a 460nm-wide photonic WIRE
  (strip) at 1.36+/-0.06 dB/cm -> ratio 0.20x (rib is the LOWER-loss,
  i.e. HIGHER-precision, geometry here -- opposite direction from slot).
  IMPORTANT CAVEAT, stated up front and not glossed over: this
  "low-loss rib" literature class is built specifically to keep the mode
  away from the fully-etched sidewall, which is most effective for WIDE
  ribs (multi-micron); the source paper's use case (straight low-loss
  sections feeding into compact wire bends) strongly suggests a wide rib,
  not this project's narrow 400-600nm single-mode rib, where the
  fully-etched sidewall is much closer to the mode's core and the
  sidewall-avoidance mechanism is weaker. No loss measurement specific to
  a narrow (<600nm) single-mode rib was found in this search. SOURCED RIB
  MULTIPLIER: reported as a BRACKET, not a single point --
    optimistic  = 0.2x  (literature ratio, if the low-loss-rib mechanism
                  transfers to narrow ribs)
    conservative = 1.0x (no benefit assumed -- floor case; no literature
                  found suggesting a narrow rib is WORSE than strip, so
                  1.0x, not >1.0x, is used as the conservative bound)
"""
import numpy as np

import fullvector_sensitivity as fv
from stress_tests import rebuild_ring_rows, rebuild_mzi_row, fisher_metrics
from atlas import load_records


SLOT_MULTIPLIER_SOURCED = 4.0          # central estimate, 3.5-5x sourced band
RIB_MULTIPLIER_OPTIMISTIC = 0.2        # sourced literature ratio (wide-rib caveat)
RIB_MULTIPLIER_CONSERVATIVE = 1.0      # no-benefit floor


def isolated_contribution(records, target_geometry, multiplier):
    """Same method as Session 4's stress_tests.part_b: strip-only baseline
    (ring+mzi, all widths) vs strip + target_geometry (ring+mzi, all widths)
    under a noise multiplier applied to target_geometry's rows only."""
    strip_records = [r for r in records if r['geometry'] == 'strip']
    target_records = [r for r in records if r['geometry'] == target_geometry]

    strip_rows = []
    for r in strip_records:
        strip_rows.extend(rebuild_ring_rows(r, fv.SIGMA_NEFF, fv.SIGMA_NG))
        strip_rows.extend(rebuild_mzi_row(r, fv.SIGMA_PHASE))
    cond0, lmin0, tinv0, n0 = fisher_metrics(strip_rows)

    sigma_neff_t = fv.SIGMA_NEFF * multiplier
    sigma_ng_t = fv.SIGMA_NG * multiplier
    rows = list(strip_rows)
    for r in target_records:
        rows.extend(rebuild_ring_rows(r, sigma_neff_t, sigma_ng_t))
        rows.extend(rebuild_mzi_row(r, fv.SIGMA_PHASE * multiplier))
    cond, lmin, tinv, n = fisher_metrics(rows)
    return dict(n0=n0, cond0=cond0, lmin0=lmin0, tinv0=tinv0,
                n=n, cond=cond, lmin=lmin, tinv=tinv,
                lmin_improvement=lmin/lmin0, cond_improvement=cond0/cond)


if __name__ == "__main__":
    records = load_records()

    print("="*78)
    print("Task 2b: rib isolated Fisher contribution under SOURCED noise multipliers")
    print("="*78)
    for label, mult in [("same as strip [Session 4 default, unjustified]", 1.0),
                         ("optimistic (0.2x, sourced wide-rib literature)", RIB_MULTIPLIER_OPTIMISTIC),
                         ("conservative (1.0x, no benefit assumed)", RIB_MULTIPLIER_CONSERVATIVE)]:
        r = isolated_contribution(records, 'rib', mult)
        print(f"\n  rib noise multiplier: {label}")
        print(f"    n_obs={r['n']}  cond(F)={r['cond']:.3e}  lambda_min={r['lmin']:.3e}")
        print(f"    lambda_min improvement over strip-only baseline: {r['lmin_improvement']:.1f}x")

    print("\n  Session 4's reported number was 571.7x, using the SAME (1.0x, unjustified)")
    print("  assumption as the first row above -- reproduced here as a consistency check.")

    print()
    print("="*78)
    print("Task 2c: resolving slot's Session-4 2-86x range with the SOURCED multiplier")
    print("="*78)
    for label, mult in [("Session 4: same as strip [optimistic, unsourced]", 1.0),
                         ("Session 4: 3x worse [unsourced round number]", 3.0),
                         ("SOURCED: 4x worse (this session, 3.5-5x band, central est.)", SLOT_MULTIPLIER_SOURCED),
                         ("Session 4: 5x worse [unsourced round number]", 5.0),
                         ("Session 4: 10x worse [unsourced round number]", 10.0)]:
        r = isolated_contribution(records, 'slot', mult)
        print(f"  slot noise multiplier: {label}")
        print(f"    lambda_min improvement over strip-only baseline: {r['lmin_improvement']:.1f}x")

    r_sourced = isolated_contribution(records, 'slot', SLOT_MULTIPLIER_SOURCED)
    print(f"\n  RESOLVED (sourced): slot gives a {r_sourced['lmin_improvement']:.1f}x lambda_min "
          f"improvement under the literature-sourced 4x noise multiplier.")
    print("  This single number replaces Session 4's open 2-86x range as the number to")
    print("  cite going forward, superseding the round-number bracket (which is retained")
    print("  above only as a consistency cross-check against the sourced result).")
