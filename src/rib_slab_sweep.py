"""
Session 5, Task 4: mechanism -- is rib's PC3-gap-filling / angle-sensitivity
behavior a property of "rib" as a family, or specifically tied to the
SLAB_H=90nm value Session 4 happened to fix it at? This matters before any
future optimizer is allowed to treat "rib" as a single fixed recipe rather
than a slab-height family with its own internal design freedom.

Sweeps slab height across a realistic partial-etch range (50-150nm, i.e.
etch depth 70-170nm out of the fixed 220nm total thickness -- all within
what's fabricable by a single additional lithography/etch step on a
standard 220nm SOI platform) at both narrow (300, 400nm -- the PC3-relevant
regime per Session 5 Task 1) and standard (500, 600nm) rib widths.

For each (slab_h, width) point, computes:
  (a) the PC3-gap projection, using EXACTLY Session 4's method: PC3
      direction is re-derived from the strip+slot-only (20-structure,
      pre-rib) fab-tolerance-scaled sensitivity set, then each rib
      variant's own fab-scaled unit vector is projected onto it.
  (b) the isolated Fisher lambda_min improvement (strip-only baseline +
      rib-at-this-slab-height, all widths, ring+mzi) -- the same method
      Task 2 used, to see whether the headline improvement number is
      robust to the slab-height choice, not just the PC3 projection.
"""
import numpy as np
import time
from functools import partial

import fullvector_sensitivity as fv
from rib_waveguide import solve_neff_rib
from atlas import load_records
from probe_db import FAB_SCALE
from stress_tests import rebuild_ring_rows, rebuild_mzi_row, fisher_metrics


def sensitivities_rib_at_slab(width_nm, slab_h_nm, thick_nm=220.0, angle_deg=0.0,
                               dw=5.0, dt=2.0, da=1.0, polarization='TE'):
    """rib_waveguide.sensitivities_rib doesn't expose slab_h_um (it was fixed
    at module-level SLAB_H in Session 4) -- bind it via a solve_fn partial
    instead of modifying that already-validated function's signature."""
    solve_fn = partial(solve_neff_rib, slab_h_um=slab_h_nm/1000)
    return fv.sensitivities(width_nm, thick_nm, angle_deg, dw, dt, da, polarization, solve_fn=solve_fn)


def pc3_direction_pre_rib():
    """Re-derive the PC3 gap direction from strip+slot only (excludes rib/
    swg), matching Session 4's stated method of projecting onto 'the gap
    BEFORE either candidate was added.'"""
    records = load_records()
    pre_rib = [r for r in records if r['geometry'] in ('strip', 'slot')]
    X = np.array([r['sens_vector_scaled'] for r in pre_rib])
    Xc = X - X.mean(axis=0)
    _, S, Vt = np.linalg.svd(Xc, full_matrices=False)
    var = (S**2) / np.sum(S**2)
    print(f"PC3 (pre-rib, strip+slot only, n={len(pre_rib)}): "
          f"variance={var[2]*100:.4f}%  loadings(w,t,a)={Vt[2].round(3)}")
    return Vt[2]


def fab_unit_vector(s):
    raw = np.array([s['dneff_dw'], s['dneff_dt'], s['dneff_da']])
    scaled = raw * FAB_SCALE
    norm = np.linalg.norm(scaled)
    return (scaled / norm) if norm > 0 else scaled


def isolated_lambda_min_improvement(records, slab_h_nm, widths_nm):
    strip_records = [r for r in records if r['geometry'] == 'strip']
    strip_rows = []
    for r in strip_records:
        strip_rows.extend(rebuild_ring_rows(r, fv.SIGMA_NEFF, fv.SIGMA_NG))
        strip_rows.extend(rebuild_mzi_row(r, fv.SIGMA_PHASE))
    _, lmin0, _, _ = fisher_metrics(strip_rows)

    rib_rows = []
    for w in widths_nm:
        s = sensitivities_rib_at_slab(w, slab_h_nm)
        r = dict(dneff_dw=s['dneff_dw'], dneff_dt=s['dneff_dt'], dneff_da=s['dneff_da'],
                  dng_dw=s['dng_dw'], dng_dt=s['dng_dt'], dng_da=s['dng_da'])
        rib_rows.extend(rebuild_ring_rows(r, fv.SIGMA_NEFF, fv.SIGMA_NG))
        rib_rows.extend(rebuild_mzi_row(r, fv.SIGMA_PHASE))
    rows = strip_rows + rib_rows
    _, lmin, _, _ = fisher_metrics(rows)
    return lmin / lmin0


if __name__ == "__main__":
    t0 = time.time()

    pc3 = pc3_direction_pre_rib()

    print("\n=== (a) PC3-gap projection vs slab height, at narrow and standard widths ===")
    slab_heights_nm = [50, 70, 90, 110, 130, 150]
    widths_nm = [300, 400, 500, 600]
    print(f"{'slab_h(nm)':>10}", end="")
    for w in widths_nm:
        print(f"  W={w}nm".rjust(12), end="")
    print()
    proj_table = {}
    for sh in slab_heights_nm:
        row = []
        for w in widths_nm:
            s = sensitivities_rib_at_slab(w, sh)
            unit = fab_unit_vector(s)
            proj = abs(np.dot(unit, pc3))
            row.append(proj)
        proj_table[sh] = row
        print(f"{sh:>10}", end="")
        for p in row:
            print(f"{p:>12.4f}", end="")
        print()

    print(f"\nSession 4's reference point was slab_h=90nm: |proj| = "
          f"{proj_table[90][widths_nm.index(400)]:.4f} @ 400nm "
          f"(Session 4 reported 0.154 for the same nominal point, from the 15-structure-")
    print("PC3 direction rather than this session's 20-structure -- pre-rib -- direction; "
          "small numeric differences from that are expected and not a discrepancy.)")

    print("\n=== (b) isolated Fisher lambda_min improvement vs slab height ===")
    print("(strip-only baseline + rib-at-this-slab-height, all 4 widths, ring+mzi, ")
    print(" strip's own SIGMA_NEFF/SIGMA_NG unchanged -- isolates slab_h's effect alone)")
    print(f"{'slab_h(nm)':>10} {'etch_depth(nm)':>15} {'lambda_min improvement':>24}")
    lmin_table = {}
    for sh in slab_heights_nm:
        imp = isolated_lambda_min_improvement(load_records(), sh, widths_nm)
        lmin_table[sh] = imp
        print(f"{sh:>10} {220-sh:>15} {imp:>23.1f}x")

    print("\n=== Task 4 answer: does slab height drive the effect, or is it robust? ===")
    proj_400 = {sh: proj_table[sh][widths_nm.index(400)] for sh in slab_heights_nm}
    proj_min, proj_max = min(proj_400.values()), max(proj_400.values())
    lmin_min, lmin_max = min(lmin_table.values()), max(lmin_table.values())
    print(f"PC3 projection @ 400nm across slab_h=50-150nm: {proj_min:.4f} to {proj_max:.4f} "
          f"({proj_max/proj_min:.2f}x spread)")
    print(f"lambda_min improvement across slab_h=50-150nm: {lmin_min:.1f}x to {lmin_max:.1f}x "
          f"({lmin_max/lmin_min:.2f}x spread)")
    best_sh = max(lmin_table, key=lmin_table.get)
    worst_sh = min(lmin_table, key=lmin_table.get)
    print(f"Best slab height tested: {best_sh}nm ({lmin_table[best_sh]:.1f}x). "
          f"Worst: {worst_sh}nm ({lmin_table[worst_sh]:.1f}x).")
    if lmin_max / lmin_min < 1.5 and proj_max / proj_min < 1.5:
        print("-> ROBUST: the effect holds across the whole realistic partial-etch range, not")
        print("   just at the specific 90nm value Session 4 happened to fix. 'Rib' can be")
        print("   treated as a family with a broad sweet spot, not a single fragile recipe.")
    else:
        print("-> SLAB-HEIGHT DEPENDENT: the effect's magnitude varies meaningfully with slab")
        print("   height -- an inverse-design phase MUST treat slab height as a free/optimized")
        print("   parameter, not fix it at 90nm by convention.")

    print(f"\nElapsed: {time.time()-t0:.1f}s")
