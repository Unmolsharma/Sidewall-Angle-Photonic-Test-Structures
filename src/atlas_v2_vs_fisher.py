"""
Session 4, Task 1 (continued): re-run Session 3's atlas-vs-Fisher
correlation using the v2 (dneff+dng, noise-weighted, 6D) atlas angle
instead of the v1 (dneff-only, 3D) angle, to check whether the
wrong-signed correlation found in Session 3 section 4 resolves.
"""
import numpy as np

from atlas import load_records
from atlas_v2 import build_matrices_v2
from atlas_vs_fisher import spearman_xy, fisher_pair_metrics

if __name__ == "__main__":
    records = load_records()
    labels, unit6, raw6, cos_sim, angle_deg, euclid = build_matrices_v2(records)
    n = len(records)

    xs_angle, ys_cond, ys_lmin, ys_det, ys_tinv = [], [], [], [], []
    for i in range(n):
        for j in range(i+1, n):
            m = fisher_pair_metrics(records[i], records[j])
            xs_angle.append(angle_deg[i, j])
            ys_cond.append(m['cond'])
            ys_lmin.append(m['lambda_min'])
            ys_det.append(m['det'])
            ys_tinv.append(m['trace_inv'])

    rho_cond = spearman_xy(xs_angle, ys_cond)
    rho_lmin = spearman_xy(xs_angle, ys_lmin)
    rho_det = spearman_xy(xs_angle, ys_det)
    rho_tinv = spearman_xy(xs_angle, ys_tinv)

    print("=== Atlas v2 (6D neff+ng) vs Fisher: pairwise correlation ===")
    print(f"n pairs = {len(xs_angle)}")
    print(f"Spearman(v2 angle, cond)       = {rho_cond:+.3f}   [Session 3 v1: +0.478]")
    print(f"Spearman(v2 angle, lambda_min) = {rho_lmin:+.3f}   [Session 3 v1: -0.369]")
    print(f"Spearman(v2 angle, det)        = {rho_det:+.3f}   [Session 3 v1: +0.558]")
    print(f"Spearman(v2 angle, trace_inv)  = {rho_tinv:+.3f}   [Session 3 v1: +0.369]")

    with open('data/atlas_v2_vs_fisher_summary.txt', 'w') as f:
        f.write(f"n pairs = {len(xs_angle)}\n")
        f.write(f"Spearman(v2 angle, cond)       = {rho_cond:+.3f}  (v1 was +0.478)\n")
        f.write(f"Spearman(v2 angle, lambda_min) = {rho_lmin:+.3f}  (v1 was -0.369)\n")
        f.write(f"Spearman(v2 angle, det)        = {rho_det:+.3f}  (v1 was +0.558)\n")
        f.write(f"Spearman(v2 angle, trace_inv)  = {rho_tinv:+.3f}  (v1 was +0.369)\n")
