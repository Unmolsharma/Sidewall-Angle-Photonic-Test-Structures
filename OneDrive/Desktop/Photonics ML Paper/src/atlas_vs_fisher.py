"""
Task 4: does the atlas explain the Fisher results?

For every pair of raw structures, compute BOTH the atlas's pairwise angle
(direction-only, geometry of the sensitivity vectors) and the resulting
2-probe Fisher metrics using each structure's 'ring' readout (same noise
model for every geometry/polarization, so the comparison is apples-to-
apples -- see progress.md for the caveat that this assumes ring readout
precision is geometry-independent, an unverified simplification carried
over from Session 2's TM-ring treatment).

Then computes the actual correlation (Spearman) between angle and each
Fisher metric across all pairs -- this is the quantitative answer to "does
the atlas predict Fisher performance," not a visual impression.

Also answers "are there unoccupied regions of parameter space" by finding
the largest empty azimuthal wedge between neighboring unit vectors in each
2D projection.
"""
import json
import numpy as np

import fullvector_sensitivity as fv
from atlas import load_records, label, build_matrices, AXIS_PAIRS, AXIS_NAMES


def spearman_xy(x, y):
    def rankdata(a):
        order = np.argsort(a)
        ranks = np.empty(len(a))
        ranks[order] = np.arange(len(a))
        return ranks
    rx, ry = rankdata(np.array(x)), rankdata(np.array(y))
    return np.corrcoef(rx, ry)[0, 1]


def fisher_pair_metrics(ra, rb):
    s_a = dict(dneff_dw=ra['dneff_dw'], dneff_dt=ra['dneff_dt'], dneff_da=ra['dneff_da'],
               dng_dw=ra['dng_dw'], dng_dt=ra['dng_dt'], dng_da=ra['dng_da'])
    s_b = dict(dneff_dw=rb['dneff_dw'], dneff_dt=rb['dneff_dt'], dneff_da=rb['dneff_da'],
               dng_dw=rb['dng_dw'], dng_dt=rb['dng_dt'], dng_da=rb['dng_da'])
    rows = []
    for s in (s_a, s_b):
        rows.append(np.array([s['dneff_dw'], s['dneff_dt'], s['dneff_da']]) / fv.SIGMA_NEFF)
        rows.append(np.array([s['dng_dw'], s['dng_dt'], s['dng_da']]) / fv.SIGMA_NG)
    J = np.array(rows)
    F = J.T @ J
    eigs = np.sort(np.linalg.eigvalsh(F))
    lambda_min = eigs[0]
    cond = eigs[-1] / eigs[0] if eigs[0] > 0 else np.inf
    det = np.linalg.det(F)
    try:
        trace_inv = np.trace(np.linalg.inv(F))
    except np.linalg.LinAlgError:
        trace_inv = np.inf
    return dict(cond=cond, lambda_min=lambda_min, det=det, trace_inv=trace_inv)


def largest_empty_wedge(unit, ax_i, ax_j):
    """Project unit vectors onto (axis_i, axis_j), compute azimuthal angle
    of each, sort, and find the largest gap between consecutive angles
    (wrapping around 360deg) -- the biggest unoccupied direction."""
    angles = np.degrees(np.arctan2(unit[:, ax_j], unit[:, ax_i])) % 360
    angles_sorted = np.sort(angles)
    gaps = np.diff(np.concatenate([angles_sorted, [angles_sorted[0] + 360]]))
    idx = np.argmax(gaps)
    start = angles_sorted[idx]
    width = gaps[idx]
    return start, width


if __name__ == "__main__":
    records = load_records()
    labels, unit, scaled, cos_sim, angle_deg, euclid = build_matrices(records)
    n = len(records)

    print("=== Atlas vs Fisher: pairwise correlation (all ring-readout pairs) ===")
    xs_angle, ys_cond, ys_lmin, ys_det, ys_tinv = [], [], [], [], []
    rows_out = []
    for i in range(n):
        for j in range(i+1, n):
            m = fisher_pair_metrics(records[i], records[j])
            xs_angle.append(angle_deg[i, j])
            ys_cond.append(m['cond'])
            ys_lmin.append(m['lambda_min'])
            ys_det.append(m['det'])
            ys_tinv.append(m['trace_inv'])
            rows_out.append(dict(pair=f"{labels[i]}<->{labels[j]}", angle_deg=angle_deg[i, j], **m))

    rho_cond = spearman_xy(xs_angle, ys_cond)
    rho_lmin = spearman_xy(xs_angle, ys_lmin)
    rho_det = spearman_xy(xs_angle, ys_det)
    rho_tinv = spearman_xy(xs_angle, ys_tinv)
    print(f"n pairs = {len(xs_angle)}")
    print(f"Spearman(angle, cond)       = {rho_cond:+.3f}  (expect negative if angle predicts better conditioning)")
    print(f"Spearman(angle, lambda_min) = {rho_lmin:+.3f}  (expect positive if angle predicts better worst-direction precision)")
    print(f"Spearman(angle, det)        = {rho_det:+.3f}")
    print(f"Spearman(angle, trace_inv)  = {rho_tinv:+.3f}  (expect negative -- larger angle -> smaller avg variance)")

    with open('data/atlas_vs_fisher_pairs.json', 'w') as f:
        json.dump(rows_out, f, indent=2)

    print("\n=== Unoccupied regions: largest empty azimuthal wedge per projection ===")
    for ax_i, ax_j, title in AXIS_PAIRS:
        start, width = largest_empty_wedge(unit, ax_i, ax_j)
        print(f"  {title}: empty from {start:.1f} to {(start+width)%360:.1f} deg "
              f"(width={width:.1f} deg of 360)")

    with open('data/atlas_vs_fisher_summary.txt', 'w') as f:
        f.write(f"n pairs = {len(xs_angle)}\n")
        f.write(f"Spearman(angle, cond)       = {rho_cond:+.3f}\n")
        f.write(f"Spearman(angle, lambda_min) = {rho_lmin:+.3f}\n")
        f.write(f"Spearman(angle, det)        = {rho_det:+.3f}\n")
        f.write(f"Spearman(angle, trace_inv)  = {rho_tinv:+.3f}\n")
