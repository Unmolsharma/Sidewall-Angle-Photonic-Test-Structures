"""
Session 4, Task 1: rebuild the atlas vector to include dng direction and
per-observable noise-weighting -- Session 3 flagged that the atlas used
only dneff direction, which section 4 of Session 3 showed was an
incomplete predictor of ring-Fisher performance.

New vector per raw structure: concatenate the ring's two ALREADY
noise-weighted Fisher rows --
    [dneff_dw, dneff_dt, dneff_da] / SIGMA_NEFF   (3 numbers)
    [dng_dw,   dng_dt,   dng_da]   / SIGMA_NG      (3 numbers)
into one 6D vector, then unit-normalize. This is literally "the direction
this structure's ring contributes to the Fisher matrix" -- not a new
modeling choice layered on top, just the existing, already-validated
ring Jacobian rows stacked instead of using dneff alone. Reuses
fv.SIGMA_NEFF / fv.SIGMA_NG unchanged from Session 2/3.

Angle/cosine/distance between two structures' v2 vectors is now the
angle/cosine/distance between their combined (neff+ng) noise-weighted
information directions in 6D, which is exactly what a 2-probe ring+ring
Fisher matrix's row space is built from.
"""
import numpy as np

import fullvector_sensitivity as fv
from atlas import load_records, label, save_matrix_csv


def build_v2_vector(r):
    neff_row = np.array([r['dneff_dw'], r['dneff_dt'], r['dneff_da']]) / fv.SIGMA_NEFF
    ng_row = np.array([r['dng_dw'], r['dng_dt'], r['dng_da']]) / fv.SIGMA_NG
    raw6 = np.concatenate([neff_row, ng_row])
    norm = np.linalg.norm(raw6)
    unit6 = raw6 / norm if norm > 0 else raw6
    return raw6, unit6


def build_matrices_v2(records):
    labels = [label(r) for r in records]
    raw6 = np.array([build_v2_vector(r)[0] for r in records])
    unit6 = np.array([build_v2_vector(r)[1] for r in records])
    n = len(records)

    cos_sim = np.clip(unit6 @ unit6.T, -1.0, 1.0)
    angle_deg = np.degrees(np.arccos(cos_sim))
    euclid = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            euclid[i, j] = np.linalg.norm(raw6[i] - raw6[j])
    return labels, unit6, raw6, cos_sim, angle_deg, euclid


def numeric_summary_v2(labels, records, angle_deg, cos_sim):
    n = len(records)
    same_class, cross_class = [], []
    for i in range(n):
        for j in range(i+1, n):
            same = (records[i]['geometry'] == records[j]['geometry'] and
                    records[i]['polarization'] == records[j]['polarization'])
            entry = (labels[i], labels[j], angle_deg[i, j], cos_sim[i, j])
            (same_class if same else cross_class).append(entry)

    def stats(group, name):
        if not group:
            return f"{name}: (none)"
        angs = [g[2] for g in group]
        return (f"{name}: n={len(group)}  angle min/mean/max = "
                f"{min(angs):.1f}/{np.mean(angs):.1f}/{max(angs):.1f} deg")

    lines = [stats(same_class, "within-class pairs (v2, 6D neff+ng direction)"),
             stats(cross_class, "cross-class pairs (v2, 6D neff+ng direction)")]
    return "\n".join(lines)


if __name__ == "__main__":
    import os
    os.makedirs('data', exist_ok=True)
    records = load_records()
    labels, unit6, raw6, cos_sim, angle_deg, euclid = build_matrices_v2(records)

    save_matrix_csv('data/atlas_v2_cosine_similarity.csv', labels, cos_sim)
    save_matrix_csv('data/atlas_v2_angle_deg.csv', labels, angle_deg)
    save_matrix_csv('data/atlas_v2_euclidean_distance.csv', labels, euclid)

    print("=== Atlas v2 (6D: noise-weighted dneff-row + dng-row) numeric summary ===")
    print(numeric_summary_v2(labels, records, angle_deg, cos_sim))

    # compare to v1 (dneff-only, 3D) side by side
    from atlas import build_matrices as build_matrices_v1
    labels1, unit1, scaled1, cos1, ang1, euc1 = build_matrices_v1(records)
    print("\n=== v1 (dneff-only) vs v2 (dneff+dng noise-weighted) pairwise angle, sample ===")
    print(f"{'pair':40} {'v1 angle':>10} {'v2 angle':>10}")
    n = len(records)
    shown = 0
    for i in range(n):
        for j in range(i+1, n):
            if shown >= 12:
                break
            print(f"{labels[i]+'<->'+labels[j]:40} {ang1[i,j]:>10.1f} {angle_deg[i,j]:>10.1f}")
            shown += 1
        if shown >= 12:
            break
