"""
Session 4, Task 2: stress-test two single-point assumptions Session 3 flagged
but didn't test.

(a) Fab-tolerance normalization (sigma_w=5nm, sigma_t=2nm, sigma_a=1deg) used
    to build the atlas's scaled/unit vectors -- sweep a range and check
    whether the atlas's qualitative conclusions (within/cross-class angle
    separation, PCA variance split, slot's angle-axis dominance) survive.

(b) Ring-readout noise model (SIGMA_NEFF, SIGMA_NG) reused identically for
    slot rings as for strip rings -- slot waveguides are well documented to
    have higher propagation loss than strip waveguides (sidewall-roughness
    scattering is proportionally worse given the field enhancement at the
    narrow gap's silicon-oxide interfaces), which would mean lower ring Q,
    broader resonance linewidth, and therefore WORSE (larger) sigma_neff/
    sigma_ng in practice. Re-run Session 3's isolated-slot-contribution
    Fisher comparison under pessimistic slot-noise assumptions (5x, 10x
    worse than strip) to see whether the 86x lambda_min improvement survives.
"""
import json
import numpy as np

import fullvector_sensitivity as fv
from atlas import load_records, label


def fab_vector(r, sigma_w, sigma_t, sigma_a):
    raw = np.array([r['dneff_dw'], r['dneff_dt'], r['dneff_da']])
    scaled = raw * np.array([sigma_w, sigma_t, sigma_a])
    norm = np.linalg.norm(scaled)
    return scaled / norm if norm > 0 else scaled


def angle_stats(records, sigma_w, sigma_t, sigma_a):
    units = [fab_vector(r, sigma_w, sigma_t, sigma_a) for r in records]
    n = len(records)
    same, cross = [], []
    for i in range(n):
        for j in range(i+1, n):
            cosang = np.clip(np.dot(units[i], units[j]), -1, 1)
            ang = np.degrees(np.arccos(cosang))
            same_class = (records[i]['geometry'] == records[j]['geometry'] and
                          records[i]['polarization'] == records[j]['polarization'])
            (same if same_class else cross).append(ang)
    return np.mean(same), np.mean(cross), units


def pca_variance(units):
    X = np.array(units)
    Xc = X - X.mean(axis=0)
    _, S, _ = np.linalg.svd(Xc, full_matrices=False)
    var = (S**2) / np.sum(S**2)
    return var


def slot_angle_rank(records, units):
    """Return (slot's mean |angle-axis component|, strip-TE's, strip-TM's)."""
    def mean_abs_a(geom, pol=None):
        idx = [i for i, r in enumerate(records)
               if r['geometry'] == geom and (pol is None or r['polarization'] == pol)]
        return np.mean([abs(units[i][2]) for i in idx])
    return mean_abs_a('slot', 'TE'), mean_abs_a('strip', 'TE'), mean_abs_a('strip', 'TM')


def part_a():
    records = load_records()
    print("=== (a) Fab-tolerance normalization sensitivity sweep ===")
    print(f"{'(sigma_w,sigma_t,sigma_a)':>28} {'within-mean':>12} {'cross-mean':>12} "
          f"{'PC1%':>7} {'PC2%':>7} {'PC3%':>7} {'slot|a|':>8} {'stTE|a|':>8} {'stTM|a|':>8}")
    configs = [
        (5, 2, 1),      # Session 3 baseline
        (10, 2, 1),     # looser width tolerance
        (2, 2, 1),      # tighter width tolerance
        (5, 5, 1),      # looser thickness tolerance
        (5, 2, 0.5),    # tighter angle tolerance
        (5, 2, 2),      # looser angle tolerance
        (1, 1, 1),      # no rescaling (raw units, extreme/control case)
    ]
    for sw, st, sa in configs:
        within, cross, units = angle_stats(records, sw, st, sa)
        var = pca_variance(units)
        slot_a, te_a, tm_a = slot_angle_rank(records, units)
        print(f"{str((sw, st, sa)):>28} {within:>12.1f} {cross:>12.1f} "
              f"{var[0]*100:>7.1f} {var[1]*100:>7.1f} {var[2]*100:>7.2f} "
              f"{slot_a:>8.3f} {te_a:>8.3f} {tm_a:>8.3f}")
    print("\nQualitative check across all configs: is slot's |angle-axis component| "
          "always the largest of the three families? ", end="")
    all_hold = True
    for sw, st, sa in configs:
        _, _, units = angle_stats(records, sw, st, sa)
        slot_a, te_a, tm_a = slot_angle_rank(records, units)
        if not (slot_a > te_a and slot_a > tm_a):
            all_hold = False
    print("YES, in every configuration tested" if all_hold else "NO -- breaks for some configs")


def rebuild_ring_rows(r, sigma_neff, sigma_ng):
    return [
        (np.array([r['dneff_dw'], r['dneff_dt'], r['dneff_da']]) / sigma_neff).tolist(),
        (np.array([r['dng_dw'], r['dng_dt'], r['dng_da']]) / sigma_ng).tolist(),
    ]


def rebuild_mzi_row(r, sigma_phase, L_um=fv.L_MZI, wl=fv.WL0):
    k = 2*np.pi/wl*L_um
    row = k*np.array([r['dneff_dw'], r['dneff_dt'], r['dneff_da']])
    return [(row/sigma_phase).tolist()]


def fisher_metrics(rows):
    J = np.array(rows)
    F = J.T @ J
    eigs = np.sort(np.linalg.eigvalsh(F))
    lambda_min = eigs[0]
    cond = eigs[-1]/eigs[0] if eigs[0] > 0 else np.inf
    trace_inv = np.trace(np.linalg.inv(F)) if eigs[0] > 1e-30 else np.inf
    return cond, lambda_min, trace_inv, J.shape[0]


def part_b():
    records = load_records()
    strip_records = [r for r in records if r['geometry'] == 'strip']
    slot_records = [r for r in records if r['geometry'] == 'slot']

    strip_rows = []
    for r in strip_records:
        strip_rows.extend(rebuild_ring_rows(r, fv.SIGMA_NEFF, fv.SIGMA_NG))
        strip_rows.extend(rebuild_mzi_row(r, fv.SIGMA_PHASE))
    cond0, lmin0, tinv0, n0 = fisher_metrics(strip_rows)

    print("\n=== (b) Slot-ring noise-model stress test ===")
    print(f"strip-only baseline: n_obs={n0}  cond={cond0:.3e}  lambda_min={lmin0:.3e}  trace_inv={tinv0:.3e}")
    print(f"{'slot noise vs strip':>22} {'sigma_neff':>12} {'sigma_ng':>10} {'n_obs':>6} "
          f"{'cond':>12} {'lambda_min':>12} {'trace_inv':>12} {'lmin improvement':>18}")

    for factor, label_ in [(1, 'same as strip'), (3, '3x worse'), (5, '5x worse'), (10, '10x worse'), (30, '30x worse (extreme)')]:
        sigma_neff_slot = fv.SIGMA_NEFF * factor
        sigma_ng_slot = fv.SIGMA_NG * factor
        rows = list(strip_rows)
        for r in slot_records:
            rows.extend(rebuild_ring_rows(r, sigma_neff_slot, sigma_ng_slot))
            rows.extend(rebuild_mzi_row(r, fv.SIGMA_PHASE * factor))
        cond, lmin, tinv, n = fisher_metrics(rows)
        print(f"{label_:>22} {sigma_neff_slot:>12.2e} {sigma_ng_slot:>10.2e} {n:>6} "
              f"{cond:>12.3e} {lmin:>12.3e} {tinv:>12.3e} {lmin/lmin0:>17.1f}x")


if __name__ == "__main__":
    part_a()
    part_b()
