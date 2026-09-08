"""
Task 3: re-evaluate Fisher metrics.

Session 2 used condition number alone and found it misleading (adding
probes that provably increase information could still make cond(F) worse).
This script recomputes every probe-pair/combo from the Session 2 landscape
using all four standard optimal-experiment-design metrics:

  - condition number  cond(F) = lambda_max/lambda_min   (numerical-stability proxy)
  - E-optimality       lambda_min(F)                     (worst-direction precision)
  - D-optimality        det(F)                            (confidence-ellipsoid volume)
  - A-optimality        trace(F^-1)                       (average parameter variance;
                                                             smaller is better)

and compares the RANKING each metric produces over the same candidate set,
to determine which one actually reflects "can we identify width, thickness,
and sidewall angle" for this problem -- rather than assuming condition
number (or any other single metric) is correct.
"""
import json
import itertools
import numpy as np

import fullvector_sensitivity as fv


def load_probe_db(path='data/probe_db.json'):
    with open(path) as f:
        records = json.load(f)
    by_id = {r['probe_id']: r for r in records}
    return records, by_id


def fisher_from_probes(probe_ids, by_id):
    rows = []
    for pid in probe_ids:
        rows.extend(by_id[pid]['fisher_rows'])
    J = np.array(rows)
    F = J.T @ J
    return F, J.shape[0]


def all_metrics(F):
    eigs = np.sort(np.linalg.eigvalsh(F))
    lambda_min, lambda_max = eigs[0], eigs[-1]
    cond = lambda_max / lambda_min if lambda_min > 0 else np.inf
    det = np.linalg.det(F)
    try:
        trace_inv = np.trace(np.linalg.inv(F))
    except np.linalg.LinAlgError:
        trace_inv = np.inf
    return dict(cond=cond, lambda_min=lambda_min, det=det, trace_inv=trace_inv, eigs=eigs.tolist())


def rank(items, key, reverse=False):
    """items: list of (label, metrics_dict). Returns list of labels sorted best-first."""
    return [lab for lab, _ in sorted(items, key=lambda x: key(x[1]), reverse=reverse)]


def spearman(rank_a, rank_b):
    """Spearman rank correlation between two label-orderings of the same items."""
    pos_a = {lab: i for i, lab in enumerate(rank_a)}
    pos_b = {lab: i for i, lab in enumerate(rank_b)}
    labels = list(pos_a.keys())
    n = len(labels)
    d2 = sum((pos_a[l] - pos_b[l])**2 for l in labels)
    return 1 - 6*d2 / (n*(n**2 - 1)) if n > 1 else float('nan')


def build_landscape(by_id, geometry_readouts, widths):
    """geometry_readouts: list of (geometry, polarization, readout) picking
    which probe family each 'slot' in the 2-probe combo comes from."""
    (g1, p1, r1), (g2, p2, r2) = geometry_readouts
    landscape = []
    for W1 in widths:
        for W2 in widths:
            id1 = f"{g1}_{p1}_{W1}_{r1}"
            id2 = f"{g2}_{p2}_{W2}_{r2}"
            if id1 == id2 or id1 not in by_id or id2 not in by_id:
                continue
            F, n = fisher_from_probes([id1, id2], by_id)
            m = all_metrics(F)
            landscape.append((f"{g1[:2]}{p1}{W1}+{g2[:2]}{p2}{W2}", m, (W1, W2)))
    return landscape


RANK_DEFICIENT_COND_THRESHOLD = 1e12  # above this, F is numerically singular (exact
                                       # rank<3 structural degeneracy, not just poorly
                                       # conditioned -- e.g. same-width ring+MZI, see
                                       # progress.md); trace(F^-1) is meaningless there
                                       # (can even go negative from roundoff) and such
                                       # points are excluded from ranking correlations,
                                       # reported separately instead of silently averaged in.


if __name__ == "__main__":
    records, by_id = load_probe_db()
    strip_widths = sorted({r['width_nm'] for r in records if r['geometry'] == 'strip'})
    slot_widths = sorted({r['width_nm'] for r in records if r['geometry'] == 'slot'})

    configs = {
        '(a) TE ring + TE ring':  [('strip', 'TE', 'ring'), ('strip', 'TE', 'ring')],
        '(b) TE ring + TM ring':  [('strip', 'TE', 'ring'), ('strip', 'TM', 'ring')],
        '(c) TE ring + MZI':      [('strip', 'TE', 'ring'), ('strip', 'TE', 'mzi')],
        '(d) TE ring + slot TE ring': [('strip', 'TE', 'ring'), ('slot', 'TE', 'ring')],
    }

    all_results = {}
    for name, geo in configs.items():
        (g1, p1, r1), (g2, p2, r2) = geo
        if g1 == g2:
            landscape = build_landscape(by_id, geo, strip_widths if g1 == 'strip' else slot_widths)
        else:
            # mixed geometry (strip ring + slot ring): widths independently indexed per side
            w1_list = sorted({r['width_nm'] for r in records if r['geometry'] == g1})
            w2_list = sorted({r['width_nm'] for r in records if r['geometry'] == g2})
            landscape = []
            for W1 in w1_list:
                for W2 in w2_list:
                    id1, id2 = f"{g1}_{p1}_{W1}_{r1}", f"{g2}_{p2}_{W2}_{r2}"
                    if id1 == id2 or id1 not in by_id or id2 not in by_id:
                        continue
                    F, n = fisher_from_probes([id1, id2], by_id)
                    m = all_metrics(F)
                    landscape.append((f"{g1[:2]}{p1}{W1}+{g2[:2]}{p2}{W2}", m, (W1, W2)))
        all_results[name] = landscape

        well_posed = [(lab, m, w) for lab, m, w in landscape if m['cond'] < RANK_DEFICIENT_COND_THRESHOLD]
        degenerate = [(lab, m, w) for lab, m, w in landscape if m['cond'] >= RANK_DEFICIENT_COND_THRESHOLD]

        conds = [m['cond'] for _, m, _ in well_posed]
        lmins = [m['lambda_min'] for _, m, _ in well_posed]
        dets = [m['det'] for _, m, _ in well_posed]
        trinvs = [m['trace_inv'] for _, m, _ in well_posed]
        print(f"\n=== {name} ({len(landscape)} combos, {len(degenerate)} exactly rank-deficient excluded from ranking) ===")
        print(f"  cond(F):      best={min(conds):.3e}  worst={max(conds):.3e}")
        print(f"  lambda_min:   best={max(lmins):.3e}  worst={min(lmins):.3e}")
        print(f"  det(F):       best={max(dets):.3e}  worst={min(dets):.3e}")
        print(f"  trace(F^-1):  best={min(trinvs):.3e}  worst={max(trinvs):.3e}")

        items = [(lab, m) for lab, m, _ in well_posed]
        rank_cond = rank(items, key=lambda m: m['cond'])                     # smaller=better
        rank_emin = rank(items, key=lambda m: m['lambda_min'], reverse=True)  # larger=better
        rank_d    = rank(items, key=lambda m: m['det'], reverse=True)         # larger=better
        rank_ainv = rank(items, key=lambda m: m['trace_inv'])                 # smaller=better

        print(f"  top-3 by cond:      {rank_cond[:3]}")
        print(f"  top-3 by lambda_min:{rank_emin[:3]}")
        print(f"  top-3 by det:       {rank_d[:3]}")
        print(f"  top-3 by trace_inv: {rank_ainv[:3]}")
        print(f"  Spearman(cond, lambda_min) = {spearman(rank_cond, rank_emin):+.3f}")
        print(f"  Spearman(cond, det)        = {spearman(rank_cond, rank_d):+.3f}")
        print(f"  Spearman(cond, trace_inv)  = {spearman(rank_cond, rank_ainv):+.3f}")
        print(f"  Spearman(det, trace_inv)   = {spearman(rank_d, rank_ainv):+.3f}")
        print(f"  Spearman(lambda_min, trace_inv) = {spearman(rank_emin, rank_ainv):+.3f}")
        if degenerate:
            print(f"  [{len(degenerate)} exactly-singular combos excluded: {[lab for lab,_,_ in degenerate][:5]}{'...' if len(degenerate)>5 else ''}]")

    # ---- everything combined: isolate slot's specific contribution ----
    print("\n=== everything combined: isolating slot's contribution ===")
    strip_only_ids = [r['probe_id'] for r in records if r['geometry'] == 'strip']
    all_ids = [r['probe_id'] for r in records]
    F_strip, n_strip = fisher_from_probes(strip_only_ids, by_id)
    F_all, n_all = fisher_from_probes(all_ids, by_id)
    m_strip = all_metrics(F_strip)
    m_all = all_metrics(F_all)
    print(f"  strip-only (TE+TM, ring+mzi, all widths), {n_strip} observables:")
    print(f"    cond={m_strip['cond']:.3e}  lambda_min={m_strip['lambda_min']:.3e}  "
          f"det={m_strip['det']:.3e}  trace_inv={m_strip['trace_inv']:.3e}")
    print(f"  + slot TE (ring+mzi, all widths), {n_all} observables:")
    print(f"    cond={m_all['cond']:.3e}  lambda_min={m_all['lambda_min']:.3e}  "
          f"det={m_all['det']:.3e}  trace_inv={m_all['trace_inv']:.3e}")
    print(f"  -> adding slot: cond {m_strip['cond']/m_all['cond']:.1f}x better, "
          f"lambda_min {m_all['lambda_min']/m_strip['lambda_min']:.1f}x better, "
          f"trace_inv {m_strip['trace_inv']/m_all['trace_inv']:.1f}x better")

    with open('data/fisher_reanalysis.json', 'w') as f:
        json.dump({name: [(lab, m) for lab, m, _ in landscape] for name, landscape in all_results.items()},
                   f, indent=2)
    print("\nSaved data/fisher_reanalysis.json")
