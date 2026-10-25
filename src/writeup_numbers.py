"""
Every number stated in WRITEUP.md, computed from a given sensitivity database
(+ rib slab-sweep derivatives), reusing the original analysis functions
(fisher_reanalysis, atlas, atlas_v2, atlas_vs_fisher, stress_tests,
noise_sourcing, angle_crb, rib_slab_sweep). Noise assumptions are the
originals (fv.SIGMA_*, noise_sourcing multipliers), unchanged.

Run on the original database first to confirm it reproduces the writeup,
then on the corrected one:

  python writeup_numbers.py --db ../data/sensitivity_db.json \
      --slab ../results_corrected/data/rib_slab_sweep_db_original.json \
      --out ../results_corrected/data/writeup_numbers_original.json
  python writeup_numbers.py --db ../results_corrected/data/sensitivity_db.json \
      --slab ../results_corrected/data/rib_slab_sweep_db_cell_exact.json \
      --out ../results_corrected/data/writeup_numbers_corrected.json

Subsets follow the writeup: the atlas angle statistics, the Spearman
correlation, the PCA split and the seven-scaling check were stated for the
first 15 structures (strip TE, strip TM, slot TE); they are also given for
all 30.
"""
import argparse
import json

import numpy as np
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import squareform

import fullvector_sensitivity as fv
from probe_db import build_probe_db, FAB_SCALE
from fisher_reanalysis import fisher_from_probes, all_metrics, build_landscape
from atlas import build_matrices, label
from atlas_v2 import build_matrices_v2
from atlas_vs_fisher import spearman_xy, fisher_pair_metrics
from stress_tests import angle_stats, pca_variance, slot_angle_rank, rebuild_ring_rows, rebuild_mzi_row, fisher_metrics
from noise_sourcing import (isolated_contribution, SLOT_MULTIPLIER_SOURCED,
                            RIB_MULTIPLIER_OPTIMISTIC, RIB_MULTIPLIER_CONSERVATIVE)
from angle_crb import rows_for, marginal_angle_crb
from information_map import pca_2d

FAMS = ['strip-TE', 'strip-TM', 'slot-TE', 'slot-TM', 'rib-TE', 'swg-TE']


def fam(r):
    return f"{r['geometry']}-{r['polarization']}"


def subset15(records):
    return [r for r in records if fam(r) in ('strip-TE', 'strip-TM', 'slot-TE')]


def atlas_stats(records):
    labels, unit, scaled, cos_sim, angle_deg, _ = build_matrices(records)
    n = len(records)
    same, cross = [], []
    for i in range(n):
        for j in range(i + 1, n):
            (same if fam(records[i]) == fam(records[j]) else cross).append(angle_deg[i, j])
    dist = 1 - cos_sim
    np.fill_diagonal(dist, 0)
    order = leaves_list(linkage(squareform(np.clip(dist, 0, None), checks=False), method='average'))
    seq = [fam(records[i]) for i in order]
    blocks = [seq[0]] + [b for a, b in zip(seq, seq[1:]) if a != b]
    return dict(within_mean=float(np.mean(same)), within_min=float(min(same)), within_max=float(max(same)),
                cross_mean=float(np.mean(cross)), cross_min=float(min(cross)), cross_max=float(max(cross)),
                n_blocks=len(blocks), n_families=len(set(seq)), clusters_cleanly=len(blocks) == len(set(seq)),
                leaf_family_order=blocks)


def spearman_block(records):
    _, _, _, _, ang1, _ = build_matrices(records)
    _, _, _, _, ang2, _ = build_matrices_v2(records)
    n = len(records)
    a1, a2, lmin = [], [], []
    for i in range(n):
        for j in range(i + 1, n):
            m = fisher_pair_metrics(records[i], records[j])
            a1.append(ang1[i, j]); a2.append(ang2[i, j]); lmin.append(m['lambda_min'])
    return dict(n_pairs=len(a1), v1_angle_vs_lambda_min=float(spearman_xy(a1, lmin)),
                v2_angle_vs_lambda_min=float(spearman_xy(a2, lmin)))


def pca_block(records):
    X = np.array([r['sens_vector_scaled'] for r in records])
    _, explained, Vt = pca_2d(X)
    return dict(n=len(records), explained_pct=(100 * explained).tolist(), pc3_loadings=Vt[2].tolist())


def seven_scalings(records):
    configs = [(5, 2, 1), (10, 2, 1), (2, 2, 1), (5, 5, 1), (5, 2, 0.5), (5, 2, 2), (1, 1, 1)]
    out = []
    for c in configs:
        within, cross, units = angle_stats(records, *c)
        var = pca_variance(units)
        sa, te, tm = slot_angle_rank(records, units)
        out.append(dict(config=list(c), within_mean=float(within), cross_mean=float(cross),
                        pc_pct=(100 * var).tolist(), slot_angle_largest=bool(sa > te and sa > tm)))
    return dict(configs=out, max_pc3_pct=max(o['pc_pct'][2] for o in out),
                slot_angle_largest_all=all(o['slot_angle_largest'] for o in out))


def pc3_pre_rib(records):
    pre = [r for r in records if r['geometry'] in ('strip', 'slot')]
    X = np.array([r['sens_vector_scaled'] for r in pre])
    _, S, Vt = np.linalg.svd(X - X.mean(axis=0), full_matrices=False)
    return Vt[2], float(100 * S[2]**2 / np.sum(S**2))


def unit_from(s):
    v = np.array([s['dneff_dw'], s['dneff_dt'], s['dneff_da']]) * FAB_SCALE
    return v / np.linalg.norm(v)


def strip_baseline_rows(records):
    rows = []
    for r in records:
        if r['geometry'] == 'strip':
            rows.extend(rebuild_ring_rows(r, fv.SIGMA_NEFF, fv.SIGMA_NG))
            rows.extend(rebuild_mzi_row(r, fv.SIGMA_PHASE))
    return rows


def main(db, slab_path):
    with open(db) as f:
        records = json.load(f)
    with open(slab_path) as f:
        slab = json.load(f)
    probes = build_probe_db(records)
    by_id = {p['probe_id']: p for p in probes}
    rec = {r['id']: r for r in records}
    N = {}

    # --- rings alone: strip TE ring + strip TE ring landscape ---
    land = build_landscape(by_id, [('strip', 'TE', 'ring'), ('strip', 'TE', 'ring')], [400, 450, 500, 550, 600])
    conds = [m['cond'] for _, m, _ in land]
    best_c = min(land, key=lambda x: x[1]['cond'])
    best_l = max(land, key=lambda x: x[1]['lambda_min'])
    N['ring_pairs'] = dict(cond_best=min(conds), cond_worst=max(conds), spread=max(conds) / min(conds),
                           best_pair_by_cond=list(best_c[2]), best_pair_by_lambda_min=list(best_l[2]),
                           best_pair_by_lambda_min_value=best_l[1]['lambda_min'])

    # --- additions to the TE 550 + 600 pair ---
    base = ['strip_TE_550_ring', 'strip_TE_600_ring']
    add = {}
    # the original writeup labelled the third TE ring "400 nm"; its 7.64 / 3.0e5 values are the 500 nm ring
    for name, extra in [('baseline', []), ('+TE ring 400', ['strip_TE_400_ring']), ('+TE ring 500', ['strip_TE_500_ring']),
                        ('+TM ring 400', ['strip_TM_400_ring']), ('+MZI 600 TE', ['strip_TE_600_mzi'])]:
        F, _ = fisher_from_probes(base + extra, by_id)
        m = all_metrics(F)
        add[name] = dict(lambda_min=m['lambda_min'], cond=m['cond'], trace_inv=m['trace_inv'])
    N['additions_550_600'] = add
    F, _ = fisher_from_probes(['strip_TE_500_ring', 'strip_TE_500_mzi'], by_id)
    N['same_width_ring_mzi_cond'] = float(np.linalg.cond(F))

    # --- atlas, Spearman, angle components ---
    r15 = subset15(records)
    N['atlas_15'] = atlas_stats(r15)
    N['atlas_30'] = atlas_stats(records)
    N['spearman_15'] = spearman_block(r15)
    N['spearman_30'] = spearman_block(records)
    comps = {}
    for fm in FAMS:
        a = [r['sens_vector_unit'][2] for r in records if fam(r) == fm]
        comps[fm] = [min(a), max(a)]
    N['angle_component_range'] = comps

    # --- PCA ---
    N['pca_15'] = pca_block(r15)
    N['pca_20_pre_rib'] = pca_block([r for r in records if r['geometry'] in ('strip', 'slot')])
    N['pca_25_with_rib'] = pca_block([r for r in records if r['geometry'] in ('strip', 'slot', 'rib')])
    N['pca_30'] = pca_block(records)
    N['seven_scalings_15'] = seven_scalings(r15)

    # --- slot ---
    N['slot_noise'] = {str(m): isolated_contribution(records, 'slot', m)['lmin_improvement']
                       for m in (1.0, 3.0, SLOT_MULTIPLIER_SOURCED, 5.0, 10.0)}
    _, unit, _, _, ang, _ = build_matrices(records)
    pairs = [(ang[i, j], label(records[i]), label(records[j]))
             for i in range(len(records)) for j in range(len(records))
             if fam(records[i]) == 'slot-TM' and fam(records[j]) == 'strip-TE']
    a, li, lj = min(pairs)
    N['slot_TM_vs_strip_TE_min_angle'] = dict(angle_deg=a, pair=f'{li} / {lj}',
                                              mean_deg=float(np.mean([p[0] for p in pairs])))

    # --- rib / SWG: PC3-gap projection and isolated contributions ---
    pc3, pc3_var = pc3_pre_rib(records)
    N['pc3_pre_rib'] = dict(loadings=pc3.tolist(), variance_pct=pc3_var)
    proj = {r['id']: float(abs(np.dot(r['sens_vector_unit'], pc3))) for r in records}
    N['pc3_projection'] = proj
    rng = {}
    for fm in FAMS:
        v = [proj[r['id']] for r in records if fam(r) == fm]
        rng[fm] = [min(v), max(v)]
    N['pc3_projection_range'] = rng
    reps = ['strip_TE_500', 'strip_TM_500', 'slot_TE_180', 'slot_TM_180'] + [f'swg_TE_{w}' for w in (400, 450, 500, 550, 600)]
    N['pc3_projection_others_writeup_set'] = [min(proj[k] for k in reps), max(proj[k] for k in reps)]
    others = [proj[r['id']] for r in records if r['geometry'] != 'rib']
    N['pc3_projection_others_all_widths'] = [min(others), max(others)]

    iso = {}
    for g, m in [('rib', RIB_MULTIPLIER_CONSERVATIVE), ('rib', RIB_MULTIPLIER_OPTIMISTIC), ('swg', 1.0), ('slot', 1.0)]:
        x = isolated_contribution(records, g, m)
        iso[f'{g} x{m}'] = dict(lmin_improvement=x['lmin_improvement'], cond_improvement=x['cond_improvement'],
                               lmin0=x['lmin0'], lmin=x['lmin'], cond0=x['cond0'], cond=x['cond'],
                               tinv0=x['tinv0'], tinv=x['tinv'], tinv_improvement=x['tinv0'] / x['tinv'])
    N['isolated'] = iso
    N['rib_vs_strip_500_raw'] = dict(
        rib=[rec['rib_TE_500'][k] for k in ('dneff_dw', 'dneff_dt', 'dneff_da')],
        strip=[rec['strip_TE_500'][k] for k in ('dneff_dw', 'dneff_dt', 'dneff_da')])

    # --- sidewall-angle CRB (angle_crb.py framings) ---
    alone = {}
    for g, m in [('strip', 1.0), ('rib', RIB_MULTIPLIER_CONSERVATIVE), ('rib', RIB_MULTIPLIER_OPTIMISTIC),
                 ('slot', SLOT_MULTIPLIER_SOURCED)]:
        alone[f'{g} x{m}'] = float(marginal_angle_crb(rows_for(records, g, 'TE', m))[3])
    N['crb_alone'] = alone
    strip_rows = rows_for(records, 'strip', 'TE', 1.0) + rows_for(records, 'strip', 'TM', 1.0)
    std0 = marginal_angle_crb(strip_rows)[3]
    crb = dict(strip_baseline=float(std0))
    for g, m in [('rib', RIB_MULTIPLIER_CONSERVATIVE), ('rib', RIB_MULTIPLIER_OPTIMISTIC), ('slot', SLOT_MULTIPLIER_SOURCED)]:
        s = marginal_angle_crb(strip_rows + rows_for(records, g, 'TE', m))[3]
        crb[f'+{g} x{m}'] = dict(std=float(s), improvement=float(std0 / s))
    N['crb_large_budget'] = crb

    def ring(g, pol, w, m=1.0):
        return rebuild_ring_rows(rec[f'{g}_{pol}_{w}'], fv.SIGMA_NEFF * m, fv.SIGMA_NG * m)
    minimal = ring('strip', 'TE', 550) + ring('strip', 'TE', 600)
    Fm, Finv, _, std_min = marginal_angle_crb(minimal)
    ev, evec = np.linalg.eigh(Fm)
    mn = dict(angle_std=float(std_min), marginal_std_wta=np.sqrt(np.diag(Finv)).tolist(),
              worst_eigvec_wta=evec[:, 0].tolist(), worst_eigval=float(ev[0]), cond=float(np.linalg.cond(Fm)))
    for name, (g, w, m) in {'+rib400 x1.0': ('rib', 400, 1.0), '+rib400 x0.2': ('rib', 400, 0.2),
                            '+rib600 x1.0': ('rib', 600, 1.0), '+slot150 x4': ('slot', 150, 4.0),
                            '+strip400 (control)': ('strip', 400, 1.0)}.items():
        s = marginal_angle_crb(minimal + ring(g, 'TE', w, m))[3]
        mn[name] = dict(std=float(s), improvement=float(std_min / s))
    best_single = {}
    for g, m in (('rib', 1.0), ('rib', 0.2), ('slot', 4.0), ('strip', 1.0)):
        ws = sorted({r['width_nm'] for r in records if r['geometry'] == g})
        imps = {w: float(std_min / marginal_angle_crb(minimal + ring(g, 'TE', w, m))[3])
                for w in ws if not (g == 'strip' and w in (550, 600))}
        best_single[f'{g} x{m}'] = imps
    mn['single_addition_all_widths'] = best_single
    N['crb_two_rings'] = mn

    # --- slab-thickness sweep (rib_slab_sweep.py method) ---
    base_rows = strip_baseline_rows(records)
    _, lmin0, _, _ = fisher_metrics(base_rows)
    sweep = {}
    for s in sorted({e['slab_nm'] for e in slab}):
        ents = [e for e in slab if e['slab_nm'] == s]
        rows = list(base_rows)
        for e in ents:
            rows.extend(rebuild_ring_rows(e, fv.SIGMA_NEFF, fv.SIGMA_NG))
            rows.extend(rebuild_mzi_row(e, fv.SIGMA_PHASE))
        _, lmin, _, _ = fisher_metrics(rows)
        sweep[str(s)] = dict(lmin_improvement=lmin / lmin0,
                             pc3_projection={str(e['width_nm']): float(abs(np.dot(unit_from(e), pc3))) for e in ents})
    N['slab_sweep'] = sweep
    return N


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--db', required=True)
    ap.add_argument('--slab', required=True)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    N = main(a.db, a.slab)
    with open(a.out, 'w') as f:
        json.dump(N, f, indent=1, default=float)
    print(json.dumps({k: v for k, v in N.items() if k not in ('pc3_projection',)}, indent=1, default=float))
