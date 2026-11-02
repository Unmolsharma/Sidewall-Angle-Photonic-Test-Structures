"""
Analytic cross-checks of femwell results, from the raw solve caches (no derived files except
results_final/data/crosscheck_fields.json from crosscheck_fields.py).

  1. Width/angle equivalence of a strip: a tilt a about the bottom edge moves the wall at height y by
     y tan(a), so the mirror-symmetric strip sees a width change -t tan(a) (t = 220 nm). Central differences
     at +-h_a give exactly -t tan(h_a)/h_a per degree (= -t pi/180 as h_a -> 0). Compared with the femwell
     ratio (dn/da)/(dn/dw) at step h for every strip (n_eff and n_g).
  2. Closed-form R = (y_c - y_pivot)/(h_wall/2) (y_c from the mode field, crosscheck_fields.py) vs femwell
     R = (dn/da)/(dn/dw) / (-h_wall tan(1 deg)) at step h (1 deg), for strips and ribs.
  3. Impact of ignoring the angle: a 1 deg strip tilt fitted as a width change (width-only weighted least
     squares over strip TE 400-600 n_eff and n_g), used to predict the rib TE 500 nm (90 nm slab); resonance
     error d lambda = lambda * (predicted - true) dn_eff / n_g(rib), lambda = 1550 nm. Noise weights: the
     placeholders (as in the original femwell_impact.py) and the sourced 1 pm level.
  4. Mesh check: derivatives at step h on the r12 and r8 meshes (results_femwell/cache/solves_{r12,r8}.jsonl),
     oxide-clad strip and rib TE 500 nm (the air-clad strip was solved on r8 only).
  5. R of the other mirror-symmetric guides (slot, SWG) and of the air-clad strip TE 500 nm (r8 mesh), same
     definition as 2.

Run from src/ (main environment):  python crosschecks.py
Writes results_final/data/crosschecks.json
"""
import json
import os

import numpy as np

import etchdepth_femwell_analysis as ed
import fullvector_sensitivity as fv
import priors_search as ps

REPO = ed.REPO
OUT = os.path.join(REPO, 'results_final', 'data')
FCACHE = os.path.join(REPO, 'results_femwell', 'cache')
T, LAM_NM, HA = 220.0, 1550.0, 1.0


def load(path):
    rows = {}
    for line in open(path, encoding='utf-8'):
        if line.strip():
            r = json.loads(line)
            rows[r['key']] = r
    return rows


def nng(n):
    return np.array([n[0], n[0] - ed.WL0 * (n[1] - n[2]) / (2 * ed.DLAM)])


def deriv(rows, geom, pol, W, S, clad, mesh, p, f=1.0):
    h = f * ed.H[p]
    k = f'{geom}|{pol}|{W}|{S}|{clad}|{mesh}|{p}|'
    return (nng(rows[k + f'{h:+.4f}']['n']) - nng(rows[k + f'{-h:+.4f}']['n'])) / (2 * h)


def main():
    os.makedirs(OUT, exist_ok=True)
    r12 = load(os.path.join(FCACHE, 'solves_r12.jsonl'))
    r8 = load(os.path.join(FCACHE, 'solves_r8.jsonl'))
    res = {}

    # 1. width/angle equivalence
    fd = T * np.tan(np.radians(HA)) / HA
    strips = {}
    for pol in ('TE', 'TM'):
        for W in (400, 450, 500, 550, 600):
            ratio = deriv(r12, 'strip', pol, W, None, 'oxide', 'r12', 'a') / deriv(r12, 'strip', pol, W, None, 'oxide', 'r12', 'w')
            strips[f'strip_{pol}_{W}'] = dict(neff_nm_per_deg=float(-ratio[0]), ng_nm_per_deg=float(-ratio[1]))
    allv = np.array([[v['neff_nm_per_deg'], v['ng_nm_per_deg']] for v in strips.values()])
    res['width_angle_equivalence'] = dict(
        analytic_t_tan_ha_over_ha_nm_per_deg=float(fd), analytic_t_pi_over_180_nm_per_deg=float(T * np.pi / 180),
        femwell=strips, femwell_min=float(allv.min()), femwell_max=float(allv.max()),
        max_rel_dev_from_analytic=float(np.max(np.abs(allv / fd - 1))),
        median_rel_dev_from_analytic=float(np.median(np.abs(allv / fd - 1))))

    # 2. closed-form R vs femwell R
    fields = json.load(open(os.path.join(OUT, 'crosscheck_fields.json')))
    comp = {}
    for uid, e in fields.items():
        g, pol, W = uid.split('_')[:3]
        S = int(uid.split('slab')[1]) if 'slab' in uid else None
        hw = (e['wall_height_um']) * 1000
        rows = r12
        ratio = deriv(rows, g, pol, int(W), S, 'oxide', 'r12', 'a') / deriv(rows, g, pol, int(W), S, 'oxide', 'r12', 'w')
        Rf = -ratio / (hw * np.tan(np.radians(HA)))
        rc = e['with_Ez']['R_closed']
        comp[uid] = dict(R_closed=rc, R_closed_without_Ez=e['without_Ez']['R_closed'], R_femwell_neff=float(Rf[0]),
                         R_femwell_ng=float(Rf[1]), rel_diff_neff=float(rc / Rf[0] - 1))
    rib = [v['rel_diff_neff'] for k, v in comp.items() if k.startswith('rib')]
    stp = [v['rel_diff_neff'] for k, v in comp.items() if k.startswith('strip')]
    res['closed_form_R'] = dict(units=comp, strip_R_exact=1.0,
                                strip_max_abs_rel_diff=float(np.max(np.abs(stp))),
                                rib_max_abs_rel_diff=float(np.max(np.abs(rib))),
                                rib_median_abs_rel_diff=float(np.median(np.abs(rib))))

    # 3. impact of ignoring the angle (1 deg)
    import noise_sourced
    imp = {}
    for label, sig in (('placeholder', (fv.SIGMA_NEFF, fv.SIGMA_NG)),
                       ('1pm', tuple(noise_sourced.sigmas(1e-12)[k] for k in ('SIGMA_NEFF', 'SIGMA_NG')))):
        Jw, Ja = [], []
        for W in (400, 450, 500, 550, 600):
            Jw.append(deriv(r12, 'strip', 'TE', W, None, 'oxide', 'r12', 'w') / np.array(sig))
            Ja.append(deriv(r12, 'strip', 'TE', W, None, 'oxide', 'r12', 'a') / np.array(sig))
        jw, ja = np.concatenate(Jw), np.concatenate(Ja)
        dw = float(jw @ ja / (jw @ jw))                              # width-only fit of a 1 deg tilt (nm)
        rw = deriv(r12, 'rib', 'TE', 500, None, 'oxide', 'r12', 'w')
        ra = deriv(r12, 'rib', 'TE', 500, None, 'oxide', 'r12', 'a')
        ng_rib = float(nng(r12['rib|TE|500|None|oxide|r12|nom|+0.0000']['n'])[1])
        err = rw * dw - ra
        imp[label] = dict(effective_dw_nm=dw, rib_dneff_error=float(err[0]), rib_dng_error=float(err[1]), n_g_rib=ng_rib,
                          resonance_error_pm=float(LAM_NM * err[0] / ng_rib * 1000))
    res['impact_1deg'] = dict(formula='d lambda = lambda * dn_eff_error / n_g(rib), lambda = 1550 nm', **imp)

    # 4. mesh check r12 vs r8
    mc = {}
    for (g, S, clad, params) in (('strip', None, 'oxide', 'wta'), ('rib', None, 'oxide', 'wta')):   # air-clad: r8 only
        for p in params:
            a, b = deriv(r12, g, 'TE', 500, S, clad, 'r12', p), deriv(r8, g, 'TE', 500, S, clad, 'r8', p)
            mc[f'{g}_{clad}_{p}'] = dict(r12=a.tolist(), r8=b.tolist(), rel_diff=(a / b - 1).tolist())
    rel = np.array([v['rel_diff'] for v in mc.values()])
    res['mesh_check'] = dict(derivatives=mc, max_abs_rel_diff=float(np.max(np.abs(rel))),
                             max_abs_rel_diff_neff=float(np.max(np.abs(rel[:, 0]))))

    # 5. R of slot, SWG, air-clad strip
    other = {}
    for g, pols, Ws in (('slot', ('TE', 'TM'), (150, 165, 180, 195, 210)), ('swg', ('TE',), (400, 450, 500, 550, 600))):
        for pol in pols:
            for W in Ws:
                ratio = deriv(r12, g, pol, W, None, 'oxide', 'r12', 'a') / deriv(r12, g, pol, W, None, 'oxide', 'r12', 'w')
                other[f'{g}_{pol}_{W}'] = (-ratio / (T * np.tan(np.radians(HA)))).tolist()
    ratio = deriv(r8, 'strip', 'TE', 500, None, 'air', 'r8', 'a') / deriv(r8, 'strip', 'TE', 500, None, 'air', 'r8', 'w')
    slot = np.array([v for k, v in other.items() if k.startswith('slot')])
    swg = np.array([v for k, v in other.items() if k.startswith('swg')])
    res['other_R'] = dict(units=other, slot_min=float(slot.min()), slot_max=float(slot.max()),
                          swg_min=float(swg.min()), swg_max=float(swg.max()),
                          strip_TE_500_air_r8_mesh=(-ratio / (T * np.tan(np.radians(HA)))).tolist())
    with open(os.path.join(OUT, 'crosschecks.json'), 'w') as f:
        json.dump(res, f, indent=1)
    w = res['width_angle_equivalence']
    print(f"1. t tan(1)/1 = {w['analytic_t_tan_ha_over_ha_nm_per_deg']:.5f} nm/deg; femwell {w['femwell_min']:.4f}-{w['femwell_max']:.4f}"
          f" (max rel dev {w['max_rel_dev_from_analytic']:.2%})")
    c = res['closed_form_R']
    print(f"2. closed-form R vs femwell: strips max {c['strip_max_abs_rel_diff']:.2%}, ribs max {c['rib_max_abs_rel_diff']:.2%},"
          f" median {c['rib_median_abs_rel_diff']:.2%}")
    for k, v in comp.items():
        print(f"   {k:20s} closed {v['R_closed']:.4f}  femwell {v['R_femwell_neff']:.4f}  ({v['rel_diff_neff']:+.2%})")
    for k, v in imp.items():
        print(f"3. impact ({k}): dw {v['effective_dw_nm']:.4f} nm, resonance error {v['resonance_error_pm']:.1f} pm")
    print(f"4. mesh r12 vs r8: max |rel diff| {res['mesh_check']['max_abs_rel_diff']:.3%}")
    o = res['other_R']
    print(f"5. slot R {o['slot_min']:.4f}-{o['slot_max']:.4f}, SWG {o['swg_min']:.4f}-{o['swg_max']:.4f}, air strip (r8) {o['strip_TE_500_air_r8_mesh']}")


if __name__ == '__main__':
    main()
