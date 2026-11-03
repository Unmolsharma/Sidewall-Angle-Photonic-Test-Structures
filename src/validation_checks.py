"""
Validation of the femwell model, from committed solves and the cited literature value:
  1. femwell n_eff and n_g at the nominal geometry vs the committed EMpy (cell-exact) database
     results_corrected/data/sensitivity_db.json, for every structure in both (EMpy is not rerun).
  2. femwell n_eff vs the published effective-index-method (EIM) values of Jiao & Pikus (rib 500 nm / 90 nm slab
     and strip 500 nm, 220 nm, n_Si 3.48, n_ox 1.44), read from results_final/data/citations.json. EIM is an
     approximation, so this is a sanity check, not a reference.
  3. Domain size: n_eff, n_g, the geometric derivatives at step h and R of the two-width headline structures on the
     1.5x domain (results_final/domain_check, src/domain_check.py) vs the standard domain, and the two-width CRB
     (stage D, realistic priors, 1 pm) with the large-domain geometric derivatives in place of the standard ones.

Run from src/ (main environment):  python validation_checks.py
Writes results_final/data/validation.json
"""
import json
import os
import re

import numpy as np

import etchdepth_femwell_analysis as ed
import montecarlo
import priors_search as ps

REPO = ed.REPO
OUT = os.path.join(REPO, 'results_final', 'data')
DOM = os.path.join(REPO, 'results_final', 'domain_check', 'cache', 'solves_r12.jsonl')
TWO = [f'strip_TE_{w}' for w in (400, 450, 500, 550, 600)] + ['rib_TE_300_slab50', 'rib_TE_500_slab50']


def nom_key(r):
    g, pol, w = r['geometry'], r['polarization'], r['width_nm']
    return f'{g}|{pol}|{w}|None|oxide|r12|nom|+0.0000'


def main():
    os.makedirs(OUT, exist_ok=True)
    rows = ed.load()
    res = {}

    # 1. femwell vs EMpy n_eff, n_g
    E = json.load(open(os.path.join(REPO, 'results_corrected', 'data', 'sensitivity_db.json')))
    cmp_ = {}
    for r in E:
        k = nom_key(r)
        if k in rows:
            n = ed.nng(rows, (r['geometry'], r['polarization'], r['width_nm'], None), 'nom', 0.0)
            cmp_[r['id']] = dict(femwell_neff=float(n[0]), empy_neff=r['neff'], rel_diff_neff=r['neff'] / n[0] - 1,
                                 femwell_ng=float(n[1]), empy_ng=r['ng'], rel_diff_ng=r['ng'] / n[1] - 1)
    by = {}
    for g in ('strip', 'rib', 'slot', 'swg'):
        v = [x for k, x in cmp_.items() if k.startswith(g)]
        if v:
            by[g] = dict(n=len(v), max_abs_rel_diff_neff=max(abs(x['rel_diff_neff']) for x in v),
                         median_abs_rel_diff_neff=float(np.median([abs(x['rel_diff_neff']) for x in v])),
                         max_abs_rel_diff_ng=max(abs(x['rel_diff_ng']) for x in v))
    res['femwell_vs_empy'] = dict(structures=cmp_, by_geometry=by)

    # 2. femwell vs published EIM
    q = json.load(open(os.path.join(OUT, 'citations.json'), encoding='utf-8'))['sources']['jiao_pikus2025']['quotes'][0]['text']
    rib_eim = float(re.search(r'neff = ([\d.]+)\) in the rib', q).group(1))
    strip_eim = float(re.search(r'strip waveguide \(neff = ([\d.]+)\)', q).group(1))
    rib_fw = float(ed.nng(rows, ('rib', 'TE', 500, 90), 'nom', 0.0)[0])
    strip_fw = float(ed.nng(rows, ('strip', 'TE', 500, None), 'nom', 0.0)[0])
    res['femwell_vs_published_eim'] = dict(
        rib_500_slab90=dict(eim=rib_eim, femwell=rib_fw, rel_diff=rib_fw / rib_eim - 1),
        strip_500=dict(eim=strip_eim, femwell=strip_fw, rel_diff=strip_fw / strip_eim - 1),
        note='EIM values are the chapter authors\' own calculation (Section 2.3); not a measurement')

    # 3. domain size
    if os.path.exists(DOM):
        dom = {}
        for line in open(DOM, encoding='utf-8'):
            if line.strip():
                r = json.loads(line)
                dom[r['key'].split('|', 2)[2]] = r
        unit = {}
        for uid in TWO:
            g, pol, W = uid.split('_')[:3]
            S = int(uid.split('slab')[1]) if 'slab' in uid else None
            u = (g, pol, int(W), S)
            base = f'{g}|{pol}|{W}|{S}|oxide|r12|'
            if base + 'nom|+0.0000' not in dom:
                continue
            n_std = ed.nng(rows, u, 'nom', 0.0)
            n_dom = ed.nng(dom, u, 'nom', 0.0)
            d = {}
            for p in ('w', 't', 'a') + (('s',) if g == 'rib' else ()):
                h = ed.H[p]
                try:
                    a = (ed.nng(rows, u, p, h) - ed.nng(rows, u, p, -h)) / (2 * h)
                    b = (ed.nng(dom, u, p, h) - ed.nng(dom, u, p, -h)) / (2 * h)
                except KeyError:
                    continue
                d[p] = dict(standard=a.tolist(), large=b.tolist(), rel_diff=(b / a - 1).tolist())
            hw = 220.0 - (S or 0)
            Rs = (np.array(d['a']['standard']) / np.array(d['w']['standard'])) / (-hw * np.tan(np.radians(1.0)))
            Rl = (np.array(d['a']['large']) / np.array(d['w']['large'])) / (-hw * np.tan(np.radians(1.0)))
            unit[uid] = dict(n_eff_standard=float(n_std[0]), n_eff_large=float(n_dom[0]), rel_diff_neff=float(n_dom[0] / n_std[0] - 1),
                             n_g_standard=float(n_std[1]), n_g_large=float(n_dom[1]), rel_diff_ng=float(n_dom[1] / n_std[1] - 1),
                             derivatives=d, R_standard=Rs.tolist(), R_large=Rl.tolist())
        res['domain'] = dict(scale=1.5, standard_um=[1.4, 1.2], large_um=[2.1, 1.8], units=unit)
        if len(unit) == len(TWO):
            alld = [abs(x) for v in unit.values() for dd in v['derivatives'].values() for x in dd['rel_diff']]
            res['domain'].update(max_abs_rel_diff_neff=max(abs(v['rel_diff_neff']) for v in unit.values()),
                                 max_abs_rel_diff_ng=max(abs(v['rel_diff_ng']) for v in unit.values()),
                                 max_abs_rel_diff_derivative=max(alld))
            # two-width CRB with the large-domain geometric derivatives (factor 1.0); index derivatives unchanged
            mixed = dict(rows)
            for k, r in dom.items():
                mixed[k] = r
            montecarlo.set_noise('1pm')
            units = [(u, ('ring', 'mzi')) for u in TWO]
            a_std = ps.evaluate(units, ps.derivatives(rows), 'D', ps.PRIORS['realistic'], factors=(1.0,))[1.0]['a']
            a_dom = ps.evaluate(units, ps.derivatives(mixed), 'D', ps.PRIORS['realistic'], factors=(1.0,))[1.0]['a']
            res['domain']['two_width_crb'] = dict(standard=a_std, large_domain=a_dom, rel_change=a_dom / a_std - 1,
                                                  noise='1pm', stage='D', priors='realistic', slab_nm=50)
    with open(os.path.join(OUT, 'validation.json'), 'w') as f:
        json.dump(res, f, indent=1)
    print('femwell vs EMpy:', json.dumps(by, indent=1))
    print('femwell vs EIM:', res['femwell_vs_published_eim'])
    if 'domain' in res:
        print('domain:', {k: v for k, v in res['domain'].items() if k != 'units'})


if __name__ == '__main__':
    main()
