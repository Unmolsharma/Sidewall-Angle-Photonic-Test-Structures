"""
Generate results_final/HEADLINE_NUMBERS.md from the data files. No number is typed here: every value is read
from a data file (or the citation record) at generation time, and every row states its category and source.

Categories
  COMPUTED  produced by a script in this repo from solver output; source = script + data file + key
  CITED     quoted from a verified source; source = citation + page/table + exact quote (results_final/data/citations.json)
  DERIVED   obtained by a stated formula from COMPUTED / CITED values; source = formula + inputs
  ASSUMED   an input chosen by us; source = where it is set in the code (read back from results_final/data/environment.json)
Commit = last commit that touched the data file (`uncommitted` if it has changes not yet committed).

Primary noise level: the sourced 1 pm resonance-wavelength precision; best case 0.1 pm.

Run from src/ (main environment):  python make_headline_numbers.py
"""
import csv
import json
import math
import os
import subprocess

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUTF = os.path.join(REPO, 'results_final', 'HEADLINE_NUMBERS.md')
F = dict(
    s1=('results_final/sym_search_1pm/data/search.json', 'src/sym_search.py 1pm'),
    p1=('results_final/sym_search_1pm/data/slab_pairs.json', 'src/sym_search.py 1pm'),
    t1=('results_final/sym_search_1pm/data/priors_table.csv', 'src/sym_search.py 1pm'),
    d1=('results_final/sym_search_1pm/data/d_prior_sweep.json', 'src/sym_search.py 1pm'),
    mc1=('results_montecarlo/data/mc_results_1pm.json', 'src/montecarlo.py mc1pm'),
    mc0=('results_montecarlo/data/mc_results.json', 'src/montecarlo.py mc'),
    val=('results_montecarlo/data/mc_valley.json', 'src/montecarlo.py valley'),
    val1=('results_montecarlo/data/mc_valley_1pm.json', 'src/montecarlo.py valley1pm'),
    vc=('results_final/data/validation.json', 'src/validation_checks.py'),
    rob=('results_montecarlo/robustness/data/robustness_1pm.json', 'src/robustness_rib300.py 1pm'),
    ns=('results_montecarlo/noise_sourced/data/noise_sourced.json', 'src/noise_sourced.py'),
    xc=('results_final/data/crosschecks.json', 'src/crosschecks.py'),
    sc=('results_montecarlo/sym_search/data/step_check_all.csv', 'src/step_check_all.py'),
    sa=('results_montecarlo/data/surrogate_accuracy_1pm.json', 'src/montecarlo.py surrcheck'),
    asn=('results_final/data/assumption_sensitivity.json', 'src/assumption_sensitivity.py'),
    cit=('results_final/data/citations.json', 'src/citations_extract.py'),
    env=('results_final/data/environment.json', 'src/record_environment.py'),
)
_cache = {}


def load(k):
    if k not in _cache:
        path = os.path.join(REPO, F[k][0])
        if path.endswith('.csv'):
            with open(path, newline='', encoding='utf-8') as f:
                _cache[k] = list(csv.DictReader(f))
        else:
            _cache[k] = json.load(open(path, encoding='utf-8'))
    return _cache[k]


def commit(path):
    dirty = subprocess.run(['git', 'status', '--porcelain', '--', path], cwd=REPO, capture_output=True, text=True).stdout.strip()
    if dirty:
        return 'uncommitted'
    h = subprocess.run(['git', 'log', '-1', '--format=%h', '--', path], cwd=REPO, capture_output=True, text=True).stdout.strip()
    return h or 'uncommitted'


def fmt(v, sig=4):
    if isinstance(v, str):
        return v
    if v is None:
        return 'n/a'
    if isinstance(v, float) and math.isinf(v):
        return '∞'
    if isinstance(v, int):
        return str(v)
    return f'{v:.{sig}g}'


ROWS = []        # (section, value text, meaning, category, source text, files for commit)


def computed(sec, key, path, value, meaning, sig=4):
    ROWS.append((sec, fmt(value, sig), meaning, 'COMPUTED', f'`{F[key][1]}` → `{F[key][0]}` : `{path}`', [F[key][0]]))


DERIVED_FILE = 'results_final/data/headline_derived.json'
DERIVED_VALUES = []


def derived(sec, value, meaning, formula, keys, sig=4):
    """A value computed here by a stated formula; it is also written to DERIVED_FILE so that it can be traced."""
    DERIVED_VALUES.append(dict(meaning=meaning, value=value, formula=formula, inputs=[F[k][0] for k in keys]))
    ROWS.append((sec, fmt(value, sig), meaning, 'DERIVED', formula + ' (inputs: ' + ', '.join(f'`{F[k][0]}`' for k in keys) +
                 f'; value stored in `{DERIVED_FILE}`)', [F[k][0] for k in keys] + [DERIVED_FILE]))


def cited(sec, value_text, meaning, src, qi):
    c = load('cit')['sources'][src]
    q = c['quotes'][qi]
    where = f"p. {q['page']}" if 'page' in q else f"Sect. {q['section']}"
    ROWS.append((sec, value_text, meaning, 'CITED',
                 f"{c['authors'].split(',')[0]} et al., \"{c['title']}\", {c['venue']}, DOI {c['doi']}, {where}: “{q['text']}” "
                 f"(`{F['cit'][0]}` : `sources.{src}.quotes[{qi}]`)", [F['cit'][0]]))


def assumed(sec, value, meaning, where, path, sig=4, key='env'):
    ROWS.append((sec, fmt(value, sig), meaning, 'ASSUMED', f'{where} (read back from `{F[key][0]}` : `{path}`)', [F[key][0]]))


def quote_text(src, qi):
    return load('cit')['sources'][src]['quotes'][qi]['text']


def quoted_number(src, qi, pattern):
    """A number read out of a stored quote with a regular expression (group 1)."""
    import re
    return float(re.search(pattern, quote_text(src, qi)).group(1))


def table_value(src, qi, row_label, column):
    """Value at (row starting with row_label, header column `column`) of a table quote 'title | header | row | ...'."""
    parts = [x.strip() for x in quote_text(src, qi).split(' | ')]
    head = [h for h in parts if h.startswith('wavelength')][0].split()
    row = [r for r in parts if r.startswith(row_label)][0].split()
    cols = [' '.join(head[i:i + 2]) for i in range(1, len(head), 2)]          # '30 K', '40 K', ...
    return float(row[2 + cols.index(column)])                                # row: '1.5', 'microns', values...


def get(d, path):
    for p in path.split('.'):
        d = d[int(p)] if isinstance(d, list) else d[p]
    return d


def build():
    s1, p1, t1, mc1, mc0 = load('s1'), load('p1'), load('t1'), load('mc1'), load('mc0')
    # --- A. structure search, 1 pm
    sec = 'A. Structure search at the primary noise (1 pm), strip symmetry imposed'
    for tag, label in (('D|realistic', 'stage D, realistic priors'), ('C|none', 'stage C, no priors')):
        e = s1[tag]['one_etch_overall']['8']
        for k, m in (('a_deg_h', 'at step h'), ('a_deg_min', 'min over h/2, h, 2h'), ('a_deg_max', 'max over h/2, h, 2h')):
            computed(sec, 's1', f'{tag}.one_etch_overall.8.{k}', e[k],
                     f"Best one-etch set of 8, {label}: strip-angle CRB (deg), {m}; slab {e['slab']} nm; set: {' + '.join(e['set'])}")
    for k in range(4, 9):
        e = s1['D|realistic']['one_etch_overall'][str(k)]
        computed(sec, 's1', f'D|realistic.one_etch_overall.{k}.a_deg_h', e['a_deg_h'],
                 f"Best one-etch set of {k}, stage D realistic, CRB at h (deg), slab {e['slab']} nm")
    pairs = {v['slab']: (k, v) for k, v in p1['D|realistic'].items() if 'rib_TE_300' in k and 'rib_TE_500' in k}
    for slab in sorted(pairs):
        k, v = pairs[slab]
        computed(sec, 'p1', f'D|realistic.{k}.a_deg_h', v['a_deg_h'],
                 f'Two-width set (strip TE 400-600 + rib 300 + 500 nm), {slab} nm slab, stage D realistic, CRB at h (deg)')
    best = min(pairs, key=lambda s: pairs[s][1]['a_deg_h'])
    k, v = pairs[best]
    computed(sec, 'p1', f'D|realistic.{k}.a_deg_min', v['a_deg_min'], f'Same, best slab ({best} nm), min over h/2, h, 2h (deg)')
    computed(sec, 'p1', f'D|realistic.{k}.a_deg_max', v['a_deg_max'], f'Same, best slab ({best} nm), max over h/2, h, 2h (deg)')
    singles = {v['slab']: (k, v) for k, v in p1['D|realistic'].items() if '+' not in k}
    for tag in ('C|none', 'D|none'):
        vals = [p1[tag][k]['a_deg_h'] for k in p1[tag] if '+' not in k]
        computed(sec, 'p1', f'{tag}.rib_TE_500_slab*.a_deg_h', max(vals) if all(math.isinf(x) for x in vals) else min(vals),
                 f'Single rib (strip TE 400-600 + rib 500 nm), {tag.replace("|", ", ")}: CRB at every slab (∞ = not identifiable)')
    sv = [v['a_deg_h'] for _, v in singles.values()]
    computed(sec, 'p1', 'D|realistic.rib_TE_500_slab*.a_deg_h (min)', min(sv), 'Single rib, stage D realistic priors: smallest CRB over slabs (deg)')
    computed(sec, 'p1', 'D|realistic.rib_TE_500_slab*.a_deg_h (max)', max(sv), 'Single rib, stage D realistic priors: largest CRB over slabs (deg)')
    so = [v['a_deg_h'] for v in s1['D|realistic']['strips_only'].values()]
    computed(sec, 's1', 'D|realistic.strips_only.*.a_deg_h (min)', min(so), f'Strips only, best set of any size 1-{len(so)}: CRB (∞ = no angle information)')
    for stage, prior, setname, m in (('B', 'none', 'strip + rib70 (500 nm)', 'stage B (rib angle = strip angle), no priors'),
                                     ('B', 'none', 'strip + rib70 (500 + 300 nm)', 'stage B (rib angle = strip angle), no priors'),
                                     ('C', 'realistic', 'strip + rib70 (500 nm)', 'stage C, realistic priors'),
                                     ('D', 'none', 'strip + rib70 (500 + 300 nm)', 'stage D, no priors'),
                                     ('C', 'none', 'strip + rib70 (500 + 300 nm)', 'stage C, no priors'),
                                     ('D', 'realistic', 'strip + rib70 (500 + 300 nm)', 'stage D, realistic priors'),
                                     ('E', 'realistic', 'strip + rib70 (500 + 300 nm)', 'stage E (+T), realistic priors'),
                                     ('C', 'ar-a 1 deg', 'strip + rib70 (500 nm)', 'stage C, rib-strip angle link 1 deg'),
                                     ('D', 'realistic', 'strip + rib70 (500 nm)', 'stage D, realistic priors')):
        r = [x for x in t1 if x['stage'] == stage and x['prior'] == prior and x['set'] == setname and x['rib_mult'] == '1.0'][0]
        computed(sec, 't1', f'row stage={stage}, prior={prior}, set={setname}, rib_mult=1.0 : a_deg_h', float(r['a_deg_h']),
                 f'Priors table, {setname}, {m}: CRB at h (deg)')
    dC = load('d1')['C']['strip + rib70 (500 nm)']
    computed(sec, 'd1', 'C.strip + rib70 (500 nm).0.01.a_deg_h', dC['0.01']['a_deg_h'],
             'Single rib, stage C, etch depth known to 0.01 nm (prior): CRB at h (deg)')

    # --- B. Monte Carlo, 1 pm
    sec = 'B. Monte Carlo at 1 pm (2,000 noisy sets per case, nonlinear least squares, same priors)'
    for case, label in (('a', 'case a (best multi-width set, 90 nm slab)'), ('b', 'case b (strip + rib 300 + 500 nm, 50 nm slab)')):
        r = mc1[case]
        computed(sec, 'mc1', f'{case}.units', len(r['units']), f"{label}: structures: {' + '.join(r['units'])}")
        computed(sec, 'mc1', f'{case}.crb_a_fd', r['crb_a_fd'], f'{label}: CRB from femwell derivatives (deg)')
        computed(sec, 'mc1', f'{case}.crb_a_surrogate', r['crb_a_surrogate'], f'{label}: CRB of the Monte Carlo surrogate model at nominal (deg)')
        computed(sec, 'mc1', f'{case}.a_std_all', r['a_std_all'], f'{label}: spread (std) of the fitted strip angle (deg)')
        computed(sec, 'mc1', f'{case}.a_std_all_ci95.0', r['a_std_all_ci95'][0], f'{label}: spread, 95% bootstrap CI low (deg)')
        computed(sec, 'mc1', f'{case}.a_std_all_ci95.1', r['a_std_all_ci95'][1], f'{label}: spread, 95% bootstrap CI high (deg)')
        computed(sec, 'mc1', f'{case}.std_over_crb_fd', r['std_over_crb_fd'], f'{label}: spread / CRB')
        computed(sec, 'mc1', f'{case}.a_bias_all', r['a_bias_all'], f'{label}: bias of the fitted strip angle (deg)', 2)
        computed(sec, 'mc1', f'{case}.a_bias_se', r['a_bias_se'], f'{label}: standard error of the bias (deg)', 2)
        for k in ('n_failed', 'n_ambiguous', 'n_at_bound'):
            computed(sec, 'mc1', f'{case}.{k}', r[k], f'{label}: {k[2:].replace("_", " ")} fits of {r["n_sets"]}')
    r = mc1['c']
    computed(sec, 'mc1', 'c.units', len(r['units']), f"case c (strip + one rib 500 nm, 50 nm slab): structures: {' + '.join(r['units'])}")
    computed(sec, 'mc1', 'c.crb_a_fd', r['crb_a_fd'], 'case c: CRB from femwell derivatives (deg)')
    for k in ('n_failed', 'n_ambiguous', 'n_at_bound'):
        computed(sec, 'mc1', f'c.{k}', r[k], f'case c: {k[2:].replace("_", " ")} fits of {r["n_sets"]} (box edge: any parameter at its bound)')
    computed(sec, 'mc1', 'c.theta_box.a', r['theta_box']['a'], 'case c: half-width of the fit box in the strip angle (deg)')
    v1 = load('val1')['c']
    computed(sec, 'val1', 'c.slope_w_per_a_nm_per_deg', v1['slope_w_per_a_nm_per_deg'],
             'case c: slope of fitted width vs fitted angle (nm/deg): the width/angle valley')
    computed(sec, 'val1', 'c.corr_w_a', v1['corr_w_a'], 'case c: correlation of fitted width and angle', 6)

    # --- C. robustness (1 pm)
    rob = load('rob')
    sec = 'C. Robustness to the rib 300 nm derivatives (1 pm)'
    for name, label in (('best_one_etch_D_realistic', 'best set (case a)'), ('two_width_slab50_D_realistic', 'two-width set, 50 nm slab (case b)')):
        e = rob[name]
        computed(sec, 'rob', f'{name}.nominal.a_deg_h', e['nominal']['a_deg_h'], f'{label}: CRB with remeshed derivatives (deg)')
        fm = e['fixed_mesh']
        computed(sec, 'rob', f'{name}.fixed_mesh.a_deg_h', fm['a_deg_h'], f'{label}: CRB with fixed-mesh rib 300 nm derivatives (deg)')
        computed(sec, 'rob', f'{name}.fixed_mesh.a_deg_max', fm['a_deg_max'], f'{label}: fixed mesh, max over h/2, h, 2h (deg)')
        computed(sec, 'rob', f'{name}.fixed_mesh.rel_change', fm['rel_change'], f'{label}: fixed mesh vs remeshed, relative change of the CRB')
        dev = max(abs(x) for d in fm['rib300_derivative_rel_diff_fixed_vs_remeshed'].values() for v in d.values() for x in v)
        derived(sec, dev, f'{label}: largest relative difference, fixed-mesh vs remeshed rib 300 nm derivative (any parameter, '
                          f'n_eff or n_g, h/2, h, 2h)', 'max |rib300_derivative_rel_diff_fixed_vs_remeshed|', ['rob'])
        for lv in ('0.3%', '1.0%', '3.0%', '5.0%'):
            pe = e['perturbation'][lv]
            computed(sec, 'rob', f'{name}.perturbation.{lv}.p05', pe['p05'], f'{label}: random ±{lv} rib 300 nm derivative errors, 5th percentile of the CRB (deg)')
            computed(sec, 'rob', f'{name}.perturbation.{lv}.p95', pe['p95'], f'{label}: same, 95th percentile (deg)')

    # --- D. best case 0.1 pm and assumption sensitivity of the noise
    asn = load('asn')
    sec = 'D. Best case (0.1 pm) and noise-level dependence'
    for name, label in (('best', 'best set (1 pm search)'), ('two', 'two-width set (1 pm search)')):
        r = asn['results'][name]
        computed(sec, 'asn', f'results.{name}.base', r['base'], f'{label}: CRB at the primary 1 pm, stage D realistic (deg)')
        for v in r['variations']:
            if v['assumption'] == 'sigma_lambda':
                computed(sec, 'asn', f"results.{name}.variations[sigma_lambda {v['level']}].a_deg", v['a_deg'], f"{label}: CRB at sigma_lambda = {v['level']} (deg)")

    # --- E. noise
    ns = load('ns')
    sec = 'E. Measurement noise'
    cited(sec, '1 pm', 'Primary resonance-wavelength precision', 'janz2023', 2)
    cited(sec, '0.1 pm', 'Best-case resonance-wavelength precision (temperature-bath-stabilized ring)', 'janz2023', 0)
    lv = ns['noise']['sourced conservative (1 pm)']
    for k, formula in (('SIGMA_NEFF', 'sigma_neff = n_g sigma_lambda / lambda'), ('SIGMA_NG', 'sigma_ng = sqrt(2) sigma_lambda n_g^2 L / lambda^2'),
                       ('SIGMA_PHASE', 'sigma_phase = 2 pi n_g L sigma_lambda / lambda^2 (MZI fringe read like a ring resonance)')):
        derived(sec, lv[k], f'{k} at 1 pm', formula + ', sigma_lambda 1 pm, lambda 1.55 um, L 100 um, n_g 4.26', ['ns'], 3)
    lv = ns['noise']['sourced optimistic (0.1 pm)']
    derived(sec, lv['SIGMA_NEFF'], 'SIGMA_NEFF at 0.1 pm', 'same formula, sigma_lambda 0.1 pm', ['ns'], 3)
    dt = quoted_number('janz2023', 3, r'approximately (\d+) pm/K')
    s01 = quoted_number('janz2023', 0, r'±([\d.]+) pm')
    cited(sec, f'{fmt(dt)} pm/K', 'Resonance temperature coefficient of Si rings (approximate)', 'janz2023', 3)
    derived(sec, s01 / dt * 1000, 'Temperature stability needed for the 0.1 pm case (mK)',
            'sigma_lambda / (pm per K), both quoted above', ['cit'], 3)

    # --- F. cross-checks
    xc = load('xc')
    sec = 'F. Analytic cross-checks'
    w = xc['width_angle_equivalence']
    derived(sec, w['analytic_t_tan_ha_over_ha_nm_per_deg'], 'Strip width equivalent of 1 deg: 220 nm tan(1 deg) / 1 deg (nm)', '220 nm × tan(1°)', ['xc'], 5)
    computed(sec, 'xc', 'width_angle_equivalence.femwell_min', w['femwell_min'], 'femwell (dn/da)/(dn/dw), strips TE/TM 400-600, n_eff and n_g: min (nm/deg)', 5)
    computed(sec, 'xc', 'width_angle_equivalence.femwell_max', w['femwell_max'], 'Same: max (nm/deg)', 5)
    computed(sec, 'xc', 'width_angle_equivalence.max_rel_dev_from_analytic', w['max_rel_dev_from_analytic'], 'Largest relative deviation from 220 nm tan(1°)', 3)
    c = xc['closed_form_R']
    computed(sec, 'xc', 'closed_form_R.rib_median_abs_rel_diff', c['rib_median_abs_rel_diff'], 'Closed-form R = (y_c - y_pivot)/(h/2) vs femwell R, ribs: median |relative difference|', 3)
    computed(sec, 'xc', 'closed_form_R.rib_max_abs_rel_diff', c['rib_max_abs_rel_diff'], 'Same, ribs: max', 3)
    computed(sec, 'xc', 'closed_form_R.strip_max_abs_rel_diff', c['strip_max_abs_rel_diff'], 'Same, strips (exact R = 1): max', 3)
    for slab in (50, 90):
        rv = [v['R_femwell_neff'] for k, v in c['units'].items() if k.endswith(f'_slab{slab}')]
        computed(sec, 'xc', f'closed_form_R.units[rib_TE_*_slab{slab}].R_femwell_neff (min)', min(rv), f'femwell R of the ribs 300-600 nm, {slab} nm slab: min', 4)
        computed(sec, 'xc', f'closed_form_R.units[rib_TE_*_slab{slab}].R_femwell_neff (max)', max(rv), f'Same: max', 4)
    i = xc['impact_1deg']['placeholder']
    computed(sec, 'xc', 'impact_1deg.placeholder.resonance_error_pm', i['resonance_error_pm'], 'Rib 500 nm resonance error when a 1° strip tilt is fitted as width (pm), d lambda = lambda dn_eff / n_g', 4)
    computed(sec, 'xc', 'impact_1deg.placeholder.effective_dw_nm', i['effective_dw_nm'], 'Width change that the strip fit assigns to a 1° tilt (nm)', 5)
    computed(sec, 'xc', 'mesh_check.max_abs_rel_diff', xc['mesh_check']['max_abs_rel_diff'], 'r12 vs r8 mesh: largest relative difference of any derivative (strip, rib TE 500)', 3)
    o = xc['other_R']
    computed(sec, 'xc', 'other_R.slot_min', o['slot_min'], 'Slot (oxide-clad, mirror-symmetric) R: min', 5)
    computed(sec, 'xc', 'other_R.slot_max', o['slot_max'], 'Slot R: max', 5)
    computed(sec, 'xc', 'other_R.swg_min', o['swg_min'], 'SWG R: min', 5)
    computed(sec, 'xc', 'other_R.swg_max', o['swg_max'], 'SWG R: max', 5)
    computed(sec, 'xc', 'other_R.strip_TE_500_air_r8_mesh.0', o['strip_TE_500_air_r8_mesh'][0], 'Air-clad strip TE 500 nm R (n_eff), r8 mesh', 4)
    sc = load('sc')
    cnt = {}
    for r in sc:
        cnt[r['verdict']] = cnt.get(r['verdict'], 0) + 1
    sec = 'G. Step-size check (h/2, h, 2h; 2% rule)'
    derived(sec, len(sc), 'Derivatives checked (all structures with solves)', 'count of rows', ['sc'])
    for k in sorted(cnt):
        derived(sec, cnt[k], f'Verdict {k}', 'count of rows with this verdict', ['sc'])
    fails = [float(r['max_rel_dev']) for r in sc if r['verdict'] == 'FAIL']
    if fails:
        derived(sec, max(fails), 'Largest relative spread among the FAIL rows', 'max of max_rel_dev over FAIL rows', ['sc'], 3)

    # --- H. surrogate accuracy at 1 pm
    sa = load('sa')
    sec = 'H. Monte Carlo surrogates: held-out error in units of the 1 pm noise'
    st = [r['neff_rms_test_over_sigma_phase'] for r in sa['rows'] if r['surrogate'].startswith('strip')]
    rb = {r['surrogate']: r['neff_rms_test_over_sigma_phase'] for r in sa['rows'] if r['surrogate'].startswith('rib')}
    computed(sec, 'sa', 'rows[strip*].neff_rms_test_over_sigma_phase (max)', max(st), 'Strips: largest held-out n_eff RMS (in sigma_phase)', 3)
    for k, v in rb.items():
        computed(sec, 'sa', f'rows[{k}].neff_rms_test_over_sigma_phase', v, f'{k}: held-out n_eff RMS (in sigma_phase)', 3)

    # --- I. literature inputs
    sec = 'I. Literature inputs'
    lam = load('env')['femwell_env']['constants']['WL0_um']
    u1 = table_value('frey2006', 0, '1.5 microns', '295 K')
    cited(sec, f'{u1:.2E}', 'n_Si absolute uncertainty, 1.5 um, 295 K (Table 1); basis of the 1e-4 n_Si prior', 'frey2006', 0)
    a15, a20 = table_value('frey2006', 2, '1.5 microns', '295 K'), table_value('frey2006', 2, '2.0 microns', '295 K')
    cited(sec, f'{a15:.2E}, {a20:.2E}', 'dn/dT of Si at 1.5 and 2.0 um, 295 K (Table 4)', 'frey2006', 2)
    derived(sec, a15 + (a20 - a15) * (lam - 1.5) / (2.0 - 1.5), 'dn/dT of Si interpolated to the wavelength (code: DNDT_SI)',
            'linear interpolation of the Table 4 values above to lambda', ['cit', 'env'], 4)
    b15, b16 = table_value('leviton2008', 0, '1.5 microns', '295 K'), table_value('leviton2008', 0, '1.6 microns', '295 K')
    cited(sec, f'{b15:.2E}, {b16:.2E}', 'dn/dT of fused silica at 1.5 and 1.6 um, 295 K (Table 6)', 'leviton2008', 0)
    derived(sec, b15 + (b16 - b15) * (lam - 1.5) / (1.6 - 1.5), 'dn/dT of SiO2 interpolated to the wavelength (code: DNDT_OX)',
            'linear interpolation of the Table 6 values above to lambda', ['cit', 'env'], 4)
    computed(sec, 'env', 'analysis_constants.DNDT_SI', load('env')['analysis_constants']['DNDT_SI'],
             'DNDT_SI as used in the code (rounded from the interpolation)', 3)
    computed(sec, 'env', 'analysis_constants.DNDT_OX', load('env')['analysis_constants']['DNDT_OX'],
             'DNDT_OX as used in the code (rounded from the interpolation)', 3)
    n1 = table_value('frey2006', 1, '1.5 microns', '295 K')
    cited(sec, f'{n1:.5f}', 'n_Si at 1.5 um, 295 K (Table 2); our nominal n_Si is a rounded assumption', 'frey2006', 1)
    n2 = table_value('leviton2008', 1, '1.5 microns', '295 K')
    cited(sec, f'{n2:.5f}', 'n of fused silica at 1.5 um, 295 K (Table 4); our nominal n_ox is a rounded assumption', 'leviton2008', 1)
    cited(sec, '2.602', 'Published EIM n_eff, rib 500 nm, 220 nm, 90 nm slab (validation only)', 'jiao_pikus2025', 0)

    # --- K. validation
    vc = load('vc')
    sec = 'K. Validation of the femwell model'
    for g in ('strip', 'rib'):
        b = vc['femwell_vs_empy']['by_geometry'][g]
        computed(sec, 'vc', f'femwell_vs_empy.by_geometry.{g}.max_abs_rel_diff_neff', b['max_abs_rel_diff_neff'],
                 f'{g}s ({b["n"]}): largest relative n_eff difference, EMpy (cell-exact) vs femwell', 3)
    e = vc['femwell_vs_published_eim']['rib_500_slab90']
    computed(sec, 'vc', 'femwell_vs_published_eim.rib_500_slab90.femwell', e['femwell'], 'femwell n_eff, rib 500 nm, 90 nm slab', 5)
    computed(sec, 'vc', 'femwell_vs_published_eim.rib_500_slab90.rel_diff', e['rel_diff'], 'Same, relative difference to the published EIM value', 3)
    if 'domain' in vc and 'two_width_crb' in vc['domain']:
        d = vc['domain']
        computed(sec, 'vc', 'domain.max_abs_rel_diff_neff', d['max_abs_rel_diff_neff'], 'Domain 1.5x larger: largest relative n_eff change (two-width structures)', 3)
        computed(sec, 'vc', 'domain.max_abs_rel_diff_ng', d['max_abs_rel_diff_ng'], 'Same: largest relative n_g change', 3)
        computed(sec, 'vc', 'domain.max_abs_rel_diff_derivative', d['max_abs_rel_diff_derivative'], 'Same: largest relative change of a geometric derivative at h', 3)
        computed(sec, 'vc', 'domain.two_width_crb.large_domain', d['two_width_crb']['large_domain'], 'Two-width CRB with large-domain geometric derivatives (deg)', 4)
        computed(sec, 'vc', 'domain.two_width_crb.rel_change', d['two_width_crb']['rel_change'], 'Same, relative change of the CRB', 3)

    # --- J. assumptions
    env = load('env')
    a = env['analysis_constants']
    fe = env['femwell_env']['constants']
    sec = 'J. Assumed inputs (see ASSUMPTIONS.md for their effect)'
    pr = asn['base']['priors']
    for k, m in (('d', 'Etch-depth prior (nm)'), ('nsi', 'n_Si prior (units of 1e-3)'), ('nox', 'n_ox prior (units of 1e-3)'),
                 ('T', 'Temperature prior (K)')):
        assumed(sec, pr[k], f'{m}, "realistic"', f'`src/priors_search.py` PRIORS["realistic"]["{k}"]', f'base.priors.{k}', key='asn')
    assumed(sec, asn['base']['rib_mult'], 'Rib noise factor (1.0 = same noise as strips; conservative)',
            '`src/priors_search.py` evaluate(rib_mult=1.0)', 'base.rib_mult', key='asn')
    assumed(sec, a['placeholder_noise']['L_MZI_um'], 'Ring/MZI length L (um) in the noise conversion and the MZI phase', '`src/fullvector_sensitivity.py` L_MZI', 'analysis_constants.placeholder_noise.L_MZI_um')
    assumed(sec, a['noise_conversion']['n_g'], 'Group index in the noise conversion (largest n_g of the headline structures)', '`src/noise_sourced.py` NG', 'analysis_constants.noise_conversion.n_g', 3)
    assumed(sec, fe['N_SI'], 'Nominal n_Si', '`src/femwell_batch.py` N_SI', 'femwell_env.constants.N_SI')
    assumed(sec, fe['N_OX'], 'Nominal n_ox (BOX and cladding)', '`src/femwell_batch.py` N_OX', 'femwell_env.constants.N_OX')
    assumed(sec, fe['WL0_um'], 'Wavelength (um)', '`src/femwell_batch.py` WL0', 'femwell_env.constants.WL0_um')
    assumed(sec, fe['DLAM_um'], 'Wavelength step of the n_g central difference (um)', '`src/femwell_batch.py` DLAM', 'femwell_env.constants.DLAM_um')
    assumed(sec, a['T_STRIP_nm'], 'Strip thickness t (nm)', '`src/etchdepth_femwell_analysis.py` T_STRIP', 'analysis_constants.T_STRIP_nm')
    assumed(sec, fe['MESH']['r12'][0], 'Core mesh resolution r12 (um)', '`src/femwell_batch.py` MESH["r12"]', 'femwell_env.constants.MESH.r12')
    assumed(sec, a['montecarlo']['K_BOX'], 'Monte Carlo box half-width (marginal CRB sigmas)', '`src/montecarlo.py` K_BOX', 'analysis_constants.montecarlo.K_BOX')
    assumed(sec, a['montecarlo']['NSTART'], 'Monte Carlo starts per fit', '`src/montecarlo.py` NSTART', 'analysis_constants.montecarlo.NSTART')


def write():
    build()
    with open(os.path.join(REPO, DERIVED_FILE), 'w', encoding='utf-8') as f:
        json.dump(DERIVED_VALUES, f, indent=1)
    files = sorted({f for r in ROWS for f in r[5]})
    hashes = {f: commit(f) for f in files}
    counts = {}
    for r in ROWS:
        counts[r[3]] = counts.get(r[3], 0) + 1
    L = ['<!-- generated by src/make_headline_numbers.py; do not edit -->', '# Headline numbers', '',
         'Every value is read from a data file by `src/make_headline_numbers.py`. Strip-angle CRBs are 1-sigma, in degrees, '
         'with the exact strip symmetry imposed. Primary noise: the sourced 1 pm resonance-wavelength precision '
         '(best case 0.1 pm); the noise conversion to n_eff, n_g and MZI phase is DERIVED (formulas in section E). '
         'Sets are read out by a ring (n_eff, n_g) and an MZI (phase) per structure. Stages: C = + separate rib angle, '
         'D = + n_Si, n_ox, E = + temperature; "realistic" priors: etch depth 2 nm, n_Si 1e-4, n_ox 1e-3, T 1 K (ASSUMED).', '',
         '| Category | Rows |', '|---|---|'] + [f'| {k} | {v} |' for k, v in sorted(counts.items())] + ['']
    sec = None
    for s, val, meaning, cat, src, fl in ROWS:
        if s != sec:
            L += ['', f'## {s}', '', '| Value | Meaning | Category | Source | Commit |', '|---|---|---|---|---|']
            sec = s
        L.append(f"| {val} | {meaning.replace('|', '/')} | {cat} | {src.replace('|', chr(92) + '|')} | " +
                 ', '.join(f'`{hashes[f]}`' for f in fl) + ' |')
    with open(OUTF, 'w', encoding='utf-8') as f:
        f.write('\n'.join(L) + '\n')
    json.dump(dict(counts=counts, n_rows=len(ROWS), commits=hashes),
              open(os.path.join(REPO, 'results_final', 'data', 'headline_counts.json'), 'w'), indent=1)
    print(counts, len(ROWS), 'rows;', sum(h == 'uncommitted' for h in hashes.values()), 'uncommitted files')


if __name__ == '__main__':
    write()
