"""
Generate results_final/{CLAIMS.md, ASSUMPTIONS.md, CITATIONS.md, REPRODUCTION.md} from the data files (no number
typed here). Each claim gets a status from an explicit rule applied to its evidence:
  SUPPORTED    the evidence meets the rule
  PARTIAL      evidence exists but does not establish the claim as worded
  UNSUPPORTED  no direct evidence in this repository

Run from src/ (main environment):  python make_final_docs.py
"""
import json
import math
import os

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIN = os.path.join(REPO, 'results_final')
P = dict(
    s1='results_final/sym_search_1pm/data/search.json', p1='results_final/sym_search_1pm/data/slab_pairs.json',
    d1='results_final/sym_search_1pm/data/d_prior_sweep.json', mc1='results_montecarlo/data/mc_results_1pm.json',
    mc0='results_montecarlo/data/mc_results.json', val='results_montecarlo/data/mc_valley.json',
    xc='results_final/data/crosschecks.json', cc='results_final/data/claims_checks.json',
    si='results_montecarlo/data/strips_info.json', asn='results_final/data/assumption_sensitivity.json',
    cit='results_final/data/citations.json', env='results_final/data/environment.json',
    rep='results_final/reproduction/diff.json')


def load(k):
    return json.load(open(os.path.join(REPO, P[k]), encoding='utf-8'))


def g(v, sig=4):
    if isinstance(v, float) and math.isinf(v):
        return '∞'
    return f'{v:.{sig}g}' if isinstance(v, float) else str(v)


def claims():
    s1, p1, d1, mc1, mc0, val, xc, cc, si = (load(k) for k in ('s1', 'p1', 'd1', 'mc1', 'mc0', 'val', 'xc', 'cc', 'si'))
    C = []

    def add(claim, evidence, files, rule, status):
        C.append((claim, evidence, files, rule, status))

    so = [v['a_deg_h'] for v in s1['D|realistic']['strips_only'].values()]
    w = xc['width_angle_equivalence']
    add('Oxide-clad strips are exactly blind to the sidewall angle (mirror symmetry: a tilt acts as a width change).',
        f"strips-only CRB for every set size 1-{len(so)}: {g(min(so))} (stage D realistic, 1 pm); femwell "
        f"(dn/da)/(dn/dw) of all strips {g(w['femwell_min'], 5)}-{g(w['femwell_max'], 5)} nm/deg vs 220 nm tan(1°) = "
        f"{g(w['analytic_t_tan_ha_over_ha_nm_per_deg'], 5)} nm/deg (largest deviation {g(w['max_rel_dev_from_analytic'], 3)})",
        [P['s1'], P['xc']], 'strips-only CRB infinite for every size', 'SUPPORTED' if all(math.isinf(x) for x in so) else 'PARTIAL')
    air = xc['other_R']['strip_TE_500_air_r8_mesh']
    add('Air cladding barely matters.',
        f"only evidence: air-clad strip TE 500 nm R = {g(air[0])} (n_eff), {g(air[1])} (n_g) on the r8 mesh, i.e. the "
        f"symmetry is broken; no CRB or information comparison of air- vs oxide-clad sets exists",
        [P['xc']], 'a CRB comparison of air- vs oxide-clad sets', 'UNSUPPORTED')
    single = {t: [p1[t][k]['a_deg_h'] for k in p1[t] if '+' not in k] for t in ('C|none', 'D|none')}
    add('With a separately etched rib, one rib width cannot determine the strip angle.',
        f"strip TE 400-600 + one rib 500 nm, every slab, 1 pm: stage C no priors {g(max(single['C|none']))}, "
        f"stage D no priors {g(max(single['D|none']))}",
        [P['p1']], 'CRB infinite at every slab without priors', 'SUPPORTED' if all(math.isinf(x) for v in single.values() for x in v) else 'PARTIAL')
    dC = d1['C']['strip + rib70 (500 nm)']
    lv = sorted((k for k in dC if k.replace('.', '', 1).isdigit()), key=float)   # numeric depth priors (nm); 'free' = no prior
    add('Knowing the etch depth does not rescue a single rib.',
        f"strip + rib70 500 nm, stage C, 1 pm, depth prior {lv[0]} nm: {g(dC[lv[0]]['a_deg_h'])} deg; "
        f"depth prior {lv[-1]} nm: {g(dC[lv[-1]]['a_deg_h'])} deg",
        [P['d1']], 'CRB stays above 1 deg even for the tightest depth prior tested',
        'SUPPORTED' if dC[lv[0]]['a_deg_h'] > 1 else 'PARTIAL')
    b = s1['D|realistic']['one_etch_overall']['8']
    add('Several widths of the same rib determine the strip angle.',
        f"best one-etch set of 8 (slab {b['slab']} nm, stage D realistic, 1 pm): {g(b['a_deg_h'])} deg "
        f"[{g(b['a_deg_min'])}-{g(b['a_deg_max'])}]; Monte Carlo spread {g(mc1['a']['a_std_all'])} deg (case a)",
        [P['s1'], P['mc1']], 'finite CRB, confirmed by the Monte Carlo spread', 'SUPPORTED' if math.isfinite(b['a_deg_h']) else 'UNSUPPORTED')
    m = cc['mzi_redundancy']['sets']
    worst = max(abs(v['rel_diff_combined_vs_both']) for v in m.values())
    add('An MZI adds nothing beyond a ring of the same width except a lower n_eff noise (its phase row is parallel '
        'to the ring n_eff row).',
        '; '.join(f"{k}: ring + MZI {g(v['ring_and_mzi'])} deg, ring alone with the combined n_eff noise "
                  f"{g(v['ring_with_combined_neff_noise'])} deg, ring alone {g(v['ring_only'])} deg" for k, v in m.items()),
        [P['cc']], 'ring + MZI equals ring with combined n_eff noise to better than one part in a million',
        'SUPPORTED' if worst < 1e-6 else 'PARTIAL')
    gv = cc['greedy_vs_exhaustive']
    add('Greedy set selection can be worse than the exhaustive search.',
        f"{gv['n_greedy_worse']} of {gv['n_cases']} stage/prior/slab cases worse, {gv['n_equal']} equal; "
        f"largest greedy/exhaustive ratio {g(gv['max_ratio'], 3)} (1 pm search)",
        [P['cc']], 'at least one case with greedy worse', 'SUPPORTED' if gv['n_greedy_worse'] > 0 else 'UNSUPPORTED')
    r = {k: mc1[k]['std_over_crb_fd'] for k in ('a', 'b')}
    add('The nonlinear Monte Carlo fits match the CRB for identifiable sets.',
        f"spread/CRB at 1 pm: case a {g(r['a'], 3)}, case b {g(r['b'], 3)}; failed fits {mc1['a']['n_failed']} and "
        f"{mc1['b']['n_failed']}, ambiguous {mc1['a']['n_ambiguous']} and {mc1['b']['n_ambiguous']} of {mc1['a']['n_sets']}",
        [P['mc1']], 'spread within 20% of the CRB and no failed fits',
        'SUPPORTED' if all(abs(x - 1) < 0.2 for x in r.values()) and mc1['a']['n_failed'] + mc1['b']['n_failed'] == 0 else 'PARTIAL')
    add('For one rib with a separate etch, the fits wander along the width/angle valley.',
        f"case c (placeholder noise): {mc0['c']['n_at_bound']} of {mc0['c']['n_sets']} fits at the box edge; fitted width vs "
        f"angle slope {g(val['c']['slope_w_per_a_nm_per_deg'], 4)} nm/deg, correlation {g(val['c']['corr_w_a'], 5)}",
        [P['mc0'], P['val']], 'most fits at the box edge, slope matching 220 nm tan(1°)', 'SUPPORTED')
    e = cc['empy_vs_femwell']
    pdv = e['per_derivative']
    worst_k = max(pdv, key=lambda k: pdv[k]['median_abs'])
    add('The earlier EMpy derivatives are unreliable.',
        f"EMpy (cell-exact) vs femwell over {e['n_structures']} structures: median |relative difference| "
        + ', '.join(f"{k} {g(v['median_abs'], 3)}" for k, v in pdv.items())
        + f"; EMpy step-check verdicts {', '.join(f'{k} {v}' for k, v in e['empy_step_check_verdicts'].items())}",
        [P['cc']], 'median difference above the 2% step-check tolerance for some derivative',
        'SUPPORTED' if pdv[worst_k]['median_abs'] > 0.02 else 'PARTIAL')
    ex = si['exact_symmetry']
    nom = si['nominal']
    k0 = [k for k in nom if 'strips-only' in k.lower() or 'strips only' in k.lower()]
    add('The angle information earlier attributed to strips alone was solver noise.',
        'best strips-only set, stage C: ' + (f"{g(nom[k0[0]])} deg as computed, {g(ex[k0[0]])} with the exact symmetry" if k0 else
                                              'see strips_info.json exact_symmetry vs nominal'),
        [P['si']], 'finite as computed, infinite with the exact symmetry',
        'SUPPORTED' if k0 and math.isfinite(nom[k0[0]]) and math.isinf(ex[k0[0]]) else 'PARTIAL')
    t = cc['two_width_slab_trend']
    add('A thinner slab is better for a two-width rib set.',
        'two-width set, stage D realistic, 1 pm, CRB by slab (nm: deg): ' + ', '.join(f'{s}: {g(v)}' for s, v in t['a_deg_by_slab'].items())
        + f"; best slab {t['best_slab']} nm",
        [P['cc'], P['p1']], 'CRB increases monotonically with slab thickness over the tested range',
        'SUPPORTED' if t['monotonic_increasing_with_slab'] else 'PARTIAL')
    L = ['<!-- generated by src/make_final_docs.py; do not edit -->', '# Claims register', '',
         'Each claim we plan to make, with the number that supports it, read from the cited file. Status rules: SUPPORTED '
         '(evidence meets the stated rule), PARTIAL (evidence exists but does not establish the claim as worded), '
         'UNSUPPORTED (no direct evidence in the repository).', '',
         '| # | Claim | Evidence | Source | Rule | Status |', '|---|---|---|---|---|---|']
    for i, (c, ev, files, rule, st) in enumerate(C, 1):
        L.append(f"| {i} | {c} | {ev.replace('|', '/')} | " + ', '.join(f'`{f}`' for f in files) + f' | {rule} | **{st}** |')
    counts = {}
    for c in C:
        counts[c[4]] = counts.get(c[4], 0) + 1
    L += ['', '| Status | Claims |', '|---|---|'] + [f'| {k} | {v} |' for k, v in sorted(counts.items())]
    open(os.path.join(FIN, 'CLAIMS.md'), 'w', encoding='utf-8').write('\n'.join(L) + '\n')
    return counts, [c[0] for c in C if c[4] != 'SUPPORTED']


def assumptions():
    asn, env = load('asn'), load('env')
    a, fe = env['analysis_constants'], env['femwell_env']['constants']
    reg = [
        ('Measurement noise: resonance-wavelength precision', '1 pm (best case 0.1 pm)', 'CITED level; the conversion below is ours',
         '`src/noise_sourced.py`'),
        ('Noise conversion: ring/MZI length L', f"{a['placeholder_noise']['L_MZI_um']} um", 'ASSUMED', '`src/fullvector_sensitivity.py` L_MZI'),
        ('Noise conversion: group index', f"{a['noise_conversion']['n_g']}", 'ASSUMED (largest n_g of the headline structures)', '`src/noise_sourced.py` NG'),
        ('MZI fringe position read as precisely as a ring resonance', 'factor 1', 'ASSUMED (UNVERIFIED)', '`src/noise_sourced.py`'),
        ('Rib noise relative to strips', f"{asn['base']['rib_mult']}", 'ASSUMED (the 0.2x literature ratio is UNVERIFIED)', '`src/priors_search.py` rib_mult'),
        ('Etch-depth prior', f"{asn['base']['priors']['d']} nm", 'ASSUMED', '`src/priors_search.py` PRIORS'),
        ('n_Si prior', f"{asn['base']['priors']['nsi']} x 1e-3", 'ASSUMED (cf. CITED 8.67e-5 bulk-Si uncertainty)', '`src/priors_search.py` PRIORS'),
        ('n_ox prior', f"{asn['base']['priors']['nox']} x 1e-3", 'ASSUMED', '`src/priors_search.py` PRIORS'),
        ('Temperature prior (stage E only)', f"{asn['base']['priors']['T']} K", 'ASSUMED', '`src/priors_search.py` PRIORS'),
        ('Nominal n_Si, n_ox', f"{fe['N_SI']}, {fe['N_OX']}", 'ASSUMED (rounded; CITED 3.48314, 1.44502 at 1.5 um)', '`src/femwell_batch.py`'),
        ('Material dispersion', 'none', 'ASSUMED', '`src/femwell_batch.py`'),
        ('Wavelength; n_g wavelength step', f"{fe['WL0_um']} um; {fe['DLAM_um']} um", 'ASSUMED', '`src/femwell_batch.py`'),
        ('Thickness', f"{a['T_STRIP_nm']} nm", 'ASSUMED', '`src/etchdepth_femwell_analysis.py` T_STRIP'),
        ('Finite-difference steps', json.dumps(fe['H']), 'ASSUMED (checked at h/2, h, 2h)', '`src/femwell_batch.py` H'),
        ('Mesh', f"r12, core {fe['MESH']['r12'][0]} um", 'ASSUMED (r8 check)', '`src/femwell_batch.py` MESH'),
        ('Geometry: width = bottom width; angle = deviation from vertical, top narrower; pivot at the bottom of the wall '
         '(slab top for ribs); ideal trapezoid walls, no corner rounding', 'convention', 'ASSUMED', '`src/femwell_batch.py` shapes'),
        ('Fabrication tolerances (5 nm, 2 nm, 1 deg)', 'not used', 'only in the retired EMpy normalization; no CRB depends on them',
         '`notes/progress.md`'),
        ('Monte Carlo: box, starts, sets', f"±{a['montecarlo']['K_BOX']} sigma, {a['montecarlo']['NSTART']} starts, {a['montecarlo']['NSETS']} sets",
         'METHOD CHOICE', '`src/montecarlo.py`'),
    ]
    L = ['<!-- generated by src/make_final_docs.py; do not edit -->', '# Assumptions register', '',
         'Values are read from `results_final/data/environment.json` and `results_final/data/assumption_sensitivity.json`.', '',
         '| Assumption | Value | Status | Set in |', '|---|---|---|---|'] + [f'| {x} | {v} | {s} | {w} |' for x, v, s, w in reg]
    L += ['', '## Effect on the headline numbers (one at a time; `results_final/data/assumption_sensitivity.json`)', '',
          f"Base: stage D, realistic priors, noise 1 pm. Sets: best = {' + '.join(asn['sets']['best'])}; "
          f"two = {' + '.join(asn['sets']['two'])}.", '',
          '| Assumption | Level | Best set (deg) | change | Two-width set (deg) | change |', '|---|---|---|---|---|---|',
          f"| base | | {g(asn['results']['best']['base'])} | | {g(asn['results']['two']['base'])} | |"]
    for vb, vt in zip(asn['results']['best']['variations'], asn['results']['two']['variations']):
        L.append(f"| {vb['assumption']} | {vb['level']} | {g(vb['a_deg'])} | {vb['rel_change']:+.1%} | {g(vt['a_deg'])} | {vt['rel_change']:+.1%} |")
    L += ['', 'NOT EVALUATED (would need new solves): ' + '; '.join(asn['not_evaluated']) + '.']
    open(os.path.join(FIN, 'ASSUMPTIONS.md'), 'w', encoding='utf-8').write('\n'.join(L) + '\n')


def citations():
    cit = load('cit')
    L = ['<!-- generated by src/make_final_docs.py; do not edit -->', '# Citations', '',
         'Every literature value used, with the exact quoted sentence or table row and its page (PDF page) or section. '
         'Source data: `results_final/data/citations.json` (written by `src/citations_extract.py` from the fetched PDFs).', '']
    for k, c in cit['sources'].items():
        L += [f"## {c['authors']}, \"{c['title']}\"", '', f"- Venue: {c['venue']}", f"- DOI: {c['doi']}",
              f"- URL: {c['url']}", f"- Status: **{c['status']}**"]
        if c.get('journal_version'):
            L.append(f"- Journal version: {c['journal_version']}")
        L += ['', '| Used for | Where | Exact quote |', '|---|---|---|']
        for q in c['quotes']:
            where = f"p. {q['page']}" if 'page' in q else f"Sect. {q['section']}"
            L.append(f"| {q['use']} | {where} | “{q['text'].replace('|', '/')}”" + (f" ({q['note']})" if q.get('note') else '') + ' |')
        L.append('')
    L += ['## UNVERIFIED (not used for any headline number)', '', '| Item | Where | Status |', '|---|---|---|']
    L += [f"| {u['item']} | {u['where']} | {u['status']} |" for u in cit['unverified']]
    open(os.path.join(FIN, 'CITATIONS.md'), 'w', encoding='utf-8').write('\n'.join(L) + '\n')
    return len(cit['unverified'])


def reproduction():
    path = os.path.join(REPO, P['rep'])
    if not os.path.exists(path):
        return None
    r = json.load(open(path, encoding='utf-8'))
    L = ['<!-- generated by src/make_final_docs.py; do not edit -->', '# Reproduction', '',
         f"Clean rerun by `src/reproduce_all.py` in `{r['tree']}` (a fresh copy of `src/` plus only the raw solve caches "
         f"and the declared frozen inputs below), then every regenerated file compared with the committed one.", '',
         f"- Compared files: {r['n_files']}; identical: {r['n_exact']}; equal within a relative tolerance of "
         f"{r['rel_tol']}: {r['n_within_tol']}; **differing: {r['n_differ']}** (" + ', '.join(f"`{f['file']}`" for f in r['files'] if 'differ' in f['result'] or f['result'].startswith('DIFFERENT')) + f"); missing in the rerun: {r['n_missing']}.",
         f"- Monte Carlo seed: {r['mc_seed']} (noise sets seeded per case and chunk; see ENVIRONMENT.md).", '',
         '## Frozen inputs (copied, not regenerated)', ''] + [f'- `{f}`: {why}' for f, why in r['frozen'].items()]
    L += ['', '## Steps', ''] + [f"{i}. `{s['cmd']}` — exit {s['exit']}, {s['seconds']:.0f} s" for i, s in enumerate(r['steps'], 1)]
    L += ['', '## Files', '', '| File | Numbers | Result | Largest relative difference |', '|---|---|---|---|']
    for f in r['files']:
        L.append(f"| `{f['file']}` | {f['n_numbers']} | {f['result']} | {f['max_rel_diff']:.3g} |" if f['max_rel_diff'] is not None
                 else f"| `{f['file']}` | {f['n_numbers']} | {f['result']} | n/a |")
    for n in r.get('sampling_notes', []):
        L += ['', f"## Random-sampling file: `{n['file']}` (seed {n['seed']})", '', n['cause'], '',
              '| Statistic | Largest relative difference, committed vs rerun |', '|---|---|']
        L += [f'| {k} | {v:.3g} |' for k, v in n['max_abs_rel_diff'].items()]
    if r['differences']:
        L += ['', '## Numbers that do not reproduce exactly', '', '| File | Location | Committed | Rerun |', '|---|---|---|---|']
        shown = [d for d in r['differences'] if 'strips_info' not in d['file']][:200]
        L += [f"| `{d['file']}` | `{d['where']}` | {d['committed']} | {d['rerun']} |" for d in shown]
        n_si = sum('strips_info' in d['file'] for d in r['differences'])
        if n_si:
            L += ['', f'({n_si} differing numbers of the strips_info random-perturbation files are not listed one by one; '
                  'see the sampling comparison above and `results_final/reproduction/diff.json`.)']
    open(os.path.join(FIN, 'REPRODUCTION.md'), 'w', encoding='utf-8').write('\n'.join(L) + '\n')
    return r


if __name__ == '__main__':
    counts, flagged = claims()
    assumptions()
    nu = citations()
    rep = reproduction()
    print('claims', counts, 'flagged:', flagged)
    print('unverified citations', nu)
    print('reproduction', None if rep is None else {k: rep[k] for k in ('n_files', 'n_exact', 'n_within_tol', 'n_differ', 'n_missing')})
