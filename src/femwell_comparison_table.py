"""
Builds results_femwell/COMPARISON.md: every number the writeup states, as
originally stated, from the original averaging, from cell-exact averaging and
from femwell, each at step h with its range over steps h/2..2h where
available. Numbers only; the verdicts are written by hand in SUMMARY.md.

Inputs:
  results_corrected/data/writeup_numbers_{original,corrected}.json, robustness_checks.json (cell_exact)
  results_femwell/data/writeup_numbers_femwell.json, robustness_checks.json (femwell)
Run from src/:  python femwell_comparison_table.py
"""
import json
import os

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
C_D = os.path.join(REPO, 'results_corrected', 'data')
F_D = os.path.join(REPO, 'results_femwell', 'data')


def load(d, n):
    return json.load(open(os.path.join(d, n)))


O = load(C_D, 'writeup_numbers_original.json')
C = load(C_D, 'writeup_numbers_corrected.json')
CR = load(C_D, 'robustness_checks.json')['cell_exact']
F = load(F_D, 'writeup_numbers_femwell.json')
FR = load(F_D, 'robustness_checks.json')['femwell']


def g(*path):
    def get(N):
        x = N
        for p in path:
            x = x[p]
        return x
    return get


def fmt(x, f):
    if x is None:
        return 'n/a'
    if isinstance(x, (bool, str)):
        return str(x)
    if f == 'e':
        return f'{x:.2e}'
    if f == 'x':
        return f'{x:,.1f}x' if abs(x) < 1e4 else f'{x:,.0f}x'
    if f == 'deg':
        return f'{x:.3g}°'
    if f == 'pct':
        return f'{x:.3g}%'
    return f'{x:.{f}f}'


def span3(get, N, R, f):
    """min..max over h/2, h, 2h"""
    try:
        v = [get(R['D_half']), get(N), get(R['D_2h'])]
    except (KeyError, TypeError):
        return 'n/a'
    if isinstance(v[0], (bool, str)):
        return ' / '.join(map(str, v))
    return f'{fmt(min(v), f)} .. {fmt(max(v), f)}'


ROWS = [
    ('Two TE rings: best condition number', '2.4e5', g('ring_pairs', 'cond_best'), 'e'),
    ('Two TE rings: worst condition number', '7.3e9', g('ring_pairs', 'cond_worst'), 'e'),
    ('Spread best/worst ring pair', '30,000x', g('ring_pairs', 'spread'), 'x'),
    ('Best pair by smallest eigenvalue', '-', g('ring_pairs', 'best_pair_by_lambda_min'), None),
    ('550+600 baseline: smallest eigenvalue', '6.33', g('additions_550_600', 'baseline', 'lambda_min'), 'e'),
    ('550+600 baseline: condition number', '2.4e5', g('additions_550_600', 'baseline', 'cond'), 'e'),
    ('550+600 baseline: trace of inverse', '-', g('additions_550_600', 'baseline', 'trace_inv'), 'e'),
    ('+ TE ring 500 nm: smallest eigenvalue', '7.64 (stated as 400 nm)', g('additions_550_600', '+TE ring 500', 'lambda_min'), 'e'),
    ('+ TE ring 500 nm: condition number', '3.0e5', g('additions_550_600', '+TE ring 500', 'cond'), 'e'),
    ('+ TE ring 400 nm: smallest eigenvalue', '-', g('additions_550_600', '+TE ring 400', 'lambda_min'), 'e'),
    ('+ TM ring 400 nm: smallest eigenvalue', '6.96', g('additions_550_600', '+TM ring 400', 'lambda_min'), 'e'),
    ('+ TM ring 400 nm: condition number', '3.6e5', g('additions_550_600', '+TM ring 400', 'cond'), 'e'),
    ('+ TM ring 400 nm: trace of inverse', '-', g('additions_550_600', '+TM ring 400', 'trace_inv'), 'e'),
    ('+ MZI 600 nm TE: smallest eigenvalue', '10.21', g('additions_550_600', '+MZI 600 TE', 'lambda_min'), 'e'),
    ('+ MZI 600 nm TE: condition number', '6.6e5', g('additions_550_600', '+MZI 600 TE', 'cond'), 'e'),
    ('+ MZI 600 nm TE: trace of inverse', '-', g('additions_550_600', '+MZI 600 TE', 'trace_inv'), 'e'),
    ('Same-width ring + MZI condition number', 'singular', g('same_width_ring_mzi_cond'), 'e'),
    ('Within-type mean angle (15 structures)', '~4°', g('atlas_15', 'within_mean'), 'deg'),
    ('Cross-type mean angle (15 structures)', '~46°', g('atlas_15', 'cross_mean'), 'deg'),
    ('Clusters cleanly by type (15)', 'yes', g('atlas_15', 'clusters_cleanly'), None),
    ('Within-type mean angle (30 structures)', '-', g('atlas_30', 'within_mean'), 'deg'),
    ('Cross-type mean angle (30 structures)', '-', g('atlas_30', 'cross_mean'), 'deg'),
    ('Spearman(angle, smallest eigenvalue), 15', '-0.37', g('spearman_15', 'v1_angle_vs_lambda_min'), 3),
    ('Spearman with group index (v2), 15', '-0.38', g('spearman_15', 'v2_angle_vs_lambda_min'), 3),
    ('Spearman, 30 structures', '-', g('spearman_30', 'v1_angle_vs_lambda_min'), 3),
    ('Slot TE angle component, most negative', '-0.62', lambda N: N['angle_component_range']['slot-TE'][0], 3),
    ('Slot TE angle component, least negative', '-0.60', lambda N: N['angle_component_range']['slot-TE'][1], 3),
    ('Strip TE angle component, most negative', '-', lambda N: N['angle_component_range']['strip-TE'][0], 3),
    ('PC1 variance, 15', '95.3%', g('pca_15', 'explained_pct', 0), 'pct'),
    ('PC2 variance, 15', '4.7%', g('pca_15', 'explained_pct', 1), 'pct'),
    ('PC3 variance, 15', '0.009%', g('pca_15', 'explained_pct', 2), 'pct'),
    ('PC3 variance, 30 (with rib)', '~0.1%', g('pca_30', 'explained_pct', 2), 'pct'),
    ('PC3 loading width (pre-rib, 15)', '0.60', lambda N: abs(N['pca_15']['pc3_loadings'][0]), 3),
    ('PC3 loading thickness (15)', '0.03', lambda N: abs(N['pca_15']['pc3_loadings'][1]), 3),
    ('PC3 loading angle (15)', '0.80', lambda N: abs(N['pca_15']['pc3_loadings'][2]), 3),
    ('Seven scalings: largest PC3 (15)', '<0.01% (really 0.014%)', g('seven_scalings_15', 'max_pc3_pct'), 'pct'),
    ('Seven scalings: slot angle component largest in all', '-', g('seven_scalings_15', 'slot_angle_largest_all'), None),
    ('Slot gain, same noise as strip', '90.8x', g('slot_noise', '1.0'), 'x'),
    ('Slot gain, 4x noise (sourced)', '7.1x', g('slot_noise', '4.0'), 'x'),
    ('Slot gain, 10x noise', '2.1x', g('slot_noise', '10.0'), 'x'),
    ('Slot TM vs nearest strip TE direction', '~2°', g('slot_TM_vs_strip_TE_min_angle', 'angle_deg'), 'deg'),
    ('Rib projection on empty direction, min', '0.087', lambda N: N['pc3_projection_range']['rib-TE'][0], 3),
    ('Rib projection on empty direction, max', '0.154', lambda N: N['pc3_projection_range']['rib-TE'][1], 3),
    ('Others projection (writeup set), max', '0.031', lambda N: N['pc3_projection_others_writeup_set'][1], 4),
    ('Others projection (all widths), max', '-', lambda N: N['pc3_projection_others_all_widths'][1], 4),
    ('Rib gain, conservative noise', '571.7x', g('isolated', 'rib x1.0', 'lmin_improvement'), 'x'),
    ('Rib gain, optimistic noise', '~4000x', g('isolated', 'rib x0.2', 'lmin_improvement'), 'x'),
    ('SWG gain', '1.1x', g('isolated', 'swg x1.0', 'lmin_improvement'), 'x'),
    ('Strip baseline angle CRB (5 widths, TE+TM, ring+MZI)', '0.0210°', g('crb_large_budget', 'strip_baseline'), 'deg'),
    ('+ rib conservative: angle CRB', '0.0009°', g('crb_large_budget', '+rib x1.0', 'std'), 'deg'),
    ('+ rib conservative: improvement', '23x', g('crb_large_budget', '+rib x1.0', 'improvement'), 'x'),
    ('+ rib optimistic: angle CRB', '0.0004°', g('crb_large_budget', '+rib x0.2', 'std'), 'deg'),
    ('+ rib optimistic: improvement', '51x', g('crb_large_budget', '+rib x0.2', 'improvement'), 'x'),
    ('+ slot sourced noise: angle CRB', '0.0081°', g('crb_large_budget', '+slot x4.0', 'std'), 'deg'),
    ('+ slot sourced noise: improvement', '2.6x', g('crb_large_budget', '+slot x4.0', 'improvement'), 'x'),
    ('Rib alone (5 widths, TE ring+MZI) angle CRB, conservative', '-', g('crb_alone', 'rib x1.0'), 'deg'),
    ('Two strip rings (550+600) angle CRB', '0.188°', g('crb_two_rings', 'angle_std'), 'deg'),
    ('Two strip rings: worst direction, width component', '-0.92', lambda N: N['crb_two_rings']['worst_eigvec_wta'][0], 3),
    ('Two strip rings: worst direction, angle component', '-0.37', lambda N: N['crb_two_rings']['worst_eigvec_wta'][2], 3),
    ('Two rings + rib 400 (cons.): improvement', '11x', g('crb_two_rings', '+rib400 x1.0', 'improvement'), 'x'),
    ('Two rings + rib 400 (cons.): angle CRB', '-', g('crb_two_rings', '+rib400 x1.0', 'std'), 'deg'),
    ('Two rings + rib 400 (opt.): improvement', '14x', g('crb_two_rings', '+rib400 x0.2', 'improvement'), 'x'),
    ('Two rings + rib 400 (opt.): angle CRB', '-', g('crb_two_rings', '+rib400 x0.2', 'std'), 'deg'),
    ('Two rings + slot 150 (4x): improvement', '1.9x', g('crb_two_rings', '+slot150 x4', 'improvement'), 'x'),
    ('Two rings + third strip ring 400: improvement', '2.0x', g('crb_two_rings', '+strip400 (control)', 'improvement'), 'x'),
] + [(f'Slab {s} nm: rib gain', st, g('slab_sweep', s, 'lmin_improvement'), 'x')
     for s, st in (('50', '789x'), ('70', '961x'), ('90', '834x'), ('110', '585x'), ('130', '352x'), ('150', '182x'))]


def val(get, N, f):
    try:
        x = get(N)
    except (KeyError, TypeError, IndexError):
        return 'n/a'
    if isinstance(x, list):
        return '+'.join(map(str, x))
    return fmt(x, f) if f else str(x)


if __name__ == '__main__':
    lines = ['# Every writeup number: original vs cell-exact vs femwell', '',
             'Generated by `src/femwell_comparison_table.py`. "h/2..2h" is the min..max over the three '
             'finite-difference steps (all derivatives taken at h/2, h or 2h). Noise assumptions unchanged. '
             'EMpy-anisotropic is not a column: the pilot showed its derivatives are not reliable, and no '
             'full anisotropic batch was run (see results_anisotropic/METHOD_COMPARISON.md).', '',
             '| Quantity | As stated (original writeup) | Original averaging | Cell-exact (h) | Cell-exact h/2..2h | femwell (h) | femwell h/2..2h |',
             '|---|---|---|---|---|---|---|']
    for name, stated, get, f in ROWS:
        lines.append(f'| {name} | {stated} | {val(get, O, f)} | {val(get, C, f)} | {span3(get, C, CR, f) if f else "n/a"} | '
                     f'{val(get, F, f)} | {span3(get, F, FR, f) if f else "n/a"} |')
    with open(os.path.join(REPO, 'results_femwell', 'COMPARISON.md'), 'w', encoding='utf-8') as fh:
        fh.write('\n'.join(lines) + '\n')
    print('\n'.join(lines))
