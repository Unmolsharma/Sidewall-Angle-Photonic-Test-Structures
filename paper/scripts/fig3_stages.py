"""Fig. 3: strip sidewall-angle CRLB across the unknown-parameter stages for strips alone, strips + one rib and
strips + two widths of the same rib (70 nm slab), without priors and with realistic priors, at 1 pm.
Unidentifiable cases (CRLB infinite) are drawn in a separate band above the axis instead of at infinity.

Data: results_final/sym_search_1pm/data/priors_table.csv (rib noise factor 1.0)."""
import csv
import math
import os

import matplotlib.pyplot as plt

import style

STAGES = ['B', 'C', 'D', 'E']
SETS = [('strip alone', 'strips only', style.C['black'], 'x', ':', -0.18),
        ('strip + rib70 (500 nm)', 'strips + one rib', style.C['rib70'], 'D', '--', 0.0),
        ('strip + rib70 (500 + 300 nm)', 'strips + two rib widths', style.C['rib70'], 'o', '-', 0.18)]


def main():
    style.setup()
    rows = list(csv.DictReader(open(os.path.join(style.REPO, 'results_final/sym_search_1pm/data/priors_table.csv'))))
    val = {(r['stage'], r['prior'], r['set']): float(r['a_deg_h']) for r in rows if r['rib_mult'] == '1.0'}
    fig, (band, ax) = plt.subplots(2, 1, figsize=(style.SINGLE, 2.75), sharex=True,
                                   gridspec_kw=dict(height_ratios=[0.5, 4], hspace=0.02))
    for setname, label, c, m, ls, dx in SETS:
        for prior, filled, dd in (('none', False, -0.05), ('realistic', True, 0.05)):
            xs, ys = [], []
            for i, st in enumerate(STAGES):
                v, x = val[(st, prior, setname)], i + dx + dd
                kw = dict(marker=m, color=c, mfc=c if filled else 'none', mew=0.8, ls='none', ms=4.5)
                if math.isinf(v):
                    band.plot(x, 0.5, **kw)
                else:
                    ax.plot(x, v, **kw)
                    xs.append(x)
                    ys.append(v)
            if len(xs) > 1:
                ax.plot(xs, ys, ls=ls, color=c, lw=0.8)
    band.set_ylim(0, 1)
    band.set_yticks([0.5])
    band.set_yticklabels([r'$\infty$'])
    band.tick_params(axis='x', which='both', bottom=False, top=False)
    band.text(0.33, 0.5, 'unidentifiable', fontsize=7, va='center', ha='center')
    ax.set_yscale('log')
    ax.set_ylim(1.5e-3, 1.5e2)
    ax.set_xlim(-0.45, 3.45)
    ax.set_xticks(range(4))
    ax.set_xticklabels(['B: + depth', 'C: + rib angle', 'D: + $n_\\mathrm{Si}$, $n_\\mathrm{ox}$', 'E: + $T$'])
    ax.set_ylabel(r'CRLB $\sigma_\theta$ (°)')
    ax.text(1.15, 45, 'strips + one rib', color=style.C['rib70'], fontsize=7, va='center')
    ax.text(1.15, 0.85, 'strips + two rib widths', color=style.C['rib70'], fontsize=7, va='center')
    ax.plot([], [], 'o', mfc='none', mec='k', mew=0.8, ms=4, label='no priors')
    ax.plot([], [], 'o', mfc='k', mec='k', ms=4, label='realistic priors')
    ax.plot([], [], 'x', color='k', ms=4, label='strips only')
    ax.legend(loc='lower right', handletextpad=0.2, borderaxespad=0.3, labelspacing=0.2)
    style.save(fig, 'fig3_stages.pdf')


if __name__ == '__main__':
    main()
