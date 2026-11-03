"""Fig. 2: angle/width ratio R = (dn_eff/da)/(dn_eff/dw)/(-h tan 1 deg) vs slab thickness for the ribs, with the
strips (R = 1 by symmetry) and the closed-form R from the mode-field centroid.

Data: results_montecarlo/sym_search/data/step_check_all.csv (raw femwell derivatives at step h, the file behind the
step-check rows of HEADLINE_NUMBERS.md) and results_final/data/crosschecks.json (closed-form R)."""
import csv
import os

import matplotlib.pyplot as plt
import numpy as np

import style

T, TAN1 = 220.0, np.tan(np.radians(1.0))


def ratios():
    rows = list(csv.DictReader(open(os.path.join(style.REPO, 'results_montecarlo/sym_search/data/step_check_all.csv'))))
    D = {}
    for r in rows:
        if r['output'] == 'neff' and r['param'] in ('w', 'a'):
            D.setdefault(r['structure'], {})[r['param']] = float(r['D_h'])
    R = {}
    for uid, d in D.items():
        if 'w' in d and 'a' in d:
            slab = int(uid.split('slab')[1]) if 'slab' in uid else 0
            R[uid] = d['a'] / d['w'] / (-(T - slab) * TAN1)
    return R


def main():
    style.setup()
    R = ratios()
    xc = style.load('results_final/data/crosschecks.json')['closed_form_R']['units']
    fig, ax = plt.subplots(figsize=(style.SINGLE, 2.45))
    widths = (300, 400, 500, 600)
    cols = plt.cm.cividis(np.linspace(0.05, 0.8, len(widths)))
    marks = ('o', 's', '^', 'D')
    for W, c, m in zip(widths, cols, marks):
        pts = sorted((int(u.split('slab')[1]), v) for u, v in R.items() if u.startswith(f'rib_TE_{W}_'))
        x, y = zip(*pts)
        ax.plot(x, y, '-', color=c, marker=m, mfc=c, mec=c, ms=3.5, lw=0.9, label=f'rib {W} nm')
        for s in (50, 90):
            k = f'rib_TE_{W}_slab{s}'
            if k in xc:
                ax.plot(s, xc[k]['R_closed'], marker=m, mfc='none', mec='k', mew=0.6, ms=6, ls='none')
    strips = [v for u, v in R.items() if u.startswith('strip')]
    ax.axhline(1.0, color=style.C['strip_TE'], lw=1.2)
    ax.text(25, 1.0 + 0.008, f'strips, TE and TM ({min(strips):.4f}–{max(strips):.4f})', color=style.C['strip_TE'],
            fontsize=7, va='bottom')
    ax.plot([], [], marker='o', mfc='none', mec='k', mew=0.6, ms=6, ls='none', label='closed form, field centroid')
    ax.legend(loc='upper right', handletextpad=0.3, borderaxespad=0.3, bbox_to_anchor=(1.0, 0.93), labelspacing=0.25)
    ax.set_xlabel('Slab thickness $s$ (nm)')
    ax.set_ylabel(r'Angle/width ratio $R$')
    ax.set_xlim(22, 158)
    ax.set_ylim(0.66, 1.06)
    ax.set_xticks([30, 50, 70, 90, 110, 130, 150])
    style.save(fig, 'fig2_ratio.pdf')


if __name__ == '__main__':
    main()
