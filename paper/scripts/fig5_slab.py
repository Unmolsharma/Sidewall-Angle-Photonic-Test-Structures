"""Fig. 5: CRLB of the two-width set (strip TE 400-600 nm + rib 300 and 500 nm) vs slab thickness, stage D,
realistic priors, 1 pm; bars span the derivative steps h/2, h, 2h.

Data: results_final/sym_search_1pm/data/slab_pairs.json ('D|realistic')."""
import matplotlib.pyplot as plt
import numpy as np

import style


def main():
    style.setup()
    P = style.load('results_final/sym_search_1pm/data/slab_pairs.json')['D|realistic']
    pts = sorted((v['slab'], v['a_deg_h'], v['a_deg_min'], v['a_deg_max']) for k, v in P.items()
                 if 'rib_TE_300' in k and 'rib_TE_500' in k)
    s, a, lo, hi = map(np.array, zip(*pts))
    fig, ax = plt.subplots(figsize=(style.SINGLE, 2.3))
    ax.errorbar(s, a, yerr=[a - lo, hi - a], fmt='-', color=style.C['rib50'], marker='^', ms=4, lw=0.9,
                elinewidth=0.7, capsize=2)
    i = int(np.argmin(a))
    ax.annotate(f'{a[i]:.2f}° at {s[i]} nm', (s[i], a[i]), xytext=(s[i] + 18, a[i] * 0.72), fontsize=7,
                arrowprops=dict(arrowstyle='-', lw=0.5, color='k'))
    ax.set_yscale('log')
    ax.set_xlabel('Slab thickness $s$ (nm)')
    ax.set_ylabel(r'CRLB $\sigma_\theta$ (°)')
    ax.set_xticks(s)
    ax.set_xlim(22, 158)
    ax.set_yticks([0.3, 0.5, 1, 2, 4])
    ax.set_yticklabels(['0.3', '0.5', '1', '2', '4'])
    ax.set_ylim(0.25, 5)
    style.save(fig, 'fig5_slab.pdf')


if __name__ == '__main__':
    main()
