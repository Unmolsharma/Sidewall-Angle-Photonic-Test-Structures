"""Fig. 8: strip-angle CRLB versus the number of copies N of each structure, for four levels of independent local
variation, with the shared-parameter bound of one copy (and 1.5x, 2x) and the shared-parameter bound at the same N.
Stage D, realistic priors, 1 pm, ring and MZI.

Data: results_review/data/feasibility.json (replication.<set>.ring_mzi)."""
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb

import style


def shade(c, f):
    """Mix the colour c with white; f = 1 is the full colour."""
    return tuple(1 - f * (1 - x) for x in to_rgb(c))


def main():
    style.setup()
    R = style.load('results_review/data/feasibility.json')
    fig, axes = plt.subplots(1, 2, figsize=(style.DOUBLE, 2.5), sharey=True)
    N = [1, 4, 16, 64, 256]
    sets = [('two', 'Two-width set (50 nm slab)', style.C['rib50']), ('best', 'Best set of 8 (90 nm slab)', style.C['rib90'])]
    marks = dict(low='o', moderate='s', high='^', worst='D')
    shades = dict(low=0.4, moderate=0.6, high=0.8, worst=1.0)
    for ax, (name, title, col), letter in zip(axes, sets, 'ab'):
        r = R['replication'][name]['ring_mzi']
        for s, d in r['scenarios'].items():
            lab = rf"{d['width_nm']:g} nm, {d['angle_deg']:g}$\degree$, {d['thickness_nm']:g} nm"
            ax.plot(N, [d['bound'][str(n)] for n in N], '-', color=shade(col, shades[s]), marker=marks[s], ms=3.5, lw=1.0, label=lab)
        ax.plot(N, [r['shared_same_N'][str(n)] for n in N], '-', color=style.C['black'], lw=0.9, marker='x', ms=3.5, label='shared parameters, same $N$')
        for f, ls, t in ((1.0, '-', r'shared, $N=1$'), (1.5, '--', r'$1.5\times$'), (2.0, ':', r'$2\times$')):
            ax.axhline(f * r['shared_N1'], color=style.C['grey'], ls=ls, lw=0.7, zorder=0)
        ax.set_xscale('log', base=2)
        ax.set_yscale('log')
        ax.set_xticks(N)
        ax.set_xticklabels([str(n) for n in N])
        ax.set_xlabel('Copies of each structure, $N$')
        ax.text(0.03, 0.05, f'({letter}) {title}', transform=ax.transAxes, fontsize=7.5, fontweight='bold', va='bottom')
        ax.legend(loc='upper right', handletextpad=0.4, borderaxespad=0.3, labelspacing=0.2, title=r'$\sigma_w,\ \sigma_\theta,\ \sigma_t$', title_fontsize=7)
    axes[0].set_ylabel(r'CRLB $\sigma_\theta$ (°)')
    axes[0].set_ylim(0.005, 300)
    style.save(fig, 'fig8_replication.pdf')


if __name__ == '__main__':
    main()
