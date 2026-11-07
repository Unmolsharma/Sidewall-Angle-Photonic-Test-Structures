"""Fig. 6: strip-angle CRLB when each structure may depart from the shared parameters (random-effects model), versus
the per-structure spread of (a) the sidewall angle, (b) the width and (c) the thickness, for the two-width set and
the best set, stage D, realistic priors, 1 pm. Dotted lines: the shared-parameter CRLB.

Data: results_review/data/review_analyses.json (variation.<set>.random_effects_crb[_w, _t], crb)."""
import matplotlib.pyplot as plt
import numpy as np

import style


def main():
    style.setup()
    V = style.load('results_review/data/review_analyses.json')['variation']
    fig, axes = plt.subplots(1, 3, figsize=(style.DOUBLE, 2.2), sharey=True)
    spec = [('random_effects_crb', 'Per-structure angle spread (°)'), ('random_effects_crb_w', 'Per-structure width spread (nm)'),
            ('random_effects_crb_t', 'Per-structure thickness spread (nm)')]
    sets = [('two', 'two rib widths (50 nm slab)', style.C['rib50'], '^', '-'), ('best', 'best set of 8 (90 nm slab)', style.C['rib90'], 'v', '--')]
    for ax, (key, xl), letter in zip(axes, spec, 'abc'):
        for name, label, c, m, ls in sets:
            d = V[name][key]
            x = np.array([float(k) for k in d])
            y = np.array(list(d.values()))
            ax.plot(x, y, ls, color=c, marker=m, ms=3.5, lw=0.9, label=label)
            ax.axhline(V[name]['crb'], color=c, ls=':', lw=0.8)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel(xl)
        ax.text(0.04, 0.96, f'({letter})', transform=ax.transAxes, fontsize=7.5, fontweight='bold', va='top')
    axes[0].set_ylabel(r'CRLB $\sigma_\theta$ (°)')
    axes[1].legend(loc='lower right', handletextpad=0.3, borderaxespad=0.3)
    axes[0].set_ylim(0.1, 20)
    style.save(fig, 'fig6_variation.pdf')


if __name__ == '__main__':
    main()
