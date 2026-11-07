"""Fig. 7: strip-angle CRLB of the two headline sets under variants of the measurement model (stage D, realistic
priors, 1 pm unless stated): correlation rho between the n_eff and n_g errors of one ring (the derived value
-1/sqrt(2) and 0.9), an unknown common wavelength-calibration offset and one known to 10 pm, temperature as a shared
unknown, an independent temperature drift per structure (10 mK), ring-only readout, and MZIs read 3x worse than rings.

Data: results_review/data/review_analyses.json (noise, mzi)."""
import matplotlib.pyplot as plt
import numpy as np

import style


def main():
    style.setup()
    R = style.load('results_review/data/review_analyses.json')
    rows = [('baseline (independent errors)', lambda n: R['noise'][n]['base']),
            (r'$\rho(n_\mathrm{eff}, n_\mathrm{g})=-0.71$ (one FSR)', lambda n: R['noise'][n]['rho']['-0.707']),
            (r'$\rho=0.9$', lambda n: R['noise'][n]['rho']['0.900']),
            ('calibration offset, known to 10 pm', lambda n: R['noise'][n]['calibration']['10 pm']),
            ('calibration offset, unknown', lambda n: R['noise'][n]['calibration']['unknown']),
            ('temperature, unknown (shared)', lambda n: R['noise'][n]['temperature_shared']['unknown']),
            ('10 mK drift per structure', lambda n: R['noise'][n]['temperature_drift_independent']['10 mK']),
            ('MZI read 3x worse than ring', lambda n: R['mzi'][n]['mzi_x3']),
            ('rings only', lambda n: R['mzi'][n]['ring_only'])]
    fig, ax = plt.subplots(figsize=(style.SINGLE, 2.7))
    y = np.arange(len(rows))[::-1]
    for name, c, m, dy in (('two', style.C['rib50'], '^', 0.13), ('best', style.C['rib90'], 'v', -0.13)):
        v = [f(name) for _, f in rows]
        ax.plot(v, y + dy, m, color=c, ms=4.5, ls='none', label={'two': 'two rib widths', 'best': 'best set of 8'}[name])
        ax.axvline(v[0], color=c, ls=':', lw=0.7)
    ax.set_yticks(y)
    ax.set_yticklabels([r for r, _ in rows], fontsize=7)
    ax.set_xscale('log')
    ax.set_xlim(0.07, 1.6)
    ax.set_xticks([0.1, 0.2, 0.5, 1])
    ax.set_xticklabels(['0.1', '0.2', '0.5', '1'])
    ax.set_xlabel(r'CRLB $\sigma_\theta$ (°)')
    ax.legend(loc='lower right', handletextpad=0.2, borderaxespad=0.3)
    style.save(fig, 'fig7_noise.pdf')


if __name__ == '__main__':
    main()
