"""Fig. 4: Monte Carlo fits at the sourced 1 pm noise (2,000 noisy sets per case, nonlinear least squares, stage D,
realistic priors).
(a) Histograms (log scale) of the fitted strip angle for strips + one rib 500 nm (50 nm slab; case c) and strips +
    rib 300 + 500 nm (50 nm slab; case b), with the fit box of case c; the text gives the slope and correlation of
    fitted width vs fitted angle for case c (mc_valley_1pm.json) and t tan(theta) per degree (crosschecks.json).
(b) Histograms of the fitted angle for case b and the best one-etch set (case a) with Gaussians of the CRLB width
    (crb_a_fd).
Data: results_montecarlo/data/{mc_results_1pm.json, mc_fits_1pm.csv, mc_valley_1pm.json}."""
import csv
import os

import matplotlib.pyplot as plt
import numpy as np

import style


def main():
    style.setup()
    R = style.load('results_montecarlo/data/mc_results_1pm.json')
    V = style.load('results_montecarlo/data/mc_valley_1pm.json')
    slope = style.load('results_final/data/crosschecks.json')['width_angle_equivalence']['analytic_t_tan_ha_over_ha_nm_per_deg']
    fits = {}
    with open(os.path.join(style.REPO, 'results_montecarlo/data/mc_fits_1pm.csv'), newline='') as f:
        for r in csv.DictReader(f):
            if r['a_hat'] and r['hat_w']:
                fits.setdefault(r['case'], []).append((float(r['a_hat']), float(r['hat_w']), r['status']))
    fig, (ax, hx) = plt.subplots(1, 2, figsize=(style.DOUBLE, 2.4), gridspec_kw=dict(width_ratios=[1.15, 1]))

    box = R['c']['theta_box']
    ba = np.linspace(-box['a'], box['a'], 61)
    for case, label, col in (('c', 'one rib', '0.45'), ('b', 'two rib widths', style.C['rib50'])):
        th = np.array([p[0] for p in fits[case]])
        ax.hist(th, bins=ba, histtype='step', color=col, lw=1.0)
    ax.set_yscale('log')
    for x0 in (-box['a'], box['a']):
        ax.axvline(x0, color='k', ls=':', lw=0.6)
    ax.text(0.07, 0.97, f"one rib (grey): {R['c']['n_at_bound']} of {R['c']['n_sets']} fits at the box edge;\n"
            f"fitted width vs angle: slope {V['c']['slope_w_per_a_nm_per_deg']:.3f} nm/°, correlation "
            f"{V['c']['corr_w_a']:.6f}\n" + r"(line $t\tan\theta$: " + f"{slope:.3f} nm/°)",
            transform=ax.transAxes, fontsize=6.5, va='top')
    ax.text(0.63, 0.6, 'two rib widths', color=style.C['rib50'], transform=ax.transAxes, fontsize=7)
    ax.set_xlabel(r'Fitted angle $\hat\theta-\theta_0$ (°)')
    ax.set_ylabel('Fits per bin')
    ax.set_xlim(-box['a'] * 1.08, box['a'] * 1.08)
    ax.set_ylim(0.8, 3e3)
    ax.text(-0.13, 1.0, '(a)', transform=ax.transAxes, fontsize=7.5, fontweight='bold', va='top')

    bins = np.linspace(-1.5, 1.5, 61)
    x = np.linspace(-1.5, 1.5, 400)
    for case, label, col, ls in (('b', 'two rib widths, 50 nm slab', style.C['rib50'], '-'),
                                 ('a', 'best set of 8, 90 nm slab', style.C['rib90'], '--')):
        th = np.array([p[0] for p in fits[case]])
        hx.hist(th, bins=bins, density=True, histtype='step', color=col, lw=1.0)
        s = R[case]['crb_a_fd']
        hx.plot(x, np.exp(-0.5 * (x / s) ** 2) / (s * np.sqrt(2 * np.pi)), ls, color='k', lw=0.7)
        peak = 1 / (s * np.sqrt(2 * np.pi))
        hx.text(s * 1.25, peak * 0.8 if case == 'a' else peak * 1.05, f"{label}\nCRLB {s:.3f}°, fits {R[case]['a_std_all']:.3f}°",
                color=col, fontsize=6.5, va='center')
    hx.set_xlabel(r'Fitted angle $\hat\theta-\theta_0$ (°)')
    hx.set_ylabel('Probability density (1/°)')
    hx.set_xlim(-1.5, 1.5)
    hx.text(-0.13, 1.0, '(b)', transform=hx.transAxes, fontsize=7.5, fontweight='bold', va='top')
    style.save(fig, 'fig4_montecarlo.pdf')


if __name__ == '__main__':
    main()
