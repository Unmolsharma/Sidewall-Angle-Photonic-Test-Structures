"""
Etch-depth nuisance study, part 2 of 2: Fisher / Cramer-Rao analysis of
sidewall-angle (SWA) precision when the rib's partial-etch depth d is an
UNKNOWN nuisance parameter, optionally with a Gaussian prior on it.

Reads (never writes) data/sensitivity_db.json for the strip records, and
reads results_etch_depth/data/rib_steps_<MESH>.json produced by
etch_depth_sims.py for the rib records. All outputs go to results_etch_depth/.

PARAMETERIZATIONS (columns of the Jacobian; units nm, nm, deg, ..., nm)
  orig3 : [w, t, a]                     existing model. Rib t-derivative taken
                                        with SLAB held fixed (as in
                                        rib_waveguide.py), i.e. slab known.
  A_dknown : [w, t, a]                  Case A with d known exactly (rib
                                        t-derivative at fixed d).
  A     : [w, t, a, d]                  one SWA shared by full and partial etch.
  B     : [w, t, a_full, a_part, d]     strip sees a_full only, rib a_part only.
  B_w   : B + separate partial-etch width bias w_part (supplementary check:
          the partial etch is a separate lithography step too).
Strip rows have zero sensitivity to d (and to a_part, w_part).

How the etch depth of ribs with DIFFERENT nominal slab thicknesses relates
(only matters for sets containing >1 slab value):
  common      : one shared additive error dd (the requested model)
  scaled      : one shared fractional etch-rate error, dd_k = (d_k/d_ref)*dd,
                d_ref = 130nm (nominal 90nm slab)
  independent : each distinct slab value is its own etch step with its own
                unknown d_k (each gets the same prior)

PRECISION METRIC
  sigma_SWA = sqrt( [F_total^-1]_SWA,SWA ), F_total = J^T J + diag(prior),
  prior = 1/sigma_prior^2 on every d column only. This is the marginal CRB
  over ALL other unknowns -- never 1/F_SWA,SWA and never a dropped-row
  sub-matrix (except "d known", where dropping d IS the exact limit
  sigma_prior -> 0). If F_total is singular, the SWA marginal is finite only
  if the null space has no SWA component (an unidentifiable parameter that
  is exactly decoupled, e.g. d for a strip-only set); otherwise it is inf.

Noise models are the existing ones, unchanged: SIGMA_NEFF, SIGMA_NG,
SIGMA_PHASE from fullvector_sensitivity.py; rib multiplier 1.0
(conservative) or 0.2 (optimistic) from noise_sourcing.py, applied to all
three sigmas of rib probes exactly as angle_crb.py does.

Run:  python etch_depth_fisher.py [--mesh M0]
"""
import os
import sys
import json
import itertools

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import fullvector_sensitivity as fv
from noise_sourcing import RIB_MULTIPLIER_CONSERVATIVE, RIB_MULTIPLIER_OPTIMISTIC

np.random.seed(0)   # nothing below is stochastic; fixed anyway for reproducibility

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, 'results_etch_depth')
MESH = sys.argv[sys.argv.index('--mesh') + 1] if '--mesh' in sys.argv else 'M0'

D_REF = 130.0                       # nm, nominal etch depth for the 90nm slab
PRIORS = ['known', 1.0, 2.0, 5.0, 10.0, None]   # sigma_prior on d [nm]; None = no prior
PRIOR_LABEL = {'known': 'd known', 1.0: '1 nm', 2.0: '2 nm', 5.0: '5 nm', 10.0: '10 nm', None: 'no prior'}
NOISE = {'conservative': RIB_MULTIPLIER_CONSERVATIVE, 'optimistic': RIB_MULTIPLIER_OPTIMISTIC}

# reference palette (dataviz skill, light mode), fixed slot order
C = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#6250d6', '#e34948']
INK, INK2, GRID = '#0b0b0b', '#52514e', '#e4e3df'


# ---------------------------------------------------------------------------
# probes: one ring (2 observables) or one MZI (1 observable) on one waveguide
# ---------------------------------------------------------------------------
def load_strip_probes():
    with open(os.path.join(REPO, 'data', 'sensitivity_db.json')) as f:
        recs = [r for r in json.load(f) if r['geometry'] == 'strip']
    probes = []
    for r in recs:
        sens = dict(neff={'w': r['dneff_dw'], 't': r['dneff_dt'], 'a': r['dneff_da']},
                    ng={'w': r['dng_dw'], 't': r['dng_dt'], 'a': r['dng_da']})
        for readout in ('ring', 'mzi'):
            probes.append(dict(id=f"strip-{r['polarization']}-{r['width_nm']}-{readout}",
                               geometry='strip', pol=r['polarization'], width=r['width_nm'],
                               slab=None, readout=readout, sens=sens))
    return probes


def load_rib_probes(mesh=MESH):
    with open(os.path.join(OUT, 'data', f'rib_steps_{mesh}.json')) as f:
        recs = json.load(f)
    probes = []
    for r in recs:
        dv = r['derivs']
        sens = dict(neff={k: dv[f'{k}|1.0']['dneff'] for k in 'wtad'},
                    ng={k: dv[f'{k}|1.0']['dng'] for k in 'wtad'})
        # existing (slab-fixed) thickness derivative: t and d move together
        for obs in ('neff', 'ng'):
            sens[obs]['t_slabfixed'] = sens[obs]['t'] + sens[obs]['d']
        for readout in ('ring', 'mzi'):
            probes.append(dict(id=f"rib-TE-{r['width_nm']}-s{r['slab_nm']}-{readout}",
                               geometry='rib', pol='TE', width=r['width_nm'],
                               slab=r['slab_nm'], d_nom=r['d_nm'], readout=readout, sens=sens))
    return probes


def param_names(case, probes, dmode):
    base = {'orig3': ['w', 't', 'a'], 'A_dknown': ['w', 't', 'a'],
            'A': ['w', 't', 'a'], 'B': ['w', 't', 'a_full', 'a_part'],
            'B_w': ['w', 't', 'a_full', 'a_part', 'w_part']}[case]
    if case in ('orig3', 'A_dknown'):
        return base
    slabs = sorted({p['slab'] for p in probes if p['geometry'] == 'rib'})
    if dmode == 'independent' and len(slabs) > 1:
        return base + [f'd@{s}' for s in slabs]
    return base + ['d']


def sens_vector(p, obs, case, names, dmode):
    """Map one probe's d(obs)/d(param) onto the case's parameter columns."""
    s = p['sens'][obs]
    v = np.zeros(len(names))
    col = {n: i for i, n in enumerate(names)}
    if p['geometry'] == 'strip':
        v[col['w']] = s['w']
        v[col['t']] = s['t']
        v[col['a_full' if case.startswith('B') else 'a']] = s['a']
        return v
    # rib
    v[col['w_part' if case == 'B_w' else 'w']] = s['w']
    v[col['t']] = s['t_slabfixed'] if case == 'orig3' else s['t']
    v[col['a_part' if case.startswith('B') else 'a']] = s['a']
    if case in ('orig3', 'A_dknown'):
        return v
    if 'd' in col:
        scale = p['d_nom'] / D_REF if dmode == 'scaled' else 1.0
        v[col['d']] = s['d'] * scale
    else:
        v[col[f"d@{p['slab']}"]] = s['d']
    return v


def jacobian(probes, case, dmode='common', rib_noise=1.0):
    names = param_names(case, probes, dmode)
    rows = []
    for p in probes:
        m = rib_noise if p['geometry'] == 'rib' else 1.0
        if p['readout'] == 'ring':
            rows.append(sens_vector(p, 'neff', case, names, dmode) / (fv.SIGMA_NEFF * m))
            rows.append(sens_vector(p, 'ng', case, names, dmode) / (fv.SIGMA_NG * m))
        else:
            k = 2 * np.pi / fv.WL0 * fv.L_MZI
            rows.append(k * sens_vector(p, 'neff', case, names, dmode) / (fv.SIGMA_PHASE * m))
    return np.array(rows), names


def fisher(probes, case, prior, dmode='common', rib_noise=1.0):
    """F_total and parameter names. prior='known' is handled by the caller
    choosing case A_dknown; here prior is None or sigma_prior in nm."""
    J, names = jacobian(probes, case, dmode, rib_noise)
    F = J.T @ J
    if prior is not None:
        for i, n in enumerate(names):
            if n.startswith('d'):
                F[i, i] += 1.0 / prior ** 2
    return F, names


def marginal_std(F, idx, rel_tol=1e-12):
    w, V = np.linalg.eigh(F)
    keep = w > rel_tol * w.max()
    null = V[:, ~keep]
    if null.size and np.max(np.abs(null[idx, :])) > 1e-8:
        return np.inf                      # target not identifiable
    Finv = (V[:, keep] / w[keep]) @ V[:, keep].T
    return float(np.sqrt(Finv[idx, idx]))


def swa_precision(probes, case, prior, dmode='common', rib_noise=1.0):
    """Returns {swa_name: sigma_deg}. 'known' -> d removed exactly."""
    if prior == 'known':
        case = {'A': 'A_dknown'}.get(case, case)
        if case.startswith('B'):
            # d known in case B: drop d columns, keep everything else
            F, names = fisher(probes, case, 1e-9, dmode, rib_noise)
            keep = [i for i, n in enumerate(names) if not n.startswith('d')]
            F, names = F[np.ix_(keep, keep)], [names[i] for i in keep]
        else:
            F, names = fisher(probes, case, None, dmode, rib_noise)
    else:
        F, names = fisher(probes, case, prior, dmode, rib_noise)
    targets = [n for n in names if n.startswith('a')]
    return {n: marginal_std(F, names.index(n)) for n in targets}


# ---------------------------------------------------------------------------
# structure sets
# ---------------------------------------------------------------------------
def pick(pool, **kw):
    out = [p for p in pool if all(p.get(k) == v for k, v in kw.items())]
    assert out, f'no probe matches {kw}'
    return out


def structure_sets(strip, rib):
    baseline = list(strip)        # strip TE+TM x 5 widths x (ring+MZI): 20 probes, 30 rows
    ring = lambda w, s: pick(rib, width=w, slab=s, readout='ring')
    sets = {
        '1 baseline (strip TE/TM rings+MZIs)': baseline,
        '2 +1 rib ring (500nm, slab 90)': baseline + ring(500, 90),
        '3a +2 rib rings, slabs 70/110': baseline + ring(500, 70) + ring(500, 110),
        '3b +3 rib rings, slabs 50/90/130': baseline + ring(500, 50) + ring(500, 90) + ring(500, 130),
        '3c +3 rib rings, slabs 70/90/110': baseline + ring(500, 70) + ring(500, 90) + ring(500, 110),
        'ctrl +3 rib rings, slab 90, w 400/500/600': baseline + ring(400, 90) + ring(500, 90) + ring(600, 90),
        'ref existing rib family (5 w, slab 90, ring+MZI)': baseline + pick(rib, slab=90),
        'min 2 strip TE rings (550, 600)': pick(strip, pol='TE', width=550, readout='ring')
                                          + pick(strip, pol='TE', width=600, readout='ring'),
        'min +1 rib ring (400nm, slab 90)': pick(strip, pol='TE', width=550, readout='ring')
                                           + pick(strip, pol='TE', width=600, readout='ring') + ring(400, 90),
    }
    return sets


# ---------------------------------------------------------------------------
# greedy budget-limited design
# ---------------------------------------------------------------------------
def design_value(F, criterion, names):
    w = np.linalg.eigvalsh(F)
    if criterion == 'E':
        return w[0]                                     # maximize lambda_min
    if w[0] <= 1e-12 * w[-1]:
        return -np.inf
    Finv = np.linalg.inv(F)
    if criterion == 'A':
        return -np.trace(Finv)                          # minimize trace(F^-1)
    idx = [i for i, n in enumerate(names) if n.startswith('a')]
    return -sum(Finv[i, i] for i in idx)                # 'SWA': minimize SWA variance(s)


def greedy(pool, case, prior, criterion, n_max=8, n_seed=3, rib_noise=1.0):
    """Exhaustive search for the best n_seed-probe set (sidesteps the
    all-singular tie at the first picks, where lambda_min = 0 for every
    candidate), then greedy one-at-a-time augmentation to n_max. Uses the
    'common' d model."""
    # Each probe's Fisher contribution. Columns that NO probe in the pool can
    # inform (d and a_part for a strip-only pool) are dropped: they are exactly
    # decoupled, and keeping them would pin lambda_min at 0 (or at the prior).
    names = param_names(case, pool, 'common')
    contrib = []
    for p in pool:
        J, _ = jacobian([p], case, 'common', rib_noise)
        contrib.append(J.T @ J)
    contrib = np.array(contrib)
    live = [i for i in range(len(names)) if np.any(contrib[:, i, i] > 0)]
    names = [names[i] for i in live]
    contrib = contrib[:, live][:, :, live]
    P0 = np.zeros((len(names), len(names)))
    if prior is not None:
        for i, n in enumerate(names):
            if n.startswith('d'):
                P0[i, i] = 1.0 / prior ** 2

    best, best_val = None, -np.inf
    for combo in itertools.combinations(range(len(pool)), n_seed):
        v = design_value(P0 + contrib[list(combo)].sum(0), criterion, names)
        if v > best_val:
            best, best_val = list(combo), v
    chosen = list(best)
    history = {n_seed: list(chosen)}
    for n in range(n_seed + 1, n_max + 1):
        cand = [(design_value(P0 + contrib[chosen + [j]].sum(0), criterion, names), j)
                for j in range(len(pool)) if j not in chosen]
        chosen.append(max(cand)[1])
        history[n] = list(chosen)
    return history


# ---------------------------------------------------------------------------
# main analysis
# ---------------------------------------------------------------------------
def fmt(x):
    return 'inf (not identifiable)' if not np.isfinite(x) else f'{x:.4f}'


def run():
    os.makedirs(os.path.join(OUT, 'figures'), exist_ok=True)
    strip = load_strip_probes()
    rib = load_rib_probes()
    sets = structure_sets(strip, rib)
    rows = []

    # --- consistency check against existing angle_crb.py (B): strip baseline
    # 0.0210 deg, + rib family (cons.) 0.0009 deg, in the original model
    base = swa_precision(sets['1 baseline (strip TE/TM rings+MZIs)'], 'orig3', None)['a']
    fam = swa_precision(sets['ref existing rib family (5 w, slab 90, ring+MZI)'], 'orig3', None)['a']
    print(f'reproduce angle_crb.py (B): baseline {base:.4f} deg (expect 0.0210), '
          f'+rib family {fam:.4f} deg (expect 0.0009)')

    # --- main table: every case x prior x set x d-model x rib noise
    for noise_name, mult in NOISE.items():
        for set_name, probes in sets.items():
            n_slabs = len({p['slab'] for p in probes if p['geometry'] == 'rib'})
            v = swa_precision(probes, 'orig3', None, rib_noise=mult)['a']
            rows.append(dict(noise=noise_name, set=set_name, case='orig3 (slab known)', prior='-',
                             dmode='-', swa='a', sigma_deg=v))
            for case in ('A', 'B', 'B_w'):
                dmodes = ['common', 'scaled', 'independent'] if n_slabs > 1 else ['common']
                for dmode in dmodes:
                    for prior in PRIORS:
                        res = swa_precision(probes, case, prior, dmode, mult)
                        for swa, val in res.items():
                            rows.append(dict(noise=noise_name, set=set_name, case=case,
                                             prior=PRIOR_LABEL[prior], dmode=dmode, swa=swa,
                                             sigma_deg=val))
    write_csv(rows, 'swa_precision_table.csv')

    # --- greedy designs
    greedy_rows = []
    pools = {'strip only': strip, 'strip + rib': strip + rib}
    for case in ('A', 'B'):
        for prior in (None, 5.0, 2.0):
            for pool_name, pool in pools.items():
                for crit in ('E', 'A', 'SWA'):
                    hist = greedy(pool, case, prior, crit)
                    for n, idx in hist.items():
                        chosen = [pool[i] for i in idx]
                        res = swa_precision(chosen, case, prior)
                        for swa, val in res.items():
                            greedy_rows.append(dict(case=case, prior=PRIOR_LABEL[prior], pool=pool_name,
                                                    criterion=crit, N=n, swa=swa, sigma_deg=val,
                                                    chosen=' | '.join(p['id'] for p in chosen)))
    write_csv(greedy_rows, 'greedy_designs.csv')

    # --- eigen-spectra for the main sets
    spectra = {}
    for set_name in list(sets)[:7]:
        for case, prior in (('A', None), ('A', 5.0), ('B', None)):
            F, names = fisher(sets[set_name], case, prior)
            w, V = np.linalg.eigh(F)      # raw units (nm, deg), as in the existing code
            spectra[(set_name, case, PRIOR_LABEL[prior])] = (w, V, names)
    with open(os.path.join(OUT, 'data', 'eigen_spectra.json'), 'w') as f:
        json.dump([dict(set=k[0], case=k[1], prior=k[2], names=v[2],
                        eigenvalues=v[0].tolist(), eigenvectors=v[1].T.tolist())
                   for k, v in spectra.items()], f, indent=1)

    make_figures(rows, greedy_rows, spectra, sets)
    return rows, greedy_rows, spectra


def write_csv(rows, name):
    import csv
    path = os.path.join(OUT, 'data', name)
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            r = dict(r)
            if isinstance(r.get('sigma_deg'), float):
                r['sigma_deg'] = 'inf' if not np.isfinite(r['sigma_deg']) else f"{r['sigma_deg']:.6g}"
            w.writerow(r)
    print(f'wrote {path} ({len(rows)} rows)')


# ---------------------------------------------------------------------------
# figures
# ---------------------------------------------------------------------------
def style(ax):
    ax.grid(True, color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    for s in ('left', 'bottom'):
        ax.spines[s].set_color(INK2)
    ax.tick_params(colors=INK2, labelsize=8)


def lookup(rows, **kw):
    out = [r for r in rows if all(r[k] == v for k, v in kw.items())]
    assert len(out) == 1, (kw, len(out))
    return out[0]['sigma_deg']


def make_figures(rows, greedy_rows, spectra, sets):
    plt.rcParams.update({'font.size': 9, 'axes.titlesize': 10, 'axes.labelcolor': INK})
    prior_order = ['1 nm', '2 nm', '5 nm', '10 nm', 'no prior']
    xs = np.arange(len(prior_order))

    # Fig 1: SWA precision vs prior on d
    plot_sets = ['1 baseline (strip TE/TM rings+MZIs)', '2 +1 rib ring (500nm, slab 90)',
                 '3b +3 rib rings, slabs 50/90/130', 'ctrl +3 rib rings, slab 90, w 400/500/600',
                 'ref existing rib family (5 w, slab 90, ring+MZI)']
    panels = [('A', 'a', 'Case A: shared SWA'), ('B', 'a_full', 'Case B: SWA of full etch (strips)'),
              ('B', 'a_part', 'Case B: SWA of partial etch (ribs)')]
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2), sharey=False)
    for ax, (case, swa, title) in zip(axes, panels):
        for i, sname in enumerate(plot_sets):
            ys = [lookup(rows, noise='conservative', set=sname, case=case, prior=pl,
                         dmode='common', swa=swa) for pl in prior_order]
            ys = np.array(ys, float)
            ok = np.isfinite(ys)
            if not ok.any():
                continue
            ax.plot(xs[ok], ys[ok], '-o', color=C[i], lw=2, ms=5, label=sname)
            known = lookup(rows, noise='conservative', set=sname, case=case, prior='d known',
                           dmode='common', swa=swa)
            if np.isfinite(known):
                ax.axhline(known, color=C[i], lw=1, ls=':', alpha=0.9)
        ax.set_yscale('log')
        ax.set_xticks(xs, prior_order)
        ax.set_xlabel('prior on etch depth d (1-sigma)')
        ax.set_title(title, color=INK)
        style(ax)
    axes[0].set_ylabel('SWA precision, 1-sigma CRB [deg]')
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, loc='lower center', ncol=3, frameon=False, fontsize=8)
    fig.text(0.5, 0.93, 'Dotted lines: same set with d known exactly.  Rib noise: conservative (1.0x strip).  '
             'Panels with a missing set: SWA not identifiable for that set.', ha='center', color=INK2, fontsize=8)
    fig.tight_layout(rect=(0, 0.12, 1, 0.92))
    fig.savefig(os.path.join(OUT, 'figures', 'fig_swa_vs_prior.png'), dpi=200)
    plt.close(fig)

    # Fig 2: SWA precision vs budget
    fig, axes = plt.subplots(2, 3, figsize=(13, 7.5), sharex=True)
    series = [('strip + rib', 'E', C[0], '-'), ('strip + rib', 'A', C[1], '-'),
              ('strip + rib', 'SWA', C[2], '-'),
              ('strip only', 'E', C[0], '--'), ('strip only', 'A', C[1], '--'),
              ('strip only', 'SWA', C[2], '--')]
    for r_i, prior in enumerate(['no prior', '2 nm']):
        for c_i, (case, swa, title) in enumerate(panels):
            ax = axes[r_i, c_i]
            for pool, crit, col, ls in series:
                pts = sorted((g['N'], g['sigma_deg']) for g in greedy_rows
                             if g['case'] == case and g['prior'] == prior and g['pool'] == pool
                             and g['criterion'] == crit and g['swa'] == swa)
                if not pts:
                    continue
                n, y = np.array(pts, float).T
                ok = np.isfinite(y)
                ax.plot(n[ok], y[ok], ls, marker='o', color=col, lw=2, ms=5,
                        label=f'{pool}, {crit}-optimal' if crit != 'SWA' else f'{pool}, SWA-targeted (ref.)')
            ax.set_yscale('log')
            ax.set_title(f'{title}\nprior on d: {prior}', color=INK)
            style(ax)
            if r_i == 1:
                ax.set_xlabel('structure budget N (rings and MZIs each count 1)')
        axes[r_i, 0].set_ylabel('SWA precision, 1-sigma CRB [deg]')
    h, l = axes[0, 0].get_legend_handles_labels()
    fig.legend(h, l, loc='lower center', ncol=3, frameon=False, fontsize=8)
    fig.text(0.5, 0.965, 'Greedy design (exhaustive best 3, then greedy). Solid: pool includes ribs; dashed: strip-only pool. '
             'Missing points: SWA not identifiable.', ha='center', color=INK2, fontsize=8)
    fig.tight_layout(rect=(0, 0.07, 1, 0.95))
    fig.savefig(os.path.join(OUT, 'figures', 'fig_swa_vs_budget.png'), dpi=200)
    plt.close(fig)

    # Fig 3: eigenvalue spectra, colored by |SWA component| of each eigenvector
    set_names = list(sets)[:7]
    cases = [('A', 'no prior'), ('A', '5 nm'), ('B', 'no prior')]
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.8), sharey=True)
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list('seq', ['#dbe8f8', '#2a78d6', '#0d2f5c'])
    floor = 1e-6
    for ax, (case, pl) in zip(axes, cases):
        for i, sname in enumerate(set_names):
            w, V, names = spectra[(sname, case, pl)]
            aidx = [j for j, n in enumerate(names) if n.startswith('a')]
            swa_w = np.sqrt((V[aidx, :] ** 2).sum(0))
            y = np.where(w > 1e-12 * w.max(), w, floor)
            sc = ax.scatter(np.full(len(w), i), y, c=swa_w, cmap=cmap, vmin=0, vmax=1, s=60,
                            edgecolors=INK2, linewidths=0.5, zorder=3)
        ax.axhline(floor, color=INK2, lw=0.8, ls=':')
        ax.set_yscale('log')
        ax.set_xticks(range(len(set_names)), [s.split(' ', 1)[1][:28] for s in set_names],
                      rotation=60, ha='right', fontsize=7)
        ax.set_title(f'Case {case}, prior on d: {pl}', color=INK)
        style(ax)
    axes[0].set_ylabel('eigenvalue of F (raw units: nm, deg)  [dotted = 0]')
    cb = fig.colorbar(sc, ax=axes, shrink=0.8, pad=0.01)
    cb.set_label('|SWA component| of eigenvector', color=INK2)
    fig.savefig(os.path.join(OUT, 'figures', 'fig_eigen_spectra.png'), dpi=200, bbox_inches='tight')
    plt.close(fig)
    print('wrote figures to', os.path.join(OUT, 'figures'))


if __name__ == '__main__':
    run()
