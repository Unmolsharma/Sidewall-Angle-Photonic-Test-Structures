"""
Every downstream analysis on a recomputed sensitivity database, with the noise
assumptions unchanged:

  1. writeup_numbers.main at step h and on the h/2 and 2h variants
     (robustness_checks.variant / slab_variant): Fisher matrices, smallest
     eigenvalues, trace of the inverse, condition numbers, ring-pair search,
     TM ring / MZI additions, direction vectors and angles, Spearman, PCA incl.
     the seven-scaling check, slot noise table, rib conservative/optimistic,
     sidewall-angle CRBs, two rings + one structure, slab sweep
  2. degeneracy ratio R = dn/da / (-h tan1deg dn/dw) for every structure
  3. the original analysis/figure scripts (fisher_reanalysis, atlas, atlas_v2,
     atlas_vs_fisher, atlas_v2_vs_fisher, information_map, stress_tests,
     noise_sourcing, angle_crb) rerun unchanged in a scratch working directory
     whose data/ holds the new database; logs kept, Figures 1-2 copied out
  4. Figure 3: step-check pass share per family, original vs cell-exact vs new

Run from src/, e.g.
  python downstream_analysis.py --data ../results_femwell/data --slab rib_slab_sweep_db_femwell.json \
      --fig ../figures/femwell --label femwell
Writes <data>/writeup_numbers_<label>.json and <data>/robustness_checks.json
"""
import argparse
import csv
import json
import os
import shutil
import subprocess
import sys
import tempfile

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import writeup_numbers as wn
from robustness_checks import variant, slab_variant, degeneracy

SRC = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(SRC)
SCRIPTS = ['fisher_reanalysis', 'atlas', 'atlas_v2', 'atlas_vs_fisher', 'atlas_v2_vs_fisher',
           'information_map', 'stress_tests', 'noise_sourcing', 'angle_crb']
DB_FILES = ('sensitivity_db.json', 'sensitivity_db.csv', 'probe_db.json', 'probe_db.csv')


def numbers(D, slabfile, label):
    db, slab = os.path.join(D, 'sensitivity_db.json'), os.path.join(D, slabfile)
    N = wn.main(db, slab)
    with open(os.path.join(D, f'writeup_numbers_{label}.json'), 'w') as f:
        json.dump(N, f, indent=1, default=float)
    records, sweep = json.load(open(db)), json.load(open(slab))
    with open(os.path.join(D, 'step_check.csv')) as f:
        steprows = {(r['structure'], r['param'], r['output']): {c: float(r[c]) for c in ('D_half', 'D_h', 'D_2h')}
                    for r in csv.DictReader(f)}
    out = {}
    tmp = tempfile.mkdtemp()
    for col in ('D_half', 'D_2h'):
        p1, p2 = os.path.join(tmp, 'db.json'), os.path.join(tmp, 'slab.json')
        json.dump(variant(records, steprows, col), open(p1, 'w'))
        json.dump(slab_variant(sweep, steprows, col), open(p2, 'w'))
        out[col] = wn.main(p1, p2)
    shutil.rmtree(tmp)
    out['degeneracy'] = degeneracy(records)
    with open(os.path.join(D, 'robustness_checks.json'), 'w') as f:
        json.dump({label: out}, f, indent=1, default=float)
    print(f"R = dn/da / (-h tan1deg dn/dw), {label}:")
    for fam in ('strip_TE', 'strip_TM', 'slot_TE', 'slot_TM', 'swg_TE', 'rib_TE'):
        v = [x for k, x in out['degeneracy'].items() if k.startswith(fam)]
        print(f"  {fam:9s} R_neff " + ' '.join(f"{x['R_neff']:.4f}" for x in v) +
              "   R_ng " + ' '.join(f"{x['R_ng']:.4f}" for x in v))


def original_scripts(D, logs, fig):
    work = tempfile.mkdtemp()
    os.makedirs(os.path.join(work, 'data'))
    os.makedirs(os.path.join(work, 'figures'))
    for fn in DB_FILES:
        shutil.copy(os.path.join(D, fn), os.path.join(work, 'data', fn))
    os.makedirs(logs, exist_ok=True)
    env = dict(os.environ, PYTHONPATH=SRC, MPLBACKEND='Agg')
    failed = {}
    for s in SCRIPTS:
        r = subprocess.run([sys.executable, os.path.join(SRC, f'{s}.py')], cwd=work, env=env,
                           capture_output=True, text=True)
        with open(os.path.join(logs, f'{s}.log'), 'w', encoding='utf-8') as f:
            f.write(r.stdout + (('\n--- stderr ---\n' + r.stderr) if r.stderr.strip() else ''))
        print(f"  {s:20s} exit {r.returncode}")
        if r.returncode:  # an original script failing is reported, not hidden; the rest still run
            failed[s] = r.stderr.strip().splitlines()[-1]
    with open(os.path.join(logs, 'failed_scripts.json'), 'w') as f:
        json.dump(failed, f, indent=1)
    os.makedirs(fig, exist_ok=True)
    for fn in ('fig1_sensitivity_atlas.png', 'fig2_information_map.png'):
        if os.path.exists(os.path.join(work, 'figures', fn)):
            shutil.copy(os.path.join(work, 'figures', fn), os.path.join(fig, fn))
    for fn in os.listdir(os.path.join(work, 'data')):
        if fn not in DB_FILES:
            shutil.copy(os.path.join(work, 'data', fn), os.path.join(D, fn))
    shutil.rmtree(work)


def fig3(D, fig, label):
    prev = os.path.join(REPO, 'results_corrected', 'data', 'step_check.csv')
    sets = (('original', prev, 'original', '#c0504d', 'Original averaging'),
            ('cell_exact', prev, 'cell_exact', '#2f6fad', 'Cell-exact averaging'),
            (label, os.path.join(D, 'step_check.csv'), label, '#4f9a55', label.capitalize()))
    fams = [('strip_TE', 'Strip TE'), ('strip_TM', 'Strip TM'), ('slot_TE', 'Slot TE'),
            ('slot_TM', 'Slot TM'), ('rib_TE', 'Rib TE'), ('swg_TE', 'SWG TE')]
    x = range(len(fams))
    fig_, ax = plt.subplots(figsize=(7, 3.4), dpi=200)
    for k, (name, path, avg, col, lab) in enumerate(sets):
        with open(path) as f:
            rows = [r for r in csv.DictReader(f) if r['averaging'] == avg and 'slab' not in r['structure']]
        share = []
        for fam, _ in fams:
            sel = [r for r in rows if r['structure'].startswith(fam)]
            share.append(100 * sum(r['verdict'] != 'FAIL' for r in sel) / len(sel) if sel else 0)
        ax.bar([i + (k - 1) * 0.27 for i in x], share, 0.26, color=col, label=lab)
    ax.set_ylim(0, 105)
    ax.set_xticks(list(x))
    ax.set_xticklabels([f[1] for f in fams], fontsize=8)
    ax.set_ylabel('Derivatives passing the step check (%)', fontsize=9)
    ax.spines[['top', 'right']].set_visible(False)
    ax.legend(frameon=False, fontsize=8, loc='lower right')
    fig_.tight_layout()
    os.makedirs(fig, exist_ok=True)
    fig_.savefig(os.path.join(fig, 'fig3_step_check.png'))
    plt.close(fig_)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--data', required=True)
    ap.add_argument('--slab', required=True)
    ap.add_argument('--fig', required=True)
    ap.add_argument('--label', required=True)
    a = ap.parse_args()
    D = os.path.abspath(a.data)
    numbers(D, a.slab, a.label)
    print("original analysis scripts on the new database:")
    original_scripts(D, os.path.join(os.path.dirname(D), 'logs', f'original_scripts_on_{a.label}_db'), os.path.abspath(a.fig))
    fig3(D, os.path.abspath(a.fig), a.label)
    print(f"figures in {os.path.abspath(a.fig)}")
