"""
Clean reproduction: copy src/ and only the raw solve caches (plus the declared frozen inputs) into a fresh tree,
rerun the whole pipeline there (derivatives -> Fisher/CRB -> search -> Monte Carlo -> robustness -> summaries),
and compare every regenerated file with the committed one.

Every script finds the repository from its own location, so the copied scripts read and write only inside the
tree. The tree (results_final/reproduction/tree) is not committed; the comparison is.

Run from src/ (main environment):  python reproduce_all.py [--compare-only]
Writes results_final/reproduction/{diff.json, steps.log}
"""
import csv
import json
import math
import os
import re
import shutil
import subprocess
import sys
import time

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BASE = os.path.join(REPO, 'results_final', 'reproduction')
TREE = os.path.join(BASE, 'tree')
VENV = 'C:/Users/unmol/.venvs/femwell/Scripts/python.exe'
REL_TOL = 1e-9
CACHES = ['results_femwell/cache/solves_r12.jsonl', 'results_femwell/cache/solves_r8.jsonl',
          'results_etchdepth_femwell/cache/solves_r12.jsonl', 'results_montecarlo/cache/solves_r12.jsonl',
          'results_montecarlo/fixedmesh/cache/solves_r12.jsonl']
FROZEN = {
    'results_priors_search/data/search.json': 'input of strips_angle_info.py; written by the pre-symmetry priors_search.py run, '
        'whose structure pool has since grown, so it cannot be regenerated as committed (its numbers are retired)',
    'results_corrected/data/sensitivity_db.json': 'EMpy database of an earlier session (EMpy is not rerun); input of claims_checks.py',
    'results_corrected/data/step_check.csv': 'EMpy step check of an earlier session; input of claims_checks.py',
    'results_femwell/data/sensitivity_db.json': 'femwell database of an earlier session (femwell_postprocess.py also reads '
        'EMpy-era files); input of claims_checks.py',
    'results_final/data/citations.json': 'quotes extracted from the fetched PDFs, which are not in the repository (copyright)',
}
M = 'python'
STEPS = [
    (M, 'step_check_all.py', {}),
    (M, 'sym_search.py', {'SYM_SEARCH_SLABS': '50,70,90,110,130,150'}),
    (M, 'sym_search.py placeholder ../results_final/sym_search_placeholder/data', {}),
    (M, 'sym_search.py 1pm ../results_final/sym_search_1pm/data', {}),
    (M, 'strips_angle_info.py', {}),
    (M, 'montecarlo.py plan', {}),
    (M, 'montecarlo.py fit', {}),
    (M, 'montecarlo.py surrcheck', {}),
    (M, 'montecarlo.py mc', {'MC_NPROC': '4'}),
    (M, 'montecarlo.py mc1pm', {'MC_NPROC': '4'}),
    (M, 'montecarlo.py valley', {}),
    (M, 'robustness_rib300.py', {}),
    (M, 'robustness_rib300.py 1pm', {}),
    (M, 'noise_sourced.py', {}),
    (VENV, 'crosscheck_fields.py', {}),
    (M, 'crosschecks.py', {}),
    (M, 'assumption_sensitivity.py', {}),
    (M, 'claims_checks.py', {}),
    (M, 'record_environment.py', {}),
    (M, 'make_headline_numbers.py', {}),
    (M, 'make_final_docs.py', {}),
]
COMPARE = [   # regenerated files compared with the committed ones
    'results_montecarlo/sym_search/data/step_check_all.csv',
    *[f'results_montecarlo/sym_search/data/{f}' for f in ('priors.json', 'priors_table.csv', 'd_prior_sweep.json', 'search.json',
                                                           'search_table.csv', 'slab_pairs.json', 'summary_tables.md')],
    *[f'results_final/sym_search_{n}/data/{f}' for n in ('placeholder', '1pm')
      for f in ('priors.json', 'priors_table.csv', 'd_prior_sweep.json', 'search.json', 'search_table.csv', 'slab_pairs.json', 'summary_tables.md')],
    'results_montecarlo/data/strips_info.json', 'results_montecarlo/data/strips_info_perturbation.csv',
    'results_montecarlo/data/mc_design.json', 'results_montecarlo/data/surrogates.json',
    'results_montecarlo/data/surrogate_accuracy.json', 'results_montecarlo/data/surrogate_accuracy.csv',
    'results_montecarlo/data/surrogate_accuracy_1pm.json',
    'results_montecarlo/data/mc_results.json', 'results_montecarlo/data/mc_fits.csv',
    'results_montecarlo/data/mc_results_1pm.json', 'results_montecarlo/data/mc_fits_1pm.csv', 'results_montecarlo/data/mc_valley.json',
    'results_montecarlo/robustness/data/robustness.json', 'results_montecarlo/robustness/data/robustness_1pm.json',
    'results_montecarlo/noise_sourced/data/noise_sourced.json',
    'results_final/data/crosscheck_fields.json', 'results_final/data/crosschecks.json',
    'results_final/data/assumption_sensitivity.json', 'results_final/data/claims_checks.json',
    'results_final/HEADLINE_NUMBERS.md', 'results_final/CLAIMS.md', 'results_final/ASSUMPTIONS.md',
]
NUM = re.compile(r'(?<![\w.])[-−+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?(?![\w])')


def setup():
    if os.path.exists(TREE):
        shutil.rmtree(TREE)
    os.makedirs(os.path.join(TREE, 'src'))
    for f in os.listdir(os.path.join(REPO, 'src')):
        if f.endswith('.py'):
            shutil.copy2(os.path.join(REPO, 'src', f), os.path.join(TREE, 'src', f))
    for f in CACHES + list(FROZEN):
        os.makedirs(os.path.dirname(os.path.join(TREE, f)), exist_ok=True)
        shutil.copy2(os.path.join(REPO, f), os.path.join(TREE, f))
    for f in COMPARE:                          # empty output folders (some scripts expect them to exist)
        os.makedirs(os.path.dirname(os.path.join(TREE, f)), exist_ok=True)


def run_steps():
    out = []
    log = open(os.path.join(BASE, 'steps.log'), 'w', encoding='utf-8')
    for py, cmd, env in STEPS:
        t0 = time.time()
        e = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', PYTHONIOENCODING='utf-8', **env)
        p = subprocess.run([py] + cmd.split(), cwd=os.path.join(TREE, 'src'), env=e, capture_output=True, text=True)
        rec = dict(cmd=(' '.join(f'{k}={v}' for k, v in env.items() if k != 'MC_NPROC') + ' ' + cmd).strip(), exit=p.returncode,
                   seconds=time.time() - t0)
        out.append(rec)
        log.write(f"### {rec['cmd']} (exit {p.returncode}, {rec['seconds']:.0f} s)\n{p.stdout[-4000:]}\n{p.stderr[-4000:]}\n")
        log.flush()
        print(rec, flush=True)
        if p.returncode:
            break
    return out


def numbers(path):
    """(location, value) pairs: JSON leaves, CSV cells, or numeric tokens of markdown (code spans removed)."""
    if path.endswith('.json'):
        out = []

        def walk(o, where):
            if isinstance(o, bool) or o is None:
                return
            if isinstance(o, (int, float)):
                out.append((where, float(o)))
            elif isinstance(o, str):
                out.extend((f'{where}#{i}', float(t.replace('−', '-'))) for i, t in enumerate(NUM.findall(o)))
            elif isinstance(o, dict):
                for k, v in o.items():
                    walk(v, f'{where}.{k}')
            elif isinstance(o, list):
                for i, v in enumerate(o):
                    walk(v, f'{where}[{i}]')
        walk(json.load(open(path, encoding='utf-8')), '')
        return out
    if path.endswith('.csv'):
        out = []
        with open(path, newline='', encoding='utf-8') as f:
            for i, row in enumerate(csv.reader(f)):
                for j, c in enumerate(row):
                    try:
                        out.append((f'r{i}c{j}', float(c)))
                    except ValueError:
                        out.extend((f'r{i}c{j}#{k}', float(t)) for k, t in enumerate(NUM.findall(c)))
        return out
    text = re.sub(r'`[^`]*`', ' ', open(path, encoding='utf-8').read())
    return [(f'#{i}', float(t.replace('−', '-'))) for i, t in enumerate(NUM.findall(text))]


def same(a, b):
    if math.isnan(a) and math.isnan(b):
        return True, 0.0
    if math.isinf(a) or math.isinf(b):
        return a == b, (0.0 if a == b else math.inf)
    if a == b:
        return True, 0.0
    d = abs(a - b) / max(abs(a), abs(b), 1e-300)
    return d <= REL_TOL, d


def compare():
    files, diffs = [], []
    for f in COMPARE:
        a, b = os.path.join(REPO, f), os.path.join(TREE, f)
        if not os.path.exists(b):
            files.append(dict(file=f, n_numbers=0, result='MISSING in rerun', max_rel_diff=None))
            continue
        na, nb = numbers(a), numbers(b)
        extra = None
        if f.endswith('.json'):                # match JSON by key path; keys present on one side only are listed
            da, db = dict(na), dict(nb)
            common = [k for k in da if k in db]
            extra = dict(only_committed=sorted(set(da) - set(db))[:20], only_rerun=sorted(set(db) - set(da))[:20])
            na, nb = [(k, da[k]) for k in common], [(k, db[k]) for k in common]
        if len(na) != len(nb):
            files.append(dict(file=f, n_numbers=len(na), result=f'DIFFERENT count ({len(na)} vs {len(nb)})', max_rel_diff=None))
            diffs.append(dict(file=f, where='(number count)', committed=len(na), rerun=len(nb)))
            continue
        worst, exact, bad = 0.0, True, 0
        for (wa, va), (wb, vb) in zip(na, nb):
            ok, d = same(va, vb)
            if d > 0:
                exact = False
            worst = max(worst, d)
            if not ok:
                bad += 1
                if bad <= 200:                 # keep the file small; the count is complete
                    diffs.append(dict(file=f, where=wa, committed=va, rerun=vb))
        res = 'identical' if exact else ('within tolerance' if bad == 0 else f'{bad} numbers differ')
        if extra and (extra['only_committed'] or extra['only_rerun']):
            res += (f"; keys only in rerun: {', '.join(extra['only_rerun'])}" if extra['only_rerun'] else '') +                    (f"; keys only in committed: {', '.join(extra['only_committed'])}" if extra['only_committed'] else '')
        files.append(dict(file=f, n_numbers=len(na), max_rel_diff=worst, result=res))
    return files, diffs


def sampling_check():
    """strips_info perturbation statistics: committed vs rerun, for the medians and 5-95% points (random draws)."""
    f = 'results_montecarlo/data/strips_info.json'
    a = json.load(open(os.path.join(REPO, f), encoding='utf-8'))['perturbation']
    b = json.load(open(os.path.join(TREE, f), encoding='utf-8'))['perturbation']
    rel = {'nominal': [], 'median': [], 'p05': [], 'p95': []}
    for lv in a:
        for name in a[lv]:
            for k in rel:
                x, y = a[lv][name][k], b[lv][name][k]
                if x:
                    rel[k].append(abs(y / x - 1))
    return dict(file=f, seed=20261027, max_abs_rel_diff={k: max(v) for k, v in rel.items()},
                cause='The perturbation draws (default_rng(20261027)) perturb every strip derivative entry in the cache. '
                      'strips_info.json was committed (3e828dd) before the strip TM index solves were added (6b84f2d), so '
                      'the rerun draws more entries and the random sequence differs. Deterministic values reproduce exactly; '
                      'the sampled statistics agree within sampling error (see max_abs_rel_diff).')


def main():
    t0 = time.time()
    if '--compare-only' not in sys.argv:
        setup()
        steps = run_steps()
    else:
        steps = json.load(open(os.path.join(BASE, 'diff.json'), encoding='utf-8'))['steps']
    files, diffs = compare()
    import montecarlo
    res = dict(tree=os.path.relpath(TREE, REPO).replace(os.sep, '/'), rel_tol=REL_TOL, mc_seed=montecarlo.SEED,
               frozen=FROZEN, caches=CACHES, steps=steps, files=files, differences=diffs,
               n_files=len(files), n_exact=sum(f['result'].startswith('identical') for f in files),
               n_within_tol=sum(f['result'].startswith('within tolerance') for f in files),
               n_differ=sum(('differ' in f['result']) or f['result'].startswith('DIFFERENT') for f in files),
               n_missing=sum(f['result'].startswith('MISSING') for f in files), minutes=(time.time() - t0) / 60,
               sampling_notes=[sampling_check()])
    with open(os.path.join(BASE, 'diff.json'), 'w', encoding='utf-8') as f:
        json.dump(res, f, indent=1, ensure_ascii=False)
    print({k: res[k] for k in ('n_files', 'n_exact', 'n_within_tol', 'n_differ', 'n_missing')})
    for f in files:
        if f['result'] not in ('identical',):
            print(f)


if __name__ == '__main__':
    main()
