"""
Etch-depth nuisance study: numerical checks on the new rib derivatives.

Reads results_etch_depth/data/{rib_steps_<MESH>,mesh_check,dlam_check}.json
(written by etch_depth_sims.py) and writes pass/fail tables as CSV next to
them. Criterion (as specified): 2% relative agreement.

  step check : D(h/2) and D(2h) vs D(h), every rib structure, every parameter
  mesh check : D on each finer mesh vs the next coarser one (M0->M1->M2->M3),
               and larger domain D1 vs M0
  dlam check : n_g-derivatives at dlam = 2.5 and 10nm vs 5nm

Relative error is ill-conditioned where a derivative passes through zero
(e.g. the existing rib dn_g/dw at 500nm is 3e-6, ~1000x below the other
n_g derivatives of that structure). Such entries are reported with their
magnitude relative to the largest same-observable derivative of the same
structure ("rel_size") so that a failure can be judged, not hidden.

Run:  python etch_depth_checks.py [--mesh M0]
"""
import os
import sys
import json
import csv

import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(REPO, 'results_etch_depth', 'data')
TOL = 0.02
PARAMS = ('w', 't', 'a', 'd')


def load(name):
    with open(os.path.join(DATA, name)) as f:
        return json.load(f)


def rel(x, ref):
    return abs(x - ref) / abs(ref) if ref != 0 else np.inf


def rel_size(rec, obs, k, f='1.0'):
    vals = [abs(rec['derivs'][f'{p}|{f}'][obs]) for p in PARAMS if f'{p}|{f}' in rec['derivs']]
    return abs(rec['derivs'][f'{k}|{f}'][obs]) / max(vals)


def step_table(mesh):
    rows = []
    for r in load(f'rib_steps_{mesh}.json'):
        for k in PARAMS:
            for obs in ('dneff', 'dng'):
                d = {f: r['derivs'][f'{k}|{f}'][obs] for f in ('0.5', '1.0', '2.0')}
                err = max(rel(d['0.5'], d['1.0']), rel(d['2.0'], d['1.0']))
                rows.append(dict(width_nm=r['width_nm'], slab_nm=r['slab_nm'], param=k, obs=obs,
                                 D_half=d['0.5'], D_h=d['1.0'], D_2h=d['2.0'],
                                 max_rel_dev=err, rel_size=rel_size(r, obs, k),
                                 passed=err < TOL))
    return rows


def mesh_table():
    recs = load('mesh_check.json')
    by = {(r['mesh'], r['width_nm'], r['slab_nm']): r for r in recs}
    structs = sorted({(r['width_nm'], r['slab_nm']) for r in recs})
    pairs = [('M0', 'M1'), ('M1', 'M2'), ('M2', 'M3'), ('M0', 'D1')]
    rows = []
    for (w, s) in structs:
        for k in PARAMS:
            for obs in ('dneff', 'dng'):
                vals = {m: by[(m, w, s)]['derivs'][f'{k}|1.0'][obs] for m in ('M0', 'M1', 'M2', 'M3', 'D1')}
                row = dict(width_nm=w, slab_nm=s, param=k, obs=obs, **{f'D_{m}': v for m, v in vals.items()},
                           rel_size_M3=rel_size(by[('M3', w, s)], obs, k))
                for a, b in pairs:
                    row[f'chg_{a}->{b}'] = rel(vals[b], vals[a])
                rows.append(row)
    return rows


def dlam_table():
    recs = load('dlam_check.json')
    by = {(r['dlam_um'], r['width_nm'], r['slab_nm']): r for r in recs}
    rows = []
    for (w, s) in sorted({(r['width_nm'], r['slab_nm']) for r in recs}):
        ref = by[(0.005, w, s)]
        for k in PARAMS:
            d = {dl: by[(dl, w, s)]['derivs'][f'{k}|1.0']['dng'] for dl in (0.0025, 0.005, 0.01)}
            err = max(rel(d[0.0025], d[0.005]), rel(d[0.01], d[0.005]))
            rows.append(dict(width_nm=w, slab_nm=s, param=k, ng=ref['ng'],
                             dng_dlam2p5=d[0.0025], dng_dlam5=d[0.005], dng_dlam10=d[0.01],
                             max_rel_dev=err, rel_size=rel_size(ref, 'dng', k), passed=err < TOL))
    return rows


def write(rows, name):
    path = os.path.join(DATA, name)
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f'wrote {path}')


if __name__ == '__main__':
    mesh = sys.argv[sys.argv.index('--mesh') + 1] if '--mesh' in sys.argv else 'M0'
    for fn, name, tab in ((step_table, f'step_check_{mesh}.csv', 'step'),
                          (mesh_table, 'mesh_check.csv', 'mesh'),
                          (dlam_table, 'dlam_check.csv', 'dlam')):
        try:
            rows = fn(mesh) if tab == 'step' else fn()
        except FileNotFoundError as e:
            print(f'[{tab}] skipped: {e.filename} not found')
            continue
        write(rows, name)
        if 'passed' in rows[0]:
            bad = [r for r in rows if not r['passed']]
            print(f'[{tab}] {len(rows) - len(bad)}/{len(rows)} pass the {TOL:.0%} criterion')
            for r in bad:
                print(f"   FAIL {r['width_nm']}nm slab{r['slab_nm']} {r['param']} "
                      f"{r.get('obs', 'dng')}: dev={r['max_rel_dev']:.3%} rel_size={r['rel_size']:.3g}")
