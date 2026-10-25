"""
Step 4: recompute every derivative, with the step-size stability check, for
BOTH averagings side by side (cellexact_averaging.SOLVERS).

Structures:
  - the 30 database structures (strip TE/TM, slot TE/TM, rib TE, SWG TE,
    5 widths each), in the same order as data/sensitivity_db.json
  - the rib slab-thickness sweep of rib_slab_sweep.py (slab 50-150 nm,
    widths 300/400/500/600 nm)

Every derivative is a central difference at steps h/2, h and 2h, with h the
original step (5 nm width, 2 nm thickness, 1 deg angle). The D(h) value is
the one used downstream (same step as the original). n_g uses the original
formula n_g = n_eff - lambda*dn_eff/dlambda, dlambda = 5 nm.

Stability criterion (decision 2):
  dev       = max(|D(h/2) - D(h)|, |D(2h) - D(h)|)
  PASS      if dev <= 2% of |D(h)|
  near-zero derivatives (|D(h)| < 10% of the largest |D(h)| among the w, t, a
  derivatives of the same structure and output) instead PASS if
  dev <= 2% of that largest |D(h)| (absolute floor). Everything else FAILS
  and is reported, not dropped.

Outputs (results_corrected/data/):
  sensitivity_db.{json,csv}, probe_db.{json,csv}    cell-exact, same schema as data/
  sensitivity_db_original_rerun.json                 original averaging, rerun here
  rib_slab_sweep_db_{cell_exact,original}.json
  step_check.csv                                     every derivative, both averagings

Run from src/:  python recompute_sensitivities.py
"""
import os
os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('MKL_NUM_THREADS', '1')
import csv
import json
import time
from multiprocessing import Pool

import numpy as np

import fullvector_sensitivity as fv
import slot_waveguide as sw
import cellexact_averaging as ce
from probe_db import save_db, build_probe_db, STRIP_GRID_SETTINGS, SLOT_GRID_SETTINGS
from extend_db import make_record

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, 'results_corrected', 'data')

H = dict(w=5.0, t=2.0, a=1.0)          # original steps
FACTORS = (0.5, 1.0, 2.0)
DLAM = 0.005
REL_TOL = 0.02
NEAR_ZERO_FRAC = 0.10

DB_STRUCTS = ([('strip', 'TE', w) for w in (400, 450, 500, 550, 600)] +
              [('strip', 'TM', w) for w in (400, 450, 500, 550, 600)] +
              [('slot', 'TE', w) for w in (150, 165, 180, 195, 210)] +
              [('slot', 'TM', w) for w in (150, 165, 180, 195, 210)] +
              [('rib', 'TE', w) for w in (400, 450, 500, 550, 600)] +
              [('swg', 'TE', w) for w in (400, 450, 500, 550, 600)])
SLAB_SWEEP = [(s, w) for s in (50, 70, 90, 110, 130, 150) for w in (300, 400, 500, 600)]


def geometries():
    yield ('nom', 0.0)
    for p in ('w', 't', 'a'):
        for f in FACTORS:
            for sgn in (+1, -1):
                yield (p, sgn * f * H[p])


def solve(task):
    avg, geom, pol, W, slab, (p, d), wl = task
    w, t, a = W / 1000, 0.22, 0.0
    if p == 'w':
        w += d / 1000
    elif p == 't':
        t += d / 1000
    elif p == 'a':
        a += d
    f = ce.solve_fn(avg, geom, None if slab is None else slab / 1000)
    return task, f(w, t, a, wl, pol)


def derivs(res, avg, geom, pol, W, slab):
    def nng(g):
        n0, npl, nmi = (res[(avg, geom, pol, W, slab, g, fv.WL0 + s * DLAM)] for s in (0, 1, -1))
        return n0, n0 - fv.WL0 * (npl - nmi) / (2 * DLAM)
    n0, g0 = nng(('nom', 0.0))
    D = {}
    for p in ('w', 't', 'a'):
        for f in FACTORS:
            h = f * H[p]
            (n1, g1), (n2, g2) = nng((p, h)), nng((p, -h))
            D[(p, 'neff', f)] = (n1 - n2) / (2 * h)
            D[(p, 'ng', f)] = (g1 - g2) / (2 * h)
    return n0, g0, D


def as_s(n0, g0, D):
    return dict(neff=n0, ng=g0,
                dneff_dw=D[('w', 'neff', 1.0)], dneff_dt=D[('t', 'neff', 1.0)], dneff_da=D[('a', 'neff', 1.0)],
                dng_dw=D[('w', 'ng', 1.0)], dng_dt=D[('t', 'ng', 1.0)], dng_da=D[('a', 'ng', 1.0)])


def step_rows(avg, sid, D):
    rows = []
    for obs in ('neff', 'ng'):
        big = max(abs(D[(p, obs, 1.0)]) for p in ('w', 't', 'a'))
        for p in ('w', 't', 'a'):
            dh, dhalf, d2h = D[(p, obs, 1.0)], D[(p, obs, 0.5)], D[(p, obs, 2.0)]
            dev = max(abs(dhalf - dh), abs(d2h - dh))
            rel = dev / abs(dh) if dh != 0 else np.inf
            near_zero = abs(dh) < NEAR_ZERO_FRAC * big
            if rel <= REL_TOL:
                verdict = 'PASS'
            elif near_zero and dev <= REL_TOL * big:
                verdict = 'PASS_NEAR_ZERO_ABS'
            else:
                verdict = 'FAIL'
            rows.append(dict(averaging=avg, structure=sid, param=p, output=obs,
                             D_half=dhalf, D_h=dh, D_2h=d2h, max_dev=dev, max_rel_dev=rel,
                             largest_D_same_output=big, frac_of_largest=abs(dh) / big,
                             near_zero=near_zero, dev_over_abs_floor=dev / (REL_TOL * big),
                             verdict=verdict))
    return rows


GRID = dict(strip=STRIP_GRID_SETTINGS, slot=SLOT_GRID_SETTINGS, rib=STRIP_GRID_SETTINGS, swg=STRIP_GRID_SETTINGS)


if __name__ == '__main__':
    t0 = time.time()
    os.makedirs(OUT, exist_ok=True)
    tasks = []
    for avg in ('cell_exact', 'original'):
        units = [(g, p, w, None) for g, p, w in DB_STRUCTS] + [('rib', 'TE', w, s) for s, w in SLAB_SWEEP]
        for geom, pol, W, slab in units:
            for g in geometries():
                for s in (0, 1, -1):
                    tasks.append((avg, geom, pol, W, slab, g, fv.WL0 + s * DLAM))
    print(f"{len(tasks)} mode solves", flush=True)
    res = {}
    with Pool(16) as pool:
        for i, (task, n) in enumerate(pool.imap_unordered(solve, tasks, chunksize=4)):
            res[task] = n
            if (i + 1) % 500 == 0:
                print(f"  {i+1}/{len(tasks)}  {time.time()-t0:.0f}s", flush=True)
    print(f"solves done in {time.time()-t0:.0f}s", flush=True)

    all_step_rows = []
    for avg in ('cell_exact', 'original'):
        records = []
        for geom, pol, W in DB_STRUCTS:
            n0, g0, D = derivs(res, avg, geom, pol, W, None)
            gap = sw.rail_gap_from_width(W / 1000)[1] * 1000 if geom == 'slot' else None
            records.append(make_record(
                geom, pol, W, as_s(n0, g0, D), gap, GRID[geom],
                convergence_status=f"{avg} averaging; recomputed {time.strftime('%Y-%m-%d')} "
                                   "(see results_corrected/data/validation.json)",
                validation_status=f"{avg} averaging; step-size check in results_corrected/data/step_check.csv"))
            all_step_rows += step_rows(avg, records[-1]['id'], D)
        sweep = []
        for s, W in SLAB_SWEEP:
            n0, g0, D = derivs(res, avg, 'rib', 'TE', W, s)
            sweep.append(dict(slab_nm=s, width_nm=W, **as_s(n0, g0, D)))
            all_step_rows += step_rows(avg, f'rib_slab{s}_TE_{W}', D)
        if avg == 'cell_exact':
            save_db(records, os.path.join(OUT, 'sensitivity_db.csv'), os.path.join(OUT, 'sensitivity_db.json'))
            probes = build_probe_db(records)
            save_db(probes, os.path.join(OUT, 'probe_db.csv'), os.path.join(OUT, 'probe_db.json'))
        else:
            with open(os.path.join(OUT, 'sensitivity_db_original_rerun.json'), 'w') as f:
                json.dump(records, f, indent=2)
        with open(os.path.join(OUT, f'rib_slab_sweep_db_{avg}.json'), 'w') as f:
            json.dump(sweep, f, indent=2)

    with open(os.path.join(OUT, 'step_check.csv'), 'w', newline='') as f:
        wr = csv.DictWriter(f, fieldnames=list(all_step_rows[0].keys()))
        wr.writeheader()
        wr.writerows(all_step_rows)

    # summary
    for avg in ('cell_exact', 'original'):
        for scope, pred in (('30 database structures', lambda r: 'slab' not in r['structure']),
                            ('slab sweep', lambda r: 'slab' in r['structure'])):
            rows = [r for r in all_step_rows if r['averaging'] == avg and pred(r)]
            c = {v: sum(r['verdict'] == v for r in rows) for v in ('PASS', 'PASS_NEAR_ZERO_ABS', 'FAIL')}
            print(f"[{avg:10s}] {scope:24s} {len(rows)} derivatives: PASS {c['PASS']}, "
                  f"PASS(near-zero abs) {c['PASS_NEAR_ZERO_ABS']}, FAIL {c['FAIL']}")
    print("\ncell-exact derivatives not passing the plain 2% relative test:")
    for r in all_step_rows:
        if r['averaging'] == 'cell_exact' and r['verdict'] != 'PASS':
            print(f"  {r['verdict']:20s} {r['structure']:22s} d{r['output']}/d{r['param']}  "
                  f"D(h/2,h,2h)=({r['D_half']:.3e}, {r['D_h']:.3e}, {r['D_2h']:.3e})  "
                  f"rel dev {r['max_rel_dev']:.1%}  |D|/largest {r['frac_of_largest']:.3f}  "
                  f"dev/abs-floor {r['dev_over_abs_floor']:.2f}")

    # original rerun vs the stored original database (should reproduce it)
    with open(os.path.join(REPO, 'data', 'sensitivity_db.json')) as f:
        stored = {r['id']: r for r in json.load(f)}
    with open(os.path.join(OUT, 'sensitivity_db_original_rerun.json')) as f:
        rerun = json.load(f)
    keys = ('neff', 'ng', 'dneff_dw', 'dneff_dt', 'dneff_da', 'dng_dw', 'dng_dt', 'dng_da')
    print("\noriginal-averaging rerun vs stored data/sensitivity_db.json (max relative difference):")
    for fam in ('strip_TE', 'strip_TM', 'slot_TE', 'slot_TM', 'rib_TE', 'swg_TE'):
        worst = max(abs(r[k] - stored[r['id']][k]) / abs(stored[r['id']][k])
                    for r in rerun if r['id'].startswith(fam) for k in keys)
        print(f"  {fam:9s} {worst:.2e}")
    print(f"\nTotal elapsed {time.time()-t0:.0f}s")
