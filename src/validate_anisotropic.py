"""
Phase 1 validation of the anisotropic (tensor) averaging, run side by side
with the original and cell-exact averagings (anisotropic_averaging.SOLVERS).

Checks (same tests and pass thresholds as validate_cellexact.py, plus n_g):
  1. strip 500x220 TE0 vs published ~2.44-2.45      (pass: |n - 2.445| <= 0.01)
     + strip FD vs EIM at 400/500/600 nm (TE and TM), information
  2. rib 500 nm / 90 nm slab vs our EIM (2.6033) and published EIM 2.602;
     FD-vs-EIM gap at 400-600 nm   (pass: FD below EIM by 0.5-5% at every width,
     the same sign/order as strip's own FD-vs-EIM gap)
  3. rib with zero slab reproduces strip             (pass: |diff| < 1e-9)
  4. mesh/domain convergence of n_eff AND n_g, same configurations as before;
     thresholds strip/slot <= 1.5%, rib <= 1.0% (the ones the original checks
     accepted), applied to n_g as well
  5. switching back reproduces the earlier numbers: every original-averaging
     and cell-exact value here vs results_corrected/data/validation.json, and
     original strip TE/TM 500 vs data/sensitivity_db.json (pass: |diff| < 1e-9)

Run from src/:  python validate_anisotropic.py
Writes results_anisotropic/data/validation.json
"""
import os
os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('MKL_NUM_THREADS', '1')
import json
import time
from multiprocessing import Pool

import numpy as np

import fullvector_sensitivity as fv
import slot_waveguide as sw
import anisotropic_averaging as an
from rib_validation import rib_eim_neff
from validate_cellexact import STRIP_MESHES, RIB_MESHES, NARROW_RIB_MESHES, SLOT_MESHES, STRIP_FINE

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, 'results_anisotropic', 'data')
AVGS = ('original', 'cell_exact', 'anisotropic')
DLAM = 0.005
WLS = (fv.WL0, fv.WL0 + DLAM, fv.WL0 - DLAM)


def job(task):
    avg, geom, w, pol, mesh, slab, wl = task
    kw = {}
    if mesh is not None:
        X, Y = (sw.make_slot_grid if geom == 'slot' else fv.make_grid)(*mesh[1:])
        kw.update(X=X, Y=Y)
    if slab is not None:
        kw.update(slab_h_um=slab)
    if avg == 'original' and geom == 'strip' and mesh is not None:
        # fv.solve_neff has no grid argument; same solve, explicit grid
        ep = fv.epsfunc_factory(w, 0.22, 0.0)
        return task, fv.solve_mode_from_eps(ep, wl, pol, 4, kw['X'], kw['Y']).neff.real
    return task, an.SOLVERS[avg][geom](w, 0.22, 0.0, wl, pol, **kw)


def units():
    u = []
    for pol in ('TE', 'TM'):
        for w in (0.4, 0.5, 0.6):
            u.append(('strip', w, pol, None, None))
    for w in (0.4, 0.45, 0.5, 0.55, 0.6):
        u.append(('rib', w, 'TE', None, None))
    u.append(('rib', 0.5, 'TE', None, 0.0))
    for m in STRIP_MESHES + STRIP_FINE:
        u.append(('strip', 0.5, 'TE', m, None))
    for m in STRIP_MESHES:
        u.append(('strip', 0.5, 'TM', m, None))
    for m in RIB_MESHES:
        u.append(('rib', 0.5, 'TE', m, None))
    for w in (0.2, 0.25, 0.3, 0.35, 0.39):
        for m in NARROW_RIB_MESHES:
            u.append(('rib', w, 'TE', m, None))
    for m in SLOT_MESHES:
        u.append(('slot', sw.RAIL0, 'TE', m, None))
    return list(dict.fromkeys(u))


if __name__ == '__main__':
    t0 = time.time()
    os.makedirs(OUT, exist_ok=True)
    tasks = [(avg,) + u + (wl,) for avg in AVGS for u in units() for wl in WLS]
    print(f"{len(tasks)} solves", flush=True)
    with Pool(16) as p:
        res = dict(p.map(job, tasks, chunksize=2))

    def get(avg, geom, w, pol='TE', mesh=None, slab=None):
        return res[(avg, geom, w, pol, mesh, slab, fv.WL0)]

    def ng(avg, geom, w, pol='TE', mesh=None, slab=None):
        n0, npl, nmi = (res[(avg, geom, w, pol, mesh, slab, wl)] for wl in WLS)
        return n0 - fv.WL0 * (npl - nmi) / (2 * DLAM)

    out = {}
    verdicts = []
    for avg in AVGS:
        o = out[avg] = {}
        gate = avg == 'anisotropic'
        print(f"\n================ averaging = {avg} ================")
        n = get(avg, 'strip', 0.5)
        ok = abs(n - 2.445) <= 0.01
        o['strip_TE_500'] = n
        o['strip_TE_500_ng'] = ng(avg, 'strip', 0.5)
        o['strip_TM_500'] = get(avg, 'strip', 0.5, 'TM')
        print(f"1. strip TE0 500x220 = {n:.5f}, n_g = {o['strip_TE_500_ng']:.4f}  "
              f"(published ~2.44-2.45; pass |n-2.445|<=0.01: {ok});  TM0 = {o['strip_TM_500']:.5f}")
        o['strip_eim'] = {}
        for pol in ('TE', 'TM'):
            for w in (0.4, 0.5, 0.6):
                fd, eim = get(avg, 'strip', w, pol), fv.eim_neff(w, fv.H0, pol)
                o['strip_eim'][f'{pol}_{int(w*1000)}'] = dict(fd=fd, eim=eim, pct=100 * (eim - fd) / fd)
                print(f"   strip {pol} {w*1000:.0f}nm  FD={fd:.4f}  EIM={eim:.4f}  EIM-FD={100*(eim-fd)/fd:+.2f}%")
        if gate:
            verdicts.append(('strip TE0 500x220 within 0.01 of 2.445', ok, f'{n:.4f}'))

        rib_eim_500 = rib_eim_neff(0.5)
        print(f"2. rib EIM (ours) 500nm/90nm slab = {rib_eim_500:.4f}  vs published 2.602 "
              f"({100*(rib_eim_500-2.602)/2.602:+.2f}%)")
        o['rib_eim_500'] = rib_eim_500
        o['rib_eim'] = {}
        gaps = []
        for w in (0.4, 0.45, 0.5, 0.55, 0.6):
            fd, eim = get(avg, 'rib', w), rib_eim_neff(w)
            gap = 100 * (eim - fd) / fd
            gaps.append(gap)
            o['rib_eim'][int(w * 1000)] = dict(fd=fd, eim=eim, pct=gap)
            print(f"   rib {w*1000:.0f}nm  FD={fd:.4f}  EIM={eim:.4f}  EIM-FD={gap:+.2f}%")
        ok = all(0.5 <= g <= 5.0 for g in gaps)
        if gate:
            verdicts.append(('rib FD below EIM by 0.5-5% at every width', ok, f'{min(gaps):.2f}..{max(gaps):.2f}%'))

        r0, s0 = get(avg, 'rib', 0.5, slab=0.0), get(avg, 'strip', 0.5)
        ok = abs(r0 - s0) < 1e-9
        o['rib_zero_slab'] = dict(rib=r0, strip=s0, diff=r0 - s0)
        print(f"3. rib(slab=0) = {r0:.10f}  strip = {s0:.10f}  diff = {r0-s0:+.2e}  (pass <1e-9: {ok})")
        if gate:
            verdicts.append(('rib with zero slab reproduces strip (|diff|<1e-9)', ok, f'{r0-s0:+.1e}'))

        print("4. mesh/domain convergence (spread as % of mean), n_eff and n_g")
        o['convergence'] = {}

        def conv(name, geom, w, pol, meshes, thr):
            for q, fn in (('neff', get), ('ng', ng)):
                vals = [fn(avg, geom, w, pol, m) for m in meshes]
                spread = 100 * (max(vals) - min(vals)) / np.mean(vals)
                o['convergence'][f'{name} {q}'] = dict(values=vals, spread_pct=spread, threshold_pct=thr)
                print(f"   {name:14s} {q:4s} {' '.join(f'{v:.5f}' for v in vals)}  spread {spread:.2f}% "
                      f"(<= {thr}%: {spread <= thr})")
                if gate:
                    verdicts.append((f'convergence {name} {q} <= {thr}%', spread <= thr, f'{spread:.2f}%'))

        conv('strip TE 500', 'strip', 0.5, 'TE', STRIP_MESHES, 1.5)
        conv('strip TM 500', 'strip', 0.5, 'TM', STRIP_MESHES, 1.5)
        conv('rib TE 500', 'rib', 0.5, 'TE', RIB_MESHES, 1.0)
        for w in (0.2, 0.25, 0.3, 0.35, 0.39):
            conv(f'rib TE {w*1000:.0f}', 'rib', w, 'TE', NARROW_RIB_MESHES, 1.0)
        conv('slot TE 180', 'slot', sw.RAIL0, 'TE', SLOT_MESHES, 1.5)
        ms = STRIP_MESHES[:1] + STRIP_FINE
        o['strip_TE_500_refinement'] = dict(neff=dict(zip(['30/18', '60/24', '90/36'], [get(avg, 'strip', 0.5, 'TE', m) for m in ms])),
                                            ng=dict(zip(['30/18', '60/24', '90/36'], [ng(avg, 'strip', 0.5, 'TE', m) for m in ms])))
        print(f"   (info) strip TE 500 core-mesh refinement 30/18, 60/24, 90/36: n_eff "
              f"{' '.join(f'{v:.5f}' for v in o['strip_TE_500_refinement']['neff'].values())}  n_g "
              f"{' '.join(f'{v:.4f}' for v in o['strip_TE_500_refinement']['ng'].values())}")

    # 5. switching back reproduces the earlier numbers
    print("\n5. reproduction of earlier runs")
    with open(os.path.join(REPO, 'results_corrected', 'data', 'validation.json')) as f:
        prev = json.load(f)
    worst = {}
    for avg in ('original', 'cell_exact'):
        p, o = prev[avg], out[avg]
        diffs = [abs(p['strip_TE_500'] - o['strip_TE_500']), abs(p['rib_zero_slab']['rib'] - o['rib_zero_slab']['rib'])]
        diffs += [abs(p['strip_eim'][k]['fd'] - o['strip_eim'][k]['fd']) for k in p['strip_eim']]
        diffs += [abs(p['rib_eim'][k]['fd'] - o['rib_eim'][int(k)]['fd']) for k in p['rib_eim']]
        for name, c in p['convergence'].items():
            diffs += [abs(a - b) for a, b in zip(c['values'], o['convergence'][f'{name} neff']['values'])]
        worst[avg] = max(diffs)
        ok = worst[avg] < 1e-9
        print(f"   {avg:10s} vs results_corrected/data/validation.json: {len(diffs)} values, max |diff| = {worst[avg]:.1e}  (pass <1e-9: {ok})")
        verdicts.append((f'{avg} averaging reproduces the previous run ({len(diffs)} n_eff values, |diff|<1e-9)', ok, f'{worst[avg]:.1e}'))
    with open(os.path.join(REPO, 'data', 'sensitivity_db.json')) as f:
        db = {r['id']: r for r in json.load(f)}
    for pol, key in (('TE', 'strip_TE_500'), ('TM', 'strip_TM_500')):
        d = abs(db[f'strip_{pol}_500']['neff'] - out['original'][key])
        print(f"   original strip {pol} 500 = {out['original'][key]:.10f} vs data/sensitivity_db.json "
              f"{db[f'strip_{pol}_500']['neff']:.10f}  |diff| = {d:.1e}")
        verdicts.append((f'original averaging reproduces stored strip {pol} 500 n_eff (|diff|<1e-9)', d < 1e-9, f'{d:.1e}'))
    out['reproduction_max_abs_diff'] = worst

    print("\n================ VERDICT ================")
    for name, ok, val in verdicts:
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {val}")
    out['verdicts'] = [dict(check=n, passed=bool(ok), value=v) for n, ok, v in verdicts]
    out['all_passed'] = all(ok for _, ok, _ in verdicts)
    print(f"  ALL PASSED: {out['all_passed']}")
    with open(os.path.join(OUT, 'validation.json'), 'w') as f:
        json.dump(out, f, indent=1)
    print(f"\nElapsed {time.time()-t0:.0f}s")
