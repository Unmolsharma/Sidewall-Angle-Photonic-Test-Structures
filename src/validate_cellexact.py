"""
Step 2: re-validate the solver with cell-exact averaging.

Repeats every validation the writeup relies on, with the corrected epsfuncs
(and, side by side, the original ones):
  1. strip 500x220 TE0 vs published ~2.44-2.45   (pass: |n - 2.445| <= 0.01)
     + strip FD vs EIM at 400/500/600 nm (TE and TM)
  2. rib 500 nm / 90 nm slab vs our EIM (2.6033) and published EIM 2.602;
     FD-vs-EIM gap at 400-600 nm   (pass: FD below EIM at every width, gap
     within 0.5-5%, i.e. the same sign/order as strip's own FD-vs-EIM gap)
  3. rib with zero slab reproduces strip   (pass: |diff| < 1e-9)
  4. mesh/domain convergence, same configurations as the original checks
     (strip, strip TM, rib 500, rib 200-390, slot); pass thresholds are the
     ones the original checks accepted: strip/slot <= 1.5%, rib <= 1.0%.
     The original 4-config narrow-rib check in rib_validation.py used < 1.0%.

Run from the repo root:  python src/validate_cellexact.py
Writes results_corrected/data/validation.json
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
import rib_waveguide as rw
import slot_waveguide as sw
import cellexact_averaging as ce
from rib_validation import rib_eim_neff

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, 'results_corrected', 'data')

STRIP_MESHES = [('baseline 1.4x1.2 30/18', 1.4, 1.2, 30, 18), ('larger domain 2.0x1.8', 2.0, 1.8, 30, 18),
                ('finer core 1.4x1.2 45/18', 1.4, 1.2, 45, 18)]
RIB_MESHES = STRIP_MESHES
NARROW_RIB_MESHES = [('baseline 1.4x1.2 30/18', 1.4, 1.2, 30, 18), ('larger domain 2.0x1.8', 2.0, 1.8, 30, 18),
                     ('finer core 1.4x1.2 45/18', 1.4, 1.2, 45, 18), ('finer core 1.4x1.2 60/24', 1.4, 1.2, 60, 24)]
SLOT_MESHES = [('baseline 1.6x1.2 50/20', 1.6, 1.2, 50, 20), ('larger domain 2.2x1.8', 2.2, 1.8, 50, 20),
               ('finer core 1.6x1.2 70/20', 1.6, 1.2, 70, 20), ('larger+finer 2.2x1.8 70/28', 2.2, 1.8, 70, 28)]
# extra (information only, not a pass/fail gate): where does the corrected strip converge to?
STRIP_FINE = [('fine 1.4x1.2 60/24', 1.4, 1.2, 60, 24), ('finest 1.4x1.2 90/36', 1.4, 1.2, 90, 36)]


def job(task):
    avg, geom, w, pol, mesh, slab = task
    if mesh is None:
        X = Y = None
    elif geom == 'slot':
        X, Y = sw.make_slot_grid(*mesh[1:])
    else:
        X, Y = fv.make_grid(*mesh[1:])
    kw = {}
    if X is not None:
        kw.update(X=X, Y=Y)
    if slab is not None:
        kw.update(slab_h_um=slab)
    if avg == 'original' and geom == 'strip' and X is not None:
        # fv.solve_neff has no grid argument; same solve, explicit grid
        ep = fv.epsfunc_factory(w, 0.22, 0.0)
        return task, fv.solve_mode_from_eps(ep, fv.WL0, pol, 4, X, Y).neff.real
    f = ce.SOLVERS[avg][geom]
    return task, f(w, 0.22, 0.0, fv.WL0, pol, **kw)


if __name__ == '__main__':
    t0 = time.time()
    os.makedirs(OUT, exist_ok=True)
    tasks = []
    for avg in ('original', 'cell_exact'):
        for pol in ('TE', 'TM'):
            for w in (0.4, 0.5, 0.6):
                tasks.append((avg, 'strip', w, pol, None, None))
        for w in (0.4, 0.45, 0.5, 0.55, 0.6):
            tasks.append((avg, 'rib', w, 'TE', None, None))
        tasks.append((avg, 'rib', 0.5, 'TE', None, 0.0))
        for m in STRIP_MESHES + STRIP_FINE:
            tasks.append((avg, 'strip', 0.5, 'TE', m, None))
        for m in STRIP_MESHES:
            tasks.append((avg, 'strip', 0.5, 'TM', m, None))
        for m in RIB_MESHES:
            tasks.append((avg, 'rib', 0.5, 'TE', m, None))
        for w in (0.2, 0.25, 0.3, 0.35, 0.39):
            for m in NARROW_RIB_MESHES:
                tasks.append((avg, 'rib', w, 'TE', m, None))
        for m in SLOT_MESHES:
            tasks.append((avg, 'slot', sw.RAIL0, 'TE', m, None))
    with Pool(16) as p:
        res = dict(p.map(job, tasks))

    def get(avg, geom, w, pol='TE', mesh=None, slab=None):
        return res[(avg, geom, w, pol, mesh, slab)]

    out = {}
    verdicts = []
    for avg in ('original', 'cell_exact'):
        o = out[avg] = {}
        print(f"\n================ averaging = {avg} ================")
        n = get(avg, 'strip', 0.5)
        ok = abs(n - 2.445) <= 0.01
        o['strip_TE_500'] = n
        print(f"1. strip TE0 500x220 = {n:.4f}  (published ~2.44-2.45; pass |n-2.445|<=0.01: {ok})")
        o['strip_eim'] = {}
        for pol in ('TE', 'TM'):
            for w in (0.4, 0.5, 0.6):
                fd, eim = get(avg, 'strip', w, pol), fv.eim_neff(w, fv.H0, pol)
                o['strip_eim'][f'{pol}_{int(w*1000)}'] = dict(fd=fd, eim=eim, pct=100 * (eim - fd) / fd)
                print(f"   strip {pol} {w*1000:.0f}nm  FD={fd:.4f}  EIM={eim:.4f}  EIM-FD={100*(eim-fd)/fd:+.2f}%")
        if avg == 'cell_exact':
            verdicts.append(('strip TE0 500x220 within 0.01 of 2.445', ok, f'{n:.4f}'))

        rib_eim_500 = rib_eim_neff(0.5)
        print(f"2. rib EIM (ours) 500nm/90nm slab = {rib_eim_500:.4f}  vs published 2.602 "
              f"({100*(rib_eim_500-2.602)/2.602:+.2f}%)")
        o['rib_eim'] = {}
        gaps = []
        for w in (0.4, 0.45, 0.5, 0.55, 0.6):
            fd, eim = get(avg, 'rib', w), rib_eim_neff(w)
            gap = 100 * (eim - fd) / fd
            gaps.append(gap)
            o['rib_eim'][int(w * 1000)] = dict(fd=fd, eim=eim, pct=gap)
            print(f"   rib {w*1000:.0f}nm  FD={fd:.4f}  EIM={eim:.4f}  EIM-FD={gap:+.2f}%")
        ok = all(0.5 <= g <= 5.0 for g in gaps)
        if avg == 'cell_exact':
            verdicts.append(('rib FD below EIM by 0.5-5% at every width', ok,
                             f'{min(gaps):.2f}..{max(gaps):.2f}%'))

        r0, s0 = get(avg, 'rib', 0.5, slab=0.0), get(avg, 'strip', 0.5)
        ok = abs(r0 - s0) < 1e-9
        o['rib_zero_slab'] = dict(rib=r0, strip=s0, diff=r0 - s0)
        print(f"3. rib(slab=0) = {r0:.10f}  strip = {s0:.10f}  diff = {r0-s0:+.2e}  (pass <1e-9: {ok})")
        if avg == 'cell_exact':
            verdicts.append(('rib with zero slab reproduces strip (|diff|<1e-9)', ok, f'{r0-s0:+.1e}'))

        print("4. mesh/domain convergence (n_eff spread as % of mean)")
        o['convergence'] = {}

        def conv(name, vals, thr):
            spread = 100 * (max(vals) - min(vals)) / np.mean(vals)
            o['convergence'][name] = dict(values=vals, spread_pct=spread, threshold_pct=thr)
            print(f"   {name:28s} {' '.join(f'{v:.5f}' for v in vals)}  spread {spread:.2f}% (<= {thr}%: {spread <= thr})")
            if avg == 'cell_exact':
                verdicts.append((f'convergence {name} <= {thr}%', spread <= thr, f'{spread:.2f}%'))

        conv('strip TE 500', [get(avg, 'strip', 0.5, 'TE', m) for m in STRIP_MESHES], 1.5)
        conv('strip TM 500', [get(avg, 'strip', 0.5, 'TM', m) for m in STRIP_MESHES], 1.5)
        conv('rib TE 500', [get(avg, 'rib', 0.5, 'TE', m) for m in RIB_MESHES], 1.0)
        for w in (0.2, 0.25, 0.3, 0.35, 0.39):
            conv(f'rib TE {w*1000:.0f}', [get(avg, 'rib', w, 'TE', m) for m in NARROW_RIB_MESHES], 1.0)
        conv('slot TE 180', [get(avg, 'slot', sw.RAIL0, 'TE', m) for m in SLOT_MESHES], 1.5)
        fine = [get(avg, 'strip', 0.5, 'TE', m) for m in STRIP_MESHES[:1] + STRIP_FINE]
        o['strip_TE_500_refinement'] = dict(zip(['30/18', '60/24', '90/36'], fine))
        print(f"   (info) strip TE 500 core-mesh refinement 30/18, 60/24, 90/36: "
              f"{' '.join(f'{v:.5f}' for v in fine)}")

    print("\n================ VERDICT (cell-exact) ================")
    for name, ok, val in verdicts:
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {val}")
    out['verdicts'] = [dict(check=n, passed=bool(ok), value=v) for n, ok, v in verdicts]
    out['all_passed'] = all(ok for _, ok, _ in verdicts)
    print(f"  ALL PASSED: {out['all_passed']}")
    with open(os.path.join(OUT, 'validation.json'), 'w') as f:
        json.dump(out, f, indent=1)
    print(f"\nElapsed {time.time()-t0:.0f}s")
