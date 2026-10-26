"""
femwell production batch (femwell is the main derivative engine).

Runs in the femwell virtual environment:
    C:/Users/unmol/.venvs/femwell/Scripts/python.exe src/femwell_batch.py meshcheck r8
    C:/Users/unmol/.venvs/femwell/Scripts/python.exe src/femwell_batch.py time r12
    C:/Users/unmol/.venvs/femwell/Scripts/python.exe src/femwell_batch.py full r12

Geometry conventions (identical to the EMpy code; constants copied from
fullvector_sensitivity.py, rib_waveguide.py, slot_waveguide.py, swg_waveguide.py):
  - lambda = 1.55 um, n_Si = 3.48, n_ox = 1.44 (no material dispersion), air = 1.0
  - nominal thickness 220 nm; width = BOTTOM width of the etched wall
  - sidewall angle a = deviation from vertical in degrees, top narrower, pivot
    at the bottom of the tilted wall (y = 0 for strip/slot/SWG, y = slab top
    for rib); wall moves inward by tan(a)(y - y_pivot)
  - rib: slab of height SLAB (default 90 nm) across the whole domain, rib
    from slab top to t; thickness t perturbed with the slab held fixed
  - slot: rails of width w, gap = max(GAP0 - (w - RAIL0), GAP_MIN), rail
    centres at +-(gap/2 + w/2); each rail tilts about its own bottom
  - SWG: zeroth-order effective medium, strip with core index
    sqrt(f n_Si^2 + (1 - f) n_ox^2), f = 0.5
  - oxide below y = 0; top cladding oxide (or air for the cladding check)
  - domain x in [-1.4, 1.4] (slot [-1.6, 1.6]), y in [-1.2, 1.2] um, as in EMpy
Steps h = 5 nm (w), 2 nm (t), 1 deg (a), central differences at h/2, h, 2h.
n_g = n_eff - lambda dn_eff/dlambda, dlambda = 5 nm, same mesh for the three
wavelengths (each perturbed geometry is remeshed once).

Mode selection: 2 modes near a guess (the EMpy cell-exact n_eff of that
structure), keep the highest-n_eff mode of the requested polarization
(te_fraction > 0.5 = TE); if none or it is > 3% from the guess, re-solve with
6 modes and no guess. Every solve is appended to a JSON-lines cache as soon
as it finishes, so a restart resumes.

Safety: before each new solve is submitted the AC-power state is checked; on
battery the runner waits and re-checks every 5 minutes (logged). The worker
count is chosen from free RAM.
"""
import ctypes
import json
import os
import sys
import time
from collections import OrderedDict
from multiprocessing import Pool

import numpy as np
from shapely.geometry import Polygon, box
from skfem import Basis, ElementTriP0
from skfem.io.meshio import from_meshio

from femwell.maxwell.waveguide import compute_modes
from femwell.mesh import mesh_from_OrderedDict

N_SI, N_OX, N_AIR, WL0 = 3.48, 1.44, 1.0, 1.55
SLAB0, RAIL0, GAP0, GAP_MIN, DUTY = 0.09, 0.180, 0.120, 0.02, 0.5
N_SWG = float(np.sqrt(DUTY * N_SI ** 2 + (1 - DUTY) * N_OX ** 2))
H = dict(w=5.0, t=2.0, a=1.0)
FACTORS = (0.5, 1.0, 2.0)
DLAM = 0.005
WLS = (WL0, WL0 + DLAM, WL0 - DLAM)
MESH = {'r12': (0.012, 0.3, 0.2), 'r8': (0.008, 0.3, 0.2)}
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, 'results_femwell')

DB_STRUCTS = ([('strip', 'TE', w) for w in (400, 450, 500, 550, 600)] +
              [('strip', 'TM', w) for w in (400, 450, 500, 550, 600)] +
              [('slot', 'TE', w) for w in (150, 165, 180, 195, 210)] +
              [('slot', 'TM', w) for w in (150, 165, 180, 195, 210)] +
              [('rib', 'TE', w) for w in (400, 450, 500, 550, 600)] +
              [('swg', 'TE', w) for w in (400, 450, 500, 550, 600)])
SLAB_SWEEP = [(s, w) for s in (50, 70, 90, 110, 130, 150) for w in (300, 400, 500, 600)]


# ---------------------------------------------------------------------------
# geometry and solve
# ---------------------------------------------------------------------------
def trapezoid(cx, w, y0, y1, a):
    hw0 = w / 2
    hwt = hw0 - np.tan(np.radians(a)) * (y1 - y0)
    return Polygon([(cx - hw0, y0), (cx + hw0, y0), (cx + hwt, y1), (cx - hwt, y1)])


def shapes(geom, w, t, a, slab):
    if geom in ('strip', 'swg'):
        return [trapezoid(0.0, w, 0.0, t, a)]
    if geom == 'slot':
        gap = max(GAP0 - (w - RAIL0), GAP_MIN)
        return [trapezoid(s * (gap / 2 + w / 2), w, 0.0, t, a) for s in (-1, 1)]
    if geom == 'rib':
        hw0, hwt = w / 2, w / 2 - np.tan(np.radians(a)) * (t - slab)
        return [Polygon([(-50, 0), (50, 0), (50, slab), (hw0, slab), (hwt, t), (-hwt, t), (-hw0, slab), (-50, slab)])]
    raise ValueError(geom)


def build(geom, w, t, a, slab, clad, mesh):
    hx, hy = (1.6, 1.2) if geom == 'slot' else (1.4, 1.2)
    dom = box(-hx, -hy, hx, hy)
    cores = [s.intersection(dom) for s in shapes(geom, w, t, a, slab)]
    union = cores[0]
    for c in cores[1:]:
        union = union.union(c)
    polys = OrderedDict((f'core{i}', c) for i, c in enumerate(cores))
    polys['box'] = box(-hx, -hy, hx, 0).difference(union)
    polys['clad'] = box(-hx, 0, hx, hy).difference(union)
    res, dist, rmax = MESH[mesh]
    m = from_meshio(mesh_from_OrderedDict(
        polys, {f'core{i}': dict(resolution=res, distance=dist) for i in range(len(cores))},
        default_resolution_max=rmax))
    basis0 = Basis(m, ElementTriP0())
    eps = basis0.zeros()
    n_core = N_SWG if geom == 'swg' else N_SI
    for i in range(len(cores)):
        eps[basis0.get_dofs(elements=f'core{i}')] = n_core ** 2
    eps[basis0.get_dofs(elements='box')] = N_OX ** 2
    eps[basis0.get_dofs(elements='clad')] = (N_OX if clad == 'oxide' else N_AIR) ** 2
    return basis0, eps


def pick(modes, pol):
    ok = [m for m in modes if (m.te_fraction > 0.5) == (pol == 'TE')]
    return max(ok, key=lambda m: m.n_eff.real) if ok else None


def solve(task):
    """task = (geom, pol, W_nm, slab_nm, clad, mesh, param, delta, guess) -> dict"""
    geom, pol, W, slab_nm, clad, mesh, p, d, guess = task
    t0 = time.time()
    w, t, a = W / 1000, 0.22, 0.0
    if p == 'w':
        w += d / 1000
    elif p == 't':
        t += d / 1000
    elif p == 'a':
        a += d
    slab = (slab_nm if slab_nm is not None else SLAB0 * 1000) / 1000
    basis0, eps = build(geom, w, t, a, slab, clad, mesh)
    ns, tes, notes = [], [], []
    for wl in WLS:
        m = pick(compute_modes(basis0, eps, wavelength=wl, num_modes=2, order=2, n_guess=guess), pol)
        if m is None or abs(m.n_eff.real - guess) > 0.03 * guess:
            notes.append(f'wl {wl}: fallback to 6 modes without guess')
            m = pick(compute_modes(basis0, eps, wavelength=wl, num_modes=6, order=2), pol)
        ns.append(float(m.n_eff.real))
        tes.append(float(m.te_fraction))
    return dict(key=key_of(task), n=ns, te_fraction=tes, n_elements=int(basis0.mesh.nelements),
                seconds=time.time() - t0, notes=notes)


def key_of(task):
    geom, pol, W, slab_nm, clad, mesh, p, d, guess = task
    return f'{geom}|{pol}|{W}|{slab_nm}|{clad}|{mesh}|{p}|{d:+.4f}'


def unit_tasks(geom, pol, W, slab_nm, clad, mesh, guess, params='wta'):
    out = [(geom, pol, W, slab_nm, clad, mesh, 'nom', 0.0, guess)]
    for p in params:
        for f in FACTORS:
            for s in (1, -1):
                out.append((geom, pol, W, slab_nm, clad, mesh, p, s * f * H[p], guess))
    return out


# ---------------------------------------------------------------------------
# machine state
# ---------------------------------------------------------------------------
class _PowerStatus(ctypes.Structure):
    _fields_ = [('ACLineStatus', ctypes.c_byte), ('BatteryFlag', ctypes.c_byte),
                ('BatteryLifePercent', ctypes.c_byte), ('SystemStatusFlag', ctypes.c_byte),
                ('BatteryLifeTime', ctypes.c_ulong), ('BatteryFullLifeTime', ctypes.c_ulong)]


def on_ac():
    s = _PowerStatus()
    ctypes.windll.kernel32.GetSystemPowerStatus(ctypes.byref(s))
    return s.ACLineStatus == 1, s.BatteryLifePercent


BATTERY_FLOOR = None  # None: pause whenever on battery; a percentage: run on battery down to it


def power_ok():
    ac, pct = on_ac()
    return (ac or (BATTERY_FLOOR is not None and pct >= BATTERY_FLOOR)), pct


class _MemStatus(ctypes.Structure):
    _fields_ = [('dwLength', ctypes.c_ulong), ('dwMemoryLoad', ctypes.c_ulong),
                ('ullTotalPhys', ctypes.c_ulonglong), ('ullAvailPhys', ctypes.c_ulonglong),
                ('ullTotalPageFile', ctypes.c_ulonglong), ('ullAvailPageFile', ctypes.c_ulonglong),
                ('ullTotalVirtual', ctypes.c_ulonglong), ('ullAvailVirtual', ctypes.c_ulonglong),
                ('ullAvailExtendedVirtual', ctypes.c_ulonglong)]


def free_gb():
    s = _MemStatus()
    s.dwLength = ctypes.sizeof(_MemStatus)
    ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(s))
    return s.ullAvailPhys / 2 ** 30


def log(msg, path):
    line = f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}"
    print(line, flush=True)
    with open(path, 'a', encoding='utf-8') as f:
        f.write(line + '\n')


# ---------------------------------------------------------------------------
# resumable runner
# ---------------------------------------------------------------------------
def load_cache(path):
    done = {}
    if os.path.exists(path):
        with open(path, encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if line:
                    r = json.loads(line)
                    done[r['key']] = r
    return done


def run(tasks, cache, logpath, workers, gb_per_worker):
    done = load_cache(cache)
    todo = [t for t in dict.fromkeys(tasks) if key_of(t) not in done]
    log(f'{len(tasks)} tasks, {len(tasks) - len(todo)} cached, {len(todo)} to run, {workers} workers', logpath)
    failed = {}
    t0 = time.time()
    n_done = 0
    tight = False
    with Pool(workers, maxtasksperchild=20) as pool:
        inflight = {}
        while todo or inflight:
            while todo and len(inflight) < (1 if tight else workers):
                ac, pct = power_ok()
                waited = False
                while not ac:
                    if not waited:
                        log(f'on battery ({pct}%{f", floor {BATTERY_FLOOR}%" if BATTERY_FLOOR is not None else ""}): '
                            f'pausing new solves, re-checking every 5 min '
                            f'({len(inflight)} in flight finish meanwhile)', logpath)
                        waited = True
                    # let in-flight solves finish and be recorded while waiting
                    for k in [k for k, r in inflight.items() if r.ready()]:
                        n_done += _collect(k, inflight.pop(k), cache, failed, logpath)
                    time.sleep(300)
                    ac, pct = power_ok()
                if waited:
                    log(f'power OK again ({pct}%), resuming', logpath)
                fg = free_gb()
                if (fg < gb_per_worker + 2) != tight:
                    tight = not tight
                    log(f'free RAM {fg:.1f} GB: ' + ('tight, dropping to 1 worker' if tight else
                                                     f'recovered, back to {workers} workers'), logpath)
                if tight and inflight:
                    break  # wait for a worker to finish before adding load
                t = todo.pop(0)
                inflight[key_of(t)] = pool.apply_async(solve, (t,))
            for k in [k for k, r in inflight.items() if r.ready()]:
                n_done += _collect(k, inflight.pop(k), cache, failed, logpath)
                if n_done % 25 == 0:
                    el = time.time() - t0
                    log(f'{n_done} solved, {len(todo)} queued, {len(inflight)} in flight, '
                        f'{el/60:.1f} min elapsed, ~{el/n_done*len(todo)/60:.0f} min left, free RAM {free_gb():.1f} GB',
                        logpath)
            time.sleep(1)
    log(f'finished: {n_done} solved this run, {len(failed)} failed', logpath)
    return failed


def _collect(k, r, cache, failed, logpath):
    try:
        res = r.get()
    except Exception as e:  # keep going; the structure is reported as skipped
        failed[k] = repr(e)
        log(f'FAILED {k}: {e!r}', logpath)
        return 0
    with open(cache, 'a', encoding='utf-8') as f:
        f.write(json.dumps(res) + '\n')
    if res['notes']:
        log(f'note {k}: {res["notes"]}', logpath)
    return 1


# ---------------------------------------------------------------------------
# task sets
# ---------------------------------------------------------------------------
def guesses():
    """EMpy cell-exact n_eff per structure, used only as the mode-search guess."""
    d = os.path.join(REPO, 'results_corrected', 'data')
    g = {(r['geometry'], r['polarization'], r['width_nm'], None): r['neff']
         for r in json.load(open(os.path.join(d, 'sensitivity_db.json')))}
    for e in json.load(open(os.path.join(d, 'rib_slab_sweep_db_cell_exact.json'))):
        g[('rib', 'TE', e['width_nm'], e['slab_nm'])] = e['neff']
    return g


def meshcheck_tasks(mesh):
    g = guesses()
    return (unit_tasks('strip', 'TE', 500, None, 'oxide', mesh, g[('strip', 'TE', 500, None)]) +
            unit_tasks('rib', 'TE', 500, None, 'oxide', mesh, g[('rib', 'TE', 500, None)]) +
            unit_tasks('strip', 'TE', 500, None, 'air', mesh, 2.39, params='wa'))


def full_tasks(mesh):
    g = guesses()
    tasks = []
    for geom, pol, W in DB_STRUCTS:
        tasks += unit_tasks(geom, pol, W, None, 'oxide', mesh, g[(geom, pol, W, None)])
    for s, W in SLAB_SWEEP:
        tasks += unit_tasks('rib', 'TE', W, s, 'oxide', mesh, g[('rib', 'TE', W, s)])
    return tasks


def choose_workers(gb_per_worker, max_workers=4, logpath=None):
    fg = free_gb()
    n = int(max(1, min(max_workers, (fg - 4) // gb_per_worker)))
    if logpath:
        log(f'free RAM {fg:.1f} GB, ~{gb_per_worker} GB per worker -> {n} workers', logpath)
    return n


if __name__ == '__main__':
    os.makedirs(os.path.join(OUT, 'cache'), exist_ok=True)
    os.makedirs(os.path.join(OUT, 'logs'), exist_ok=True)
    what, mesh = sys.argv[1], sys.argv[2]
    gbw = 3.5 if mesh == 'r8' else 2.0
    if what == 'meshcheck':
        logpath = os.path.join(OUT, 'logs', f'meshcheck_{mesh}.log')
        run(meshcheck_tasks(mesh), os.path.join(OUT, 'cache', f'solves_{mesh}.jsonl'), logpath,
            choose_workers(gbw, max_workers=3, logpath=logpath), gbw)
    elif what == 'time':
        logpath = os.path.join(OUT, 'logs', f'timing_{mesh}.log')
        g = guesses()
        sample = [('strip', 'TE', 500, None), ('strip', 'TM', 500, None), ('slot', 'TE', 180, None),
                  ('slot', 'TM', 180, None), ('rib', 'TE', 500, None), ('swg', 'TE', 500, None)]
        tasks = [(geom, pol, W, s, 'oxide', mesh, 'nom', 0.0, g[(geom, pol, W, s)]) for geom, pol, W, s in sample]
        tasks += [(geom, pol, W, s, 'oxide', mesh, 'w', 5.0, g[(geom, pol, W, s)]) for geom, pol, W, s in sample]
        run(tasks, os.path.join(OUT, 'cache', f'solves_{mesh}.jsonl'), logpath, choose_workers(gbw, logpath=logpath), gbw)
        rows = load_cache(os.path.join(OUT, 'cache', f'solves_{mesh}.jsonl'))
        for geom, pol, W, s in sample:
            rs = [rows[key_of(t)] for t in tasks if t[:4] == (geom, pol, W, s) and key_of(t) in rows]
            log(f'timing {geom} {pol} {W}: {np.mean([r["seconds"] for r in rs]):.0f} s per geometry '
                f'(3 wavelengths, {rs[0]["n_elements"]} elements), n_eff {rs[0]["n"][0]:.5f}, '
                f'te_fraction {rs[0]["te_fraction"][0]:.3f}', logpath)
    elif what == 'full':
        # optional third argument caps the worker count (e.g. 2 to keep memory low when resuming)
        cap = int(sys.argv[3]) if len(sys.argv) > 3 else 4
        # optional fourth argument: run on battery while charge >= this percent
        if len(sys.argv) > 4:
            BATTERY_FLOOR = int(sys.argv[4])
        logpath = os.path.join(OUT, 'logs', f'full_{mesh}.log')
        failed = run(full_tasks(mesh), os.path.join(OUT, 'cache', f'solves_{mesh}.jsonl'), logpath,
                     choose_workers(gbw, max_workers=cap, logpath=logpath), gbw)
        with open(os.path.join(OUT, 'logs', f'full_{mesh}_failed.json'), 'w') as f:
            json.dump(failed, f, indent=1)
