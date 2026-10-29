"""
femwell solves at the Monte Carlo surrogate design points (results_montecarlo/data/mc_design.json, written
by montecarlo.py plan).

Runs in the femwell virtual environment, from src/:
    C:/Users/unmol/.venvs/femwell/Scripts/python.exe montecarlo_batch.py [max_workers]

Each task is one geometry at absolute values: width (bottom width of the wall, nm), thickness (nm), wall
angle (deg; 0 for strips, whose angle enters only through the mid-height width under the symmetry), slab
(nm, ribs), n_Si, n_ox. Same mesh (r12), wavelengths and mode selection as femwell_batch.py. Every solve is
appended to results_montecarlo/cache/solves_r12.jsonl; a restart resumes. The runner (power and RAM checks)
is femwell_batch.run.
"""
import json
import os
import sys
import time

import femwell_batch as fb
from femwell.maxwell.waveguide import compute_modes

OUT = os.path.join(fb.REPO, 'results_montecarlo')
CACHE = os.path.join(OUT, 'cache', 'solves_r12.jsonl')


def key_of(task):
    geom, pol, W, S, w, t, a, slab, nsi, nox, guess = task
    # montecarlo.mc_key builds the same key from the design file
    return f'mc|{geom}|{pol}|{W}|{S}|{w:.4f}|{t:.4f}|{a:.5f}|{slab}|{nsi:.7f}|{nox:.7f}'


def solve(task):
    geom, pol, W, S, w, t, a, slab, nsi, nox, guess = task
    t0 = time.time()
    basis0, eps = fb.build(geom, w / 1000, t / 1000, a, slab / 1000 if slab is not None else fb.SLAB0, 'oxide', 'r12', nsi, nox)
    ns, tes, notes = [], [], []
    for wl in fb.WLS:
        m = fb.pick(compute_modes(basis0, eps, wavelength=wl, num_modes=2, order=2, n_guess=guess), pol)
        if m is None or abs(m.n_eff.real - guess) > 0.03 * guess:
            notes.append(f'wl {wl}: fallback to 6 modes without guess')
            m = fb.pick(compute_modes(basis0, eps, wavelength=wl, num_modes=6, order=2), pol)
        ns.append(float(m.n_eff.real))
        tes.append(float(m.te_fraction))
    return dict(key=key_of(task), n=ns, te_fraction=tes, n_elements=int(basis0.mesh.nelements),
                seconds=time.time() - t0, notes=notes)


def tasks():
    design = json.load(open(os.path.join(OUT, 'data', 'mc_design.json')))
    out = []
    for d in design['designs']:
        for p in d['points']:
            out.append((d['geom'], d['pol'], d['W'], d['S'], p['w'], p['t'], p['a'], p['slab'], p['nsi'], p['nox'], d['guess']))
    return out


if __name__ == '__main__':
    os.makedirs(os.path.join(OUT, 'cache'), exist_ok=True)
    os.makedirs(os.path.join(OUT, 'logs'), exist_ok=True)
    cap = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    if len(sys.argv) > 2:   # optional: run on battery while charge >= this percent
        fb.BATTERY_FLOOR = int(sys.argv[2])
    logpath = os.path.join(OUT, 'logs', 'mc_solves_r12.log')
    failed = fb.run(tasks(), CACHE, logpath, fb.choose_workers(2.0, max_workers=cap, logpath=logpath), 2.0,
                    solver=solve, keyer=key_of)
    with open(os.path.join(OUT, 'logs', 'mc_solves_r12_failed.json'), 'w') as f:
        json.dump(failed, f, indent=1)
