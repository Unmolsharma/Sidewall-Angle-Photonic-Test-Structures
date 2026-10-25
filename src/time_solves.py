"""
Phase 1.3: time one anisotropic mode solve per structure type, serially and
in parallel, and estimate the Phase 2 runtime.

  serial    : one process, BLAS threads at the library default, and again with 1 thread
  parallel  : multiprocessing.Pool of P workers (1 BLAS thread each) solving a
              batch of N_BATCH solves of one type; wall seconds per solve =
              batch wall time / N_BATCH, per-solve latency = that x P

Each configuration runs in a fresh subprocess (`python time_solves.py worker ...`)
so the BLAS thread setting takes effect before numpy is imported.

Run from src/:  python time_solves.py
Writes results_anisotropic/data/timing.json
"""
import json
import os
import subprocess
import sys
import time

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TYPES = [('strip', 'TE', 0.5), ('strip', 'TM', 0.5), ('slot', 'TE', 0.18), ('slot', 'TM', 0.18),
         ('rib', 'TE', 0.5), ('swg', 'TE', 0.5)]
N_SERIAL = 3
N_BATCH = 32
POOLS = (4, 8, 16)


def one(t):
    import anisotropic_averaging as an
    import fullvector_sensitivity as fv
    g, pol, w, k = t
    # tiny width offset per task so no two solves are identical
    return an.SOLVERS['anisotropic'][g](w + 1e-4 * k, 0.22, 0.0, fv.WL0, pol)


def worker(threads, mode, n, procs=1):
    out = {}
    if mode == 'serial':
        for g, pol, w in TYPES:
            one((g, pol, w, 0))  # warm-up (imports, caches)
            ts = []
            for k in range(n):
                t0 = time.perf_counter()
                one((g, pol, w, k))
                ts.append(time.perf_counter() - t0)
            out[f'{g}_{pol}'] = ts
    else:
        from multiprocessing import Pool
        with Pool(procs) as pool:
            pool.map(one, [(g, pol, w, 0) for g, pol, w in TYPES] * 2)  # warm-up
            for g, pol, w in TYPES:
                t0 = time.perf_counter()
                pool.map(one, [(g, pol, w, k) for k in range(n)], chunksize=1)
                out[f'{g}_{pol}'] = time.perf_counter() - t0
    print(json.dumps(out))


def run(threads, mode, n, procs=1):
    env = dict(os.environ)
    for v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
        if threads == 'default':
            env.pop(v, None)
        else:
            env[v] = threads
    r = subprocess.run([sys.executable, os.path.abspath(__file__), 'worker', threads, mode, str(n), str(procs)],
                       capture_output=True, text=True, env=env, cwd=os.path.dirname(os.path.abspath(__file__)),
                       timeout=3600)
    if r.returncode:
        raise RuntimeError(r.stderr)
    return json.loads(r.stdout.strip().splitlines()[-1])


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == 'worker':
        worker(sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5]))
        sys.exit()
    res = {}
    for threads in ('default', '1'):
        s = run(threads, 'serial', N_SERIAL)
        res[f'serial_threads_{threads}'] = {k: sum(v) / len(v) for k, v in s.items()}
        print(f"serial, BLAS threads={threads}: " + ', '.join(f"{k} {sum(v)/len(v):.2f}s" for k, v in s.items()), flush=True)
    for P in POOLS:
        p = run('1', 'pool', N_BATCH, P)
        res[f'pool_{P}'] = {k: v / N_BATCH for k, v in p.items()}
        res[f'pool_{P}_latency'] = {k: v * P / N_BATCH for k, v in p.items()}
        print(f"Pool({P}), {N_BATCH} solves each: " +
              ', '.join(f"{k} {v/N_BATCH:.3f}s/solve (latency {v*P/N_BATCH:.2f}s)" for k, v in p.items()), flush=True)

    # Phase 2 core solve counts (anisotropic only; original and cell-exact already exist)
    per_struct = 19 * 3   # (nominal + 3 params x 3 steps x 2 signs) x 3 wavelengths
    counts = dict(strip_TE=5 * per_struct, strip_TM=5 * per_struct, slot_TE=5 * per_struct, slot_TM=5 * per_struct,
                  rib_TE=(5 + 24) * per_struct,   # database + slab sweep (6 slabs x 4 widths)
                  swg_TE=5 * per_struct)
    res['phase2_core_solve_counts'] = counts
    modes = [k for k in res if k.startswith('serial') or (k.startswith('pool') and not k.endswith('latency'))]
    est = {m: sum(counts[k] * res[m][k] for k in counts) / 60 for m in modes}
    res['phase2_core_estimate_minutes'] = est
    print(f"\nPhase 2 core recompute ({sum(counts.values())} solves: 30 structures + slab sweep, all steps, n_eff and n_g):")
    for k, v in est.items():
        print(f"  {k:24s} {v:.1f} min")
    res['fastest'] = min(est, key=est.get)
    print(f"fastest: {res['fastest']}")
    os.makedirs(os.path.join(REPO, 'results_anisotropic', 'data'), exist_ok=True)
    with open(os.path.join(REPO, 'results_anisotropic', 'data', 'timing.json'), 'w') as f:
        json.dump(res, f, indent=1)
