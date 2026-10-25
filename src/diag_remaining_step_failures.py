"""
DIAGNOSTIC ONLY (not used by any result): why do some cell-exact derivatives
still fail the h/2, h, 2h check?

Hypothesis: cell-exact averaging conserves silicon area, but it still uses the
ARITHMETIC mean of eps in partially filled cells. That is the right average
for the field component parallel to an interface but not for the component
normal to it (which wants the harmonic mean). TM's dominant E is normal to the
strip top surface, and slot TE's dominant E is normal to the gap walls, which
is where the failures concentrate.

Test on strip TM 500 nm, dn_eff/dt and dn_eff/da (both fail on the production grid):
  (1) cell-exact on the production grid and on two finer core meshes: a pure
      discretization effect should shrink as cells shrink;
  (2) cell-exact fills with a crude anisotropic tensor on the production grid:
      eps_yy harmonic in cells cut only by a horizontal interface, eps_xx
      harmonic in cells cut only by a vertical one, arithmetic otherwise.
Run from src/:  python diag_remaining_step_failures.py
"""
import os
os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
import json
from multiprocessing import Pool

import numpy as np

import fullvector_sensitivity as fv
import cellexact_averaging as ce

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
H = dict(t=2.0, a=1.0)
MESHES = {'M0 30/18 (production)': (1.4, 1.2, 30, 18), 'M2 60/24': (1.4, 1.2, 60, 24), 'M3 90/36': (1.4, 1.2, 90, 36)}


def eps_tensor(X, Y, w, t, a):
    fill = ce.fill_strip(X, Y, w, t, a)
    e_si, e_ox = fv.N_SI**2, fv.N_OX**2
    arith = fill * e_si + (1 - fill) * e_ox
    harm = 1 / (fill / e_si + (1 - fill) / e_ox)
    xc = (X[:-1] + X[1:]) / 2
    yc = (Y[:-1] + Y[1:]) / 2
    xx, yy = np.meshgrid(xc, yc, indexing='ij')
    part = (fill > 1e-12) & (fill < 1 - 1e-12)
    # cell cut by the top/bottom face only: its centre is well inside the core laterally
    hw = w / 2 - np.tan(np.radians(a)) * np.clip(yy, 0, t)
    dxc = np.diff(X)[:, None] * np.ones_like(yy)
    horiz = part & (abs(xx) + dxc / 2 <= hw - 1e-12)
    vert = part & ~horiz & (yy - np.diff(Y)[None, :] / 2 >= 0) & (yy + np.diff(Y)[None, :] / 2 <= t)
    exx, eyy = arith.copy(), arith.copy()
    eyy[horiz] = harm[horiz]
    exx[vert] = harm[vert]
    out = np.zeros(fill.shape + (5,))
    out[..., 0], out[..., 3], out[..., 4] = exx, eyy, arith
    return lambda x_, y_: out


def solve(task):
    mode, mesh, p, d = task
    X, Y = fv.make_grid(*MESHES[mesh])
    w, t, a = 0.5, 0.22 + (d / 1000 if p == 't' else 0), (d if p == 'a' else 0.0)
    ep = ce.epsfunc_strip_ce(X, Y, w, t, a) if mode == 'cell_exact' else eps_tensor(X, Y, w, t, a)
    return task, fv.solve_mode_from_eps(ep, fv.WL0, 'TM', 4, X, Y).neff.real


if __name__ == '__main__':
    tasks = []
    for mode, meshes in (('cell_exact', list(MESHES)), ('tensor', ['M0 30/18 (production)'])):
        for mesh in meshes:
            for p in ('t', 'a'):
                for f in (0.5, 1, 2):
                    for s in (1, -1):
                        tasks.append((mode, mesh, p, s * f * H[p]))
    with Pool(12) as pool:
        res = dict(pool.map(solve, tasks))
    out = {}
    print("strip TM 500 nm: D(h/2), D(h), D(2h) and max relative deviation")
    for mode, mesh, p in sorted({(m, me, p) for m, me, p, _ in tasks}):
        D = [(res[(mode, mesh, p, f * H[p])] - res[(mode, mesh, p, -f * H[p])]) / (2 * f * H[p]) for f in (0.5, 1, 2)]
        dev = max(abs(D[0] - D[1]), abs(D[2] - D[1])) / abs(D[1])
        out[f'{mode}|{mesh}|dneff/d{p}'] = dict(D=D, max_rel_dev=dev)
        print(f"  {mode:10s} {mesh:24s} dneff/d{p}: {D[0]:.4e} {D[1]:.4e} {D[2]:.4e}  dev {dev:.1%}")
    with open(os.path.join(REPO, 'results_corrected', 'data', 'diag_remaining_step_failures.json'), 'w') as f:
        json.dump(out, f, indent=1)
