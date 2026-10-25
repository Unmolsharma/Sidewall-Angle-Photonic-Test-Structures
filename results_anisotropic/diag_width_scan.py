"""DIAGNOSTIC (pilot follow-up): strip TE n_eff vs width, 480-520 nm in 0.5 nm steps,
production mesh, cell-exact vs anisotropic. Residual from a cubic fit shows any
grid-periodic structure (a wall crosses a cell edge every 26.7 nm of wall motion,
i.e. every 53.3 nm of width, since both walls move)."""
import os, sys, json
os.environ.setdefault('OMP_NUM_THREADS', '1')
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'src'))
from multiprocessing import Pool
import numpy as np
import anisotropic_averaging as an, fullvector_sensitivity as fv
WS = np.round(np.arange(480, 520.01, 0.5), 2)
def f(t):
    avg, w = t
    return t, an.SOLVERS[avg]['strip'](w / 1000, 0.22, 0.0, fv.WL0, 'TE')
if __name__ == '__main__':
    with Pool(4) as p:
        r = dict(p.map(f, [(a, w) for a in ('cell_exact', 'anisotropic') for w in WS]))
    out = {}
    for a in ('cell_exact', 'anisotropic'):
        n = np.array([r[(a, w)] for w in WS])
        res = n - np.polyval(np.polyfit(WS, n, 3), WS)
        d = np.diff(n) / 0.5
        out[a] = dict(w=WS.tolist(), neff=n.tolist())
        print(f"{a:12s} fit residual rms {res.std():.2e}, max {abs(res).max():.2e};  local slope dn/dw (per nm): "
              f"min {d.min():.3e} max {d.max():.3e} mean {d.mean():.3e} (spread {100*(d.max()-d.min())/abs(d.mean()):.0f}% of mean)")
        print("   slope by 2.5 nm window: " + ' '.join(f"{x:.2e}" for x in (n[5::5] - n[:-5:5]) / 2.5))
    json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'diag_width_scan.json'), 'w'))
