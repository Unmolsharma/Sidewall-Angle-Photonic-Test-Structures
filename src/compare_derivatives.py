"""
Compare every derivative method against the femwell reference for strip TE 500
and rib TE 500 (90 nm slab):
  femwell        : conforming FEM mesh, FD in geometry (femwell_crosscheck.py)
  perturbation   : boundary-perturbation integral from EMpy fields on a
                   boundary-aligned grid (perturbation_derivatives.py), refine 1 and 2
  aligned FD     : EMpy FD with grids whose lines follow the moving boundary
                   (same cell counts, lines morph with w or t), refine 1-3;
                   width and thickness only (a tilted wall cannot be grid-aligned)
  grid-shift     : FD averaged over N grid offsets (gridshift_validate.py), N = 1, 4, 8
                   (N = 1 is plain FD on the production mesh)

Pass mark for a method: within 5% of femwell (n_eff derivatives), and the
angle/width ratio within 3%.

Run from src/ after femwell_crosscheck.py and gridshift_validate.py:
    python compare_derivatives.py
Writes results_anisotropic/data/derivative_method_comparison.json
"""
import json
import os
from multiprocessing import Pool

import numpy as np

import fullvector_sensitivity as fv
import cellexact_averaging as ce
import perturbation_derivatives as pd

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D = os.path.join(REPO, 'results_anisotropic', 'data')
H = dict(w=5.0, t=2.0, a=1.0)
WALL = dict(strip=0.22, rib=0.13)


# --- aligned ("morphed") FD: fixed cell counts, lines follow the boundaries ---
def _axis(breaks, core, half, counts, ncl):
    c0, c1 = core
    pts = sorted({c0, c1, *[b for b in breaks if c0 < b < c1]})
    ax = [np.linspace(-half, c0, ncl, endpoint=False)]
    for (a, b), n in zip(zip(pts[:-1], pts[1:]), counts):
        ax.append(np.linspace(a, b, n, endpoint=False))
    ax.append(np.linspace(c1, half, ncl))
    return np.concatenate(ax)


def _counts(breaks, core, step):
    c0, c1 = core
    pts = sorted({c0, c1, *[b for b in breaks if c0 < b < c1]})
    return [max(int(np.ceil((b - a) / step - 1e-9)), 1) for a, b in zip(pts[:-1], pts[1:])]


def aligned_solve(task):
    geom, w, t, r, wl = task
    ncl = int(round(18 * min(r, 2)))
    yb = [0.0, t] + ([0.09] if geom == 'rib' else [])
    yb0 = [0.0, 0.22] + ([0.09] if geom == 'rib' else [])
    cx = _counts([-0.25, 0.25], (-0.4, 0.4), 0.8 / 30 / r)
    cy = _counts(yb0, (-0.1, 0.5), 0.02 / r)
    X = _axis([-w / 2, w / 2], (-0.4, 0.4), 1.4, cx, ncl)
    Y = _axis(yb, (-0.1, 0.5), 1.2, cy, ncl)
    f = ce.fill_rib(X, Y, w, t, 0.0, 0.09) if geom == 'rib' else ce.fill_strip(X, Y, w, t, 0.0)
    assert ((f > 1e-9) & (f < 1 - 1e-9)).sum() == 0
    fn = ce.solve_neff_rib_ce if geom == 'rib' else ce.solve_neff_strip_ce
    return task, fn(w, t, 0.0, wl, 'TE', X=X, Y=Y)


def aligned_fd():
    wls = (fv.WL0, fv.WL0 + 0.005, fv.WL0 - 0.005)
    tasks = []
    for geom in ('strip', 'rib'):
        for r in (1, 2, 3):
            for w, t in ((0.5, 0.22), (0.505, 0.22), (0.495, 0.22), (0.5, 0.222), (0.5, 0.218),
                         (0.5025, 0.22), (0.4975, 0.22), (0.51, 0.22), (0.49, 0.22), (0.5, 0.221), (0.5, 0.219),
                         (0.5, 0.224), (0.5, 0.216)):
                tasks += [(geom, w, t, r, wl) for wl in wls]
    with Pool(4) as p:
        res = dict(p.map(aligned_solve, tasks, chunksize=1))

    def ng(geom, w, t, r):
        n0, a, b = (res[(geom, w, t, r, wl)] for wl in wls)
        return n0, n0 - fv.WL0 * (a - b) / 0.01
    out = {}
    for geom in ('strip', 'rib'):
        for r in (1, 2, 3):
            o = out[f'{geom}|refine{r}'] = dict(neff=ng(geom, 0.5, 0.22, r)[0], ng=ng(geom, 0.5, 0.22, r)[1])
            for p, pairs in (('w', {0.5: ((0.5025, 0.22), (0.4975, 0.22)), 1.0: ((0.505, 0.22), (0.495, 0.22)),
                                    2.0: ((0.51, 0.22), (0.49, 0.22))}),
                             ('t', {0.5: ((0.5, 0.221), (0.5, 0.219)), 1.0: ((0.5, 0.222), (0.5, 0.218)),
                                    2.0: ((0.5, 0.224), (0.5, 0.216))})):
                for f, (a, b) in pairs.items():
                    (n1, g1), (n2, g2) = ng(geom, *a, r), ng(geom, *b, r)
                    o[f'{p}|neff|{f}'] = (n1 - n2) / (2 * f * H[p])
                    o[f'{p}|ng|{f}'] = (g1 - g2) / (2 * f * H[p])
    return out


def pct(x, ref):
    return 100 * (x - ref) / abs(ref)


if __name__ == '__main__':
    fem = json.load(open(os.path.join(D, 'femwell_derivatives_r12.json')))
    p8 = os.path.join(D, 'femwell_derivatives_r8.json')
    fem8 = json.load(open(p8)) if os.path.exists(p8) else fem   # r8 mesh check optional
    gs = json.load(open(os.path.join(D, 'gridshift_validation.json')))
    out = dict(perturbation={}, aligned_fd=aligned_fd())
    for geom in ('strip', 'rib'):
        for r in (1.0, 2.0):
            out['perturbation'][f'{geom}|refine{r:g}'] = pd.derivatives(geom, 500, refine=r)

    rows = []
    for geom in ('strip', 'rib'):
        ref = fem[f'{geom}_TE_500_oxide']
        ref8 = fem8[f'{geom}_TE_500_oxide']
        print(f"\n===== {geom} TE 500 (oxide clad) =====   femwell n_eff {ref['neff']:.5f}, n_g {ref['ng']:.4f}")
        print("  femwell step check: " + ', '.join(f"{k} {'PASS' if v['passed'] else 'FAIL'} ({100*v['max_rel_dev']:.2f}%)"
                                                  for k, v in ref['step_check'].items()))
        for obs in ('neff', 'ng'):
            for p in ('w', 't', 'a'):
                k = f'{p}|{obs}|1.0'
                fr = ref['D'][k]
                line = dict(structure=geom, derivative=f'd{obs}/d{p}', femwell=fr,
                            femwell_r8_pct=pct(ref8['D'][k], fr))
                for r in (1, 2):
                    line[f'perturbation_r{r}_pct'] = pct(out['perturbation'][f'{geom}|refine{r}'][f'd{obs}_d{p}'], fr)
                if p != 'a':
                    for r in (1, 2, 3):
                        line[f'aligned_fd_r{r}_pct'] = pct(out['aligned_fd'][f'{geom}|refine{r}'][k], fr)
                for avg in ('cell_exact', 'anisotropic'):
                    g = gs[f'{avg}|{geom}_TE_500']
                    for n in (1, 4, 8):
                        line[f'{avg}_N{n}_pct'] = pct(g['30/18'][f'N{n}'][k], fr)
                    line[f'{avg}_N8_45/18_pct'] = pct(g['45/18']['N8'][k], fr)
                rows.append(line)
                print(f"  d{obs:4s}/d{p}  femwell {fr:+.4e} (r8 {line['femwell_r8_pct']:+.2f}%)  |  "
                      f"pert r1/r2 {line['perturbation_r1_pct']:+6.1f}/{line['perturbation_r2_pct']:+6.1f}%  |  "
                      + (f"aligned r1/r2/r3 {line['aligned_fd_r1_pct']:+5.1f}/{line['aligned_fd_r2_pct']:+5.1f}/{line['aligned_fd_r3_pct']:+5.1f}%  |  "
                         if p != 'a' else ' ' * 43 + '|  ')
                      + '  '.join(f"{a[:4]} N1/N4/N8 {line[f'{a}_N1_pct']:+5.1f}/{line[f'{a}_N4_pct']:+5.1f}/{line[f'{a}_N8_pct']:+5.1f}% (45/18 {line[f'{a}_N8_45/18_pct']:+5.1f}%)"
                                  for a in ('cell_exact', 'anisotropic')))
        # angle/width ratio
        R = {}
        for obs in ('neff', 'ng'):
            fr = ref['D'][f'a|{obs}|1.0'] / ref['D'][f'w|{obs}|1.0']
            R[obs] = dict(femwell=fr, femwell_r8=ref8['D'][f'a|{obs}|1.0'] / ref8['D'][f'w|{obs}|1.0'])
            for r in (1, 2):
                q = out['perturbation'][f'{geom}|refine{r}']
                R[obs][f'perturbation_r{r}'] = q[f'd{obs}_da'] / q[f'd{obs}_dw']
            for avg in ('cell_exact', 'anisotropic'):
                g = gs[f'{avg}|{geom}_TE_500']
                for n in (1, 4, 8):
                    R[obs][f'{avg}_N{n}'] = g['30/18'][f'N{n}'][f'a|{obs}|1.0'] / g['30/18'][f'N{n}'][f'w|{obs}|1.0']
                R[obs][f'{avg}_N8_45/18'] = g['45/18']['N8'][f'a|{obs}|1.0'] / g['45/18']['N8'][f'w|{obs}|1.0']
            R[obs]['R_normalized_femwell'] = fr / (-WALL[geom] * 1000 * np.tan(np.radians(1)))
            print(f"  ratio (d{obs}/da)/(d{obs}/dw) [nm/deg]: femwell {fr:+.4f} (R = {R[obs]['R_normalized_femwell']:.4f}; r8 {R[obs]['femwell_r8']:+.4f})  "
                  + '  '.join(f"{k} {v:+.4f}" for k, v in R[obs].items() if k not in ('femwell', 'femwell_r8', 'R_normalized_femwell')))
        out[f'ratio_{geom}'] = R
    out['rows'] = rows
    with open(os.path.join(D, 'derivative_method_comparison.json'), 'w') as f:
        json.dump(out, f, indent=1, default=float)
