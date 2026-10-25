"""
Step 1 check for cellexact_averaging.py: is silicon area conserved?

(a) Exactness: total modelled silicon area vs the analytic cross-section
    area, and per-cell fills vs brute-force supersampling of the boundary cells.
(b) Conservation under boundary motion: for every parameter (w, t, a, slab,
    slot rail width), sweep the nominal position finely across the grid and
    compute
        rate = d(modelled Si area)/d(param) / d(true Si area)/d(param)
    (ideal 1.000 everywhere) with the ORIGINAL and the CELL-EXACT averaging.
    For the original factories the modelled area is recovered from eps at
    the cell centres: fill = (eps - n_ox^2)/(n_si^2 - n_ox^2).

Run from the repo root:  python src/verify_cellexact.py
Writes results_corrected/data/area_conservation.json
"""
import os
import json

import numpy as np

import fullvector_sensitivity as fv
import rib_waveguide as rw
import slot_waveguide as sw
import cellexact_averaging as ce

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, 'results_corrected', 'data')
D2R = np.pi / 180


def area_from_eps(eps, X, Y):
    fill = (eps - fv.N_OX**2) / (fv.N_SI**2 - fv.N_OX**2)
    return (fill * ce.cell_areas(X, Y)).sum()


def centres(X, Y):
    return (X[:-1] + X[1:]) / 2, (Y[:-1] + Y[1:]) / 2


# modelled area as a function of one parameter, for each geometry/averaging
def model_area(geom, averaging, p):
    if geom == 'slot':
        X, Y = sw.X_SLOT, sw.Y_SLOT
    else:
        X, Y = fv.X, fv.Y
    xc, yc = centres(X, Y)
    if geom == 'strip':
        if averaging == 'original':
            return area_from_eps(fv.epsfunc_factory(p['w'], p['t'], p['a'])(xc, yc), X, Y)
        return (ce.fill_strip(X, Y, p['w'], p['t'], p['a']) * ce.cell_areas(X, Y)).sum()
    if geom == 'rib':
        if averaging == 'original':
            return area_from_eps(rw.epsfunc_factory_rib(p['w'], p['t'], p['a'], p['s'])(xc, yc), X, Y)
        return (ce.fill_rib(X, Y, p['w'], p['t'], p['a'], p['s']) * ce.cell_areas(X, Y)).sum()
    if geom == 'slot':
        rail, gap = sw.rail_gap_from_width(p['w'])
        if averaging == 'original':
            return area_from_eps(sw.epsfunc_factory_slot(rail, gap, p['t'], p['a'])(xc, yc), X, Y)
        return (ce.fill_slot(X, Y, rail, gap, p['t'], p['a']) * ce.cell_areas(X, Y)).sum()


def true_area(geom, p):
    ta = np.tan(p['a'] * D2R)
    if geom == 'strip':
        return p['t'] * p['w'] - ta * p['t']**2
    if geom == 'rib':
        h = p['t'] - p['s']
        return p['s'] * (fv.X[-1] - fv.X[0]) + h * p['w'] - ta * h**2
    if geom == 'slot':
        return 2 * (p['t'] * p['w'] - ta * p['t']**2)


def rate_sweep(geom, key, centre, span, step, base, dh):
    """rate(position) = [A(pos+dh)-A(pos-dh)]_model / [A(pos+dh)-A(pos-dh)]_true"""
    out = {}
    positions = np.arange(centre - span, centre + span + 1e-12, step)
    for avg in ('original', 'cell_exact'):
        rates = []
        for v in positions:
            pp, pm = dict(base), dict(base)
            pp[key], pm[key] = v + dh, v - dh
            dm = model_area(geom, avg, pp) - model_area(geom, avg, pm)
            dt = true_area(geom, pp) - true_area(geom, pm)
            rates.append(dm / dt)
        rates = np.array(rates)
        nominal = rates[np.argmin(abs(positions - centre))]
        out[avg] = dict(min=float(rates.min()), max=float(rates.max()), at_nominal=float(nominal))
    return out


def supersample_check(n_sub=400):
    """Per-cell fills vs brute-force n_sub x n_sub point sampling, on a tapered
    strip and a tapered rib, over every partially filled cell."""
    X, Y = fv.X, fv.Y
    res = {}
    for name, fill, inside in [
        ('strip w500 t220 a3', ce.fill_strip(X, Y, 0.5, 0.22, 3.0),
         lambda x, y: (y >= 0) & (y <= 0.22) & (abs(x) <= 0.25 - np.tan(3 * D2R) * y)),
        ('rib w450 t221 a-1 slab91', ce.fill_rib(X, Y, 0.45, 0.221, -1.0, 0.091),
         lambda x, y: ((y >= 0) & (y <= 0.091)) |
                      ((y > 0.091) & (y <= 0.221) & (abs(x) <= 0.225 - np.tan(-1 * D2R) * (y - 0.091)))),
    ]:
        part = np.argwhere((fill > 1e-9) & (fill < 1 - 1e-9))
        u = (np.arange(n_sub) + 0.5) / n_sub
        errs = []
        for i, j in part:
            xs = X[i] + u * (X[i + 1] - X[i])
            ys = Y[j] + u * (Y[j + 1] - Y[j])
            xx, yy = np.meshgrid(xs, ys, indexing='ij')
            errs.append(abs(inside(xx, yy).mean() - fill[i, j]))
        res[name] = dict(n_partial_cells=len(part), max_abs_fill_error=float(max(errs)))
    return res


if __name__ == '__main__':
    os.makedirs(OUT, exist_ok=True)
    report = {}

    # the window the original factories actually use: mean spacing of the cell
    # centres EMpy passes in
    xc0, yc0 = centres(fv.X, fv.Y)
    dxm, dym = np.diff(xc0).mean(), np.diff(yc0).mean()
    core_dx = np.diff(fv.X)[(fv.X[:-1] >= -0.4) & (fv.X[:-1] < 0.4)].mean()
    core_dy = np.diff(fv.Y)[(fv.Y[:-1] >= -0.1) & (fv.Y[:-1] < 0.5)].mean()
    report['grid'] = dict(original_window_nm=[dxm * 1000, dym * 1000],
                          core_cell_nm=[core_dx * 1000, core_dy * 1000])
    print(f"production grid: original averaging window = {dxm*1000:.1f} x {dym*1000:.1f} nm "
          f"(mean spacing); core cells = {core_dx*1000:.1f} x {core_dy*1000:.1f} nm")

    print('\n(a) total silicon area, modelled vs analytic (um^2)')
    report['total_area'] = {}
    cases = [('strip', dict(w=0.5, t=0.22, a=0.0)), ('strip', dict(w=0.4137, t=0.2213, a=1.7)),
             ('strip', dict(w=0.6, t=0.219, a=-1.0)),
             ('rib', dict(w=0.5, t=0.22, a=0.0, s=0.09)), ('rib', dict(w=0.4321, t=0.2207, a=2.0, s=0.0613)),
             ('slot', dict(w=0.18, t=0.22, a=0.0)), ('slot', dict(w=0.1953, t=0.218, a=1.0))]
    for geom, p in cases:
        a_true = true_area(geom, p)
        a_ce = model_area(geom, 'cell_exact', p)
        a_or = model_area(geom, 'original', p)
        key = f"{geom} {p}"
        report['total_area'][key] = dict(true=a_true, cell_exact=a_ce, original=a_or,
                                         cell_exact_rel_err=(a_ce - a_true) / a_true,
                                         original_rel_err=(a_or - a_true) / a_true)
        print(f"  {key:52s} true={a_true:.10f}  cell-exact err={(a_ce-a_true)/a_true:+.2e}  "
              f"original err={(a_or-a_true)/a_true:+.2e}")

    print('\n    per-cell fill vs 400x400 brute-force sampling (boundary cells only):')
    report['supersample'] = supersample_check()
    for k, v in report['supersample'].items():
        print(f"  {k:28s} {v['n_partial_cells']} partial cells, max |fill error| = "
              f"{v['max_abs_fill_error']:.2e} (sampling resolution 1/400 = 2.5e-3)")

    print('\n(b) rate of modelled Si change per unit of true Si change (ideal 1.000),')
    print('    swept over boundary positions spanning several cells, dh = 0.5 nm (1 deg for a)')
    sweeps = [
        ('strip', 'w', 0.50, 0.03, 0.0005, dict(w=0.5, t=0.22, a=0.0), 0.0005, 'strip width (vertical sidewall)'),
        ('strip', 't', 0.22, 0.03, 0.0005, dict(w=0.5, t=0.22, a=0.0), 0.0005, 'strip thickness (horizontal top)'),
        ('strip', 'a', 0.0, 0.0, 1.0, dict(w=0.5, t=0.22, a=0.0), 1.0, 'strip sidewall angle (+-1 deg at 0)'),
        ('rib', 's', 0.09, 0.03, 0.0005, dict(w=0.5, t=0.22, a=0.0, s=0.09), 0.0005, 'rib slab height'),
        ('rib', 'w', 0.50, 0.03, 0.0005, dict(w=0.5, t=0.22, a=0.0, s=0.09), 0.0005, 'rib width'),
        ('rib', 'a', 0.0, 0.0, 1.0, dict(w=0.5, t=0.22, a=0.0, s=0.09), 1.0, 'rib sidewall angle (+-1 deg at 0)'),
        ('slot', 'w', 0.18, 0.02, 0.0005, dict(w=0.18, t=0.22, a=0.0), 0.0005, 'slot rail width (gap co-varies)'),
    ]
    report['rates'] = {}
    for geom, key, c, span, step, base, dh, name in sweeps:
        r = rate_sweep(geom, key, c, span, step, base, dh)
        report['rates'][name] = r
        print(f"  {name:36s} original: {r['original']['min']:.3f}..{r['original']['max']:.3f} "
              f"(nominal {r['original']['at_nominal']:.3f})   cell-exact: "
              f"{r['cell_exact']['min']:.6f}..{r['cell_exact']['max']:.6f}")

    with open(os.path.join(OUT, 'area_conservation.json'), 'w') as f:
        json.dump(report, f, indent=1)
    print(f"\nSaved {os.path.join(OUT, 'area_conservation.json')}")
