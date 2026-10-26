"""
Turn the femwell solve cache (results_femwell/cache/solves_<mesh>.jsonl, written
by femwell_batch.py in the femwell venv) into derivatives, step checks and a
sensitivity database with the same schema as data/sensitivity_db.json.

  python femwell_postprocess.py meshcheck            r12 vs r8 for strip TE 500, rib TE 500, air-clad strip TE 500
  python femwell_postprocess.py database r12         30 structures + slab sweep -> results_femwell/data/

Same derivative definitions and step criterion as recompute_sensitivities.py
(2% relative; near-zero entries, |D(h)| < 10% of the largest |D(h)| for the same
structure and output, pass if the step disagreement is within 2% of that largest
value). Run from src/.
"""
import csv
import json
import os
import sys

import numpy as np

from extend_db import make_record
from probe_db import save_db, build_probe_db
from recompute_sensitivities import H, FACTORS, DB_STRUCTS, SLAB_SWEEP, step_rows, as_s

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, 'results_femwell')
WL0, DLAM = 1.55, 0.005
WALL = dict(strip=0.22, swg=0.22, slot=0.22, rib=0.13)
KEYS = {('w', 'neff'): 'dneff_dw', ('t', 'neff'): 'dneff_dt', ('a', 'neff'): 'dneff_da',
        ('w', 'ng'): 'dng_dw', ('t', 'ng'): 'dng_dt', ('a', 'ng'): 'dng_da'}


def load(mesh):
    rows = {}
    path = os.path.join(OUT, 'cache', f'solves_{mesh}.jsonl')
    if not os.path.exists(path):
        return rows
    with open(path, encoding='utf-8') as f:
        for line in f:
            if line.strip():
                r = json.loads(line)
                rows[r['key']] = r
    return rows


def key(geom, pol, W, slab, clad, mesh, p, d):
    return f'{geom}|{pol}|{W}|{slab}|{clad}|{mesh}|{p}|{d:+.4f}'


def nng(rows, k):
    n0, npl, nmi = rows[k]['n']
    return n0, n0 - WL0 * (npl - nmi) / (2 * DLAM)


def unit(rows, geom, pol, W, slab, clad, mesh, params='wta'):
    """n_eff, n_g, D[(p, obs, f)] or None if any solve is missing."""
    try:
        n0, g0 = nng(rows, key(geom, pol, W, slab, clad, mesh, 'nom', 0.0))
        D = {}
        for p in params:
            for f in FACTORS:
                h = f * H[p]
                (n1, g1), (n2, g2) = nng(rows, key(geom, pol, W, slab, clad, mesh, p, h)), \
                    nng(rows, key(geom, pol, W, slab, clad, mesh, p, -h))
                D[(p, 'neff', f)] = (n1 - n2) / (2 * h)
                D[(p, 'ng', f)] = (g1 - g2) / (2 * h)
    except KeyError as e:
        return None, str(e)
    return (n0, g0, D), None


def ratio(geom, D, obs='neff'):
    r = D[('a', obs, 1.0)] / D[('w', obs, 1.0)]
    return r, r / (-WALL[geom] * 1000 * np.tan(np.radians(1.0)))


def meshcheck():
    r12, r8 = load('r12'), load('r8')
    # the r12 strip/rib/air reference was computed by femwell_crosscheck.py (same mesh settings,
    # 1-mode search); use the cached batch values when present, else that file
    ref = json.load(open(os.path.join(REPO, 'results_anisotropic', 'data', 'femwell_derivatives_r12.json')))
    out, lines = {}, []
    for geom, clad, params in (('strip', 'oxide', 'wta'), ('rib', 'oxide', 'wta'), ('strip', 'air', 'wa')):
        u8, err = unit(r8, geom, 'TE', 500, None, clad, 'r8', params)
        if u8 is None:
            lines.append(f'{geom} {clad}: r8 incomplete ({err})')
            continue
        u12, _ = unit(r12, geom, 'TE', 500, None, clad, 'r12', params) if r12 else (None, None)
        if u12 is None:
            R = ref[f'{geom}_TE_500_{clad}']
            D12 = {tuple(k.split('|')[:2]) + (float(k.split('|')[2]),): v for k, v in R['D'].items()}
            n12, g12, src = R['neff'], R['ng'], 'femwell_crosscheck r12'
        else:
            n12, g12, D12 = u12
            src = 'batch r12'
        n8, g8, D8 = u8
        o = out[f'{geom}_TE_500_{clad}'] = dict(r12_source=src, neff=[n12, n8], ng=[g12, g8], derivs={})
        lines.append(f'\n{geom} TE 500 {clad}: n_eff r12 {n12:.6f} r8 {n8:.6f}; n_g r12 {g12:.5f} r8 {g8:.5f}  (r12 from {src})')
        for obs in ('neff', 'ng'):
            big = max(abs(D8[(p, obs, 1.0)]) for p in params)
            for p in params:
                a, b = D12[(p, obs, 1.0)], D8[(p, obs, 1.0)]
                rel = (a - b) / abs(b)
                near = abs(b) < 0.1 * big
                ok = abs(rel) <= 0.02 or (near and abs(a - b) <= 0.02 * big)
                o['derivs'][f'd{obs}/d{p}'] = dict(r12=a, r8=b, rel_pct=100 * rel, near_zero=near, within_2pct=bool(ok))
                lines.append(f'  d{obs:4s}/d{p}  r12 {a:+.5e}  r8 {b:+.5e}  diff {100*rel:+.2f}%'
                             f'{"  (near-zero: " + f"{100*abs(a-b)/big:.2f}% of largest)" if near else ""}  {"OK" if ok else "DIFFERS"}')
            # step check on r8 too
        if 'w' in params and 'a' in params:
            o['ratio'] = {m: ratio(geom, D)[0] for m, D in (('r12', D12), ('r8', D8))}
            o['R'] = {m: ratio(geom, D)[1] for m, D in (('r12', D12), ('r8', D8))}
            o['ratio_ng'] = {m: ratio(geom, D, 'ng')[0] for m, D in (('r12', D12), ('r8', D8))}
            lines.append(f'  (dn_eff/da)/(dn_eff/dw): r12 {o["ratio"]["r12"]:+.4f}  r8 {o["ratio"]["r8"]:+.4f} nm/deg;  '
                         f'R r12 {o["R"]["r12"]:.4f}  r8 {o["R"]["r8"]:.4f}')
        rows8 = step_rows('femwell_r8', f'{geom}_{clad}', D8) if params == 'wta' else []
        o['r8_step_fail'] = [f"d{r['output']}/d{r['param']}" for r in rows8 if r['verdict'] == 'FAIL']
        lines.append(f'  r8 step check failures: {o["r8_step_fail"] or "none"}')
    out['all_within_2pct'] = all(d['within_2pct'] for k, v in out.items() if isinstance(v, dict)
                                 for d in v.get('derivs', {}).values())
    lines.append(f'\nall r12 vs r8 derivatives within 2% (near-zero rule applied): {out["all_within_2pct"]}')
    print('\n'.join(lines))
    os.makedirs(os.path.join(OUT, 'data'), exist_ok=True)
    with open(os.path.join(OUT, 'data', 'meshcheck_r12_vs_r8.json'), 'w') as f:
        json.dump(out, f, indent=1)


def database(mesh):
    rows = load(mesh)
    d = os.path.join(OUT, 'data')
    os.makedirs(d, exist_ok=True)
    settings = dict(solver='femwell 0.1.12 (FEM, order 2, conforming mesh)', mesh=mesh,
                    core_resolution_um={'r12': 0.012, 'r8': 0.008}[mesh], distance_um=0.3, max_element_um=0.2,
                    half_x=1.4, half_y=1.2)
    records, sweep, steps, skipped = [], [], [], {}
    for geom, pol, W in DB_STRUCTS:
        u, err = unit(rows, geom, pol, W, None, 'oxide', mesh)
        if u is None:
            skipped[f'{geom}_{pol}_{W}'] = err
            continue
        n0, g0, D = u
        gap = (max(120 - (W - 180), 20)) if geom == 'slot' else None
        s = dict(settings, half_x=1.6) if geom == 'slot' else settings
        records.append(make_record(geom, pol, W, as_s(n0, g0, D), gap, s,
                                   convergence_status=f'femwell {mesh}; mesh check in results_femwell/data/meshcheck_r12_vs_r8.json',
                                   validation_status='femwell; step-size check in results_femwell/data/step_check.csv'))
        steps += step_rows('femwell', records[-1]['id'], D)
    for s_nm, W in SLAB_SWEEP:
        u, err = unit(rows, 'rib', 'TE', W, s_nm, 'oxide', mesh)
        if u is None:
            skipped[f'rib_slab{s_nm}_TE_{W}'] = err
            continue
        n0, g0, D = u
        sweep.append(dict(slab_nm=s_nm, width_nm=W, **as_s(n0, g0, D)))
        steps += step_rows('femwell', f'rib_slab{s_nm}_TE_{W}', D)
    # a slab whose four widths are not all complete would give a gain from fewer rib widths
    # than the other slabs; drop it entirely rather than compare unlike sets
    widths = {}
    for e in sweep:
        widths.setdefault(e['slab_nm'], []).append(e['width_nm'])
    partial = [s for s, ws in widths.items() if len(ws) < len({W for _, W in SLAB_SWEEP})]
    for s in partial:
        skipped[f'rib_slab{s}_TE (whole slab)'] = f'only widths {widths[s]} complete; slab dropped'
    sweep = [e for e in sweep if e['slab_nm'] not in partial]
    steps = [r for r in steps if not any(r['structure'].startswith(f'rib_slab{s}_') for s in partial)]
    save_db(records, os.path.join(d, 'sensitivity_db.csv'), os.path.join(d, 'sensitivity_db.json'))
    save_db(build_probe_db(records), os.path.join(d, 'probe_db.csv'), os.path.join(d, 'probe_db.json'))
    with open(os.path.join(d, 'rib_slab_sweep_db_femwell.json'), 'w') as f:
        json.dump(sweep, f, indent=2)
    with open(os.path.join(d, 'step_check.csv'), 'w', newline='') as f:
        wr = csv.DictWriter(f, fieldnames=list(steps[0].keys()))
        wr.writeheader()
        wr.writerows(steps)
    with open(os.path.join(d, 'skipped.json'), 'w') as f:
        json.dump(skipped, f, indent=1)
    notes = [r for r in rows.values() if r.get('notes')]
    print(f'{len(records)} database structures, {len(sweep)} slab-sweep structures, skipped {len(skipped)}: {skipped}')
    print(f'{len(notes)} solves needed the 6-mode fallback')
    for scope, pred in (('30 database structures', lambda r: 'slab' not in r['structure']),
                        ('slab sweep', lambda r: 'slab' in r['structure'])):
        sel = [r for r in steps if pred(r)]
        c = {v: sum(r['verdict'] == v for r in sel) for v in ('PASS', 'PASS_NEAR_ZERO_ABS', 'FAIL')}
        print(f'[femwell] {scope:24s} {len(sel)} derivatives: PASS {c["PASS"]}, near-zero PASS {c["PASS_NEAR_ZERO_ABS"]}, FAIL {c["FAIL"]}')
    for r in steps:
        if r['verdict'] != 'PASS':
            print(f"  {r['verdict']:20s} {r['structure']:22s} d{r['output']}/d{r['param']}  "
                  f"D(h/2,h,2h)=({r['D_half']:.4e}, {r['D_h']:.4e}, {r['D_2h']:.4e})  rel dev {r['max_rel_dev']:.2%}  "
                  f"|D|/largest {r['frac_of_largest']:.3f}")


if __name__ == '__main__':
    if sys.argv[1] == 'meshcheck':
        meshcheck()
    else:
        database(sys.argv[2])
