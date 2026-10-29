"""
Extra femwell solves for the etch-depth / nuisance-parameter study.

Runs in the femwell virtual environment, from src/:
    C:/Users/unmol/.venvs/femwell/Scripts/python.exe etchdepth_femwell_batch.py time
    C:/Users/unmol/.venvs/femwell/Scripts/python.exe etchdepth_femwell_batch.py full [max_workers]
    (also: striptm, ribslab — the later strip TM and extra rib-slab units)

Structures: strip TE 400-600 nm, strip TM 500 nm, rib TE 300-600 nm at 70 and 90 nm slab.
New derivatives (central differences at h/2, h, 2h, same mesh r12 and conventions as femwell_batch.py):
  rib:   s   (slab thickness, rib top fixed; h = 2 nm), nsi, nox (h = 0.002)
  strip: nsi, nox
The w, t, a solves already exist in results_femwell/cache/solves_r12.jsonl and are reused.
Every solve is appended to results_etchdepth_femwell/cache/solves_r12.jsonl; a restart resumes.
"""
import json
import os
import sys

import femwell_batch as fb

OUT = os.path.join(fb.REPO, 'results_etchdepth_femwell')
CACHE = os.path.join(OUT, 'cache', 'solves_r12.jsonl')
STRIPS = [('strip', 'TE', w, None) for w in (400, 450, 500, 550, 600)] + [('strip', 'TM', 500, None)]
RIBS = [('rib', 'TE', w, s) for s in (70, 90) for w in (300, 400, 500, 600)]
# added later: index derivatives for the other strip TM widths (structure search pool)
STRIPS_TM_EXTRA = [('strip', 'TM', w, None) for w in (400, 450, 550, 600)]
# added later: slab / index derivatives for two rib widths at the other slabs (best slab for a two-width set)
RIBS_SLAB_EXTRA = [('rib', 'TE', w, s) for s in (50, 110, 130, 150) for w in (300, 500)]


def guesses(units):
    rows = fb.load_cache(os.path.join(fb.OUT, 'cache', 'solves_r12.jsonl'))
    return {u: rows[f'{u[0]}|{u[1]}|{u[2]}|{u[3]}|oxide|r12|nom|+0.0000']['n'][0] for u in units}


def tasks(units=STRIPS + RIBS):
    g = guesses(units)
    out = []
    for u in units:
        params = ['nsi', 'nox'] + (['s'] if u[0] == 'rib' else [])
        out += [t for t in fb.unit_tasks(*u, 'oxide', 'r12', g[u], params=params) if t[6] != 'nom']
    return out


if __name__ == '__main__':
    os.makedirs(os.path.join(OUT, 'cache'), exist_ok=True)
    os.makedirs(os.path.join(OUT, 'logs'), exist_ok=True)
    what = sys.argv[1]
    if len(sys.argv) > 3:   # optional: run on battery while charge >= this percent (striptm, ribslab, full)
        fb.BATTERY_FLOOR = int(sys.argv[3])
    if what == 'time':
        logpath = os.path.join(OUT, 'logs', 'timing_r12.log')
        sample = [t for t in tasks() if t[7] == fb.H[t[6]] and
                  (t[:4] in (('strip', 'TE', 500, None), ('rib', 'TE', 500, 90)))]
        fb.run(sample, CACHE, logpath, fb.choose_workers(2.0, max_workers=2, logpath=logpath), 2.0)
        rows = fb.load_cache(CACHE)
        for t in sample:
            r = rows[fb.key_of(t)]
            fb.log(f'timing {fb.key_of(t)}: {r["seconds"]:.1f} s (3 wavelengths, {r["n_elements"]} elements), '
                   f'n_eff {r["n"][0]:.6f}', logpath)
    elif what == 'full':
        cap = int(sys.argv[2]) if len(sys.argv) > 2 else 4
        logpath = os.path.join(OUT, 'logs', 'full_r12.log')
        failed = fb.run(tasks(), CACHE, logpath, fb.choose_workers(2.0, max_workers=cap, logpath=logpath), 2.0)
        with open(os.path.join(OUT, 'logs', 'full_r12_failed.json'), 'w') as f:
            json.dump(failed, f, indent=1)
    elif what == 'striptm':
        cap = int(sys.argv[2]) if len(sys.argv) > 2 else 4
        logpath = os.path.join(OUT, 'logs', 'striptm_r12.log')
        failed = fb.run(tasks(STRIPS_TM_EXTRA), CACHE, logpath, fb.choose_workers(2.0, max_workers=cap, logpath=logpath), 2.0)
        with open(os.path.join(OUT, 'logs', 'striptm_r12_failed.json'), 'w') as f:
            json.dump(failed, f, indent=1)
    elif what == 'ribslab':
        cap = int(sys.argv[2]) if len(sys.argv) > 2 else 4
        logpath = os.path.join(OUT, 'logs', 'ribslab_r12.log')
        failed = fb.run(tasks(RIBS_SLAB_EXTRA), CACHE, logpath, fb.choose_workers(2.0, max_workers=cap, logpath=logpath), 2.0)
        with open(os.path.join(OUT, 'logs', 'ribslab_r12_failed.json'), 'w') as f:
            json.dump(failed, f, indent=1)
