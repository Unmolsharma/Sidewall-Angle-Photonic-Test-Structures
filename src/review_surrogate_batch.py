"""
femwell solves for the review surrogate check (results_review/data/surrogate_check_design.json, written by
review_surrogate_check.py design). Same solver as montecarlo_batch.py.

Runs in the femwell virtual environment, from src/:
    C:/Users/unmol/.venvs/femwell/Scripts/python.exe review_surrogate_batch.py [max_workers] [battery_floor_percent]
Writes results_review/cache/solves_r12.jsonl and results_review/logs/surrogate_solves.log
"""
import json
import os
import sys

import femwell_batch as fb
import montecarlo_batch as mb

OUT = os.path.join(fb.REPO, 'results_review')
CACHE = os.path.join(OUT, 'cache', 'solves_r12.jsonl')


def key_of(task):
    geom, pol, W, S, w, t, a, slab, nsi, nox, guess = task
    return f'review|{geom}|{pol}|{W}|{S}|{w:.4f}|{t:.4f}|{a:.5f}|{slab}|{nsi:.7f}|{nox:.7f}'


def solve(task):
    r = mb.solve(task)
    r['key'] = key_of(task)
    return r


def tasks():
    D = json.load(open(os.path.join(OUT, 'data', 'surrogate_check_design.json')))
    return [(p['geom'], p['pol'], p['W'], p['S'], p['w'], p['t'], p['a'], p['slab'], p['nsi'], p['nox'], p['guess'])
            for c in D['cases'].values() for p in c['points']]


if __name__ == '__main__':
    os.makedirs(os.path.join(OUT, 'cache'), exist_ok=True)
    os.makedirs(os.path.join(OUT, 'logs'), exist_ok=True)
    cap = int(sys.argv[1]) if len(sys.argv) > 1 else 2
    if len(sys.argv) > 2:          # optional: run on battery while charge >= this percent
        fb.BATTERY_FLOOR = int(sys.argv[2])
    logpath = os.path.join(OUT, 'logs', 'surrogate_solves.log')
    failed = fb.run(tasks(), CACHE, logpath, max(cap, fb.choose_workers(2.0, max_workers=cap, logpath=logpath)), 2.0,
                    solver=solve, keyer=key_of)
    with open(os.path.join(OUT, 'logs', 'surrogate_solves_failed.json'), 'w') as f:
        json.dump(failed, f, indent=1)
