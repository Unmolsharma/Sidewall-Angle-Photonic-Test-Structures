"""
Provenance audit: trace every number in results_montecarlo/HEADLINE_NUMBERS.md and in every
results*/SUMMARY.md to a data file, and list citations and assumptions with their verification status.

Run from src/ (main environment):  python provenance_audit.py
Writes results_montecarlo/provenance/{AUDIT.md, audit.json}

Method
- Numbers audited: tokens with a decimal point or an exponent (e.g. 0.162, 6.5e-6, 1,009.5); integers
  are skipped (they are mostly labels: widths, slabs, set sizes, counts, dates). Text inside `code
  spans` (file paths) is skipped.
- A number is TRACED if some value in its source pool rounds to it at the precision shown (also
  accepting |value|, value x 100 for percentages and value - 1 or 1 - value as a relative change).
  * HEADLINE_NUMBERS.md: the pool of a table row is the files cited in that row (backticked paths);
    a row without its own citation uses the paths cited in the section text above the table.
  * SUMMARY.md: the pool of a section is the files cited in it and in its parent sections; a section
    with none uses every data/ and logs/ folder of that results folder.
  Files holding more than MAX_VALUES values are left out, and the chance-match rate of the pools (share
  of random 3-digit numbers that would 'trace') is reported with the results.
  A match by coincidence is possible, most of all in large pools, so TRACED means "a matching value
  exists in the pool", not proof; UNTRACED means no file in the pool holds it (derived ratios, values
  quoted from other folders, hand calculations, typos).
- Citations: DOIs and arXiv identifiers found in the audited files and in src/*.py are listed, each
  VERIFIED only if it is in the VERIFIED table below (checked against the source text in a session).
- Assumptions: the registry ASSUMPTIONS below (unsourced modelling inputs).
"""
import csv
import glob
import json
import os
import re

import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, 'results_montecarlo', 'provenance')
HEADLINE = os.path.join(REPO, 'results_montecarlo', 'HEADLINE_NUMBERS.md')
NUM = re.compile(r'(?<![\w.])[-−+]?(?:\d{1,3}(?:,\d{3})+|\d+)(?:\.\d+)?(?:[eE][-−+]?\d+)?(?![\w])')
PATH = re.compile(r'`([^`]+)`')
SKIP_DIRS = ('cache', 'provenance', '__pycache__')

VERIFIED = {
    '10.48550/arXiv.2306.15643': 'Janz et al. 2023 (arXiv:2306.15643 v1): full text read on 2026-10-29; '
                                 '"±0.1 pm repeatability of the laser" (Fig. 6 discussion) and "Wavelength '
                                 'changes down to 1 pm or less can be routinely measured" (Introduction)',
    'arXiv:2306.15643': 'same as 10.48550/arXiv.2306.15643',
}

ASSUMPTIONS = [
    ('sigma_neff 6.5e-6, sigma_ng 1e-3, sigma_phase 1e-3 rad (placeholder noise)', 'all CRBs unless stated',
     'UNSOURCED placeholders; the sourced 1 pm level gives 2.75e-6 / 1.07e-3 / 1.11e-3 (noise_sourced.json)'),
    ('Conversion of sigma_lambda to sigma_neff, sigma_ng, sigma_phase: n_g = 4.26, L = 100 um, one FSR',
     'noise_sourced.json', 'ASSUMPTION (our conversion; L and n_g chosen, not measured)'),
    ('MZI fringe read as precisely as a ring resonance', 'noise_sourced.json', 'UNVERIFIED'),
    ('Rib noise multiplier 1.0x (conservative) / 0.2x (optimistic)', 'all rib rows',
     '0.2x from a wide-rib loss ratio quoted in src/noise_sourcing.py (citation UNVERIFIED)'),
    ('Etch-depth prior 2 nm', 'realistic priors', 'ASSUMPTION (no verified source)'),
    ('n_ox prior 1e-3', 'realistic priors', 'ASSUMPTION'),
    ('n_Si prior 1e-4', 'realistic priors', 'quoted from Frey et al. (arXiv:physics/0606168) in '
                                            'results_priors_search/SUMMARY.md; citation UNVERIFIED in this audit'),
    ('Temperature prior 1 K', 'stage E realistic', 'ASSUMPTION'),
    ('Rib-strip angle link 1, 2, 5 deg', 'priors table', 'ASSUMPTION (scenarios)'),
    ('n_Si 3.48, n_ox 1.44, no dispersion; lambda 1.55 um; t 220 nm', 'all solves', 'ASSUMPTION (nominal values)'),
    ('Strip symmetry: exact R = 1 for oxide-clad strips', 'all CRBs with strips', 'DERIVED (mirror symmetry; '
     'femwell residuals agree to <= 8.5e-3, results_montecarlo/data/strips_info.json)'),
    ('Monte Carlo truth = nominal; box +-6 CRB sigma; 5 starts', 'mc_results.json', 'METHOD CHOICE'),
]


def fmt_step(tok):
    """Half a unit in the last shown digit of a number token."""
    t = tok.replace(',', '').replace('−', '-').lstrip('+-')
    mant, _, exp = t.lower().partition('e')
    dec = len(mant.split('.')[1]) if '.' in mant else 0
    return 0.5 * 10 ** (-dec) * (10 ** int(exp) if exp else 1)


def to_float(tok):
    return float(tok.replace(',', '').replace('−', '-'))


def numbers_in_text(s):
    return [float(m.replace(',', '').replace('−', '-')) for m in NUM.findall(s) if m.strip('+-−')]


def file_values(path, _cache={}):
    if path in _cache:
        return _cache[path]
    vals = []
    try:
        if path.endswith('.json'):
            def walk(o):
                if isinstance(o, bool):
                    return
                if isinstance(o, (int, float)):
                    vals.append(float(o))
                elif isinstance(o, dict):
                    for k, v in o.items():
                        walk(v)
                        try:
                            vals.append(float(k))
                        except ValueError:
                            pass
                elif isinstance(o, list):
                    for v in o:
                        walk(v)
            walk(json.load(open(path, encoding='utf-8')))
        elif path.endswith('.csv'):
            with open(path, newline='', encoding='utf-8') as f:
                for row in csv.reader(f):
                    for c in row:
                        try:
                            vals.append(float(c))
                        except ValueError:
                            vals += numbers_in_text(c)
        else:
            vals = numbers_in_text(open(path, encoding='utf-8', errors='replace').read())
    except (OSError, ValueError, json.JSONDecodeError):
        pass
    a = np.array([v for v in vals if np.isfinite(v)])
    _cache[path] = a
    return a


MAX_VALUES = 5000      # files holding more values than this (designs, per-draw tables) are left out of pools


def rel(p):
    return os.path.relpath(p, REPO).replace(os.sep, '/')


def traced(tok, pool, percent=False):
    """Pool file holding a value that rounds to tok at its shown precision (|v|; for a % token also 100 v
    and 100 |v - 1|), or None."""
    x, step = abs(to_float(tok)), fmt_step(tok)
    tol = step * (1 + 1e-9) + 1e-15
    for p in pool:
        v = np.abs(file_values(p))
        if v.size == 0:
            continue
        cands = (v, v * 100, np.abs(v - 1) * 100) if percent else (v,)
        if any(np.any(np.abs(c - x) <= tol) for c in cands):
            return rel(p)
    return None


def strip_code(s):
    return PATH.sub(lambda m: ' ' * len(m.group(0)), s)


def resolve(paths, base):
    out = []
    for p in paths:
        p = p.strip().split(' ')[0].rstrip(',;:')
        if '{' in p and '}' in p:             # e.g. data/surrogate_accuracy.{json,csv}
            pre, rest = p.split('{', 1)
            alts, post = rest.split('}', 1)
            out += [x for alt in alts.split(',') for x in resolve([pre + alt + post], base)]
            continue
        for root in (base, REPO):
            full = os.path.normpath(os.path.join(root, p))
            if os.path.isfile(full):
                out.append(full)
                break
            if os.path.isdir(full) and p not in ('.', '/', ''):
                out += [f for f in glob.glob(os.path.join(full, '**', '*'), recursive=True) if os.path.isfile(f)]
                break
    return [f for f in dict.fromkeys(out) if f.endswith(('.json', '.csv', '.log', '.md', '.txt'))
            and not any(f'{os.sep}{d}{os.sep}' in f for d in SKIP_DIRS)
            and os.path.basename(f) not in ('SUMMARY.md', 'HEADLINE_NUMBERS.md')
            and file_values(f).size <= MAX_VALUES]


def folder_pool(folder):
    """Fallback pool for a section that cites no file: every data/ and logs/ folder under the results folder."""
    dirs = [os.path.relpath(d, folder) for d in glob.glob(os.path.join(folder, '**', ''), recursive=True)
            if os.path.basename(os.path.dirname(d)) in ('data', 'logs')]
    return resolve(dirs, folder)


def tokens(line):
    """(token, is_percent) for every audited number in a line (decimal point or exponent; code spans skipped)."""
    s = strip_code(line)
    out = []
    for m in NUM.finditer(s):
        tok = m.group(0)
        if re.search(r'[.eE]', tok.replace(',', '')):
            out.append((tok, bool(re.match(r'(?:\s?[–-]\s?[-−+]?[\d.,]+)?\s?%', s[m.end():]))))
    return out


def sections(text):
    """Split markdown into (start_line, lines) at headings."""
    out, cur, start = [], [], 1
    for i, line in enumerate(text.splitlines(), 1):
        if line.startswith('#') and cur:
            out.append((start, cur))
            cur, start = [], i
        cur.append(line)
    out.append((start, cur))
    return out


def false_match_rate(pool, n=200, seed=1):
    """Share of random 3-significant-digit numbers in 0.01-5 that would 'trace' in this pool."""
    if not pool:
        return 0.0
    rng = np.random.default_rng(seed)
    return float(np.mean([traced(f'{x:.3g}', pool) is not None for x in rng.uniform(0.01, 5, n)]))


def audit_summary(path):
    text = open(path, encoding='utf-8').read()
    folder = os.path.dirname(path)
    rows, rates = [], []
    stack = []                                # (heading level, cited paths) of the enclosing sections
    for start, lines in sections(text):
        level = len(lines[0]) - len(lines[0].lstrip('#')) if lines[0].startswith('#') else 0
        stack = [(lv, q) for lv, q in stack if lv < level] + [(level, [q for line in lines for q in PATH.findall(line)])]
        cited = resolve([q for _, qs in stack for q in qs], folder)   # a subsection inherits its parents' citations
        pool = cited or folder_pool(folder)
        toks = [(start + k, line, t, pc) for k, line in enumerate(lines) for t, pc in tokens(line)]
        if toks:
            rates.append((len(toks), false_match_rate(pool, n=50)))
        for i, line, t, pc in toks:
            rows.append(dict(line=i, number=t, source=traced(t, pool, pc), pool='section' if cited else 'folder',
                             context=line.strip()[:140]))
    fm = sum(n * r for n, r in rates) / max(1, sum(n for n, _ in rates))
    return rows, fm


def audit_headline(path):
    text = open(path, encoding='utf-8').read()
    base = os.path.dirname(path)
    rows, section_paths, top, rates = [], [], '', []
    for i, line in enumerate(text.splitlines(), 1):
        if line.startswith('## '):
            top = line
        if line.startswith('#'):
            section_paths = []
        own = [q for q in PATH.findall(line) if '/' in q or q.endswith(('.json', '.csv', '.md'))]
        if not line.startswith('|'):
            section_paths += own
            continue
        if set(line.replace('|', '').strip()) <= {'-', ' '}:
            continue
        pool = resolve(own or section_paths, base)
        toks = tokens(line)
        if toks and 'RETIRED' not in top:
            rates.append((len(toks), false_match_rate(pool, n=50)))
        for t, pc in toks:
            rows.append(dict(line=i, number=t, source=traced(t, pool, pc) if pool else None,
                             cited=[rel(q) for q in pool][:4], retired='RETIRED' in top, context=line.strip()[:140]))
    fm = sum(n * r for n, r in rates) / max(1, sum(n for n, _ in rates))
    return rows, fm


def citations():
    found = {}
    files = glob.glob(os.path.join(REPO, 'results*', '*.md')) + glob.glob(os.path.join(REPO, 'src', '*.py'))
    pat = re.compile(r'(10\.\d{4,9}/[^\s,;)\]"\']+|arXiv:\s?[\w./-]+\d)')
    for f in files:
        for m in pat.findall(open(f, encoding='utf-8', errors='replace').read()):
            m = m.rstrip('.').replace('arXiv: ', 'arXiv:')
            found.setdefault(m, set()).add(os.path.relpath(f, REPO).replace('\\', '/'))
    named = []                                 # citations given by name without DOI (from src/noise_sourcing.py)
    for f in files:
        for line in open(f, encoding='utf-8', errors='replace'):
            if re.search(r'et al\.|Opt\. Express|Proc\. SPIE|Photonics Research|J\. Lightwave', line):
                named.append((os.path.relpath(f, REPO).replace('\\', '/'), line.strip()[:150]))
    return found, named


def main():
    os.makedirs(OUT, exist_ok=True)
    head, head_fm = audit_headline(HEADLINE)
    summaries, fms = {}, {}
    for p in sorted(glob.glob(os.path.join(REPO, 'results*', 'SUMMARY.md'))):
        summaries[rel(p)], fms[rel(p)] = audit_summary(p)
    found, named = citations()
    cites = [dict(id=k, files=sorted(v), status='VERIFIED' if k in VERIFIED else 'UNVERIFIED',
                  note=VERIFIED.get(k, '')) for k, v in sorted(found.items())]
    res = dict(false_match_rate=dict(headline=head_fm, **fms), headline=head, summaries=summaries, citations=cites,
               named_citations=[dict(file=f, line=l) for f, l in named],
               assumptions=[dict(item=a, used_in=b, status=c) for a, b, c in ASSUMPTIONS])
    with open(os.path.join(OUT, 'audit.json'), 'w', encoding='utf-8') as f:
        json.dump(res, f, indent=1, ensure_ascii=False)

    L = ['<!-- generated by src/provenance_audit.py; do not edit -->', '# Provenance audit', '']
    cur = [r for r in head if not r['retired']]
    bad = [r for r in cur if not r['source']]
    L.append(f'## HEADLINE_NUMBERS.md: {len(cur)} numbers in current/conditional rows, {len(cur) - len(bad)} '
             f'traced to the cited file, **{len(bad)} untraced**')
    L.append('(Numbers in the RETIRED table are not audited: they are listed as retired, with the file '
             'where they appeared.) Chance-match rate of the cited files (share of random 3-digit numbers '
             f'that would "trace"): {head_fm:.0%}.\n')
    if bad:
        L += ['| line | number | cited files | row |', '|---|---|---|---|']
        L += [f"| {r['line']} | {r['number']} | {', '.join(r['cited']) or 'none'} | {r['context'].replace('|', '/')} |"
              for r in bad]
    for name, rows in summaries.items():
        miss = [r for r in rows if not r['source']]
        L.append(f'\n## {name}: {len(rows)} numbers, {len(rows) - len(miss)} traced, **{len(miss)} untraced**\n')
        L.append(f'Chance-match rate of the pools: {fms[name]:.0%} (each section is checked against the files it '
                 'cites; a section citing no file uses the folder data/ and logs/).\n')
        if miss:
            L += ['| line | number | context |', '|---|---|---|']
            L += [f"| {r['line']} | {r['number']} | {r['context'].replace('|', '/')} |" for r in miss]
    L.append('\n## Citations (DOI / arXiv identifiers)\n')
    L += ['| identifier | status | where | note |', '|---|---|---|---|']
    L += [f"| {c['id']} | {c['status']} | {', '.join(c['files'][:3])} | {c['note'][:120]} |" for c in cites]
    L.append('\n## Citations given by name only (all UNVERIFIED unless listed above)\n')
    L += [f"- `{f}`: {l}" for f, l in dict.fromkeys(named)]
    L.append('\n## Unsourced or assumed inputs\n')
    L += ['| item | used in | status |', '|---|---|---|']
    L += [f'| {a} | {b} | {c} |' for a, b, c in ASSUMPTIONS]
    with open(os.path.join(OUT, 'AUDIT.md'), 'w', encoding='utf-8') as f:
        f.write('\n'.join(L) + '\n')
    print(f'headline: {len(cur)} audited, {len(bad)} untraced, chance-match rate {head_fm:.0%}')
    for name, rows in summaries.items():
        print(f'{name}: {len(rows)} numbers, {sum(not r["source"] for r in rows)} untraced, '
              f'chance-match rate {fms[name]:.0%}')
    print(f'citations: {len(cites)} identifiers, {sum(c["status"] == "VERIFIED" for c in cites)} verified; '
          f'{len(named)} named-only lines')


if __name__ == '__main__':
    main()
