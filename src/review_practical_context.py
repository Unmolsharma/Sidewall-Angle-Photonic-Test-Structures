"""
Literature context for the practical-value comparison (review, step 6): published wafer-scale and die-scale resonance
and linewidth variation, read from the verified abstracts in results_review/data/literature_records.json, with
frequency spreads converted to wavelength by d lambda = lambda^2 d f / c at lambda = 1.55 um. Every number is parsed
from the stored abstract text (regular expression), not typed here.

Run from src/ (main environment):  python review_practical_context.py
Writes results_review/data/practical_context.json
"""
import json
import os
import re

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LIT = os.path.join(REPO, 'results_review', 'data', 'literature_records.json')
OUT = os.path.join(REPO, 'results_review', 'data', 'practical_context.json')
LAM, C0 = 1.55e-6, 299792458.0


def main():
    r = json.load(open(LIT, encoding='utf-8'))
    z = r['zortman2010']['abstract']
    m = re.search(r'within (\d+(?:\.\d+)?) ?THz across a wafer and (\d+(?:\.\d+)?) ?GHz within a single die', z)
    wafer_THz, die_GHz = float(m.group(1)), float(m.group(2))
    j = r['jayatilleka2021']['abstract']
    trim_pm = float(re.search(r'within ±(\d+) pm', j).group(1))
    s = r['selvaraja2010']['abstract']
    within_chip_nm = float(re.search(r'linewidth control better than (\d+(?:\.\d+)?) nm within chip', s).group(1))
    chip_to_chip_nm = float(re.search(r'better than (\d+(?:\.\d+)?) nm chip-to-chip', s).group(1))
    out = dict(
        zortman2010=dict(quote=m.group(0), wafer_THz=wafer_THz, die_GHz=die_GHz,
                         wafer_nm_at_1550=LAM ** 2 * wafer_THz * 1e12 / C0 * 1e9,
                         die_pm_at_1550=LAM ** 2 * die_GHz * 1e9 / C0 * 1e12,
                         note='microdisks, 6-inch SOI 0.35 um process; frequency spread converted with d lambda = lambda^2 d f / c'),
        jayatilleka2021=dict(quote=re.search(r'within ±\d+ pm of a target wavelength of \d+ ?nm', j).group(0), trimmed_to_pm=trim_pm,
                             note='after post-fabrication trimming, 300 mm SOI wafer'),
        selvaraja2010=dict(quote=re.search(r'linewidth control better than .*?chip-to-chip', s).group(0),
                           within_chip_nm=within_chip_nm, chip_to_chip_nm=chip_to_chip_nm))
    json.dump(out, open(OUT, 'w', encoding='utf-8'), indent=1, ensure_ascii=False)
    print(json.dumps(out, indent=1, ensure_ascii=False))


if __name__ == '__main__':
    main()
