"""
Verify the literature candidates against Crossref (and arXiv for preprints): title, authors, venue, volume, pages, year
and abstract where the registry has one. Entries given by title are found with a Crossref bibliographic query and
accepted only if the returned title matches (normalized string similarity >= 0.9).

Run from src/ (main environment):  python literature_verify.py
Writes results_review/data/literature_records.json
"""
import difflib
import json
import os
import re
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, 'results_review', 'data', 'literature_records.json')

C = [  # (key, doi or None, title for a title search when the DOI is not known, topic)
    ('xing2018', '10.1364/PRJ.6.001008', None, 'extraction'),
    ('liu2025', '10.1364/OE.558406', None, 'extraction'),
    ('butt2024', '10.1063/5.0173914', None, 'extraction'),
    ('dwivedi2015', '10.1109/JLT.2015.2476603', None, 'extraction'),
    ('atsumi2012', '10.1364/OE.20.026969', None, 'extraction'),
    ('zortman2010', '10.1364/OE.18.023598', None, 'extraction'),
    ('zhang2024', '10.1364/OE.534015', None, 'extraction'),
    ('lu2017', '10.1364/OE.25.009712', None, 'variability'),
    ('xing2023', '10.1021/acsphotonics.2c01194', None, 'variability'),
    ('bogaerts2018', '10.1002/lpor.201700237', None, 'variability'),
    ('selvaraja2010', '10.1109/JSTQE.2009.2026550', None, 'variability'),
    ('robinson2025', '10.1364/OE.577064', None, 'variability'),
    ('chrostowski2014', None, 'Impact of fabrication non-uniformity on chip-scale silicon photonic integrated circuits', 'variability'),
    ('jayatilleka2021', None, 'Post-fabrication trimming of silicon photonic ring resonators at wafer-scale', 'variability'),
    ('mdpi2025', '10.3390/photonics12030234', None, 'variability'),
    ('aim2023', None, 'Manufacturability and performance of microdisk resonators from the AIM Photonics foundry', 'variability'),
    ('chrostowski2015', '10.1017/CBO9781316084168', None, 'variability'),
    ('bogaerts2012', '10.1002/lpor.201100017', None, 'variability'),
    ('orji2018', '10.1038/s41928-018-0150-9', None, 'sidewall'),
    ('raymond2005', '10.1063/1.2062993', None, 'sidewall'),
    ('mdpi2021', '10.3390/app11052366', None, 'sidewall'),
    ('biberman2012', '10.1364/OL.37.004236', None, 'sidewall'),
    ('morichetti2010', '10.1103/PhysRevLett.104.033902', None, 'sidewall'),
    ('dai2011', '10.1364/OE.19.010940', None, 'sidewall'),
    ('padmaraju2014', '10.1515/nanoph-2013-0013', None, 'sidewall'),
    ('dai2015', '10.1364/OE.23.032452', None, 'sidewall'),
    ('yoo2022', '10.1186/s42649-022-00070-5', None, 'sidewall'),
    ('resist2022', 'arXiv:2204.03591', None, 'sidewall'),
    ('kay1993', None, None, 'design'),
    ('klein2025', '10.1117/12.3051113', None, 'design'),
    ('klein2025arxiv', 'arXiv:2504.12133', None, 'design'),
    ('raue2009', '10.1093/bioinformatics/btp358', None, 'design'),
    ('bellman1970', '10.1016/0025-5564(70)90132-X', None, 'design'),
    ('pukelsheim2006', '10.1137/1.9780898719109', None, 'design'),
    ('silver2007', None, 'Fundamental limits of optical critical dimension metrology: a simulation study', 'design'),
    ('alassaad2007', None, 'Error analysis in inverse scatterometry. I. Modeling', 'design'),
    ('henn2012', None, 'A maximum likelihood approach to the inverse problem of scatterometry', 'design'),
    ('bouchet2021', None, 'Maximum information states for coherent scattering measurements', 'design'),
    ('johnson2002', '10.1103/PhysRevE.65.066611', None, 'solver'),
    ('farjadpour2006', '10.1364/OL.31.002972', None, 'solver'),
    ('kottke2008', '10.1103/PhysRevE.77.036611', None, 'solver'),
    ('fallahkhair2008', '10.1109/JLT.2008.923643', None, 'solver'),
    ('koshiba1992', None, 'Simple and efficient finite-element analysis of microwave and optical waveguides', 'solver'),
    ('gustafsson2020', '10.21105/joss.02369', None, 'solver'),
    ('geuzaine2009', '10.1002/nme.2579', None, 'solver'),
    ('janz2023', 'arXiv:2306.15643', None, 'noise'),
    ('frey2006', '10.1117/12.672850', None, 'noise'),
    ('leviton2008', 'arXiv:0805.0091', None, 'noise'),
    ('jiao2025', '10.5772/intechopen.1011647', None, 'noise'),
]


def get(url):
    req = urllib.request.Request(url, headers={'User-Agent': 'literature-verify/1.0'})
    for k in range(4):
        try:
            with urllib.request.urlopen(req, timeout=30) as r:
                return r.read()
        except Exception as e:
            err = e
            time.sleep(2 + 3 * k)
    raise err


def clean(s):
    return re.sub(r'\s+', ' ', re.sub(r'<[^>]+>', ' ', s or '')).strip()


def crossref(doi):
    m = json.loads(get('https://api.crossref.org/works/' + urllib.parse.quote(doi)))['message']
    return dict(source='crossref', doi=m.get('DOI'), title=clean(' '.join(m.get('title', []))),
                authors=[(a.get('given', '') + ' ' + a.get('family', '')).strip() or a.get('name', '') for a in m.get('author', [])],
                venue=clean(' '.join(m.get('container-title', []))), volume=m.get('volume'), issue=m.get('issue'),
                page=m.get('page') or m.get('article-number'), year=(m.get('issued', {}).get('date-parts') or [[None]])[0][0],
                type=m.get('type'), publisher=m.get('publisher'), abstract=clean(m.get('abstract')))


def arxiv(aid):
    x = get('http://export.arxiv.org/api/query?id_list=' + aid.split(':', 1)[1])
    ns = {'a': 'http://www.w3.org/2005/Atom', 'arxiv': 'http://arxiv.org/schemas/atom'}
    e = ET.fromstring(x).find('a:entry', ns)
    jr = e.find('arxiv:journal_ref', ns)
    dd = e.find('arxiv:doi', ns)
    return dict(source='arxiv', doi=dd.text if dd is not None else '10.48550/arXiv.' + aid.split(':', 1)[1], arxiv=aid,
                title=clean(e.find('a:title', ns).text), authors=[clean(a.find('a:name', ns).text) for a in e.findall('a:author', ns)],
                venue=f'arXiv:{aid.split(":", 1)[1]}' + (f' (journal-ref: {jr.text})' if jr is not None else ''),
                year=int(e.find('a:published', ns).text[:4]), abstract=clean(e.find('a:summary', ns).text))


def search(title):
    q = urllib.parse.quote(title)
    items = json.loads(get(f'https://api.crossref.org/works?query.bibliographic={q}&rows=5'))['message']['items']
    best = None
    for it in items:
        t = clean(' '.join(it.get('title', [])))
        sim = difflib.SequenceMatcher(None, t.lower(), title.lower()).ratio()
        if best is None or sim > best[0]:
            best = (sim, it.get('DOI'))
    return best


def main():
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    recs = {}
    for key, ident, title, topic in C:
        try:
            if ident is None and title is None:
                recs[key] = dict(topic=topic, status='NOT IN A REGISTRY (book without DOI; verified separately)')
                continue
            if ident is None:
                sim, doi = search(title)
                if sim < 0.9:
                    recs[key] = dict(topic=topic, status=f'UNVERIFIED (best Crossref title match {sim:.2f})', query=title)
                    print(key, 'UNVERIFIED', round(sim, 2))
                    continue
                ident = doi
            r = arxiv(ident) if ident.startswith('arXiv:') else crossref(ident)
            r.update(topic=topic, status='VERIFIED (registry record)')
            recs[key] = r
            print(key, r['year'], r['title'][:80], '| abstract' if r.get('abstract') else '| no abstract')
        except Exception as e:
            recs[key] = dict(topic=topic, status=f'UNVERIFIED (lookup failed: {e!r})', ident=ident, query=title)
            print(key, 'FAILED', e)
        time.sleep(1.0)
    json.dump(recs, open(OUT, 'w', encoding='utf-8'), indent=1, ensure_ascii=False)


if __name__ == '__main__' and not os.environ.get('FILL_ONLY'):
    main()
    fill_missing()


def openalex_abstract(doi):
    j = json.loads(get('https://api.openalex.org/works/doi:' + urllib.parse.quote(doi)))
    inv = j.get('abstract_inverted_index')
    if not inv:
        return None
    pos = sorted((p, w) for w, ps in inv.items() for p in ps)
    return ' '.join(w for _, w in pos)


def fill_missing():
    """Second pass: abstracts from OpenAlex for registry records without one; retry failed arXiv lookups."""
    recs = json.load(open(OUT, encoding='utf-8'))
    ids = {k: i for k, i, t, tp in C}
    for key, r in recs.items():
        try:
            if r.get('status', '').startswith('UNVERIFIED (lookup failed') and str(ids.get(key, '')).startswith('arXiv:'):
                time.sleep(5)
                topic = r['topic']
                r = arxiv(ids[key])
                r.update(topic=topic, status='VERIFIED (registry record)')
                recs[key] = r
                print(key, 'arXiv ok')
            if r.get('doi') and not r.get('abstract') and r.get('source') == 'crossref':
                a = openalex_abstract(r['doi'])
                if a:
                    r['abstract'], r['abstract_source'] = a, 'OpenAlex'
                    print(key, 'OpenAlex abstract')
                else:
                    print(key, 'no abstract anywhere')
            elif r.get('abstract'):
                r.setdefault('abstract_source', r.get('source'))
        except Exception as e:
            print(key, 'second pass failed', e)
        time.sleep(1.0)
    json.dump(recs, open(OUT, 'w', encoding='utf-8'), indent=1, ensure_ascii=False)


if __name__ == '__main__' and os.environ.get('FILL_ONLY'):
    fill_missing()
