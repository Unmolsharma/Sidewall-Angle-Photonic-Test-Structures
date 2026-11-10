"""
Verbatim quotes (with PDF page numbers) for every literature value, written to results_final/data/citations.json.

The PDFs are not in the repository (copyright). To rerun: download each PDF from the URL in `cites`, extract
each page with `pdftotext -layout -f P -l P <name>.pdf <name>_pP.txt` into one folder, and pass that folder:
    python citations_extract.py <folder>
Quotes are located by their opening and closing words, so a mismatch raises instead of inventing text. The
IntechOpen chapter is HTML; its quotes were read from the page by a fetch agent (2026-11-01) and are stored as text.
"""
import json
import os
import re
import sys

C = sys.argv[1] if len(sys.argv) > 1 else 'cites'
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'results_final', 'data', 'citations.json')


def page(f, p):
    return open(os.path.join(C, f'{f}_p{p}.txt'), encoding='utf-8', errors='replace').read()


MISSING = '<page file missing>'


def have(f, p):
    """Sources whose page files are not in the folder keep their entry from the existing citations.json."""
    return os.path.exists(os.path.join(C, f'{f}_p{p}.txt'))


def quote(f, p, start, end):
    """Verbatim text from the first occurrence of `start` through `end` on that page, whitespace collapsed."""
    if not have(f, p):
        return MISSING
    t = re.sub(r'\s+', ' ', page(f, p)).replace('\ufffd', '±')
    i = t.index(start)
    j = t.index(end, i) + len(end)
    return t[i:j]


def row(f, p, label):
    """A table row: the line starting with `label` (whitespace collapsed)."""
    if not have(f, p):
        return MISSING
    for line in page(f, p).splitlines():
        if line.strip().startswith(label):
            return re.sub(r'\s+', ' ', line.strip()).replace('�', '–')
    raise KeyError(label)


def header(f, p, startswith):
    if not have(f, p):
        return MISSING
    for line in page(f, p).splitlines():
        if startswith in line:
            return re.sub(r'\s+', ' ', line.strip()).replace('�', '–')
    raise KeyError(startswith)


cites = {
    'janz2023': dict(
        authors='S. Janz, S. Dedyulin, D.-X. Xu, M. Vachon, S. Wang, R. Cheriton, J. Weber',
        title='Measurement Accuracy in Silicon Photonic Ring Resonator Thermometers: Identifying and Mitigating Intrinsic Impairments',
        venue='arXiv:2306.15643 [physics.optics], v1, submitted 27 Jun 2023', doi='10.48550/arXiv.2306.15643',
        journal_version='UNVERIFIED (none listed on the arXiv page)', url='https://arxiv.org/pdf/2306.15643',
        status='VERIFIED (full text read)',
        quotes=[
            dict(use='sigma_lambda = 0.1 pm (best case)', page=9,
                 text=quote('janz2023', 9, 'so the resonance wavelength uncertainty in Fig. 6', 'used to acquire the spectra.')),
            dict(use='context of the 0.1 pm figure (temperature bath)', page=9,
                 text=quote('janz2023', 9, 'The temperature bath stability', 'bench-top stage used for the Ring A and B')),
            dict(use='sigma_lambda = 1 pm (primary)', page=1,
                 text=quote('janz2023', 1, 'Wavelength changes', 'optical communication component testing')),
            dict(use='resonance temperature coefficient ~80 pm/K', page=1,
                 text=quote('janz2023', 1, 'In silicon waveguides this temperature dependence is approximately 80 pm/K',
                            'operating wavelength.')),
            dict(use='laser run-to-run reproducibility', page=5,
                 text=quote('janz2023', 5, 'reproducibility from run to run was within', '0.1 pm')),
            dict(use='ring length L = 100 um used in the paper', page=1,
                 text=quote('janz2023', 1, 'For a ring with L = 100', 'cavity length'),
                 note='the micro sign of "100 um" is lost in the PDF text extraction'),
        ]),
    'frey2006': dict(
        authors='Bradley J. Frey, Douglas B. Leviton, Timothy J. Madison',
        title='Temperature dependent refractive index of silicon and germanium',
        venue='arXiv:physics/0606168 (to be published in Proc. SPIE 6273, per the arXiv comments)',
        doi='10.48550/arXiv.physics/0606168; related DOI 10.1117/12.672850', url='https://arxiv.org/pdf/physics/0606168',
        status='VERIFIED (full text read; page numbers are PDF pages)',
        quotes=[
            dict(use='n_Si absolute uncertainty at 1.5 um, 295 K (basis of the 1e-4 n_Si prior)', page=3,
                 text=header('frey2006', 3, 'Table 1:') + ' | ' + header('frey2006', 3, 'wavelength 30 K') + ' | ' +
                 row('frey2006', 3, '1.5 microns 8.61E-05')),
            dict(use='n_Si at 1.5 um, 295 K (nominal 3.48 is a rounded assumption)', page=3,
                 text=header('frey2006', 3, 'Table 2:') + ' | ' + header('frey2006', 3, 'wavelength 30 K 40 K') + ' | ' + row('frey2006', 3, '1.5 microns 3.45309')),
            dict(use='dn/dT of Si at 1.5 and 2.0 um, 295 K (DNDT_SI = 1.86e-4 by linear interpolation to 1.55 um)', page=5,
                 text=header('frey2006', 5, 'Table 4:') + ' | ' + header('frey2006', 5, 'wavelength 30 K 40 K') + ' | ' + row('frey2006', 5, '1.5 microns 5.86E-06') + ' | ' +
                 row('frey2006', 5, '2.0 microns 5.39E-06')),
        ]),
    'leviton2008': dict(
        authors='Douglas B. Leviton, Bradley J. Frey',
        title='Temperature-dependent absolute refractive index measurements of synthetic fused silica',
        venue='arXiv:0805.0091 [physics.optics], submitted 1 May 2008', doi='10.48550/arXiv.0805.0091',
        journal_version='UNVERIFIED (none listed)', url='https://arxiv.org/pdf/0805.0091',
        status='VERIFIED (full text read; page numbers are PDF pages)',
        quotes=[
            dict(use='dn/dT of fused silica at 1.5 and 1.6 um, 295 K (DNDT_OX = 8.4e-6 by interpolation to 1.55 um, rounded)',
                 page=7, text=header('leviton2008', 7, 'thermo-optic coefficient, dn/dT, of Corning') + ' | ' + header('leviton2008', 7, 'wavelength      30 K') +
                 ' | ' + row('leviton2008', 7, '1.5 microns   9.92E-07') + ' | ' + row('leviton2008', 7, '1.6 microns   9.31E-07')),
            dict(use='n of fused silica at 1.5 um, 295 K (nominal 1.44 is a rounded assumption; deposited cladding differs)',
                 page=7, text=header('leviton2008', 7, 'absolute refractive index, n, of Corning') + ' | ' + header('leviton2008', 7, 'wavelength     30 K') + ' | ' + row('leviton2008', 7, '1.5 microns   1.44369')),
            dict(use='absolute index uncertainty of bulk fused silica (abstract)', page=1,
                 text=quote('leviton2008', 1, 'with an absolute', 'x 10-5.')),
        ]),
    'johnson2002': dict(
        authors='Steven G. Johnson, M. Ibanescu, M. A. Skorobogatiy, O. Weisberg, J. D. Joannopoulos, Y. Fink',
        title="Perturbation theory for Maxwell's equations with shifting material boundaries",
        venue='Phys. Rev. E 65, 066611 (2002)', doi='10.1103/PhysRevE.65.066611',
        url='https://journals.aps.org/pre/abstract/10.1103/PhysRevE.65.066611; https://math.mit.edu/~stevenj/papers/JohnsonIb02.pdf',
        status='VERIFIED bibliographically (APS page and PDF page 1); used only as the method behind the '
               'field-centroid cross-check, no number taken from it',
        quotes=[dict(use='title page', page=1, text=quote('johnson2002', 1, 'PHYSICAL REVIEW E, VOLUME 65, 066611',
                                                          'shifting material boundaries'))]),
    'jiao_pikus2025': dict(
        authors='Lihong (Heidi) Jiao, Fedor Pikus',
        title='Silicon Waveguides and Temperature-Dependent Performance',
        venue='in: Advanced Optoelectronics and Photonic Technologies - Fundamentals, Devices, and Emerging Applications '
              '[Working Title], eds. D. G. Popescu, M. S. Nazir, IntechOpen (online first, published 26 August 2025)',
        doi='10.5772/intechopen.1011647', url='https://www.intechopen.com/online-first/1224995',
        status='VERIFIED (HTML page read by a fetch agent, Section 2.3 quoted verbatim); values are the authors\' own '
               'EIM calculation, not a measurement',
        quotes=[
            dict(use='published EIM n_eff: rib 2.602, strip 2.494 (validation only, not a headline input)',
                 section='2.3', text='The effective index (neff = 2.602) in the rib waveguide is higher than that in the strip '
                 'waveguide (neff = 2.494) due to reduced confinement in the rib waveguide with identical core width.'),
            dict(use='geometry of the 2.602 value', section='2.3',
                 text='shows the effective index for a rib waveguide with the following parameters: total height (H) = 220 nm, '
                 'slab height (h) = 90 nm, width (W) = 500 nm.'),
            dict(use='material indices of the EIM calculation', section='2.3',
                 text='For these calculations, the material indices used are n1 = 3.48, n2 = n3 = 1.44 at wavelength of 1.55 μm.'),
        ]),
    # Full text supplied by the author (2026-11-06), kept outside git in literature_pdfs/xing2018_photonres.pdf. Pages are
    # PDF pages; journal page = PDF page + 1007 (the article spans pp. 1008-1020). Text pages extracted with
    # `pdftotext -enc UTF-8 -f P -l P` (reading order); xing2018L_p6.txt is the same page with -layout, for Table 3.
    'xing2018': dict(
        authors='Yufei Xing, Jiaxing Dong, Sarvagya Dwivedi, Umar Khan, Wim Bogaerts',
        title='Accurate extraction of fabricated geometry using optical measurement',
        venue='Photonics Research 6(11), 1008-1020 (2018)', doi='10.1364/PRJ.6.001008',
        url='https://doi.org/10.1364/PRJ.6.001008',
        status='VERIFIED (full text read, all 13 pages; page numbers are PDF pages, journal page = PDF page + 1007)',
        quotes=[
            dict(use='sidewall angle assumed from the foundry handbook, not measured (journal p. 1010)', page=3,
                 text=quote('xing2018', 3, 'According to the IMEC technology handbook', 'bottom width of the trapezoid.')),
            dict(use='error from the unknown wall shape is not compensated (journal p. 1010)', page=3,
                 text=quote('xing2018', 3, 'Using a simulated model can introduce a simulation error', 'some parameters will be only relative.')),
            dict(use='85 deg restated for the measured 450 nm waveguides (journal p. 1017)', page=10,
                 text=quote('xing2018', 10, 'According to the technology handbook', 'sidewall angle of 85°.')),
            dict(use='meaning of the quoted fitting errors: 1.96 sigma (journal p. 1013)', page=6,
                 text=quote('xing2018', 6, 'The fitting error we presented', 'mately 95%.')),
            dict(use='Table 3, order-150 MZI: n_g error and propagated w, t errors (nm), 1.96 sigma (journal p. 1013)', page=6,
                 text=row('xing2018L', 6, '150 ')),
            dict(use='simulated measurement noise of the fit study (journal p. 1012)', page=5,
                 text=quote('xing2018', 5, '0.2 dBm to the spectrum', 'typical measurement noise.'),
                 note='a symbol before "0.2" (probably ±) is lost in the PDF text extraction'),
            dict(use='measured fitting errors of n_eff and n_g and their width/thickness equivalents (journal p. 1018)', page=11,
                 text=quote('xing2018', 11, 'The average fitting error for neff is', '0.18 nm in thickness.')),
            dict(use='total extraction error = fitting + mapping (journal p. 1018)', page=11,
                 text=quote('xing2018', 11, 'Errorw', '0.26 nm'),
                 note='the "=" and "+" signs are lost in the PDF text extraction; the equations read '
                      'Error_w = 0.31 + 0.06 = 0.37 nm and Error_t = 0.18 + 0.08 = 0.26 nm (the text above says 0.33 nm)'),
            dict(use='per-die width and thickness spread over 44 MZI pairs (journal p. 1018)', page=11,
                 text=quote('xing2018', 11, 'The standard deviations are 1.26 nm', 'respectively (Table 4).')),
            dict(use='local width and thickness variation between nearby devices (journal p. 1019)', page=12,
                 text=quote('xing2018', 12, 'We did observe a systematic trend', 'local variation range we set in the extraction.')),
            dict(use='arms of one MZI assumed to share n_eff (journal p. 1011)', page=4,
                 text=quote('xing2018', 4, 'since the distance between waveguides', 'device-to-device variation.')),
            dict(use='measurement temperature (journal p. 1018)', page=11,
                 text=quote('xing2018', 11, 'The optical measurement was conducted', 'controlled at 20°C.')),
            dict(use='wavelength calibration in the procedure (journal p. 1019)', page=12,
                 text=quote('xing2018', 12, '1. Calibrate the laser', 'stability of the laser.')),
            dict(use='reference [11], Dwivedi et al. 2015 (journal p. 1020)', page=13,
                 text=quote('xing2018', 13, '11. S. Dwivedi', '4471–4477 (2015).')),
            dict(use='definition of local variation: devices less than a few hundred microns apart (journal p. 1015)', page=8,
                 text=quote('xing2018', 8, 'On the other hand, the local', 'few hundred microns apart).')),
            dict(use='location-dependent thickness difference on the two dies (journal p. 1019)', page=12,
                 text=quote('xing2018', 12, 'dependent difference in thickness', '1.5 nm.')),
        ]),
    # NIST certificate fetched 2026-11-09 from the URL below; PDF pages extracted with `pdftotext -enc UTF-8 -f P -l P`
    'nist2519a': dict(
        authors='Sarah L. Gilbert, William C. Swann, Chih-Ming Wang',
        title='Hydrogen cyanide H13C14N absorption reference for 1530 nm to 1565 nm wavelength calibration - SRM 2519a',
        venue='NIST Special Publication 260-137 (2005); Certificate of Analysis, SRM 2519a', doi='10.6028/NIST.SP.260-137',
        url='https://tsapps.nist.gov/srmext/certificates/archives/2519a.pdf',
        status='VERIFIED (NIST certificate and SP 260-137 read; page numbers are certificate pages)',
        quotes=[
            dict(use='certified line-centre uncertainty and coverage factor', page=1,
                 text=quote('nist2519a', 1, 'The center wavelengths of the lines listed', 'coverage factor k = 2')),
            dict(use='suggested procedure, 3 pm to 30 pm', page=4,
                 text=header('nist2519a', 4, 'Suggested Procedure for Moderate-Accuracy')),
            dict(use='suggested procedure, below 3 pm', page=4,
                 text=header('nist2519a', 4, 'Suggested Procedure for High-Accuracy')),
        ]),
}
unverified = [
    dict(item='published strip n_eff ~2.44-2.45 (500 x 220 nm, TE)', where='src/fullvector_sensitivity.py, src/validate_*.py, '
         'notes/progress.md', status='UNVERIFIED: no source is named anywhere in the repo; not used for any headline number'),
    dict(item='Baehr-Jones et al., Opt. Express 19, 11529 (2011) and "Low-loss silicon slot waveguides and couplers ..." '
         'Opt. Express 19, 26275 (2011): slot/strip loss ratio', where='src/noise_sourcing.py',
         status='UNVERIFIED (not fetched); only used for slot results, which are dropped'),
    dict(item='"Record low-loss hybrid rib/wire waveguides for silicon photonic circuits" (rib/strip loss ratio 0.2x)',
         where='src/noise_sourcing.py', status='UNVERIFIED (no full citation in the repo; not fetched); the 0.2x optimistic '
         'rib noise factor is therefore treated as an ASSUMPTION'),
    dict(item='Dwivedi et al., J. Lightwave Technol. 33(21), 4471 (2015), DOI 10.1109/JLT.2015.2476603', where='searched 2026-10-29',
         status='UNVERIFIED (full text not readable: HTTP 403); no value taken'),
    dict(item='Farjadpour et al., Opt. Lett. 31, 2972 (2006); Kottke et al., PRE 77, 036611 (2008); Fallahkhair & Li (2008)',
         where='src/anisotropic_averaging.py, notes', status='UNVERIFIED (method references of the retired EMpy work; not fetched)'),
]
os.makedirs(os.path.dirname(OUT), exist_ok=True)
old = json.load(open(OUT, encoding='utf-8'))['sources'] if os.path.exists(OUT) else {}
for k, v in cites.items():
    if any(MISSING in (q.get('text') or '') for q in v['quotes']):
        if k not in old:
            raise SystemExit(f'{k}: page files missing from {C} and no earlier entry in {OUT}')
        cites[k] = old[k]
        print(k, 'page files not in the folder: kept from the existing citations.json')
json.dump(dict(sources=cites, unverified=unverified), open(OUT, 'w', encoding='utf-8'), indent=1, ensure_ascii=False)
for k, v in cites.items():
    for q in v['quotes']:
        print(k, q.get('page', q.get('section')), '|', q['text'][:160])
