"""
Build paper/LITERATURE.md (annotated bibliography) and paper/refs.bib from the verified records in
results_review/data/literature_records.json plus a few entries verified by hand (book without DOI, software, arXiv
preprints read in full earlier, a patent). Bibliographic fields come from the records; the annotations below are ours
and are written from each record's abstract (or, where none exists, from the title and the text we read).

Run from src/ (main environment):  python make_literature.py
"""
import html
import json
import os

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REC = os.path.join(REPO, 'results_review', 'data', 'literature_records.json')
LIT = os.path.join(REPO, 'paper', 'LITERATURE.md')
BIB = os.path.join(REPO, 'paper', 'refs.bib')

MANUAL = {
    'kay1993': dict(type='book', authors=['Steven M. Kay'], title='Fundamentals of Statistical Signal Processing, Volume I: Estimation Theory',
                    venue='Prentice Hall, Upper Saddle River, NJ', year=1993, doi=None, isbn='0133457117',
                    status='VERIFIED (library catalogue and publisher listing, ISBN 0133457117)', topic='design'),
    'janz2023': dict(type='preprint', authors=['S. Janz', 'S. Dedyulin', 'D.-X. Xu', 'M. Vachon', 'S. Wang', 'R. Cheriton', 'J. Weber'],
                     title='Measurement Accuracy in Silicon Photonic Ring Resonator Thermometers: Identifying and Mitigating Intrinsic Impairments',
                     venue='arXiv:2306.15643', year=2023, doi='10.48550/arXiv.2306.15643',
                     status='VERIFIED (full text read; quotes in results_final/CITATIONS.md)', topic='noise'),
    'leviton2008': dict(type='preprint', authors=['Douglas B. Leviton', 'Bradley J. Frey'],
                        title='Temperature-dependent absolute refractive index measurements of synthetic fused silica',
                        venue='arXiv:0805.0091', year=2008, doi='10.48550/arXiv.0805.0091',
                        status='VERIFIED (full text read; quotes in results_final/CITATIONS.md)', topic='noise'),
    'nist2519a': dict(type='report', authors=['Sarah L. Gilbert', 'William C. Swann', 'Chih-Ming Wang'],
                      title='Hydrogen cyanide H13C14N absorption reference for 1530 nm to 1565 nm wavelength calibration - SRM 2519a',
                      venue='NIST Special Publication 260-137, Gaithersburg, MD', year=2005, doi='10.6028/NIST.SP.260-137',
                      status='VERIFIED (NIST certificate and SP 260-137 read; quotes in results_final/CITATIONS.md)', topic='noise'),
    'femwell': dict(type='software', authors=['Helge Gehring'], title='femwell: mode solver for photonic and electric waveguides based on FEM, version 0.1.12',
                    venue='https://github.com/HelgeGehring/femwell', year=None, doi=None, status='VERIFIED (installed package metadata)', topic='solver'),
    'empy': dict(type='software', authors=[], title='EMpy (ElectroMagneticPython), version 2.2.3', venue='https://github.com/lbolla/EMpy',
                 year=None, doi=None, status='VERIFIED (installed package metadata)', topic='solver'),
}

# registry records upgraded after the full text was read (quotes in results_final/CITATIONS.md)
FULLTEXT = {
    'xing2018': dict(page='1008-1020', abstract_source=None,
                     status='VERIFIED (full text read, pp. 1008-1020; quotes with pages in results_final/CITATIONS.md)'),
}

# (summary, cite for, relation to our work)
A = {
    'xing2018': ('Extracts the linewidth and thickness of SOI strip waveguides from n_eff and n_g fitted to a low-order (m ~ 15) and a '
                 'high-order (m ~ 150) MZI, through a third-order polynomial map from (w, t) to (n_eff, n_g), with total errors of 0.37 nm '
                 'and 0.26 nm on an imec MPW die. The sidewall angle is fixed at 85 deg, taken from the foundry handbook, and the authors '
                 'state that the error from the unknown wall shape cannot be compensated, so some extracted parameters are only relative.',
                 'state of the art of optical geometry extraction; the angle held at a handbook value; measured fitting errors '
                 '(n_eff 1.1e-5, n_g 1.0e-3); local width variation up to 3 nm between nearby devices',
                 'Extracts width and thickness with the angle fixed at an assumed 85 deg; we ask whether the angle itself can be '
                 'extracted. Its measured local width variation is the scale against which our co-location requirement must be read.'),
    'liu2025': ('Extracts the dimensions of SOI strip waveguides (imec iSiPP50G) from the n_eff and n_g of a single high-order MZI, with a '
                'mapping model from geometry to indices, and discusses the feasibility and the extraction errors.',
                'single-MZI extraction of n_eff and n_g, and its error analysis',
                'Again width and thickness; motivates our question of which readouts add independent information (we find the MZI phase is parallel to the ring n_eff).'),
    'butt2024': ('Extracts core thickness, width offset and core and cladding indices of silicon nitride waveguides with unbalanced MZIs '
                 'across a 300 mm wafer, with relative standard deviations of 0.3% (2.6 nm) for the width offset and 0.5% (1.1 nm) for the '
                 'thickness, and verifies the maps by predicting other devices.',
                 'wafer-scale multi-parameter extraction and its use for prediction',
                 'Closest in spirit to our parameter set (geometry plus indices); four parameters, no sidewall angle, SiN rather than SOI.'),
    'dwivedi2015': ('Measures n_eff and its thermo-optic coefficient over the C-band with three MZIs, accounting for fabrication variability, '
                    'with wafer-scale maps for 450, 600 and 800 nm wide waveguides.',
                    'MZI-based index and dispersion extraction',
                    'Uses several interferometers whose shared geometry is assumed; our Monte Carlo and variation analysis test that kind of assumption.'),
    'atsumi2012': ('For ring biosensors read in TE and TM, determines the fabricated waveguide size from group indices obtained from the FSR, '
                   'to +-2 nm in width and +-1 nm in height, and shows that ignoring cross-section deviations biases the retrieved quantities.',
                   'group-index (FSR) based size determination; bias from unmodelled geometry',
                   'Same mechanism as our "cost of ignoring the angle": an unmodelled geometric parameter biases what is extracted.'),
    'zortman2010': ('Quantifies process-induced resonance variation of microdisks across dies, wafers and lots, and introduces a TE/TM method '
                    'to extract thickness and width variation; resonances are held within 1 THz across a wafer and 105 GHz within a die.',
                    'wafer- and die-scale resonance spread; TE/TM extraction of width and thickness',
                    'Gives the scale against which a 711 pm/deg angle error can be compared; extraction again without the angle.'),
    'zhang2024': ('Extracts n_eff and n_g of SOI rib waveguides (450 x 200 nm, 70 nm slab) from two racetrack rings of different perimeter, '
                  'with ~1e-3 uncertainty, and estimates the cross-section, which agrees with STEM.',
                  'rib-waveguide index extraction with two racetracks; cross-section check against STEM',
                  'A rib with a slab, as in our proposal, but the angle is not extracted.'),
    'lu2017': ('Extracts waveguide width and height with sub-nanometre accuracy from racetrack resonators, maps their wafer statistics, and '
               'feeds spatially correlated variations into a layout-aware Monte Carlo circuit simulation.',
               'racetrack-based width/height extraction; correlated variability models',
               'One of the closest: extraction from resonators assuming the structures share the extracted parameters. We show the angle is invisible to such strips and quantify how the shared-parameter assumption limits an angle extraction.'),
    'nrc2026': ('Patent application: determines waveguide geometry from multiple consecutive resonances of one ring of known cavity length by '
                'constrained regression. Its background states that methods using several devices must assume they share the same width '
                'and thickness variations, which "is recognized as often not satisfied", and it treats a non-negligible sidewall angle as '
                'an input to the index model.',
                'single-ring geometry extraction; the shared-variation assumption stated as a known weakness',
                'Closest to our finding on local variation: it avoids multi-device fitting for the same reason our variation analysis exposes. It does not extract the angle.'),
    'xing2023': ('Extracts a granular map of line width and thickness on a silicon photonics wafer and separates layout- and '
                 'location-dependent systematic variation from random variation with a hierarchical model; systematic contributions dominate '
                 'and die-scale line-width variation correlates with pattern density.',
                 'spatial structure of width/thickness variation; pattern-density dependence',
                 'Supports our concern that structures of different width or density see different fabrication offsets.'),
    'bogaerts2018': ('Reviews the silicon photonic circuit design flow and its links to electronic design automation, including variability '
                     'analysis and compact models.', 'context: variability-aware design', 'Background for why process parameters matter to circuit yield.'),
    'selvaraja2010': ('Reports linewidth uniformity of silicon photonic devices made with 193 and 248 nm lithography, assessed with rings, MZIs '
                      'and AWGs: better than 0.6 nm within chip and 2 nm chip-to-chip.',
                      'measured within-chip and chip-to-chip linewidth control',
                      'Sets a realistic scale for per-structure width differences (sub-nm within a chip), far above the ~0.04 nm our two-width set tolerates.'),
    'robinson2025': ('Statistics of more than 4,600 microdisks, racetracks and microrings on 300 mm wafers, used to evaluate DWDM link '
                     'initialization; microrings vary most.', 'recent wafer-scale resonator statistics', 'Context for the size of resonance variation.'),
    'chrostowski2014': ('A study of 371 identical resonators on a 16 x 9 mm chip shows that the wavelength mismatch between devices grows '
                        'linearly with their separation.', 'distance-dependent (correlated) local variation',
                        'Directly relevant to co-locating the test structures: the closer they are, the better the shared-parameter model.'),
    'jayatilleka2021': ('Post-fabrication trimming of silicon rings with a Ge-implanted section and on-chip annealing, setting resonances across '
                        'a 300 mm wafer to within +-32 pm of target.', 'scale of post-trimming resonance accuracy',
                        'A reference point for how small a resonance error matters in practice.'),
    'mdpi2025': ('Wafer-scale study of ring coupling and loss from 2242 resonators on 59 dies of a 200 mm wafer, with intra-wafer '
                 'dependence of the extracted properties.', 'wafer-scale ring test-structure statistics', 'Context for test-structure based process monitoring.'),
    'aim2023': ('Passive performance statistics of many AIM Photonics microdisk resonators, including insertion loss, FSR, Q and '
                'temperature dependence.', 'foundry device statistics', 'Context.'),
    'chrostowski2015': ('Textbook on silicon photonic design, simulation, fabrication and test.', 'general background on SOI waveguides and test structures', 'Background.'),
    'bogaerts2012': ('Review of silicon microring resonators: theory, sensitivity to fabrication and perturbations, and applications.',
                     'ring resonator basics (resonance, FSR, group index)', 'Background for the ring observables we use.'),
    'orji2018': ('Review of dimensional metrology for integrated circuits (SEM, AFM, scatterometry, TEM): advantages, limitations and open '
                 'problems; its summary table lists cross-section TEM as destructive and scatterometry as model-dependent with non-unique '
                 'inverse solutions.', 'limitations of physical dimensional metrology', 'Positions optical test structures as a complement.'),
    'raymond2005': ('Review of scatterometry for high-volume silicon manufacturing, describing the forward (signature acquisition) and '
                    'inverse (parameter extraction) problems.', 'scatterometry is model-based inverse metrology', 'Our problem is an inverse problem of the same kind, posed with guided-mode observables.'),
    'mdpi2021': ('Studies mode conversion between TM0 and higher-order TE modes in tapered SOI ridge waveguides caused by vertical asymmetry, '
                 'including the effect of angled sidewalls and asymmetric cladding.', 'optical effect of sidewall angle (mode hybridization)',
                 'Shows the angle matters through vertical asymmetry, the same symmetry that hides it from strips in our analysis.'),
    'dai2015': ('Analyses mode hybridization in tapered SOI nanowires with angled sidewalls: the TM fundamental mode converts to '
                'higher-order TE modes along a taper, with a conversion ratio that can approach 100% for an adiabatic taper.',
                'the optical effect of a non-vertical sidewall (vertical asymmetry, mode hybridization)',
                'Shows that the angle acts through the vertical asymmetry it creates; a symmetric strip has none, which is why strips cannot see it to first order.'),
    'yoo2022': ('Measures sidewall roughness of vertical and undercut structures with 3D-AFM, with die-level reproducibility better than 2% '
                'relative standard deviation over 13 dies.',
                'physical sidewall metrology (3D-AFM) and its throughput',
                'The physical metrology our optical structures would complement: precise but slow and local.'),
    'biberman2012': ('Silicon rings with intrinsic Q of 2.2e7 (2.7 dB/m), with losses limited by bend loss.', 'loss/Q context', 'Context for ring readout precision.'),
    'morichetti2010': ('Direct observation of backscattering from sidewall roughness in silicon wires, its scaling with waveguide parameters, and '
                       'the role of roughness in polarization rotation.', 'sidewall effects on loss and polarization', 'Context: sidewall geometry beyond angle (roughness) is a separate limitation of our ideal-wall model.'),
    'dai2011': ('Polarization splitter-rotator using an adiabatic taper and an asymmetric directional coupler, based on TM0-TE1 conversion in '
                'a taper.', 'polarization conversion through mode hybridization', 'Context for angle-sensitive devices.'),
    'padmaraju2014': ('Review of thermal sensitivity of silicon microring devices and approaches to stabilize them.',
                      'temperature sensitivity of ring measurements', 'Supports our temperature-drift analysis.'),
    'kay1993': ('Textbook on estimation theory, including the Cramer-Rao lower bound and maximum-likelihood and least-squares estimation.',
                'the CRLB and Fisher information', 'Our analysis tool.'),
    'klein2025': ('Uses the Fisher information matrix to choose wavelength and incidence angles that maximize EUV scatterometry sensitivity on a '
                  'damascene sample and break correlations between parameters.', 'Fisher-information experiment design in dimensional metrology',
                  'Closest methodological precedent: FIM-based design of the measurement; we design the test structures themselves.'),
    'raue2009': ('Profile-likelihood approach to detect structural and practical non-identifiability in partially observed models, with '
                 'confidence intervals and experiment planning.', 'structural vs practical identifiability', 'Our strip result is a structural non-identifiability; the local-variation result is practical.'),
    'bellman1970': ('Introduces structural identifiability of model parameters.', 'definition of structural identifiability', 'Terminology.'),
    'pukelsheim2006': ('Monograph on optimal design of experiments.', 'optimal design criteria', 'Background for the exhaustive design search.'),
    'silver2007': ('SEMATECH study of the limits of optical critical-dimension scatterometry: model comparison and sensitivity across process stacks and '
                   'measurement modalities.', 'sensitivity limits of OCD metrology', 'Precedent for simulation-based limits of an optical metrology method.'),
    'alassaad2007': ('Models random and offset errors in inverse scatterometry and propagates them to profile-parameter estimates.',
                     'systematic (offset) errors in inverse metrology', 'Analogous to our calibration-offset and bias analyses.'),
    'henn2012': ('Maximum-likelihood estimation for scatterometry that estimates the error model from the data and removes systematic errors seen '
                 'with least squares.', 'error-model estimation in inverse metrology', 'Background for treating noise parameters as unknowns.'),
    'bouchet2021': ('Maximum-information states for coherent scattering measurements (title-level; abstract not available in the registries).',
                    'Fisher information in optical measurement', 'Background.'),
    'johnson2002': ('Correct first-order perturbation theory for shifted dielectric interfaces, using surface integrals of the unperturbed fields.',
                    'the boundary-perturbation formula behind R', 'Our closed-form R uses this result.'),
    'kottke2008': ('First-order perturbation theory for shifted anisotropic interfaces, leading to a subpixel smoothing scheme.',
                   'perturbation theory and subpixel smoothing', 'Explains why grid-based derivatives need careful averaging.'),
    'farjadpour2006': ('Subpixel smoothing of the dielectric function for FDTD, designed from perturbation theory, with quadratic convergence.',
                       'accuracy of discretized interfaces', 'Context for the EMpy derivative errors.'),
    'fallahkhair2008': ('Full-vector finite-difference mode solver for anisotropic waveguides.', 'the finite-difference solver family used by EMpy', 'Methods comparison.'),
    'koshiba1992': ('Edge-element finite-element formulation for microwave and optical waveguides that avoids spurious modes.',
                    'vector FEM for waveguide modes', 'Background for the femwell formulation.'),
    'gustafsson2020': ('scikit-fem, a Python finite-element assembly library.', 'software used', 'Methods.'),
    'geuzaine2009': ('Gmsh, an open-source finite-element mesh generator.', 'software used', 'Methods.'),
    'femwell': ('FEM mode solver for photonic waveguides (software).', 'software used', 'Methods.'),
    'empy': ('Python electromagnetics package with a finite-difference mode solver (software).', 'software of the earlier comparison', 'Methods comparison.'),
    'nist2519a': ('NIST standard reference gas cell for the 1530-1565 nm band: 54 H13C14N line centres certified with expanded (k = 2) '
                  'uncertainties of 0.04-0.24 pm, with suggested calibration procedures for accuracies below 3 pm and from 3 to 30 pm.',
                  'accuracy attainable for an absolute wavelength reference',
                  'Sets the scale of the absolute reference that the two-width set needs; our reference sweep (0.1, 1 and 5 pm) brackets the certified lines and the two suggested procedures.'),
    'janz2023': ('Accuracy of silicon ring thermometers: resonance wavelengths located to +-0.1 pm with a calibrated laser in a stable bath; '
                 'changes of 1 pm or less routinely measured.', 'wavelength precision for the noise model', 'Source of our 1 pm and 0.1 pm levels.'),
    'frey2006': ('Absolute refractive index and dn/dT of silicon and germanium from cryogenic to room temperature.', 'n_Si, dn/dT of Si', 'Material inputs.'),
    'leviton2008': ('Absolute refractive index and dn/dT of Corning 7980 fused silica.', 'n and dn/dT of SiO2', 'Material inputs.'),
    'jiao2025': ('Book chapter on silicon waveguides and their temperature dependence, with effective-index-method values for strip and rib.',
                 'EIM n_eff of the 500 nm rib (2.602)', 'Validation only.'),
}
CLOSEST = ['xing2018', 'lu2017', 'butt2024', 'nrc2026', 'zhang2024']
TOPICS = [('extraction', 'Optical extraction of waveguide geometry'), ('variability', 'Fabrication variability and its effect on devices'),
          ('sidewall', 'Sidewall angle: measurement and optical effects'), ('design', 'Experiment design, Fisher information and identifiability'),
          ('solver', 'Mode solvers, perturbation theory and derivative accuracy'), ('noise', 'Measurement noise and material inputs')]


def cite(r):
    au = r.get('authors') or []
    a = ', '.join(au[:6]) + (' et al.' if len(au) > 6 else '')
    parts = [a, f"\"{html.unescape(r.get('title', ''))}\"", html.unescape(r.get('venue') or '')]
    for k, lab in (('volume', 'vol. '), ('issue', 'no. '), ('page', '')):
        if r.get(k):
            parts.append(f'{lab}{r[k]}')
    if r.get('year'):
        parts.append(str(r['year']))
    s = ', '.join(p for p in parts if p)
    if r.get('doi'):
        s += f", doi:{r['doi']}"
    if r.get('url'):
        s += f", {r['url']}"
    return s


def tex(s):
    s = html.unescape(s or '')
    for a, b in (('&', r'\&'), ('%', r'\%'), ('_', r'\_'), ('#', r'\#')):
        s = s.replace(a, b)
    return s


def bib(key, r):
    au = ' and '.join(tex(x) for x in (r.get('authors') or []))
    t = r.get('type') or ''
    f = dict(author=au, title='{' + tex(r.get('title')) + '}', year=r.get('year'))
    if t in ('journal-article',):
        kind = 'article'
        f.update(journal=tex(r.get('venue')), volume=r.get('volume'), number=r.get('issue'), pages=(r.get('page') or '').replace('-', '--'))
    elif t in ('proceedings-article',):
        kind = 'inproceedings'
        f.update(booktitle=tex(r.get('venue')), volume=r.get('volume'), pages=(r.get('page') or '').replace('-', '--'))
    elif t in ('book-chapter',):
        kind = 'incollection'
        f.update(booktitle=tex(r.get('venue')), publisher=tex(r.get('publisher')))
    elif t in ('book', 'monograph', 'edited-book'):
        kind = 'book'
        f.update(publisher=tex(r.get('publisher') or r.get('venue')))
    elif t == 'preprint' or r.get('source') == 'arxiv':
        kind = 'article'
        f.update(journal='arXiv preprint ' + tex(r.get('venue')))
    elif t == 'software':
        kind = 'misc'
        f.update(howpublished=r'\url{' + r['venue'] + '}')
    elif r.get('source') == 'patent':
        kind = 'misc'
        f.update(howpublished=tex(r.get('venue')), note=r'\url{' + r['url'] + '}')
    else:
        kind = 'misc'
        f.update(howpublished=tex(r.get('venue')))
    if r.get('doi'):
        f['doi'] = r['doi']
    if r.get('isbn'):
        f['isbn'] = r['isbn']
    body = ',\n'.join(f'  {k:9s}= {{{v}}}' if not str(v).startswith('{') else f'  {k:9s}= {v}' for k, v in f.items() if v not in (None, ''))
    return f'@{kind}{{{key},\n{body}\n}}\n'


def main():
    R = json.load(open(REC, encoding='utf-8'))
    for k, v in MANUAL.items():
        R[k] = v
    for k, v in FULLTEXT.items():
        R[k] = {**R[k], **v}
    ok = {k: v for k, v in R.items() if str(v.get('status', '')).startswith('VERIFIED') and k in A}
    L = ['<!-- generated by src/make_literature.py from results_review/data/literature_records.json; do not edit -->',
         '# Annotated bibliography', '',
         f'{len(ok)} verified entries. Each was checked against its registry record (Crossref, arXiv, OpenAlex for the abstract), a '
         'publisher or patent page, or (marked) the full text; summaries are written from the abstract or the text read. Entries that '
         'could not be verified are listed at the end and are not cited in the paper.', '',
         '## The five papers closest to our contribution', '']
    for k in CLOSEST:
        if k in ok:
            L.append(f'- **{k}**: {cite(ok[k])}. *How ours differs:* {A[k][2]}')
    L += ['', 'All five extract width, thickness (and, in Butt et al., material indices) from rings or interferometers and hold the '
          'sidewall angle at an assumed value or leave it out. None asks whether the angle is identifiable, which structures make it '
          'identifiable, or how the shared-parameter assumption limits it; those are the questions of this paper.', '']
    for t, title in TOPICS:
        L += [f'## {title}', '']
        for k, v in ok.items():
            if v.get('topic') == t:
                s, c, rel = A[k]
                L += [f'### {k}', '', f'- **Citation:** {cite(v)}', f'- **Verified:** {v["status"]}' +
                      (f" (abstract from {v.get('abstract_source')})" if v.get('abstract_source') else ''),
                      f'- **Summary:** {s}', f'- **Cite for:** {c}', f'- **Relation to our work:** {rel}', '']
    bad = {k: v for k, v in R.items() if not str(v.get('status', '')).startswith('VERIFIED')}
    L += ['## Not verified (not cited)', ''] + [f"- {k}: {v.get('status')}" for k, v in bad.items()]
    open(LIT, 'w', encoding='utf-8').write('\n'.join(L) + '\n')
    with open(BIB, 'w', encoding='utf-8') as fh:
        fh.write('% generated by src/make_literature.py from verified records; do not edit by hand\n\n')
        for k, v in ok.items():
            fh.write(bib(k, v) + '\n')
    print(len(ok), 'verified entries;', len(bad), 'not verified:', list(bad))


if __name__ == '__main__':
    main()
