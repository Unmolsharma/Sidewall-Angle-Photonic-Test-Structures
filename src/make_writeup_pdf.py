"""
Build WRITEUP.pdf from WRITEUP.md (repo root).

Also draws Figure 3 (step-size check, original vs cell-exact averaging) from
results_etch_depth/data/diag_averaging.json, since it isn't in figures/.
Figures 1 and 2 are the existing figures/ PNGs, used as they are.

Handles only the markdown the writeup uses: #/##/### headings, paragraphs,
**bold**, `code`, pipe tables and ![caption](path) images.

Run:  python make_writeup_pdf.py
"""
import os
import re
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from reportlab.lib import colors
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, Image,
                                Table, TableStyle, KeepTogether)

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MD = os.path.join(REPO, 'WRITEUP.md')
PDF = os.path.join(REPO, 'WRITEUP.pdf')
FIG3 = os.path.join(REPO, 'results_etch_depth', 'figures', 'step_check_diag.png')

# Times New Roman from the Windows font folder
FONT_DIR = os.path.join(os.environ.get('WINDIR', r'C:\Windows'), 'Fonts')
for name, fn in (('TNR', 'times.ttf'), ('TNR-Bold', 'timesbd.ttf'),
                 ('TNR-Italic', 'timesi.ttf'), ('TNR-BoldItalic', 'timesbi.ttf')):
    pdfmetrics.registerFont(TTFont(name, os.path.join(FONT_DIR, fn)))
pdfmetrics.registerFontFamily('TNR', normal='TNR', bold='TNR-Bold',
                              italic='TNR-Italic', boldItalic='TNR-BoldItalic')
plt.rcParams['font.family'] = 'Times New Roman'


def make_fig3():
    with open(os.path.join(REPO, 'results_etch_depth', 'data', 'diag_averaging.json')) as f:
        d = json.load(f)
    rows = [('w', 'neff', 'Width\n(n_eff)'), ('t', 'neff', 'Thickness\n(n_eff)'),
            ('a', 'neff', 'Sidewall angle\n(n_eff)'), ('d', 'neff', 'Etch depth\n(n_eff)'),
            ('d', 'ng', 'Etch depth\n(n_g)')]
    orig = [100 * d[f'existing|{k}|{o}']['max_rel_dev'] for k, o, _ in rows]
    fixed = [100 * d[f'cell_exact|{k}|{o}']['max_rel_dev'] for k, o, _ in rows]
    x = range(len(rows))
    fig, ax = plt.subplots(figsize=(7, 3.4), dpi=200)
    ax.bar([i - 0.2 for i in x], orig, 0.38, color='#c0504d', label='Original averaging')
    ax.bar([i + 0.2 for i in x], fixed, 0.38, color='#2f6fad', label='Cell-exact averaging')
    ax.axhline(2, color='#444444', ls='--', lw=1)
    ax.text(len(rows) - 0.5, 2.8, '2% limit', ha='right', fontsize=8, color='#444444')
    ax.set_xticks(list(x))
    ax.set_xticklabels([r[2] for r in rows], fontsize=8)
    ax.set_ylabel('Largest step disagreement (%)', fontsize=9)
    ax.spines[['top', 'right']].set_visible(False)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    os.makedirs(os.path.dirname(FIG3), exist_ok=True)
    fig.savefig(FIG3)
    plt.close(fig)


def inline(text):
    text = text.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')
    text = re.sub(r'\*\*(.+?)\*\*', r'<b>\1</b>', text)
    text = re.sub(r'`(.+?)`', r'<font face="Courier">\1</font>', text)
    return text


def build():
    ss = getSampleStyleSheet()
    body = ParagraphStyle('body', parent=ss['BodyText'], fontName='TNR', fontSize=11.5, leading=15, spaceAfter=7)
    cap = ParagraphStyle('cap', parent=body, fontSize=9, leading=12, textColor=colors.HexColor('#444444'))
    cell = ParagraphStyle('cell', parent=body, fontSize=10, leading=12, spaceAfter=0)
    h = {1: ParagraphStyle('h1', parent=ss['Title'], fontName='TNR-Bold', fontSize=20, spaceAfter=14),
         2: ParagraphStyle('h2', parent=ss['Heading2'], fontName='TNR-Bold', fontSize=15, spaceBefore=12, spaceAfter=6),
         3: ParagraphStyle('h3', parent=ss['Heading3'], fontName='TNR-BoldItalic', fontSize=12.5, spaceBefore=8, spaceAfter=4)}
    width = letter[0] - 2 * inch

    story, para, table = [], [], []

    def flush_para():
        if para:
            story.append(Paragraph(inline(' '.join(para)), body))
            para.clear()

    def flush_table():
        if table:
            rows = [[Paragraph(inline(c.strip()), cell) for c in r.strip().strip('|').split('|')]
                    for r in table if not re.match(r'^\|[\s\-|]+\|$', r.strip())]
            t = Table(rows, colWidths=[width / len(rows[0])] * len(rows[0]), hAlign='LEFT')
            t.setStyle(TableStyle([
                ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#e8edf3')),
                ('LINEBELOW', (0, 0), (-1, 0), 0.8, colors.HexColor('#7a8899')),
                ('LINEBELOW', (0, -1), (-1, -1), 0.5, colors.HexColor('#7a8899')),
                ('VALIGN', (0, 0), (-1, -1), 'TOP'),
                ('TOPPADDING', (0, 0), (-1, -1), 3), ('BOTTOMPADDING', (0, 0), (-1, -1), 3)]))
            story.extend([t, Spacer(1, 9)])
            table.clear()

    with open(MD, encoding='utf-8') as f:
        lines = f.read().splitlines()
    for line in lines:
        img = re.match(r'!\[(.*)\]\((.*)\)', line.strip())
        head = re.match(r'(#{1,3}) (.*)', line)
        if line.strip().startswith('|'):
            flush_para()
            table.append(line)
            continue
        flush_table()
        if not line.strip():
            flush_para()
        elif head:
            flush_para()
            story.append(Paragraph(inline(head.group(2)), h[len(head.group(1))]))
        elif img:
            flush_para()
            path = os.path.join(REPO, img.group(2))
            iw, ih = Image(path).imageWidth, Image(path).imageHeight
            w = width if iw / ih > 1.2 else width * 0.75
            story.append(KeepTogether([Image(path, w, w * ih / iw), Spacer(1, 4),
                                       Paragraph(inline(img.group(1)), cap)]))
            story.append(Spacer(1, 8))
        else:
            para.append(line.strip())
    flush_para()
    flush_table()

    doc = SimpleDocTemplate(PDF, pagesize=letter, leftMargin=inch, rightMargin=inch,
                            topMargin=0.9 * inch, bottomMargin=0.9 * inch,
                            title='Measuring sidewall angle optically with test structures')
    doc.build(story, onLaterPages=lambda c, d: (c.setFont('TNR', 9),
                                                c.drawCentredString(letter[0] / 2, 0.5 * inch, str(d.page))))
    print(f'wrote {PDF}')


if __name__ == '__main__':
    make_fig3()
    build()
