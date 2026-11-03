"""Render figure PDFs to PNG (for looking at them): python preview.py <out_dir> [names...]"""
import os
import sys

import pypdfium2 as pdfium

import style

out = sys.argv[1]
os.makedirs(out, exist_ok=True)
names = sys.argv[2:] or sorted(f[:-4] for f in os.listdir(style.FIG) if f.endswith('.pdf'))
for n in names:
    pdf = pdfium.PdfDocument(os.path.join(style.FIG, n + '.pdf'))
    page = pdf[0]
    w, h = page.get_size()
    page.render(scale=250 / 72).to_pil().save(os.path.join(out, n + '.png'))
    print(n, f'{w / 72:.2f} x {h / 72:.2f} in')
