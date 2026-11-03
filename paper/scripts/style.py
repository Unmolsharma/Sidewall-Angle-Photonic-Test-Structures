"""Shared figure style and data access for the paper figures."""
import json
import os

import matplotlib as mpl

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
FIG = os.path.join(REPO, 'paper', 'figures')
SINGLE, DOUBLE = 3.5, 7.16          # column widths (in)

# Okabe-Ito; one colour per structure family, used in every figure
C = dict(strip_TE='#0072B2', strip_TM='#56B4E9', rib50='#D55E00', rib90='#E69F00', rib70='#CC79A7',
         other='#009E73', black='#000000', grey='#7F7F7F')
MARK = dict(strip_TE='o', strip_TM='s', rib50='^', rib90='v', rib70='D', other='P')


def setup():
    mpl.rcParams.update({
        'font.family': 'serif', 'font.serif': ['Times New Roman', 'STIXGeneral', 'DejaVu Serif'],
        'mathtext.fontset': 'stix', 'font.size': 8, 'axes.labelsize': 8, 'axes.titlesize': 8,
        'xtick.labelsize': 7.5, 'ytick.labelsize': 7.5, 'legend.fontsize': 7.5,
        'axes.linewidth': 0.6, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
        'xtick.minor.width': 0.4, 'ytick.minor.width': 0.4, 'lines.linewidth': 1.0, 'lines.markersize': 4,
        'xtick.direction': 'in', 'ytick.direction': 'in', 'xtick.top': True, 'ytick.right': True,
        'savefig.dpi': 300, 'figure.constrained_layout.use': True,
        'figure.constrained_layout.h_pad': 0.02, 'figure.constrained_layout.w_pad': 0.02, 'pdf.fonttype': 42, 'figure.dpi': 150, 'legend.frameon': False})


def load(rel):
    return json.load(open(os.path.join(REPO, rel), encoding='utf-8'))


def panel(ax, letter, x=-0.02, y=1.0):
    ax.text(x, y, f'({letter})', transform=ax.transAxes, ha='right', va='top', fontsize=8, fontweight='bold')


def save(fig, name):
    os.makedirs(FIG, exist_ok=True)
    fig.savefig(os.path.join(FIG, name))          # exact figure size (constrained layout), no cropping
    print('wrote', name)
