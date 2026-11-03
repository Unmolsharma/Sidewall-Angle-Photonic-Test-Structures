"""Fig. 1: (a) strip and (b) rib cross-sections with the TE mode intensity and the wall tilt; (c) the first-order
wall-perturbation weight along the right wall vs height above the pivot, with its centroid y_c.

Runs in the femwell virtual environment (it re-solves the two nominal modes with the production settings):
    C:/Users/unmol/.venvs/femwell/Scripts/python.exe paper/scripts/fig1_sections.py
The weight and centroid follow src/crosscheck_fields.py exactly; the script checks that its centroids equal the
stored ones in results_final/data/crosscheck_fields.json (to 1e-9 um) and that R_closed matches.
"""
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection
from matplotlib.patches import Polygon as MPoly

import style

sys.path.insert(0, os.path.join(style.REPO, 'src'))
import crosscheck_fields as cf          # noqa: E402
import femwell_batch as fb              # noqa: E402
from femwell.maxwell.waveguide import compute_modes   # noqa: E402

UNITS = [('strip', 'TE', 500, None, 'strip_TE_500'), ('rib', 'TE', 500, 90, 'rib_TE_500_slab90')]
TILT_SHOWN = 8.0                       # deg, exaggerated for visibility in (a, b) only


def mode(geom, W, S):
    w, t = W / 1000, 0.22
    slab = (S or 90) / 1000
    basis0, eps = fb.build(geom, w, t, 0.0, slab, 'oxide', 'r12')
    guess = {'strip': 2.4, 'rib': 2.5}[geom]
    m = fb.pick(compute_modes(basis0, eps, wavelength=fb.WL0, num_modes=2, order=2, n_guess=guess), 'TE')
    return basis0, m, w, t, slab


def weight_profile(basis0, m, w, t, y0):
    (et, et_b), (ez, ez_b) = m.basis.split(m.E)
    Et, Ez = et_b.interpolate(et).value, ez_b.interpolate(ez).value
    X, dx = et_b.global_coordinates().value, et_b.dx
    core = np.zeros(m.basis.mesh.nelements, bool)
    core[basis0.get_dofs(elements='core0')] = True
    sel = core[:, None] & (X[0] >= w / 2 - cf.DELTA) & (X[0] <= w / 2) & (X[1] >= y0) & (X[1] <= t)
    e1, e2 = fb.N_SI ** 2, fb.N_OX ** 2
    wt = (e1 - e2) * (np.abs(Et[1]) ** 2 + np.abs(Ez) ** 2) + (e1 - e2) * (e1 / e2) * np.abs(Et[0]) ** 2
    ww, yy = (wt * dx)[sel], X[1][sel]
    yc = float(np.sum(ww * yy) / np.sum(ww))
    edges = np.linspace(y0, t, 11)
    num, _ = np.histogram(yy, bins=edges, weights=ww)            # weight per bin ...
    den, _ = np.histogram(yy, bins=edges, weights=dx[sel])       # ... per unit area of the sampled band
    prof = num / den
    return yc, 0.5 * (edges[1:] + edges[:-1]), prof / prof.max()


def intensity(m):
    (et, et_b), _ = m.basis.split(m.E)
    v = et_b.interpolate(et).value
    I = (np.abs(v[0]) ** 2 + np.abs(v[1]) ** 2).mean(axis=1)
    return I / I.max()


def outline(ax, geom, w, t, slab, a, **kw):
    if geom == 'strip':
        y0, hw0 = 0.0, w / 2
    else:
        y0, hw0 = slab, w / 2
    hwt = hw0 - np.tan(np.radians(a)) * (t - y0)
    pts = [(-hw0, y0), (hw0, y0), (hwt, t), (-hwt, t)]
    ax.add_patch(MPoly(pts, closed=True, fill=False, **kw))
    if geom == 'rib':
        ax.plot([-0.9, -hw0], [slab, slab], **{k: v for k, v in kw.items() if k in ('color', 'lw', 'ls')})
        ax.plot([hw0, 0.9], [slab, slab], **{k: v for k, v in kw.items() if k in ('color', 'lw', 'ls')})
        ax.plot([-0.9, 0.9], [0, 0], color=kw.get('color'), lw=kw.get('lw', 0.8))


def main():
    style.setup()
    stored = style.load('results_final/data/crosscheck_fields.json')
    xc = style.load('results_final/data/crosschecks.json')['closed_form_R']['units']
    fig = plt.figure(figsize=(style.DOUBLE, 1.95))
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 0.62])
    axes = [fig.add_subplot(gs[0, i]) for i in range(3)]
    prof_ax = axes[2]
    for k, (geom, pol, W, S, uid) in enumerate(UNITS):
        basis0, m, w, t, slab = mode(geom, W, S)
        y0 = 0.0 if geom == 'strip' else slab
        yc, yb, prof = weight_profile(basis0, m, w, t, y0)
        ref = stored[uid]['with_Ez']
        assert abs(yc - ref['y_c_um']) < 1e-9, (uid, yc, ref['y_c_um'])
        Rc = (yc - y0) / ((t - y0) / 2)
        assert abs(Rc - ref['R_closed']) < 1e-9
        ax = axes[k]
        mesh = m.basis.mesh
        tris = mesh.p[:, mesh.t].T                      # (nel, 3, 2)
        pc = PolyCollection(tris, array=intensity(m), cmap='cividis', edgecolors='none', antialiased=False, rasterized=True)
        ax.add_collection(pc)
        outline(ax, geom, w, t, slab, 0.0, color='w', lw=0.8)
        outline(ax, geom, w, t, slab, TILT_SHOWN, color='w', lw=0.8, ls=(0, (2, 1.5)))
        ax.plot(w / 2, y0, 'o', color='w', ms=2.5)
        ax.set_xlim(-0.62, 0.62)
        ax.set_ylim(-0.18, 0.40)
        ax.set_aspect('equal')
        ax.set_xlabel(r'$x$ (µm)')
        if k == 0:
            ax.set_ylabel(r'$y$ (µm)')
        else:
            ax.set_yticklabels([])
        name = '(a) strip, 500 nm' if geom == 'strip' else '(b) rib, 500 nm, 90 nm slab'
        ax.text(-0.59, 0.375, name, color='w', fontsize=7.5, va='top', fontweight='bold')
        Rf = xc[uid]['R_femwell_neff']
        ax.text(0.59, -0.155, rf'$R={Rf:.3f}$', color='w', fontsize=7.5, ha='right', va='bottom')
        c = style.C['strip_TE'] if geom == 'strip' else style.C['rib90']
        hgt = (yb - y0) / (t - y0)
        prof_ax.plot(prof, hgt, '-', color=c, marker='os'[k], ms=3, lw=0.9, label='strip' if k == 0 else 'rib')
        prof_ax.axhline((yc - y0) / (t - y0), color=c, ls='--', lw=0.8)
        prof_ax.text(0.03, (yc - y0) / (t - y0) + (0.03 if k == 0 else -0.03), 'strip $y_c$' if k == 0 else 'rib $y_c$',
                     color=c, fontsize=7, va='bottom' if k == 0 else 'top')
    cb = fig.colorbar(pc, ax=axes[:2], shrink=0.85, pad=0.01, fraction=0.025, aspect=15)
    cb.set_label(r'$|\mathbf{E}_t|^2$ (normalized)')
    prof_ax.set_xlim(0, 1.05)
    prof_ax.set_ylim(0, 1)
    prof_ax.set_xlabel('Wall weight (normalized)')
    prof_ax.set_ylabel(r'$(y-y_\mathrm{pivot})/h$')
    prof_ax.yaxis.set_label_position('right')
    prof_ax.yaxis.tick_right()
    prof_ax.text(0.04, 0.97, '(c)', transform=prof_ax.transAxes, fontsize=7.5, fontweight='bold', va='top')
    style.save(fig, 'fig1_sections.pdf')


if __name__ == '__main__':
    main()
