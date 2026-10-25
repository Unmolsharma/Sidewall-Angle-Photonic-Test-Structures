"""
Etch-depth nuisance study, part 1 of 2: rib sensitivities in the
(width w, thickness t, sidewall angle a, etch depth d) parameterization.

WHY THIS EXISTS
The existing rib model (rib_waveguide.py) fixes the slab at SLAB_H=90nm and,
when it perturbs thickness t, holds the SLAB fixed -- so etch depth silently
moves with t and is never itself an unknown. Here the partial-etch depth d is
an explicit fabrication parameter:

    slab thickness = t - d        (t = total SOI Si thickness at the rib)

so perturbing t at fixed d raises the slab and the rib top together (a
thicker wafer, same timed etch), and perturbing d at fixed t trades rib
height for slab height. Strip waveguides do not see d at all.

Nothing in the existing modules is modified. The geometry/eps averaging is
rib_waveguide.epsfunc_factory_rib (called with slab = t - d) and the mode
solve/TE classification is fullvector_sensitivity.solve_mode_from_eps, both
unchanged. n_g uses the same formula as fullvector_sensitivity.neff_ng,
n_g = n_eff - lambda * dn_eff/dlambda (central difference, dlam = 5nm); this
is asserted to agree bit-for-bit with the existing function (check_ng_method).

Units: w, t, d in nm; a in degrees of deviation from vertical (existing
convention: 0 = vertical sidewall, positive = top narrower). SWA = 90 - a,
so d/dSWA = -d/da; precisions (std devs) are identical in magnitude.

CLI (all outputs go to results_etch_depth/data/, never to data/):
    python etch_depth_sims.py ngcheck   # n_g method identity check
    python etch_depth_sims.py mesh      # mesh/domain convergence of derivatives
    python etch_depth_sims.py steps [--mesh M0]   # h/2, h, 2h for every rib
    python etch_depth_sims.py dlam      # wavelength-step stability of n_g derivs
"""
import os
# one BLAS thread per worker process; parallelism comes from the process pool
os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('MKL_NUM_THREADS', '1')

import sys
import json
import time
from functools import partial
from multiprocessing import Pool

import numpy as np

import fullvector_sensitivity as fv
from rib_waveguide import epsfunc_factory_rib, make_rib_grid

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT_DATA = os.path.join(REPO, 'results_etch_depth', 'data')

PARAMS = ('w', 't', 'a', 'd')
# w, t, a steps are the existing project's (fullvector_sensitivity.sensitivities
# defaults). d has no precedent; 2nm is chosen to match the t step since both
# are vertical dimensions resolved by the same y-grid.
H0 = dict(w=5.0, t=2.0, a=1.0, d=2.0)
STEP_FACTORS = (0.5, 1.0, 2.0)

T_NOM = 220.0        # nm, existing nominal
A_NOM = 0.0          # deg from vertical, existing nominal (SWA = 90 deg)
DLAM_UM = 0.005      # existing n_g wavelength step

RIB_WIDTHS = [400, 450, 500, 550, 600]          # existing rib widths
RIB_SLABS = [50, 70, 90, 110, 130, 150]         # existing slab sweep range (rib_slab_sweep.py)

# Grids. M0 is the production grid used for every existing strip/rib number.
MESHES = {
    'M0': dict(half_x=1.4, half_y=1.2, n_core=30, n_clad=18),
    'M1': dict(half_x=1.4, half_y=1.2, n_core=45, n_clad=18),
    'M2': dict(half_x=1.4, half_y=1.2, n_core=60, n_clad=24),
    'M3': dict(half_x=1.4, half_y=1.2, n_core=90, n_clad=36),
    'D1': dict(half_x=2.4, half_y=1.8, n_core=30, n_clad=30),   # larger domain, M0 core density
}
MESH_CHECK_STRUCTS = [(400, 50), (400, 150), (500, 90), (600, 50), (600, 150)]


def nominal(width_nm, slab_nm):
    return dict(w=float(width_nm), t=T_NOM, a=A_NOM, d=T_NOM - slab_nm)


def solve_neff(p, wl, pol, X, Y):
    """n_eff of the requested-polarization fundamental rib mode at parameter
    point p = {w, t, a, d}. Slab thickness is t - d."""
    slab_um = (p['t'] - p['d']) / 1000
    ep = epsfunc_factory_rib(p['w'] / 1000, p['t'] / 1000, p['a'], slab_um)
    return fv.solve_mode_from_eps(ep, wl, pol, 4, X, Y).neff.real


def neff_ng(p, pol, X, Y, dlam=DLAM_UM):
    """Same formula as fullvector_sensitivity.neff_ng:
    n_g = n_eff - lambda * (n(lambda+dlam) - n(lambda-dlam)) / (2 dlam).
    (Material dispersion is not modeled -- constant N_SI, N_OX -- exactly as
    in the existing code, so this is the waveguide-dispersion group index.)"""
    n0 = solve_neff(p, fv.WL0, pol, X, Y)
    n_plus = solve_neff(p, fv.WL0 + dlam, pol, X, Y)
    n_minus = solve_neff(p, fv.WL0 - dlam, pol, X, Y)
    return n0, n0 - fv.WL0 * (n_plus - n_minus) / (2 * dlam)


def shifted(p, k, delta):
    q = dict(p)
    q[k] += delta
    return q


def structure_job(width_nm, slab_nm, mesh='M0', pol='TE', factors=STEP_FACTORS,
                  params=PARAMS, dlam=DLAM_UM):
    """Central-difference derivatives of (n_eff, n_g) w.r.t. each parameter, at
    each step factor (step = factor * H0[k]). Returns a JSON-able dict."""
    X, Y = make_rib_grid(**MESHES[mesh])
    p0 = nominal(width_nm, slab_nm)
    t0 = time.time()
    n0, g0 = neff_ng(p0, pol, X, Y, dlam)
    derivs = {}
    for k in params:
        for f in factors:
            h = f * H0[k]
            n_p, g_p = neff_ng(shifted(p0, k, +h), pol, X, Y, dlam)
            n_m, g_m = neff_ng(shifted(p0, k, -h), pol, X, Y, dlam)
            derivs[f'{k}|{f}'] = dict(h=h, dneff=(n_p - n_m) / (2 * h), dng=(g_p - g_m) / (2 * h))
    return dict(width_nm=width_nm, slab_nm=slab_nm, d_nm=p0['d'], mesh=mesh, pol=pol,
                dlam_um=dlam, neff=n0, ng=g0, derivs=derivs, seconds=time.time() - t0)


def _run_job(kwargs):
    return structure_job(**kwargs)


def run_pool(jobs, n_proc=14):
    with Pool(min(n_proc, len(jobs))) as pool:
        return pool.map(_run_job, jobs, chunksize=1)


def save(obj, name):
    os.makedirs(OUT_DATA, exist_ok=True)
    path = os.path.join(OUT_DATA, name)
    with open(path, 'w') as f:
        json.dump(obj, f, indent=1)
    print(f'wrote {path}')


def check_ng_method():
    """n_g here must be the existing method, not a reimplementation that
    happens to be close: compare against fullvector_sensitivity.neff_ng with
    an injected solve_fn at the same point. Expect agreement to solver
    repeatability: the ARPACK eigensolve uses a random start vector, so two
    identical calls differ at the ~1e-14 level (measured), hence 1e-10."""
    X, Y = make_rib_grid(**MESHES['M0'])
    p = nominal(500, 90)

    def solve_fn(width_um, thick_um, angle_deg, wl, polarization):
        q = dict(w=width_um * 1000, t=thick_um * 1000, a=angle_deg, d=p['d'])
        return solve_neff(q, wl, polarization, X, Y)

    n_ref, g_ref = fv.neff_ng(p['w'] / 1000, p['t'] / 1000, p['a'], 'TE', solve_fn=solve_fn)
    n_new, g_new = neff_ng(p, 'TE', X, Y)
    print(f'existing fv.neff_ng : neff={n_ref:.10f} ng={g_ref:.10f}')
    print(f'this module         : neff={n_new:.10f} ng={g_new:.10f}')
    diff = max(abs(n_ref - n_new), abs(g_ref - g_new))
    assert diff < 1e-10, f'n_g method differs from existing code by {diff:.2e}'
    print(f'n_g method matches existing code: OK (max abs diff {diff:.1e})')
    return dict(neff_existing=n_ref, ng_existing=g_ref, neff_new=n_new, ng_new=g_new, max_abs_diff=diff)


if __name__ == '__main__':
    cmd = sys.argv[1] if len(sys.argv) > 1 else ''
    t_start = time.time()
    if cmd == 'ngcheck':
        save(check_ng_method(), 'ng_method_check.json')
    elif cmd == 'mesh':
        jobs = [dict(width_nm=w, slab_nm=s, mesh=m, factors=(1.0,))
                for m in MESHES for (w, s) in MESH_CHECK_STRUCTS]
        save(run_pool(jobs), 'mesh_check.json')
    elif cmd == 'steps':
        mesh = sys.argv[sys.argv.index('--mesh') + 1] if '--mesh' in sys.argv else 'M0'
        jobs = [dict(width_nm=w, slab_nm=s, mesh=mesh) for s in RIB_SLABS for w in RIB_WIDTHS]
        save(run_pool(jobs), f'rib_steps_{mesh}.json')
    elif cmd == 'dlam':
        jobs = [dict(width_nm=w, slab_nm=s, mesh='M0', factors=(1.0,), dlam=dl)
                for dl in (0.0025, 0.005, 0.01) for (w, s) in MESH_CHECK_STRUCTS]
        save(run_pool(jobs), 'dlam_check.json')
    else:
        print(__doc__)
        sys.exit(1)
    print(f'elapsed {time.time() - t_start:.0f}s')
