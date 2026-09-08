"""
Task 1: reusable sensitivity database.

Every simulated structure (a specific geometry family + polarization +
width, at nominal thickness/angle/wavelength) becomes one RAW record with
its full sensitivity vector, a normalized (fab-tolerance-scaled) direction
vector, convergence/validation status, and simulation settings. Each raw
record is then wrapped by one or more READOUTS (ring, mzi) that turn its
neff/ng sensitivities into noise-weighted Fisher-information rows -- the
readout is a measurement choice layered on top of the physical structure,
not a new physical simulation, so it's kept as a separate, derived table
rather than re-simulating anything.

Two tables are written to data/:
  - sensitivity_db.csv / .json   : one row per (geometry, polarization, width)
                                    raw simulation
  - probe_db.csv / .json         : one row per (geometry, polarization,
                                    width, readout) -- the unit Fisher
                                    analyses (fisher_reanalysis.py) operate on

Normalization convention (documented, not silently assumed): each raw
sensitivity vector (dneff/dw, dneff/dt, dneff/da) is scaled by assumed
typical fabrication tolerances (sigma_w=5nm, sigma_t=2nm, sigma_a=1deg --
the same step sizes already used for the finite-difference derivatives,
for consistency) before unit-normalizing, so that the three axes are
comparable in the atlas despite having different physical units. This is a
modeling choice; its sensitivity is not exhaustively tested here.
"""
import json
import csv
import numpy as np

import fullvector_sensitivity as fv
import slot_waveguide as sw

SIGMA_W_NM  = 5.0    # assumed typical fab tolerance, width, nm
SIGMA_T_NM  = 2.0    # assumed typical fab tolerance, thickness, nm
SIGMA_A_DEG = 1.0    # assumed typical fab tolerance, sidewall angle, deg
FAB_SCALE = np.array([SIGMA_W_NM, SIGMA_T_NM, SIGMA_A_DEG])

STRIP_WIDTHS = [400, 450, 500, 550, 600]   # nm
SLOT_WIDTHS  = [150, 165, 180, 195, 210]   # nm, rail width (gap co-varies, see slot_waveguide.py)

STRIP_GRID_SETTINGS = dict(half_x=1.4, half_y=1.2, n_core=30, n_clad=18, n_modes=4)
SLOT_GRID_SETTINGS  = dict(half_x=1.6, half_y=1.2, n_core=50, n_clad=20, n_modes=6)

STRIP_VALIDATION = {
    'TE': "validated: matches published neff~2.44-2.45 @ 500x220nm (Session 1); "
          "EIM cross-check +1.3 to +4.0% (expected EIM high-index-contrast bias)",
    'TM': "EIM cross-check +0.35 to +3.1% (same sign/magnitude as TE's EIM bias, "
          "supporting evidence the TM solve is correct); no directly-accessible "
          "published TM0@500x220nm reference found (Session 2, search attempted, "
          "sources paywalled) -- flagged as not independently literature-confirmed",
}
SLOT_VALIDATION = {
    'TE': "large-gap decoupling check: converges to independent (already-validated) "
          "strip TE0 solver to within 0.13% at 1000nm gap (Session 3); no "
          "directly-accessible published neff reference found for this exact "
          "geometry (search attempted, sources paywalled) -- flagged as not "
          "independently literature-confirmed",
}
STRIP_CONVERGENCE = "checked: <=1.5% neff spread across domain half-width 1.4-2.0um, " \
                     "core mesh 30-45pts (Session 1/2)"
SLOT_CONVERGENCE  = "checked: 1.40% neff spread across domain half-width 1.6-2.2um, " \
                     "core mesh 50-70pts (Session 3); dneff_dw, dneff_da at rail=210nm " \
                     "re-verified on a finer grid (70x24 vs 50x20) to rule out a " \
                     "discretization artifact for the width>=195nm nonlinear-sensitivity regime"


def _fab_vector(s):
    raw = np.array([s['dneff_dw'], s['dneff_dt'], s['dneff_da']])
    scaled = raw * FAB_SCALE
    norm = np.linalg.norm(scaled)
    unit = (scaled / norm) if norm > 0 else scaled
    return raw, scaled, unit


def build_sensitivity_db():
    """Raw physical records: one per (geometry, polarization, width)."""
    records = []

    for pol in ('TE', 'TM'):
        for W in STRIP_WIDTHS:
            s = fv.sensitivities(W, polarization=pol)  # solve_fn defaults to strip
            raw, scaled, unit = _fab_vector(s)
            records.append(dict(
                geometry='strip', polarization=pol, width_nm=W, gap_nm=None,
                thickness_nm=220.0, angle_deg=0.0, wavelength_um=fv.WL0,
                neff=s['neff'], ng=s['ng'],
                dneff_dw=s['dneff_dw'], dneff_dt=s['dneff_dt'], dneff_da=s['dneff_da'],
                dng_dw=s['dng_dw'], dng_dt=s['dng_dt'], dng_da=s['dng_da'],
                sens_vector_raw=raw.tolist(), sens_vector_scaled=scaled.tolist(),
                sens_vector_unit=unit.tolist(),
                sim_settings=dict(fd_steps=dict(dw_nm=5.0, dt_nm=2.0, da_deg=1.0), **STRIP_GRID_SETTINGS),
                convergence_status=STRIP_CONVERGENCE,
                validation_status=STRIP_VALIDATION[pol],
            ))

    for pol in ('TE',):  # slot TM mode exists (verified) but not simulated this session -- see progress.md
        for W in SLOT_WIDTHS:
            s = sw.sensitivities_slot(W, polarization=pol)
            _, gap_um = sw.rail_gap_from_width(W/1000)
            raw, scaled, unit = _fab_vector(s)
            records.append(dict(
                geometry='slot', polarization=pol, width_nm=W, gap_nm=gap_um*1000,
                thickness_nm=220.0, angle_deg=0.0, wavelength_um=fv.WL0,
                neff=s['neff'], ng=s['ng'],
                dneff_dw=s['dneff_dw'], dneff_dt=s['dneff_dt'], dneff_da=s['dneff_da'],
                dng_dw=s['dng_dw'], dng_dt=s['dng_dt'], dng_da=s['dng_da'],
                sens_vector_raw=raw.tolist(), sens_vector_scaled=scaled.tolist(),
                sens_vector_unit=unit.tolist(),
                sim_settings=dict(fd_steps=dict(dw_nm=5.0, dt_nm=2.0, da_deg=1.0), **SLOT_GRID_SETTINGS),
                convergence_status=SLOT_CONVERGENCE,
                validation_status=SLOT_VALIDATION[pol],
            ))

    for i, r in enumerate(records):
        r['id'] = f"{r['geometry']}_{r['polarization']}_{r['width_nm']}"
    return records


def build_probe_db(sensitivity_records):
    """Derived table: wrap each raw record with ring and mzi readouts,
    producing noise-weighted Fisher rows. Reuses Session 2's noise model
    (SIGMA_NEFF, SIGMA_NG, SIGMA_PHASE) and MZI phase formula unchanged."""
    probes = []
    for r in sensitivity_records:
        s = dict(dneff_dw=r['dneff_dw'], dneff_dt=r['dneff_dt'], dneff_da=r['dneff_da'],
                 dng_dw=r['dng_dw'], dng_dt=r['dng_dt'], dng_da=r['dng_da'])

        ring_rows = [
            (np.array([s['dneff_dw'], s['dneff_dt'], s['dneff_da']]) / fv.SIGMA_NEFF).tolist(),
            (np.array([s['dng_dw'], s['dng_dt'], s['dng_da']]) / fv.SIGMA_NG).tolist(),
        ]
        probes.append(dict(
            probe_id=f"{r['id']}_ring", raw_id=r['id'],
            geometry=r['geometry'], polarization=r['polarization'], width_nm=r['width_nm'],
            readout='ring', n_observables=2, fisher_rows=ring_rows,
            noise_model=dict(sigma_neff=fv.SIGMA_NEFF, sigma_ng=fv.SIGMA_NG),
        ))

        m = fv.mzi_phase_sensitivities(s)
        mzi_rows = [(np.array([m['dphase_dw'], m['dphase_dt'], m['dphase_da']]) / fv.SIGMA_PHASE).tolist()]
        probes.append(dict(
            probe_id=f"{r['id']}_mzi", raw_id=r['id'],
            geometry=r['geometry'], polarization=r['polarization'], width_nm=r['width_nm'],
            readout='mzi', n_observables=1, fisher_rows=mzi_rows,
            noise_model=dict(sigma_phase=fv.SIGMA_PHASE, L_um=fv.L_MZI),
        ))
    return probes


def save_db(records, csv_path, json_path):
    with open(json_path, 'w') as f:
        json.dump(records, f, indent=2)
    # flatten nested fields to JSON strings for the CSV view
    flat_rows = []
    for r in records:
        flat = {}
        for k, v in r.items():
            flat[k] = json.dumps(v) if isinstance(v, (list, dict)) else v
        flat_rows.append(flat)
    fieldnames = list(flat_rows[0].keys()) if flat_rows else []
    with open(csv_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(flat_rows)


if __name__ == "__main__":
    import time
    import os
    t0 = time.time()
    os.makedirs('data', exist_ok=True)

    print("Building sensitivity_db (raw physical records)...")
    sens_records = build_sensitivity_db()
    save_db(sens_records, 'data/sensitivity_db.csv', 'data/sensitivity_db.json')
    print(f"  {len(sens_records)} records -> data/sensitivity_db.{{csv,json}}")

    print("Building probe_db (ring/mzi readout wrappers + Fisher rows)...")
    probe_records = build_probe_db(sens_records)
    save_db(probe_records, 'data/probe_db.csv', 'data/probe_db.json')
    print(f"  {len(probe_records)} records -> data/probe_db.{{csv,json}}")

    print(f"\nElapsed: {time.time()-t0:.1f}s")
