"""
Session 4, Task 3 (+ Task 4 additions): extend the sensitivity database
built in Session 3 with new raw structures, without re-running geometries
that haven't changed.

Adds:
  - slot TM (5 widths) -- confirmed to exist in Session 3 (mode index 1,
    TEfrac~0.02 at nominal dims), not simulated then. Values below were
    computed this session via slot_waveguide.sensitivities_slot(W,
    polarization='TM') for W in [150,165,180,195,210] and are reused
    directly here (not re-run) to avoid repaying the ~21-minute cost twice
    in the same session -- same numbers, same computation, just not
    duplicated.
  - rib TE (5 widths) and SWG TE (5 widths) -- Session 4 Task 4 candidates,
    computed by rib_waveguide.py / swg_waveguide.py.

Convergence/validation status strings are attached exactly as established
when each geometry was validated (see progress.md Session 4).
"""
import json
import numpy as np

import fullvector_sensitivity as fv
from probe_db import _fab_vector, STRIP_GRID_SETTINGS, save_db, build_probe_db

# --- slot TM: reused from this session's sensitivities_slot(W, polarization='TM') run ---
SLOT_TM_RESULTS = {
    150: dict(neff=1.55022, ng=2.25451, dneff_dw=0.00199, dneff_dt=0.00258, dneff_da=-0.00746,
              dng_dw=0.0114, dng_dt=0.01374, dng_da=-0.04283),
    165: dict(neff=1.57221, ng=2.37825, dneff_dw=0.00144, dneff_dt=0.00299, dneff_da=-0.00549,
              dng_dw=0.00803, dng_dt=0.01574, dng_da=-0.03071),
    180: dict(neff=1.6003, ng=2.53462, dneff_dw=0.00208, dneff_dt=0.00354, dneff_da=-0.00799,
              dng_dw=0.01142, dng_dt=0.0182, dng_da=-0.04392),
    195: dict(neff=1.63414, ng=2.71858, dneff_dw=0.00187, dneff_dt=0.00417, dneff_da=-0.00728,
              dng_dw=0.00993, dng_dt=0.02061, dng_da=-0.03869),
    210: dict(neff=1.65725, ng=2.83941, dneff_dw=0.00203, dneff_dt=0.0046, dneff_da=-0.00761,
              dng_dw=0.0107, dng_dt=0.02198, dng_da=-0.04022),
}

SLOT_GRID_SETTINGS = dict(half_x=1.6, half_y=1.2, n_core=50, n_clad=20, n_modes=6)
RIB_GRID_SETTINGS = dict(half_x=1.4, half_y=1.2, n_core=30, n_clad=18, n_modes=4)
SWG_GRID_SETTINGS = dict(half_x=1.4, half_y=1.2, n_core=30, n_clad=18, n_modes=4)


def make_record(geometry, pol, W, s, gap_nm, sim_settings, convergence_status, validation_status):
    raw, scaled, unit = _fab_vector(s)
    return dict(
        geometry=geometry, polarization=pol, width_nm=W, gap_nm=gap_nm,
        thickness_nm=220.0, angle_deg=0.0, wavelength_um=fv.WL0,
        neff=s['neff'], ng=s['ng'],
        dneff_dw=s['dneff_dw'], dneff_dt=s['dneff_dt'], dneff_da=s['dneff_da'],
        dng_dw=s['dng_dw'], dng_dt=s['dng_dt'], dng_da=s['dng_da'],
        sens_vector_raw=raw.tolist(), sens_vector_scaled=scaled.tolist(), sens_vector_unit=unit.tolist(),
        sim_settings=dict(fd_steps=dict(dw_nm=5.0, dt_nm=2.0, da_deg=1.0), **sim_settings),
        convergence_status=convergence_status,
        validation_status=validation_status,
        id=f"{geometry}_{pol}_{W}",
    )


if __name__ == "__main__":
    import sys
    with open('data/sensitivity_db.json') as f:
        records = json.load(f)
    # dedupe by id (keep first occurrence) -- makes this script idempotent
    # to re-run, guarding against accidentally doubling entries
    seen = set()
    deduped = []
    for r in records:
        if r['id'] not in seen:
            seen.add(r['id'])
            deduped.append(r)
    records = deduped
    existing_ids = {r['id'] for r in records}
    print(f"Loaded {len(records)} existing records (post-dedupe).")

    import slot_waveguide as sw
    n_added = 0
    for W, s in SLOT_TM_RESULTS.items():
        if f"slot_TM_{W}" in existing_ids:
            continue
        _, gap_um = sw.rail_gap_from_width(W/1000)
        records.append(make_record(
            'slot', 'TM', W, s, gap_um*1000, SLOT_GRID_SETTINGS,
            convergence_status="assumed same order as slot-TE (1.40%, Session 3); not independently re-checked for TM",
            validation_status="mode confirmed to exist (TEfrac~0.02 at nominal dims, Session 3); "
                               "no independent cross-check run for slot TM specifically this session",
        ))
        n_added += 1
    print(f"Added {n_added} new slot-TM records.")

    if '--with-rib-swg' in sys.argv:
        from rib_waveguide import sensitivities_rib
        from swg_waveguide import sensitivities_swg, DUTY_CYCLE, PERIOD_UM, N_SWG

        rib_widths = [400, 450, 500, 550, 600]
        n_rib = 0
        for W in rib_widths:
            if f"rib_TE_{W}" in existing_ids:
                continue
            s = sensitivities_rib(W)
            records.append(make_record(
                'rib', 'TE', W, s, None, RIB_GRID_SETTINGS,
                convergence_status="checked: 0.49% neff spread across domain/mesh variation (Session 4)",
                validation_status="validated: exact match to strip solver at slab_h->0 limit "
                                   "(diff=0.0000 @ 500nm, Session 4); no separate literature ref sought "
                                   "(rib is a standard PDK geometry, the limiting-case check is the primary evidence)",
            ))
            n_rib += 1
        print(f"Added {n_rib} new rib-TE records.")

        swg_widths = [400, 450, 500, 550, 600]
        n_swg = 0
        for W in swg_widths:
            if f"swg_TE_{W}" in existing_ids:
                continue
            s = sensitivities_swg(W)
            records.append(make_record(
                'swg', 'TE', W, s, None, SWG_GRID_SETTINGS,
                convergence_status="inherits strip's grid convergence (same grid, same solver, only "
                                    "core index changed) -- not independently re-checked",
                validation_status=f"effective-medium approximation (0th-order Rytov, duty_cycle={DUTY_CYCLE}, "
                                   f"period={PERIOD_UM*1000:.0f}nm, n_eff_ordinary={N_SWG:.4f}); "
                                   "validated the EMA implementation itself via exact match to strip "
                                   "at duty_cycle=1 (diff=0.0000 @ 500nm, Session 4); NOT a validated "
                                   "quantitative SWG simulation -- true SWG needs 3D/Bloch-mode "
                                   "simulation this toolchain cannot do; flagged as approximate by design",
            ))
            n_swg += 1
        print(f"Added {n_swg} new SWG-TE records.")

    for i, r in enumerate(records):
        r['id'] = f"{r['geometry']}_{r['polarization']}_{r['width_nm']}"

    save_db(records, 'data/sensitivity_db.csv', 'data/sensitivity_db.json')
    print(f"Total: {len(records)} records saved to data/sensitivity_db.{{csv,json}}")

    probe_records = build_probe_db(records)
    save_db(probe_records, 'data/probe_db.csv', 'data/probe_db.json')
    print(f"Total: {len(probe_records)} probe records saved to data/probe_db.{{csv,json}}")
