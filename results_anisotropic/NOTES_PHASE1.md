# Phase 1 notes: anisotropic averaging validation

Source: `src/validate_anisotropic.py`, output `data/validation.json`, log `logs/phase1_validation.log`.

All n_eff checks passed (strip TE 500x220 = 2.4446; rib 500/90 FD 1.57% below our EIM 2.6033;
rib with zero slab = strip to 1.8e-14; n_eff mesh spreads 0.07-0.53%; original and cell-exact
averaging reproduce the previous run to <= 2.8e-13).

Accepted with known exceptions (n_g mesh spread; the n_g check was new in this run and uses the
n_eff thresholds):

| Structure | Original | Cell-exact | Anisotropic | Threshold |
|---|---|---|---|---|
| slot TE 180 n_g | 2.23% | 1.55% | 1.57% | 1.5% |
| rib TE 200 n_g | 1.10% | 0.33% | 1.23% | 1.0% |
| rib TE 250 n_g | 1.23% | 0.12% | 1.07% | 1.0% |

- Slot TE n_g is not mesh-converged to 1.5% with any averaging (anisotropic values 2.287, 2.275,
  2.262, 2.251 over the four meshes; monotone with core refinement).
- For narrow ribs (200-250 nm, not part of the 30-structure database) cell-exact n_g converged
  better than anisotropic. The anisotropic spread comes from core-mesh refinement (rib 200:
  3.152 -> 3.176 -> 3.182), not domain size.
