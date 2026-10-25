# Sidewall-Angle Photonic Test Structures

Silicon photonic test structures (rings, MZIs, waveguides) can recover waveguide
width and thickness optically, but sidewall angle is much harder to measure.
This project treats test-structure selection as an experiment-design problem:
a full-vector finite-difference mode solver gives the sensitivity of effective
and group index to width, thickness and sidewall angle, and Fisher information
and Cramér-Rao bounds show how well a given set of structures can determine
each parameter. The main finding is that cross-sections that are symmetric top
to bottom (strip, slot, subwavelength grating) respond to a sidewall tilt almost
exactly as they do to a width change, so they cannot separate the two, while a
rib waveguide's slab breaks that symmetry and makes sidewall angle measurable.

**Status:** work in progress. Results are being recomputed with improved
sub-pixel material averaging and cross-checked against an independent solver.

## Repository layout

- `src/` — mode-solver geometry, sensitivity calculations, Fisher / Cramér-Rao
  analyses, validation and plotting scripts
- `data/` — sensitivity databases and analysis outputs from the original runs
- `results/`, `results_corrected/`, `results_etch_depth/` — outputs and reports
  for each stage (original vs. corrected averaging, etch-depth follow-up)
- `figures/` — figures, one folder per averaging method
- `notes/` — project brief and progress notes
- `WRITEUP.md` / `WRITEUP.pdf` — current write-up of methods and findings
