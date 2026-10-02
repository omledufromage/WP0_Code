# WP0 simulation code

R simulation studies for the WP0 paper (`../Paper/location_scale/`).

## Coordination with the Paper session

Read `../COORDINATION.md` at the start of every session. It covers folder
ownership, how figures and numbers are handed to the paper, the table of
simulation settings the paper quotes, and open requests between the two sides.
Update it when any of those change.

## Layout

- `R/`: shared helpers. `R/utils.R` has `run_simulation()` (runs a scenario grid,
  seed `seed + i` for scenario `i`), `save_results()` and `load_latest_results()`.
- `simulations/NN_<study>.R`: one script per study, copied from
  `simulations/00_template.R`. Run them from the project root, e.g.
  `Rscript simulations/01_<study>.R`.
- `results/<study>/`: raw results as timestamped `.rds` files that also store
  the settings and session info. This folder is git-ignored.
- Final figures for the paper go to `../Paper/Figures/ls_<study>_<name>.pdf`,
  written by the study's script.
