# multihostSIR 0.1.0

First release.

## Building and simulating communities

* `run_multihost_sir()` assembles random host communities across a range of
  species richness levels, simulates an epidemic in each, and returns tidy
  per-community results together with mean trajectories and per-level timing.
* `build_community()`, `transmission_matrix()`, `community_R0()` and
  `simulate_epidemic()` expose the individual steps for use on their own.
* `allometry()`, `community_spec()` and `trophic_coefficients()` hold the
  assumptions behind a run; every argument has a default and can be overridden.
* Community rates are derived from body mass through quarter-power allometric
  scaling; between-species transmission is asymmetric, with a unidirectional
  prey-to-predator pathway; and community R0 is the dominant eigenvalue of the
  next-generation matrix.

## Bringing your own data

* `read_community()` reads a species table from a data frame, CSV or Excel file,
  filling any missing rates from body mass and recording what was supplied.
* `read_transmission_matrix()` reads and validates a transmission matrix.
* `run_fixed_community()` solves a single, fully specified community once and
  deterministically; when a community is complete and a matrix is supplied, no
  random iteration is used.
* Every result reports which values were supplied and which came from the
  built-in defaults.

## Output

* A set of plotting functions covers richness comparisons, macro-driver
  relationships, predator structure, the transmission matrix and the allometric
  relationships; `plot_gallery()` produces them together.
* `export_results()` writes results to CSV and, when openxlsx is installed,
  Excel, for both replicate studies and single deterministic runs.

## Documentation

* Two vignettes, an introduction and a fuller treatment of the model and the
  simulation study, plus a plain-text guide and a runnable demo.
