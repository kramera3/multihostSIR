# ============================================================================
# multihostSIR: getting started
# ============================================================================
# Run with:  demo(getting_started, package = "multihostSIR")
# Or:        source(system.file("demo", "getting_started.R",
#                              package = "multihostSIR"))
#
# Plain R throughout. No knitr, no rmarkdown.
# ============================================================================

library(multihostSIR)
library(ggplot2)

set.seed(42)

section <- function(x) cat("\n\n", strrep("=", 70), "\n", x, "\n",
                           strrep("=", 70), "\n\n", sep = "")

# ---------------------------------------------------------------------------
section("0. Are all dependencies in place?")
# ---------------------------------------------------------------------------

check_dependencies()

# ---------------------------------------------------------------------------
section("1. The three building blocks")
# ---------------------------------------------------------------------------

# Every simulation is assembled from three objects, each with sensible
# defaults, and each of which prints itself so you can see what a run will
# actually do before starting it.

allo <- allometry()
print(allo)

spec <- community_spec()
print(spec)

coeffs <- trophic_coefficients()
str(coeffs)

# allometry() holds the body-mass power laws Y = a * M^b. The defaults encode
# quarter-power scaling: physiological rates go as M^-0.25, so large animals
# live longer and recover more slowly, while density goes as M^-0.75 following
# Damuth's law, so large animals are rarer.

print(plot_allometry(allo))

# ---------------------------------------------------------------------------
section("2. One community")
# ---------------------------------------------------------------------------

comm <- build_community(6, spec = spec, allo = allo)
print(comm[, c("species", "type", "mass", "mu", "gamma", "N", "competence")])

# The transmission matrix. Diagonal entries are within-species rates
# beta_ii = competence_i / N_i. Off-diagonal entries are the geometric mean
# sqrt(beta_ii * beta_jj) scaled by a trophic coefficient.
# Rows receive infection; columns are the source.

beta <- transmission_matrix(comm)
print(round(beta, 6))

# Under the defaults, predators contract infection from prey but never
# transmit back, which shows up as a block of zeros:

cat("\nPredator -> prey entries (should all be zero):\n")
print(beta[comm$type == "Prey", comm$type == "Predator", drop = FALSE])

print(plot_transmission_matrix(comm))

# Community R0 is the dominant eigenvalue of the next-generation matrix
# K[i,j] = beta[i,j] * N[i] * D[j], with D[j] = 1 / (gamma_j + mu_j + alpha_j).

cat(sprintf("\nCommunity R0 = %.4f\n", community_R0(comm, beta)))

# ---------------------------------------------------------------------------
section("3. One epidemic")
# ---------------------------------------------------------------------------

sim <- simulate_epidemic(comm, beta = beta, duration = 365)
print(species_metrics(sim))

# The epidemic is seeded in the most competent host by default. Pass an integer
# to seed_species to choose a different one.

cat(sprintf("\nSeeded in species %d (%s, %g kg)\n",
            sim$seed_species,
            comm$type[sim$seed_species],
            comm$mass[sim$seed_species]))

# ---------------------------------------------------------------------------
section("4. Replicate communities across richness levels")
# ---------------------------------------------------------------------------

# run_multihost_sir() is the main entry point. It draws n_iterations
# independent random communities at each richness level and records
# epidemiological outcomes, mean trajectories, and per-level wall-clock time.
#
# n_iterations is kept small here so the demo runs quickly. Use hundreds or
# thousands for real work.

res <- run_multihost_sir(
  richness_levels = 3:7,
  n_iterations    = 1000,
  duration        = 365,
  spec            = spec,
  allo            = allo,
  coefficients    = coeffs,
  seed            = 42
)

print(res)
print(summary(res))

print(plot_r0_richness(res))
print(plot_trajectories(res, xlim = c(0, 120)))
print(plot_peak_prevalence(res))

# ---------------------------------------------------------------------------
section("5. Run-time scaling")
# ---------------------------------------------------------------------------

# Each richness level is timed separately. Fitting log(time) ~ log(richness)
# gives an approximate polynomial order p such that run time grows as
# richness^p. The n x n matrix-vector product in the force-of-infection term
# and the n x n eigendecomposition used for R0 make an exponent between 2 and 3
# the expectation at large n.

cx <- timing_complexity(res)
cat(sprintf("time ~ richness^%.2f   (R2 = %.3f)\n",
            cx$exponent_richness, cx$r_squared_richness))
cat(sprintf("time ~ state_dim^%.2f  (R2 = %.3f)\n",
            cx$exponent_state_dim, cx$r_squared_state_dim))

print(plot_timing(res))
print(plot_timing(res, log_scale = TRUE))

# A short untimed warm-up batch runs before the timed loop so that one-off
# costs, such as JIT compilation and cluster start-up, do not inflate the first
# level and distort the fitted curve.

# ---------------------------------------------------------------------------
section("6. Changing the assumptions")
# ---------------------------------------------------------------------------

# Every default is an argument. Here: predators do transmit back to prey,
# recovery is slower, and predators are less rare.

custom <- run_multihost_sir(
  richness_levels = 3:5,
  n_iterations    = 1000,
  duration        = 120,
  spec = community_spec(predator_density_factor = 0.4),
  allo = allometry(a_recovery = 0.05),
  coefficients = trophic_coefficients(pred_to_prey = 0.1),
  seed = 1
)

print(summary(custom)$by_richness)

# ---------------------------------------------------------------------------
section("7. Exporting results")
# ---------------------------------------------------------------------------

out_dir <- file.path(tempdir(), "multihostSIR_demo")
paths <- export_results(res, dir = out_dir, formats = "csv")
cat("Written:\n")
cat(paste0("  ", paths, collapse = "\n"), "\n")

# Pass formats = c("csv", "xlsx") to also write an Excel workbook. That needs
# the optional openxlsx package.

# ---------------------------------------------------------------------------
section("Done")
# ---------------------------------------------------------------------------

cat("For larger runs, enable parallel execution:\n\n")
cat('  big <- run_multihost_sir(\n')
cat('    richness_levels = 3:30,\n')
cat('    n_iterations    = 1000,\n')
cat('    seed            = 42,\n')
cat('    parallel        = TRUE,\n')
cat('    n_cores         = parallel::detectCores() - 1,\n')
cat('    progress        = TRUE\n')
cat('  )\n\n')
cat("See GUIDE.md in the repository for the full walkthrough.\n")
