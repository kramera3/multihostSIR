# ============================================================================
# multihostSIR – an R package for simulating multi-host pathogen dynamics in 
# ecological communities
# Mayank Gangwar and Andrew Kramer
# ============================================================================
#
# Reproduces the complete simulation study. Plain R throughout: no knitr, no
# rmarkdown, no pandoc. Figures are written as PNG, tables as CSV, and a plain
# text log records everything printed to the console.
#
# Run from the repository root:
#
#   Rscript analysis/reproduce_study.R
#
# or from an R session:
#
#   source("analysis/reproduce_study.R")
#
# All model code lives in R/. This script only sets parameters, runs the study
# and writes results.
#
# Note on the previous version of this analysis: it was an .Rmd with
# `error = TRUE` set globally, which silently swallowed genuine failures. This
# script fails loudly instead.
# ============================================================================

# ---------------------------------------------------------------------------
# 0. Setup
# ---------------------------------------------------------------------------

if (!requireNamespace("multihostSIR", quietly = TRUE)) {
  stop("multihostSIR is not installed. Run:\n",
       "  source(\"install.R\")\n",
       "from the repository root, or:\n",
       "  remotes::install_github(\"kramera3/multihostSIR\")",
       call. = FALSE)
}

library(multihostSIR)
library(ggplot2)

check_dependencies()

start_time <- Sys.time()

# ---------------------------------------------------------------------------
# 1. Study settings
#
# Every input is collected here. Nothing below this block needs editing to
# change the study design.
# ---------------------------------------------------------------------------

# -- Simulation size --------------------------------------------------------
# Set n_iterations to ______ and richness_levels to 3:n for the full study.
# The values below keep a run to a couple of minutes.
n_iterations    <- 200
richness_levels <- 3:8
sim_duration    <- 365          # days

# -- Execution --------------------------------------------------------------
use_parallel <- TRUE
n_cores      <- max(1L, parallel::detectCores() - 1L)
random_seed  <- 42

# -- Output -----------------------------------------------------------------
output_dir <- "analysis/output"
figure_dir <- "analysis/figures"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

# Mirror everything printed below into a log file as well as the console.
# on.exit() only works inside functions, so cleanup happens explicitly at the
# end of the script. If the run stops early, call sink() once at the prompt to
# restore normal console output.
log_path <- file.path(output_dir, "run_log.txt")
log_con <- file(log_path, open = "wt")
sink(log_con, split = TRUE)

# -- Community structure ----------------------------------------------------
spec <- community_spec(
  prey_prob               = 0.8,
  prey_masses             = c(0.02, 0.1, 2, 10, 50),  # e.g. mouse -> rodent -> deer
  predator_masses         = c(1, 5, 25),              # e.g. weasel -> wolf
  predator_density_factor = 0.1,
  competence_log_scale    = TRUE
)

# -- Allometric scaling -----------------------------------------------------
# Worked example, a 0.02 kg host against a 50 kg host:
#   mu    : 0.005 * 0.02^-0.25 = 0.0133 /day    vs 0.0019 /day
#   gamma : 0.100 * 0.02^-0.25 = 0.266  /day    vs 0.0380 /day
#   N     : 500   * 0.02^-0.75 = 9403 ind/km2   vs 26.6 ind/km2
#   alpha : 0.010 * 0.02^-0.25 = 0.0266 /day    vs 0.0038 /day
allo <- allometry(
  a_mort      = 0.005, b_mort      = -0.25,
  a_recovery  = 0.100, b_recovery  = -0.25,
  a_density   = 500,   b_density   = -0.75,
  a_virulence = 0.010, b_virulence = -0.25
)

# -- Trophic transmission routes --------------------------------------------
coeffs <- trophic_coefficients(
  prey_prey    = 0.50,
  pred_pred    = 0.25,
  prey_to_pred = 0.25,
  pred_to_prey = 0.00   # blocked
)

cat("\n=== Study settings ===\n\n")
print(spec)
cat("\n")
print(allo)

# Design note ---------------------------------------------------------------
# Each richness level is split into into every Low / Medium / High
# competence composition and simulated each composition n_iterations times. To
# record the time required to run n_iterations iterations PER species
# community, each richness level now runs n_iterations iterations in which the
# competence composition is sampled at random per iteration. Every iteration
# still records its realized composition, so all downstream summaries remain
# valid.

# ---------------------------------------------------------------------------
# 2. Run the study
# ---------------------------------------------------------------------------

cat("\n=== Running simulation ===\n\n")

results <- run_multihost_sir(
  richness_levels = richness_levels,
  n_iterations    = n_iterations,
  duration        = sim_duration,
  spec            = spec,
  allo            = allo,
  coefficients    = coeffs,
  seed            = random_seed,
  parallel        = use_parallel,
  n_cores         = n_cores,
  progress        = TRUE
)

cat("\n")
print(results)

# ---------------------------------------------------------------------------
# 3. Summary tables
# ---------------------------------------------------------------------------

cat("\n=== Summary by species richness ===\n\n")
res_summary <- summary(results)
print(res_summary)

cat("\n=== Timing per richness level ===\n\n")
print(results$timing, row.names = FALSE)

cx <- timing_complexity(results)
cat(sprintf("\ntime ~ richness^%.2f   (R2 = %.3f)\n",
            cx$exponent_richness, cx$r_squared_richness))
cat(sprintf("time ~ state_dim^%.2f  (R2 = %.3f)\n",
            cx$exponent_state_dim, cx$r_squared_state_dim))
cat(sprintf("Total simulation time: %.1f s\n",
            sum(results$timing$elapsed_sec)))

# ---------------------------------------------------------------------------
# 4. Figures
# ---------------------------------------------------------------------------

cat("\n=== Writing figures ===\n\n")

save_fig <- function(name, plot, width = 7, height = 5, dpi = 300) {
  path <- file.path(figure_dir, paste0(name, ".png"))
  ggsave(path, plot, width = width, height = height, dpi = dpi)
  cat(sprintf("  %s\n", path))
  invisible(path)
}

# -- 4.1 Epidemic outcomes by richness --------------------------------------
save_fig("r0_by_richness",      plot_r0_richness(results))
save_fig("peak_by_richness",    plot_peak_prevalence(results))
save_fig("epidemic_trajectory", plot_trajectories(results, xlim = c(0, 120)))

# -- 4.2 Model assumptions --------------------------------------------------
save_fig("allometric_scaling", plot_allometry(allo), height = 5.5)

set.seed(random_seed)
demo_comm <- build_community(6, spec = spec, allo = allo)
save_fig("transmission_matrix",
         plot_transmission_matrix(demo_comm, show_values = TRUE),
         width = 7.5, height = 5.5)

# -- 4.3 Run-time scaling ---------------------------------------------------
save_fig("timing",         plot_timing(results))
save_fig("timing_loglog",  plot_timing(results, log_scale = TRUE))

# -- 4.4 Competence, density and body size ----------------------------------
# These use the raw results table directly, so any further figure can be built
# the same way.

d <- results$results
threshold_red <- "#B4322A"

scatter_vs_r0 <- function(xvar, title, xlab, log_x = TRUE) {
  # Build a tiny frame with fixed column names rather than relying on tidy
  # evaluation, so this works in a plain script with no rlang attached.
  dd <- data.frame(x = d[[xvar]], R0 = d$R0,
                   richness = factor(d$richness))
  p <- ggplot(dd, aes(x, R0, colour = richness)) +
    geom_point(alpha = 0.5, size = 1.5) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 1) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = threshold_red) +
    scale_colour_viridis_d(option = "D", begin = 0.15, end = 0.9,
                           name = "Richness") +
    scale_y_log10() +
    labs(title = title, x = xlab,
         y = expression(Community ~ R[0] ~ "(log)")) +
    theme_multihost()
  if (log_x) p <- p + scale_x_log10()
  p
}

save_fig("r0_by_competence",
         scatter_vs_r0("avg_competence",
                       expression(bold(R[0] ~ "by average competence")),
                       expression("Mean species competence," ~ bar(beta)),
                       log_x = FALSE))

save_fig("r0_by_density",
         scatter_vs_r0("total_density",
                       expression(bold(R[0] ~ "by host density")),
                       expression("Total host density (ind" ~ km^-2 * ", log)")))

save_fig("r0_by_body_size",
         scatter_vs_r0("avg_body_mass",
                       expression(bold(R[0] ~ "by body size")),
                       "Mean community body mass (kg, log)"))

# -- 4.5 Trophic structure --------------------------------------------------
d$trophic_bin <- cut(d$avg_trophic_level,
                     breaks = c(-Inf, 1.05, 1.15, 1.25, 1.35, Inf),
                     labels = c("1.0", "1.1", "1.2", "1.3", ">1.35"))

p_trophic <- ggplot(d, aes(trophic_bin, R0, fill = trophic_bin)) +
  geom_boxplot(alpha = 0.92, outlier.alpha = 0.15, colour = "grey25") +
  geom_hline(yintercept = 1, linetype = "dashed", colour = threshold_red) +
  scale_fill_viridis_d(option = "D", begin = 0.2, end = 0.85, guide = "none") +
  scale_y_log10() +
  labs(title = expression(bold(R[0] ~ "by trophic level")),
       subtitle = "1.0 = all prey; higher = more predators",
       x = "Mean community trophic level",
       y = expression(Community ~ R[0] ~ "(log)")) +
  theme_multihost()
save_fig("r0_by_trophic_level", p_trophic)

# Which trophic group actually carries more infection at the peak is checked
# against the data rather than assumed.
prey_rate <- mean(d$prey_infected_at_peak / pmax(d$prey_density, 1))
pred_rate <- mean(d$pred_infected_at_peak / pmax(d$predator_density, 1))
cat(sprintf(
  "\nMean attack rate at peak -- prey: %.4f, predators: %.4f  =>  %s carry more.\n",
  prey_rate, pred_rate, if (prey_rate >= pred_rate) "prey" else "predators"))

# -- 4.6 Competence sensitivity ---------------------------------------------
cs <- competence_summary(results)

p_comp <- ggplot(cs, aes(prop_high_competence, mean_R0,
                         colour = factor(richness))) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  scale_colour_viridis_d(option = "D", name = "Richness") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(title = "Competence sensitivity",
       x = "Proportion of high-competence (reservoir) species",
       y = expression("Mean community" ~ R[0])) +
  theme_multihost()
save_fig("competence_sensitivity", p_comp)

# ---------------------------------------------------------------------------
# 5. Export tables
# ---------------------------------------------------------------------------

cat("\n=== Writing tables ===\n\n")

formats <- if (requireNamespace("openxlsx", quietly = TRUE)) {
  c("csv", "xlsx")
} else {
  cat("  (openxlsx not installed; writing CSV only)\n")
  "csv"
}

paths <- export_results(results, dir = output_dir, formats = formats)
cat(paste0("  ", paths, collapse = "\n"), "\n")

# ---------------------------------------------------------------------------
# 6. Session information
# ---------------------------------------------------------------------------

cat("\n=== Session information ===\n\n")
print(sessionInfo())

cat(sprintf("\nTotal elapsed: %s\n", format(Sys.time() - start_time)))
cat(sprintf("Figures : %s\n", normalizePath(figure_dir)))
cat(sprintf("Tables  : %s\n", normalizePath(output_dir)))
cat(sprintf("Log     : %s\n", normalizePath(log_path)))

# Restore console output.
sink()
close(log_con)
