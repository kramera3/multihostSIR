#' A clean theme for multihostSIR figures
#'
#' A light wrapper around [ggplot2::theme_bw()] with a muted grid, a stronger
#' panel border and a bold title. Every plotting function in the package applies
#' it by default.
#'
#' @param base_size Base font size in points.
#' @return A ggplot2 theme object.
#' @examples
#' library(ggplot2)
#' ggplot(mtcars, aes(wt, mpg)) + geom_point() + theme_multihost()
#' @export
theme_multihost <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = "grey88"),
          panel.border = element_rect(colour = "grey45", linewidth = 0.8),
          plot.title = element_text(face = "bold",
                                    size = base_size * 1.2,
                                    margin = margin(b = 4)),
          plot.subtitle = element_text(size = base_size * 0.9,
                                       colour = "grey30",
                                       margin = margin(b = 6)),
          legend.key.height = unit(0.42, "cm"),
          plot.margin = margin(8, 10, 6, 8))
}

# Colours used for threshold lines and highlights.
.mh_red <- "#B4322A"

#' Choose a legible subset of richness levels
#' @noRd
display_levels <- function(levels, max_levels = 8L) {
  levels <- sort(unique(levels))
  if (length(levels) <= max_levels) return(levels)
  sort(unique(round(seq(min(levels), max(levels), length.out = max_levels))))
}

#' @noRd
check_result_object <- function(x) {
  if (!inherits(x, "multihost_sir")) {
    stop("`x` must be created by run_multihost_sir().", call. = FALSE)
  }
  invisible(TRUE)
}

#' Community R0 by species richness
#'
#' Boxplots of community \eqn{R_0} against host species richness, on a log
#' \eqn{y} axis, with the epidemic threshold \eqn{R_0 = 1} marked.
#'
#' @param x A `"multihost_sir"` object from [run_multihost_sir()].
#' @param max_levels Maximum number of richness levels to display. When more
#'   levels were simulated, an evenly spaced subset is shown so the axis stays
#'   readable. All data remain in `x$results`.
#' @return A ggplot object.
#' @examples
#' res <- run_multihost_sir(3:6, n_iterations = 5, duration = 120, seed = 1)
#' plot_r0_richness(res)
#' @export
plot_r0_richness <- function(x, max_levels = 8L) {
  check_result_object(x)
  keep <- display_levels(x$results$richness, max_levels)
  d <- x$results[x$results$richness %in% keep, , drop = FALSE]
  d$richness_f <- factor(d$richness, levels = sort(keep))

  ggplot(d, aes(x = .data$richness_f, y = .data$R0, fill = .data$richness_f)) +
    geom_boxplot(alpha = 0.92, outlier.alpha = 0.15, outlier.size = 0.8,
                 linewidth = 0.5, colour = "grey25") +
    geom_hline(yintercept = 1, linetype = "dashed",
               colour = .mh_red, linewidth = 0.9) +
    scale_fill_viridis_d(option = "D", begin = 0.15, end = 0.9,
                         guide = "none") +
    scale_y_log10(labels = label_number(accuracy = 0.1)) +
    labs(title = expression(bold(R[0] ~ "by species richness")),
         subtitle = "Dashed line marks the epidemic threshold",
         x = "Host species richness",
         y = expression(Community ~ R[0] ~ "(log scale)")) +
    theme_multihost()
}

#' Peak prevalence by species richness
#'
#' @inheritParams plot_r0_richness
#' @return A ggplot object.
#' @examples
#' res <- run_multihost_sir(3:6, n_iterations = 5, duration = 120, seed = 1)
#' plot_peak_prevalence(res)
#' @export
plot_peak_prevalence <- function(x, max_levels = 8L) {
  check_result_object(x)
  keep <- display_levels(x$results$richness, max_levels)
  d <- x$results[x$results$richness %in% keep, , drop = FALSE]
  d$richness_f <- factor(d$richness, levels = sort(keep))

  ggplot(d, aes(x = .data$richness_f, y = .data$peak_prevalence,
                fill = .data$richness_f)) +
    geom_boxplot(alpha = 0.92, outlier.alpha = 0.15, outlier.size = 0.8,
                 linewidth = 0.5, colour = "grey25") +
    scale_fill_viridis_d(option = "D", begin = 0.15, end = 0.9,
                         guide = "none") +
    scale_y_log10(labels = label_number(scale_cut = scales::cut_short_scale())) +
    labs(title = "Peak prevalence by species richness",
         x = "Host species richness",
         y = "Peak infected individuals (log scale)") +
    theme_multihost()
}

#' Mean epidemic trajectories
#'
#' Mean daily community-wide prevalence for each richness level.
#'
#' @inheritParams plot_r0_richness
#' @param xlim Optional length-2 numeric vector limiting the time axis. Defaults
#'   to the full horizon.
#' @return A ggplot object.
#' @examples
#' res <- run_multihost_sir(3:6, n_iterations = 5, duration = 150, seed = 1)
#' plot_trajectories(res, xlim = c(0, 120))
#' @export
plot_trajectories <- function(x, max_levels = 7L, xlim = NULL) {
  check_result_object(x)
  keep <- display_levels(x$trajectories$richness, max_levels)
  d <- x$trajectories[x$trajectories$richness %in% keep, , drop = FALSE]
  d$richness_f <- factor(d$richness, levels = sort(keep))

  p <- ggplot(d, aes(x = .data$time, y = .data$mean_prevalence,
                     colour = .data$richness_f)) +
    geom_line(linewidth = 1.1) +
    scale_colour_viridis_d(option = "D", begin = 0.15, end = 0.9,
                           name = "Richness") +
    labs(title = "Mean epidemic trajectories",
         subtitle = "Community-wide infected individuals, averaged over iterations",
         x = "Day", y = "Infected individuals") +
    theme_multihost()

  if (!is.null(xlim)) {
    if (!is.numeric(xlim) || length(xlim) != 2L) {
      stop("`xlim` must be a length-2 numeric vector.", call. = FALSE)
    }
    p <- p + ggplot2::coord_cartesian(xlim = xlim)
  }
  p
}

#' Computation time by community size
#'
#' Wall-clock time to run all iterations at each richness level, with the fitted
#' log-log complexity exponent reported in the subtitle.
#'
#' @inheritParams plot_r0_richness
#' @param log_scale If `TRUE`, draw both axes on a log scale with a linear fit,
#'   which makes the complexity exponent the slope of the line.
#' @return A ggplot object.
#' @examples
#' res <- run_multihost_sir(3:7, n_iterations = 3, duration = 100, seed = 1)
#' plot_timing(res)
#' @export
plot_timing <- function(x, log_scale = FALSE) {
  check_result_object(x)
  check_flag(log_scale, "log_scale")
  cx <- timing_complexity(x)

  subtitle <- if (is.na(cx$exponent_richness)) {
    sprintf("%d iterations per level", x$settings$n_iterations)
  } else {
    sprintf("%d iterations per level; time ~ richness^%.2f",
            x$settings$n_iterations, cx$exponent_richness)
  }

  p <- ggplot(x$timing, aes(x = .data$richness, y = .data$elapsed_sec)) +
    labs(title = "Computation time by community size",
         subtitle = subtitle,
         x = "Host species richness", y = "Elapsed time (s)") +
    theme_multihost()

  if (log_scale) {
    p <- p + geom_point(colour = "#24450F", size = 2.6) +
      geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                  colour = .mh_red, linewidth = 1) +
      scale_x_log10() + scale_y_log10()
  } else {
    p <- p + geom_line(colour = "#4E8F00", linewidth = 1.2) +
      geom_point(colour = "#24450F", size = 2.6)
  }
  p
}

#' Transmission matrix heat map
#'
#' Visualises the asymmetric between-species transmission matrix. Rows are
#' recipients of infection, columns are sources, so the blocked
#' predator-to-prey corner is immediately visible under the default trophic
#' coefficients.
#'
#' @param community A community data frame from [build_community()].
#' @param beta Optional transmission matrix. Defaults to
#'   `transmission_matrix(community)`.
#' @param show_values If `TRUE` (default), print the actual transmission rate
#'   inside each tile. Turn off for large communities where the numbers crowd.
#' @param relative If `FALSE` (default), colour tiles by their actual value, so
#'   the legend reads in transmission-rate units. Set `TRUE` to colour by value
#'   relative to the matrix maximum instead, which can help when rates span
#'   orders of magnitude.
#' @return A ggplot object.
#' @examples
#' set.seed(42)
#' comm <- build_community(6)
#' plot_transmission_matrix(comm)                    # actual values, labelled
#' plot_transmission_matrix(comm, show_values = FALSE)
#' @export
plot_transmission_matrix <- function(community,
                                     beta = transmission_matrix(community),
                                     show_values = TRUE,
                                     relative = FALSE) {
  check_community(community)
  check_flag(show_values, "show_values")
  check_flag(relative, "relative")

  n <- nrow(community)
  labels <- sprintf("%d: %s\n(%g kg)", seq_len(n),
                    as.character(community$type), community$mass)

  grid <- expand.grid(recipient = seq_len(n), source = seq_len(n))
  grid$value <- as.numeric(beta[cbind(grid$recipient, grid$source)])
  grid$fill_value <- if (relative && max(grid$value) > 0) {
    grid$value / max(grid$value)
  } else {
    grid$value
  }
  grid$recipient_f <- factor(labels[grid$recipient], levels = rev(labels))
  grid$source_f <- factor(labels[grid$source], levels = labels)

  p <- ggplot(grid, aes(x = .data$source_f, y = .data$recipient_f,
                        fill = .data$fill_value)) +
    geom_tile(colour = "white", linewidth = 0.7) +
    scale_fill_gradientn(
      colours = c("#FFFFFF", "#FEE0D2", "#FC9272", "#DE2D26", "#A50F15"),
      name = if (relative) "Relative rate" else expression(beta[ij]),
      labels = if (relative) scales::percent_format(accuracy = 1) else
        scales::label_number()) +
    labs(title = "Transmission matrix",
         subtitle = "Rows receive infection; columns are the source",
         x = "Source of infection", y = "Recipient") +
    theme_multihost() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 30, hjust = 1))

  if (show_values) {
    p <- p + ggplot2::geom_text(aes(label = format(signif(.data$value, 3),
                                                   scientific = FALSE)),
                                size = 3, colour = "grey10")
  }
  p
}

#' Allometric scaling reference curves
#'
#' Draws the four body-mass power laws encoded in an [allometry()] object on
#' log-log axes.
#'
#' @param allo An allometry object from [allometry()].
#' @param mass_range Length-2 numeric vector giving the body-mass range (kg).
#' @return A ggplot object.
#' @examples
#' plot_allometry(allometry())
#' @export
plot_allometry <- function(allo = allometry(), mass_range = c(0.01, 100)) {
  if (!inherits(allo, "multihost_allometry")) {
    stop("`allo` must be created by allometry().", call. = FALSE)
  }
  if (!is.numeric(mass_range) || length(mass_range) != 2L ||
      any(mass_range <= 0) || mass_range[1L] >= mass_range[2L]) {
    stop("`mass_range` must be two increasing positive numbers.", call. = FALSE)
  }

  m <- 10^seq(log10(mass_range[1L]), log10(mass_range[2L]), length.out = 250)
  params <- c("Density (N)", "Recovery rate", "Mortality rate", "Virulence")
  d <- data.frame(
    mass = rep(m, times = 4L),
    parameter = factor(rep(params, each = length(m)), levels = params),
    value = c(allo$a_density   * m^allo$b_density,
              allo$a_recovery  * m^allo$b_recovery,
              allo$a_mort      * m^allo$b_mort,
              allo$a_virulence * m^allo$b_virulence),
    stringsAsFactors = FALSE)

  ggplot(d, aes(x = .data$mass, y = .data$value, colour = .data$parameter)) +
    geom_line(linewidth = 1.2) +
    scale_x_log10(labels = scales::label_log()) +
    scale_y_log10(labels = scales::label_log()) +
    scale_colour_viridis_d(option = "D", begin = 0.1, end = 0.85, name = NULL) +
    labs(title = "Allometric scaling",
         subtitle = sprintf("Rates scale as M^%.2f; density as M^%.2f",
                            allo$b_mort, allo$b_density),
         x = "Body mass (kg, log scale)",
         y = "Parameter value (log scale)") +
    theme_multihost() +
    theme(legend.position = "bottom")
}

#' Plot a multi-host simulation
#'
#' Convenience wrapper dispatching to the individual plotting functions.
#'
#' @param x A `"multihost_sir"` object from [run_multihost_sir()].
#' @param which One of `"r0"`, `"peak"`, `"trajectories"` or `"timing"`.
#' @param ... Passed to the underlying plotting function.
#' @return A ggplot object.
#' @examples
#' res <- run_multihost_sir(3:5, n_iterations = 5, duration = 120, seed = 1)
#' plot(res, which = "r0")
#' plot(res, which = "trajectories")
#' @export
plot.multihost_sir <- function(x, which = c("r0", "peak", "trajectories",
                                            "timing"), ...) {
  which <- match.arg(which)
  switch(which,
         r0           = plot_r0_richness(x, ...),
         peak         = plot_peak_prevalence(x, ...),
         trajectories = plot_trajectories(x, ...),
         timing       = plot_timing(x, ...))
}

#' Epidemic trajectory for a fixed community
#'
#' Plots the community-wide number infected over time for a deterministic run
#' from [run_fixed_community()]. With `by_species = TRUE`, one line is drawn per
#' species instead of the community total.
#'
#' @param x A `"multihost_fixed"` object from [run_fixed_community()].
#' @param by_species If `TRUE`, draw one curve per species; if `FALSE` (default),
#'   draw the community total.
#' @return A ggplot object.
#' @examples
#' species <- data.frame(
#'   type = c("Prey", "Prey", "Predator"),
#'   mass = c(0.02, 50, 5), mu = c(0.013, 0.002, 0.003),
#'   gamma = c(0.27, 0.04, 0.07), alpha = c(0.027, 0.004, 0.007),
#'   N = c(9000, 27, 6), competence = c(0.5, 0.1, 0.2))
#' comm <- read_community(species)
#' beta <- read_transmission_matrix(matrix(c(
#'   0.05, 0.01, 0, 0.01, 0.03, 0, 0.02, 0.02, 0.04), nrow = 3, byrow = TRUE),
#'   community = comm)
#' fit <- run_fixed_community(comm, beta, duration = 150)
#' plot_fixed_trajectory(fit)
#' plot_fixed_trajectory(fit, by_species = TRUE)
#' @seealso [run_fixed_community()], [plot.multihost_fixed()]
#' @export
plot_fixed_trajectory <- function(x, by_species = FALSE) {
  if (!inherits(x, "multihost_fixed")) {
    stop("`x` must be created by run_fixed_community().", call. = FALSE)
  }
  check_flag(by_species, "by_species")

  if (!by_species) {
    d <- x$trajectory
    return(
      ggplot(d, aes(x = .data$time, y = .data$prevalence)) +
        geom_line(linewidth = 1.1, colour = "#2C6E9B") +
        labs(title = "Epidemic trajectory",
             subtitle = sprintf("Fixed community of %d species; R0 = %.2f",
                                nrow(x$community), x$R0),
             x = "Day", y = "Infected individuals") +
        theme_multihost())
  }

  wide <- x$species_curves
  sp <- setdiff(names(wide), "time")
  long <- data.frame(
    time    = rep(wide$time, times = length(sp)),
    species = factor(rep(sp, each = nrow(wide)), levels = sp),
    infected = unlist(wide[sp], use.names = FALSE),
    stringsAsFactors = FALSE)

  ggplot(long, aes(x = .data$time, y = .data$infected,
                   colour = .data$species)) +
    geom_line(linewidth = 1) +
    scale_colour_viridis_d(option = "D", begin = 0.1, end = 0.9,
                           name = "Species") +
    labs(title = "Epidemic trajectory by species",
         subtitle = sprintf("Fixed community; R0 = %.2f", x$R0),
         x = "Day", y = "Infected individuals") +
    theme_multihost()
}

#' Plot a fixed-community run
#'
#' Convenience wrapper for the plots that apply to a single deterministic
#' community from [run_fixed_community()]. The richness-comparison plots used for
#' [run_multihost_sir()] do not apply here, because a fixed run is one community
#' rather than many.
#'
#' @param x A `"multihost_fixed"` object from [run_fixed_community()].
#' @param which One of `"trajectory"` (community total, the default),
#'   `"species"` (one curve per species), or `"matrix"` (the transmission-matrix
#'   heatmap).
#' @param ... Passed to the underlying plotting function.
#' @return A ggplot object.
#' @examples
#' species <- data.frame(
#'   type = c("Prey", "Prey", "Predator"),
#'   mass = c(0.02, 50, 5), mu = c(0.013, 0.002, 0.003),
#'   gamma = c(0.27, 0.04, 0.07), alpha = c(0.027, 0.004, 0.007),
#'   N = c(9000, 27, 6), competence = c(0.5, 0.1, 0.2))
#' comm <- read_community(species)
#' beta <- read_transmission_matrix(matrix(c(
#'   0.05, 0.01, 0, 0.01, 0.03, 0, 0.02, 0.02, 0.04), nrow = 3, byrow = TRUE),
#'   community = comm)
#' fit <- run_fixed_community(comm, beta, duration = 150)
#' plot(fit)
#' plot(fit, which = "species")
#' plot(fit, which = "matrix")
#' @seealso [run_fixed_community()], [plot_fixed_trajectory()]
#' @export
plot.multihost_fixed <- function(x, which = c("trajectory", "species",
                                              "matrix"), ...) {
  which <- match.arg(which)
  switch(which,
         trajectory = plot_fixed_trajectory(x, by_species = FALSE),
         species    = plot_fixed_trajectory(x, by_species = TRUE),
         matrix     = plot_transmission_matrix(x$community, x$beta, ...))
}


## ---------------------------------------------------------------------------
##
## Three families:
##   * two more richness-comparison boxplots (cumulative infected, days to peak)
##   * macro-driver scatterplots (an outcome against a community-average driver)
##   * predator-structure plots (an outcome against predator proportion)
##
## The scatter and predator plots share one internal engine, mh_scatter(), so
## the many named wrappers stay thin and consistent.
## ---------------------------------------------------------------------------

#' Cumulative proportion infected by species richness
#'
#' @param x A `"multihost_sir"` object from [run_multihost_sir()].
#' @param max_levels Largest number of richness levels to show.
#' @return A ggplot object.
#' @examples
#' res <- run_multihost_sir(3:6, n_iterations = 8, duration = 150, seed = 1)
#' plot_cumulative_richness(res)
#' @export
plot_cumulative_richness <- function(x, max_levels = 8L) {
  check_result_object(x)
  keep <- display_levels(x$results$richness, max_levels)
  d <- x$results[x$results$richness %in% keep, , drop = FALSE]
  d$richness_f <- factor(d$richness, levels = sort(keep))

  ggplot(d, aes(x = .data$richness_f, y = .data$cum_infected_prop,
                fill = .data$richness_f)) +
    geom_boxplot(alpha = 0.92, outlier.alpha = 0.15, outlier.size = 0.8,
                 linewidth = 0.5, colour = "grey25") +
    scale_fill_viridis_d(option = "D", begin = 0.15, end = 0.9,
                         guide = "none") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(title = "Cumulative infections by species richness",
         x = "Host species richness",
         y = "Cumulative proportion infected") +
    theme_multihost()
}

#' Days to epidemic peak by species richness
#'
#' @param x A `"multihost_sir"` object from [run_multihost_sir()].
#' @param max_levels Largest number of richness levels to show.
#' @return A ggplot object.
#' @examples
#' res <- run_multihost_sir(3:6, n_iterations = 8, duration = 150, seed = 1)
#' plot_days_to_peak_richness(res)
#' @export
plot_days_to_peak_richness <- function(x, max_levels = 8L) {
  check_result_object(x)
  keep <- display_levels(x$results$richness, max_levels)
  d <- x$results[x$results$richness %in% keep, , drop = FALSE]
  d$richness_f <- factor(d$richness, levels = sort(keep))

  ggplot(d, aes(x = .data$richness_f, y = .data$days_to_peak,
                fill = .data$richness_f)) +
    geom_boxplot(alpha = 0.92, outlier.alpha = 0.15, outlier.size = 0.8,
                 linewidth = 0.5, colour = "grey25") +
    scale_fill_viridis_d(option = "D", begin = 0.15, end = 0.9,
                         guide = "none") +
    labs(title = "Time to epidemic peak by species richness",
         x = "Host species richness",
         y = "Days to peak") +
    theme_multihost()
}

#' Scatter of one outcome against one driver, coloured by richness
#'
#' Shared engine for the macro-driver and predator-structure plots. Each point
#' is one simulated community.
#'
#' @param x A `"multihost_sir"` object.
#' @param xvar,yvar Column names in `x$results`.
#' @param xlab,ylab,title Axis labels and title.
#' @param threshold If not `NULL`, draw a dashed horizontal reference line at
#'   this value (used to mark R0 = 1).
#' @param xpercent,ypercent Format that axis as a percentage.
#' @param smooth Add a smoothed trend line.
#' @param logy Put the y axis on a log scale.
#' @return A ggplot object.
#' @noRd
mh_scatter <- function(x, xvar, yvar, xlab, ylab, title,
                       threshold = NULL, xpercent = FALSE, ypercent = FALSE,
                       smooth = TRUE, logy = FALSE) {
  check_result_object(x)
  d <- x$results
  d$richness_f <- factor(d$richness, levels = sort(unique(d$richness)))

  p <- ggplot(d, aes(x = .data[[xvar]], y = .data[[yvar]])) +
    geom_point(aes(colour = .data$richness_f), alpha = 0.55, size = 1.6)

  if (smooth) {
    p <- p + geom_smooth(method = "loess", se = TRUE, colour = "grey20",
                         linewidth = 0.8, formula = y ~ x)
  }
  if (!is.null(threshold)) {
    p <- p + geom_hline(yintercept = threshold, linetype = "dashed",
                        colour = .mh_red, linewidth = 0.9)
  }

  p <- p + scale_colour_viridis_d(option = "D", begin = 0.15, end = 0.9,
                                  name = "Richness")

  if (xpercent) {
    p <- p + scale_x_continuous(labels = scales::percent_format(accuracy = 1))
  }
  if (ypercent) {
    p <- p + scale_y_continuous(labels = scales::percent_format(accuracy = 1))
  } else if (logy) {
    p <- p + scale_y_log10(labels = label_number(scale_cut = scales::cut_short_scale()))
  }

  p + labs(title = title, x = xlab, y = ylab) + theme_multihost()
}

# --- Macro-driver plots (7) ------------------------------------------------

#' Community R0 against average trophic level
#' @param x A `"multihost_sir"` object from [run_multihost_sir()].
#' @return A ggplot object.
#' @examples
#' res <- run_multihost_sir(3:7, n_iterations = 10, duration = 150, seed = 1)
#' plot_r0_trophic(res)
#' @export
plot_r0_trophic <- function(x) {
  mh_scatter(x, "avg_trophic_level", "R0",
             "Average trophic level", expression(Community ~ R[0]),
             expression(bold(R[0] ~ "vs average trophic level")),
             threshold = 1)
}

#' Community R0 against mean community body mass
#' @inheritParams plot_r0_trophic
#' @return A ggplot object.
#' @examples
#' res <- run_multihost_sir(3:7, n_iterations = 10, duration = 150, seed = 1)
#' plot_r0_body_mass(res)
#' @export
plot_r0_body_mass <- function(x) {
  mh_scatter(x, "avg_body_mass", "R0",
             "Mean community body mass (kg)", expression(Community ~ R[0]),
             expression(bold(R[0] ~ "vs mean body mass")),
             threshold = 1)
}

#' Community R0 against total host density
#' @inheritParams plot_r0_trophic
#' @return A ggplot object.
#' @examples
#' res <- run_multihost_sir(3:7, n_iterations = 10, duration = 150, seed = 1)
#' plot_r0_density(res)
#' @export
plot_r0_density <- function(x) {
  mh_scatter(x, "total_density", "R0",
             "Total host density (ind km^-2)", expression(Community ~ R[0]),
             expression(bold(R[0] ~ "vs total host density")),
             threshold = 1)
}

#' Community R0 against average species competence
#' @inheritParams plot_r0_trophic
#' @return A ggplot object.
#' @examples
#' res <- run_multihost_sir(3:7, n_iterations = 10, duration = 150, seed = 1)
#' plot_r0_competence(res)
#' @export
plot_r0_competence <- function(x) {
  mh_scatter(x, "avg_competence", "R0",
             "Average species competence", expression(Community ~ R[0]),
             expression(bold(R[0] ~ "vs average competence")),
             threshold = 1)
}

#' Cumulative infections against average competence
#' @inheritParams plot_r0_trophic
#' @return A ggplot object.
#' @examples
#' res <- run_multihost_sir(3:7, n_iterations = 10, duration = 150, seed = 1)
#' plot_cumulative_competence(res)
#' @export
plot_cumulative_competence <- function(x) {
  mh_scatter(x, "avg_competence", "cum_infected_prop",
             "Average species competence", "Cumulative proportion infected",
             "Cumulative infections vs average competence",
             ypercent = TRUE)
}

#' Peak prevalence against total system density
#' @inheritParams plot_r0_trophic
#' @return A ggplot object.
#' @examples
#' res <- run_multihost_sir(3:7, n_iterations = 10, duration = 150, seed = 1)
#' plot_peak_density(res)
#' @export
plot_peak_density <- function(x) {
  mh_scatter(x, "total_density", "peak_prevalence",
             "Total host density (ind km^-2)",
             "Peak infected individuals (log scale)",
             "Peak prevalence vs total system density",
             logy = TRUE)
}

#' Days to peak against community R0
#' @inheritParams plot_r0_trophic
#' @return A ggplot object.
#' @examples
#' res <- run_multihost_sir(3:7, n_iterations = 10, duration = 150, seed = 1)
#' plot_days_to_peak_r0(res)
#' @export
plot_days_to_peak_r0 <- function(x) {
  mh_scatter(x, "R0", "days_to_peak",
             expression(Community ~ R[0]), "Days to peak",
             expression(bold("Days to peak vs" ~ R[0])))
}

# --- Predator-structure plots (3) ------------------------------------------

#' Peak prevalence against predator proportion
#' @inheritParams plot_r0_trophic
#' @return A ggplot object.
#' @examples
#' res <- run_multihost_sir(3:7, n_iterations = 10, duration = 150, seed = 1)
#' plot_peak_predator(res)
#' @export
plot_peak_predator <- function(x) {
  mh_scatter(x, "predator_proportion", "peak_prevalence",
             "Predator proportion", "Peak infected individuals (log scale)",
             "Peak prevalence vs predator proportion",
             xpercent = TRUE, logy = TRUE)
}

#' Community R0 against predator proportion
#' @inheritParams plot_r0_trophic
#' @return A ggplot object.
#' @examples
#' res <- run_multihost_sir(3:7, n_iterations = 10, duration = 150, seed = 1)
#' plot_r0_predator(res)
#' @export
plot_r0_predator <- function(x) {
  mh_scatter(x, "predator_proportion", "R0",
             "Predator proportion", expression(Community ~ R[0]),
             expression(bold(R[0] ~ "vs predator proportion")),
             threshold = 1, xpercent = TRUE)
}

#' Trophic infection partitioning at peak against predator proportion
#'
#' Shows how the infected individuals present at the epidemic peak split between
#' prey and predators, as predator proportion rises. Prey and predator shares
#' are drawn as separate smoothed series.
#'
#' @inheritParams plot_r0_trophic
#' @return A ggplot object.
#' @examples
#' res <- run_multihost_sir(3:7, n_iterations = 10, duration = 150, seed = 1)
#' plot_trophic_partition(res)
#' @export

plot_gallery <- function(x) {
  check_result_object(x)
  list(
    # richness comparisons
    r0_richness           = plot_r0_richness(x),
    peak_richness         = plot_peak_prevalence(x),
    cumulative_richness   = plot_cumulative_richness(x),
    days_to_peak_richness = plot_days_to_peak_richness(x),
    trajectories          = plot_trajectories(x),
    # macro-drivers
    r0_trophic            = plot_r0_trophic(x),
    r0_body_mass          = plot_r0_body_mass(x),
    r0_density            = plot_r0_density(x),
    r0_competence         = plot_r0_competence(x),
    cumulative_competence = plot_cumulative_competence(x),
    peak_density          = plot_peak_density(x),
    days_to_peak_r0       = plot_days_to_peak_r0(x),
    # predator structure
    peak_predator         = plot_peak_predator(x),
    r0_predator           = plot_r0_predator(x),
    # methods and performance
    timing                = plot_timing(x),
    allometry             = plot_allometry())
}
