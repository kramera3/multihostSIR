#' Summarise a multi-host simulation by richness
#'
#' @param object A `"multihost_sir"` object from [run_multihost_sir()].
#' @param digits Number of digits used when rounding the table.
#' @param ... Ignored.
#'
#' @return An object of class `"summary.multihost_sir"`, a list containing the
#'   per-richness table `by_richness`, the timing table, and the fitted
#'   complexity exponents.
#'
#' @examples
#' res <- run_multihost_sir(2:4, n_iterations = 5, duration = 120, seed = 1)
#' summary(res)
#'
#' @export
summary.multihost_sir <- function(object, digits = 3, ...) {
  r <- object$results
  split_by <- factor(r$richness, levels = sort(unique(r$richness)))

  agg <- function(f, col) {
    vapply(split(r[[col]], split_by), f, numeric(1L))
  }

  tab <- data.frame(
    richness        = as.integer(levels(split_by)),
    n_sims          = as.integer(table(split_by)),
    median_R0       = round(agg(median, "R0"), digits),
    IQR_R0          = round(agg(IQR, "R0"), digits),
    pct_R0_above_1  = round(100 * agg(function(z) mean(z > 1), "R0"), 1),
    median_peak     = round(agg(median, "peak_prevalence"), 1),
    median_days     = round(agg(median, "days_to_peak"), 0),
    median_cum_prop = round(agg(median, "cum_infected_prop"), 4),
    stringsAsFactors = FALSE,
    row.names = NULL)

  structure(
    list(by_richness = tab,
         timing = object$timing,
         complexity = timing_complexity(object),
         settings = object$settings),
    class = "summary.multihost_sir")
}

#' @param x A `"summary.multihost_sir"` object.
#' @param ... Ignored.
#' @rdname summary.multihost_sir
#' @export
print.summary.multihost_sir <- function(x, ...) {
  cat("Multi-host SIR simulation summary\n")
  cat(sprintf("%d iterations per richness level, %d-day horizon\n\n",
              x$settings$n_iterations, x$settings$duration))
  print(x$by_richness, row.names = FALSE)
  cat("\nRun-time scaling\n")
  cat(sprintf("  time ~ richness^%.2f   (R2 = %.3f)\n",
              x$complexity$exponent_richness, x$complexity$r_squared_richness))
  cat(sprintf("  time ~ state_dim^%.2f  (R2 = %.3f)\n",
              x$complexity$exponent_state_dim, x$complexity$r_squared_state_dim))
  cat(sprintf("  total elapsed: %.1f s\n", sum(x$timing$elapsed_sec)))
  invisible(x)
}

#' Empirical run-time complexity
#'
#' Fits `log(elapsed) ~ log(richness)` and `log(elapsed) ~ log(state_dim)` to the
#' timing table, giving approximate polynomial orders \eqn{p} such that
#' run time scales as \eqn{richness^p}. The \eqn{n \times n} matrix-vector
#' product in the force-of-infection term and the \eqn{n \times n}
#' eigendecomposition used for \eqn{R_0} make an exponent between 2 and 3 the
#' expectation at large \eqn{n}.
#'
#' @param x A `"multihost_sir"` object from [run_multihost_sir()].
#'
#' @return A list with the two exponents, their \eqn{R^2} values, and the two
#'   fitted `lm` objects. Exponents are `NA` when fewer than two richness levels
#'   were simulated.
#'
#' @examples
#' res <- run_multihost_sir(2:5, n_iterations = 3, duration = 100, seed = 1)
#' timing_complexity(res)$exponent_richness
#'
#' @export
timing_complexity <- function(x) {
  if (!inherits(x, "multihost_sir")) {
    stop("`x` must be created by run_multihost_sir().", call. = FALSE)
  }
  tm <- x$timing

  # A log-log fit needs at least two distinct, strictly positive timings.
  usable <- nrow(tm) >= 2L && all(tm$elapsed_sec > 0) &&
    length(unique(tm$richness)) >= 2L

  if (!usable) {
    return(list(exponent_richness = NA_real_, exponent_state_dim = NA_real_,
                r_squared_richness = NA_real_, r_squared_state_dim = NA_real_,
                fit_richness = NULL, fit_state_dim = NULL))
  }

  fit_r <- lm(log(elapsed_sec) ~ log(richness), data = tm)
  fit_s <- lm(log(elapsed_sec) ~ log(state_dim), data = tm)

  list(exponent_richness  = unname(coef(fit_r)[2L]),
       exponent_state_dim = unname(coef(fit_s)[2L]),
       r_squared_richness  = summary(fit_r)$r.squared,
       r_squared_state_dim = summary(fit_s)$r.squared,
       fit_richness = fit_r,
       fit_state_dim = fit_s)
}

#' Mean outcomes by predator proportion
#'
#' @param x A `"multihost_sir"` object.
#' @return A data frame of mean \eqn{R_0}, peak prevalence and cumulative
#'   infected proportion for each realised predator proportion.
#' @examples
#' res <- run_multihost_sir(3:5, n_iterations = 5, duration = 120, seed = 1)
#' predator_summary(res)
#' @export
predator_summary <- function(x) {
  if (!inherits(x, "multihost_sir")) {
    stop("`x` must be created by run_multihost_sir().", call. = FALSE)
  }
  r <- x$results
  key <- round(r$predator_proportion, 3)
  out <- aggregate(r[, c("R0", "peak_prevalence", "cum_infected_prop")],
                   by = list(predator_proportion = key), FUN = mean)
  names(out)[-1L] <- c("mean_R0", "mean_peak", "mean_cum_prop")
  out[order(out$predator_proportion), , drop = FALSE]
}

#' Mean outcomes by richness and high-competence proportion
#'
#' @param x A `"multihost_sir"` object.
#' @return A data frame with mean and standard deviation of \eqn{R_0}, and mean
#'   peak prevalence, for each combination of richness and the proportion of
#'   high-competence hosts.
#' @examples
#' res <- run_multihost_sir(3:5, n_iterations = 5, duration = 120, seed = 1)
#' head(competence_summary(res))
#' @export
competence_summary <- function(x) {
  if (!inherits(x, "multihost_sir")) {
    stop("`x` must be created by run_multihost_sir().", call. = FALSE)
  }
  r <- x$results
  grp <- list(richness = r$richness,
              prop_high_competence = round(r$prop_high_competence, 3))

  m_r0   <- aggregate(r["R0"], by = grp, FUN = mean)
  s_r0   <- aggregate(r["R0"], by = grp, FUN = sd)
  m_peak <- aggregate(r["peak_prevalence"], by = grp, FUN = mean)

  out <- data.frame(
    richness = m_r0$richness,
    prop_high_competence = m_r0$prop_high_competence,
    mean_R0 = m_r0$R0,
    sd_R0 = s_r0$R0,
    mean_peak = m_peak$peak_prevalence,
    stringsAsFactors = FALSE)
  out[order(out$richness, out$prop_high_competence), , drop = FALSE]
}
