#' Describe the assumptions behind a run
#'
#' Builds a short, human-readable record of which quantities were supplied by the
#' user and which fell back to the package's built-in assumptions. Printed
#' alongside results so the reader can see, at a glance, what was assumed.
#'
#' @param allo An allometry object.
#' @param spec A community specification, or `NULL` when a fixed community was
#'   supplied.
#' @param coefficients Trophic coefficients, or `NULL`.
#' @return A character vector, one line per assumption.
#' @noRd
describe_assumptions <- function(allo = NULL, spec = NULL,
                                 coefficients = NULL) {
  lines <- character(0)

  if (!is.null(allo)) {
    lines <- c(lines,
      "Allometric scaling (rate = a * mass^b), per day, density per km^2:",
      sprintf("  mortality  mu    = %g * mass^%g", allo$a_mort, allo$b_mort),
      sprintf("  recovery   gamma = %g * mass^%g", allo$a_recovery, allo$b_recovery),
      sprintf("  density    N     = %g * mass^%g", allo$a_density, allo$b_density),
      sprintf("  virulence  alpha = %g * mass^%g", allo$a_virulence, allo$b_virulence))
  }

  if (!is.null(spec)) {
    lines <- c(lines,
      "Community assembly:",
      sprintf("  P(prey) = %.2f; predator density factor = %.3g",
              spec$prey_prob, spec$predator_density_factor),
      sprintf("  competence sampled on the %s scale",
              if (spec$competence_log_scale) "log" else "linear"))
  }

  if (!is.null(coefficients)) {
    lines <- c(lines,
      "Trophic transmission coefficients (recipient <- source):",
      sprintf("  prey<-prey %.2f   pred<-pred %.2f   pred<-prey %.2f   prey<-pred %.2f",
              coefficients$prey_prey, coefficients$pred_pred,
              coefficients$prey_to_pred, coefficients$pred_to_prey))
  }

  lines
}

#' Print the assumptions for a fixed-community run
#' @noRd
print_assumptions <- function(assumptions, x) {
  cat("Assumptions used\n")
  if (isTRUE(assumptions$deterministic)) {
    cat("  Fully specified: every rate and the transmission matrix were\n")
    cat("  supplied, so the result is exact (no random iteration).\n")
  } else {
    cat("  Partly filled from built-in assumptions:\n")
    if (!assumptions$beta_supplied) {
      cat("    - transmission matrix built from competence and trophic rules\n")
    }
    if (length(assumptions$filled_columns) > 0L) {
      cat(sprintf("    - filled from allometry: %s\n",
                  paste(assumptions$filled_columns, collapse = ", ")))
    }
    cat("  Because some values were generated, treat this as one draw; use\n")
    cat("  run_multihost_sir() to average over many.\n")
  }
  invisible(NULL)
}
