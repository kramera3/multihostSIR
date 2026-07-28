#' multihostSIR: Multi-Host SIR Models with Allometric Scaling
#'
#' Construct and simulate Susceptible-Infected-Recovered models for an arbitrary
#' number of interacting host species. Many wildlife and zoonotic pathogens --
#' including SARS-CoV-2, highly pathogenic avian influenza, rabies and bovine
#' tuberculosis -- infect several species at once. Single-host SIR models cannot
#' represent the asymmetric transmission, predator-prey contact structure, and
#' species-specific susceptibility and virulence that govern such systems.
#'
#' @section Main entry points:
#' \describe{
#'   \item{[allometry()]}{Coefficients of the body-mass power laws.}
#'   \item{[community_spec()]}{Rules used to assemble random host communities.}
#'   \item{[build_community()]}{Draw one host community.}
#'   \item{[transmission_matrix()]}{Asymmetric between-species transmission matrix.}
#'   \item{[community_R0()]}{Dominant eigenvalue of the next-generation matrix.}
#'   \item{[simulate_epidemic()]}{Solve the multi-species SIR system once.}
#'   \item{[run_multihost_sir()]}{Replicate communities across richness levels,
#'     with per-level timing.}
#' }
#'
#' @section Model:
#' For species \eqn{i} in a community of \eqn{n} species,
#' \deqn{dS_i/dt = b_i - \mu_i S_i - \Lambda_i S_i}
#' \deqn{dI_i/dt = \Lambda_i S_i - (\mu_i + \gamma_i + \alpha_i) I_i}
#' \deqn{dR_i/dt = \gamma_i I_i - \mu_i R_i}
#' \deqn{dC_i/dt = \Lambda_i S_i}
#' where \eqn{\Lambda_i = \sum_j \beta_{ij} I_j} is the force of infection and
#' \eqn{C_i} is a cumulative-incidence accumulator that does not feed back into
#' the dynamics. Births \eqn{b_i = \mu_i N_i} are constant, so the disease-free
#' equilibrium is \eqn{S_i = N_i}. Total population size is therefore *not*
#' conserved when virulence \eqn{\alpha_i > 0}.
#'
#' @references
#' Damuth, J. (1981) Population density and body size in mammals.
#' *Nature* **290**, 699-700.
#'
#' Diekmann, O., Heesterbeek, J.A.P. and Metz, J.A.J. (1990) On the definition
#' and the computation of the basic reproduction ratio R0. *Journal of
#' Mathematical Biology* **28**, 365-382.
#'
#' Kleiber, M. (1932) Body size and metabolism. *Hilgardia* **6**, 315-353.
#'
#' Peters, R.H. (1983) *The Ecological Implications of Body Size*. Cambridge
#' University Press.
#'
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_hline geom_line geom_point
#' @importFrom ggplot2 geom_smooth geom_tile geom_text labs
#' @importFrom ggplot2 scale_colour_viridis_c scale_colour_viridis_d
#' @importFrom ggplot2 scale_fill_gradientn scale_fill_viridis_d
#' @importFrom ggplot2 scale_x_continuous scale_y_continuous scale_x_log10 scale_y_log10 coord_cartesian
#' @importFrom ggplot2 theme theme_bw element_blank element_line element_rect
#' @importFrom ggplot2 element_text margin unit
#' @importFrom scales label_number label_log percent_format cut_short_scale
#' @importFrom rlang .data
#' @importFrom stats aggregate coef lm median sd setNames IQR runif
#' @importFrom utils write.csv packageVersion install.packages compareVersion
#' @importFrom deSolve ode
#' @importFrom parallel makePSOCKcluster stopCluster clusterEvalQ clusterExport parLapply clusterSetRNGStream
#' @keywords internal
"_PACKAGE"

## ---------------------------------------------------------------------------
## Internal validation helpers (not exported)
## ---------------------------------------------------------------------------

#' Validate a single finite number
#' @noRd
check_scalar_num <- function(x, name, lo = -Inf, hi = Inf, integer = FALSE) {
  if (is.null(x) || length(x) != 1L || !is.numeric(x) || !is.finite(x)) {
    stop(sprintf("`%s` must be a single finite number (got length %d).",
                 name, length(x)), call. = FALSE)
  }
  if (integer && x != round(x)) {
    stop(sprintf("`%s` must be a whole number (got %s).", name, format(x)),
         call. = FALSE)
  }
  if (x < lo || x > hi) {
    stop(sprintf("`%s` must lie in [%s, %s] (got %s).",
                 name, format(lo), format(hi), format(x)), call. = FALSE)
  }
  invisible(as.numeric(x))
}

#' Validate a numeric vector of strictly positive values
#' @noRd
check_positive_vector <- function(x, name) {
  if (!is.numeric(x) || length(x) < 1L || anyNA(x) || any(!is.finite(x)) ||
      any(x <= 0)) {
    stop(sprintf("`%s` must be a numeric vector of positive, finite values.",
                 name), call. = FALSE)
  }
  invisible(as.numeric(x))
}

#' Validate a single TRUE/FALSE
#' @noRd
check_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(sprintf("`%s` must be a single TRUE or FALSE.", name), call. = FALSE)
  }
  invisible(x)
}

#' Validate a community data frame
#' @noRd
check_community <- function(community, name = "community") {
  required <- c("type", "mass", "mu", "gamma", "alpha", "N", "beta_ii")
  if (!is.data.frame(community)) {
    stop(sprintf("`%s` must be a data frame (see build_community()).", name),
         call. = FALSE)
  }
  missing_cols <- setdiff(required, names(community))
  if (length(missing_cols) > 0L) {
    stop(sprintf("`%s` is missing required column(s): %s.",
                 name, paste(missing_cols, collapse = ", ")), call. = FALSE)
  }
  if (nrow(community) < 1L) {
    stop(sprintf("`%s` must contain at least one species.", name), call. = FALSE)
  }
  if (!all(community$type %in% c("Prey", "Predator"))) {
    stop(sprintf("`%s$type` must contain only \"Prey\" or \"Predator\".", name),
         call. = FALSE)
  }
  invisible(TRUE)
}

#' Null-coalescing helper
#' @noRd
`%||%` <- function(x, y) if (is.null(x)) y else x

#' Sample with replacement from a pool, safely
#'
#' `sample(x, k)` interprets a length-1 numeric `x` as `1:x`, so a single-value
#' pool such as `c(50)` would silently be treated as `1:50`. Indexing through
#' `sample.int()` avoids that trap.
#'
#' @param pool Vector to draw from.
#' @param size Number of draws.
#' @return A vector of length `size`.
#' @noRd
sample_pool <- function(pool, size) {
  pool[sample.int(length(pool), size, replace = TRUE)]
}
