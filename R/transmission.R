#' Trophic transmission coefficients
#'
#' Scaling factors applied to the geometric mean of two species within-species
#' transmission rates, according to the trophic types of the recipient and the
#' source. The defaults encode a system in which prey mix freely with one
#' another, predators contract infection from their prey, and prey are not
#' infected by predators.
#'
#' @param prey_prey Coefficient for prey infecting prey.
#' @param pred_pred Coefficient for predators infecting predators.
#' @param prey_to_pred Coefficient for prey infecting predators, for example
#'   through consumption of infected prey.
#' @param pred_to_prey Coefficient for predators infecting prey. Zero by default,
#'   blocking this route entirely.
#'
#' @return An object of class `"multihost_trophic"`.
#'
#' @examples
#' trophic_coefficients()
#'
#' # Allow a weak predator-to-prey route
#' trophic_coefficients(pred_to_prey = 0.05)
#'
#' @seealso [transmission_matrix()]
#' @export
trophic_coefficients <- function(prey_prey = 0.50, pred_pred = 0.25,
                                 prey_to_pred = 0.25, pred_to_prey = 0.00) {
  check_scalar_num(prey_prey,    "prey_prey",    lo = 0)
  check_scalar_num(pred_pred,    "pred_pred",    lo = 0)
  check_scalar_num(prey_to_pred, "prey_to_pred", lo = 0)
  check_scalar_num(pred_to_prey, "pred_to_prey", lo = 0)

  structure(list(prey_prey = prey_prey, pred_pred = pred_pred,
                 prey_to_pred = prey_to_pred, pred_to_prey = pred_to_prey),
            class = "multihost_trophic")
}

#' Between-species transmission matrix
#'
#' Builds the \eqn{n \times n} matrix \eqn{\beta} whose element
#' \eqn{\beta_{ij}} is the rate at which one infected individual of species
#' \eqn{j} infects susceptibles of species \eqn{i}.
#'
#' @details
#' Diagonal entries are the within-species transmission rates
#' \eqn{\beta_{ii} = \mathrm{competence}_i / N_i}. Off-diagonal entries are the
#' geometric mean of the two within-species rates, scaled by a trophic
#' coefficient:
#'
#' \deqn{\beta_{ij} = c(\mathrm{type}_i, \mathrm{type}_j)
#'        \sqrt{\beta_{ii}\,\beta_{jj}}, \quad i \neq j}
#'
#' The matrix is asymmetric whenever `prey_to_pred` differs from
#' `pred_to_prey`. Row \eqn{i} is the *recipient* of infection and column
#' \eqn{j} is the *source*.
#'
#' @param community A community data frame, typically from [build_community()].
#'   Must contain columns `type` and `beta_ii`.
#' @param coefficients Trophic coefficients from [trophic_coefficients()].
#'
#' @return An \eqn{n \times n} numeric matrix with dimnames taken from
#'   `community$species` where available.
#'
#' @examples
#' set.seed(42)
#' comm <- build_community(4)
#' beta <- transmission_matrix(comm)
#' round(beta, 6)
#'
#' # Prey are never infected by predators under the defaults
#' beta[comm$type == "Prey", comm$type == "Predator", drop = FALSE]
#'
#' @seealso [community_R0()], [plot_transmission_matrix()]
#' @export
transmission_matrix <- function(community,
                                coefficients = trophic_coefficients()) {
  check_community(community)
  if (!inherits(coefficients, "multihost_trophic")) {
    stop("`coefficients` must be created by trophic_coefficients().",
         call. = FALSE)
  }

  b <- community$beta_ii
  n <- length(b)
  is_pred <- community$type == "Predator"

  # geometric mean of within-species rates: gm[i, j] = sqrt(b_i * b_j)
  gm <- sqrt(outer(b, b))

  # Rows are recipients, columns are sources.
  coef_mat <- matrix(coefficients$prey_prey, nrow = n, ncol = n)
  coef_mat[is_pred,  is_pred]  <- coefficients$pred_pred
  coef_mat[is_pred,  !is_pred] <- coefficients$prey_to_pred   # prey infect predators
  coef_mat[!is_pred, is_pred]  <- coefficients$pred_to_prey   # predators infect prey

  beta <- coef_mat * gm
  diag(beta) <- b

  labels <- if (!is.null(community$species)) {
    as.character(community$species)
  } else {
    as.character(seq_len(n))
  }
  dimnames(beta) <- list(labels, labels)
  beta
}
