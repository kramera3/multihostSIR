#' Next-generation matrix
#'
#' Constructs the next-generation matrix \eqn{K} for a multi-host community,
#' \deqn{K_{ij} = \beta_{ij} \, N_i \, D_j, \qquad
#'       D_j = \frac{1}{\gamma_j + \mu_j + \alpha_j}}
#' where \eqn{D_j} is the mean infectious period of species \eqn{j} and
#' \eqn{N_i} is the disease-free density of species \eqn{i}. Element
#' \eqn{K_{ij}} is the expected number of species-\eqn{i} infections generated
#' by one infected individual of species \eqn{j} in an otherwise susceptible
#' community.
#'
#' @param community A community data frame with columns `N`, `gamma`, `mu` and
#'   `alpha`, typically from [build_community()].
#' @param beta An \eqn{n \times n} transmission matrix. Defaults to
#'   `transmission_matrix(community)`.
#'
#' @return An \eqn{n \times n} numeric matrix.
#'
#' @examples
#' set.seed(7)
#' comm <- build_community(3)
#' next_generation_matrix(comm)
#'
#' @seealso [community_R0()]
#' @export
next_generation_matrix <- function(community,
                                   beta = transmission_matrix(community)) {
  check_community(community)
  n <- nrow(community)
  if (!is.matrix(beta) || nrow(beta) != n || ncol(beta) != n) {
    stop(sprintf("`beta` must be a %d x %d matrix.", n, n), call. = FALSE)
  }

  denom <- community$gamma + community$mu + community$alpha
  if (any(denom <= 0)) {
    stop("gamma + mu + alpha must be strictly positive for every species.",
         call. = FALSE)
  }
  duration <- 1 / denom

  # K[i, j] = beta[i, j] * N[i] * duration[j]
  beta * outer(community$N, duration)
}

#' Community basic reproduction number
#'
#' Computes \eqn{R_0} for the whole community as the dominant eigenvalue of the
#' next-generation matrix. Values above one indicate that a pathogen introduced
#' into the community will spread.
#'
#' @inheritParams next_generation_matrix
#'
#' @return A single numeric value, the dominant real eigenvalue of \eqn{K}.
#'
#' @details
#' Because \eqn{K} has non-negative entries, the Perron-Frobenius theorem
#' guarantees that the spectral radius is attained by a real, non-negative
#' eigenvalue. The function returns the largest real part among the eigenvalues,
#' which coincides with that value.
#'
#' @examples
#' set.seed(42)
#' comm <- build_community(5)
#' community_R0(comm)
#'
#' # R0 rises with host competence
#' more_competent <- comm
#' more_competent$competence <- more_competent$competence * 2
#' more_competent$beta_ii <- more_competent$competence / more_competent$N
#' community_R0(more_competent)
#'
#' # But R0 is invariant when every density is scaled by the same factor:
#' # beta_ii = competence / N falls as 1/N while K carries a factor of N, so
#' # the two cancel under this frequency-dependent parameterisation.
#' denser <- comm
#' denser$N <- denser$N * 10
#' denser$beta_ii <- denser$competence / denser$N
#' isTRUE(all.equal(community_R0(denser), community_R0(comm)))
#'
#' @references
#' Diekmann, O., Heesterbeek, J.A.P. and Metz, J.A.J. (1990) On the definition
#' and the computation of the basic reproduction ratio R0 in models for
#' infectious diseases in heterogeneous populations. *Journal of Mathematical
#' Biology* **28**, 365-382.
#'
#' @seealso [next_generation_matrix()], [transmission_matrix()]
#' @export
community_R0 <- function(community, beta = transmission_matrix(community)) {
  K <- next_generation_matrix(community, beta)
  max(Re(eigen(K, only.values = TRUE)$values))
}
