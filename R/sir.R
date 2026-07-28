#' Right-hand side of the multi-species SIR system
#'
#' State vector layout, of length `4 * n`:
#' `[S_1..S_n | I_1..I_n | R_1..R_n | C_1..C_n]`. `C` accumulates cumulative
#' incidence and does not feed back into the dynamics.
#'
#' @param time Current time (unused; required by [deSolve::ode()]).
#' @param state Numeric state vector.
#' @param parameters List with elements `n`, `beta`, `mu`, `gamma`, `alpha`
#'   and `births`.
#' @return A list whose single element is the vector of derivatives.
#' @noRd
multi_species_SIR <- function(time, state, parameters) {
  n <- parameters$n
  S <- state[seq_len(n)]
  I <- state[n + seq_len(n)]
  R <- state[2L * n + seq_len(n)]

  lambda <- as.numeric(parameters$beta %*% I)   # force of infection

  dS <- parameters$births - parameters$mu * S - lambda * S
  dI <- lambda * S - (parameters$mu + parameters$gamma + parameters$alpha) * I
  dR <- parameters$gamma * I - parameters$mu * R
  dC <- lambda * S

  list(c(dS, dI, dR, dC))
}

#' Discrete extinction event
#'
#' Numerical integration can leave vanishingly small residual infections (for
#' example 1e-15 individuals) that sustain phantom transmission forever. This
#' event snaps any compartment holding fewer than `extinction_threshold`
#' infected individuals to exactly zero.
#'
#' @inheritParams multi_species_SIR
#' @return The modified state vector.
#' @noRd
die_out_event <- function(time, state, parameters) {
  n <- parameters$n
  idx <- n + seq_len(n)
  below <- state[idx] < parameters$extinction_threshold
  state[idx][below] <- 0
  state
}

#' Simulate one multi-species epidemic
#'
#' Solves the multi-species SIR system for a single community over a fixed
#' horizon, seeding the epidemic in one species.
#'
#' @param community A community data frame from [build_community()].
#' @param beta Transmission matrix. Defaults to `transmission_matrix(community)`.
#' @param duration Length of the simulation in days.
#' @param seed_species Index of the species in which the epidemic starts, or
#'   `"max_competence"` (default) to seed the most competent host.
#' @param seed_size Number of infected individuals introduced at time zero.
#' @param extinction_threshold Infected compartments falling below this many
#'   individuals are set to zero at each integer time step. Set to `0` to
#'   disable.
#' @param ... Further arguments passed to [deSolve::ode()], for example
#'   `method = "lsoda"`.
#'
#' @return A list of class `"multihost_epidemic"` with elements:
#' \describe{
#'   \item{out}{Data frame of the solution. Columns are `time`, then `S1..Sn`,
#'     `I1..In`, `R1..Rn`, `C1..Cn`.}
#'   \item{community}{The community that was simulated.}
#'   \item{beta}{The transmission matrix used.}
#'   \item{R0}{Community \eqn{R_0}.}
#'   \item{seed_species}{Index of the seeded species.}
#' }
#'
#' @examples
#' set.seed(42)
#' comm <- build_community(4)
#' sim <- simulate_epidemic(comm, duration = 200)
#' str(sim$out[1:3, 1:5])
#' sim$R0
#'
#' # Per-species outcomes
#' species_metrics(sim)
#'
#' @seealso [species_metrics()], [run_multihost_sir()]
#' @export
simulate_epidemic <- function(community,
                              beta = transmission_matrix(community),
                              duration = 365,
                              seed_species = "max_competence",
                              seed_size = 1,
                              extinction_threshold = 1,
                              ...) {

  check_community(community)
  check_scalar_num(duration, "duration", lo = 1, integer = TRUE)
  check_scalar_num(seed_size, "seed_size", lo = 0)
  check_scalar_num(extinction_threshold, "extinction_threshold", lo = 0)

  n <- nrow(community)

  if (identical(seed_species, "max_competence")) {
    if (is.null(community$competence)) {
      stop(paste0("`seed_species = \"max_competence\"` needs a `competence` ",
                  "column; supply an integer index instead."), call. = FALSE)
    }
    seed_idx <- which.max(community$competence)
  } else {
    check_scalar_num(seed_species, "seed_species", lo = 1, hi = n,
                     integer = TRUE)
    seed_idx <- as.integer(seed_species)
  }

  S0 <- community$N
  I0 <- rep(0, n)
  I0[seed_idx] <- seed_size
  R0_state <- rep(0, n)
  C0 <- I0   # seed individuals count as the first cumulative infections

  state <- setNames(
    c(S0, I0, R0_state, C0),
    c(paste0("S", seq_len(n)), paste0("I", seq_len(n)),
      paste0("R", seq_len(n)), paste0("C", seq_len(n))))

  parameters <- list(
    n = n, beta = beta,
    mu = community$mu, gamma = community$gamma, alpha = community$alpha,
    births = community$mu * community$N,
    extinction_threshold = extinction_threshold)

  times <- 0:as.integer(duration)

  event_arg <- if (extinction_threshold > 0) {
    list(func = die_out_event, time = seq_len(as.integer(duration)))
  } else {
    NULL
  }

  out <- deSolve::ode(y = state, times = times, func = multi_species_SIR,
                      parms = parameters, events = event_arg, ...)
  out <- as.data.frame(out)

  structure(
    list(out = out,
         community = community,
         beta = beta,
         R0 = community_R0(community, beta),
         seed_species = seed_idx),
    class = "multihost_epidemic")
}
