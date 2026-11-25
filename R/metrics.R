#' Per-species epidemic metrics
#'
#' Summarises the outcome of a single simulated epidemic for each host species.
#'
#' @param sim An object from [simulate_epidemic()].
#'
#' @return A data frame with one row per species and columns:
#' \describe{
#'   \item{species}{Species index.}
#'   \item{type}{Trophic type.}
#'   \item{mass}{Body mass (kg).}
#'   \item{peak_infected}{Largest number infected simultaneously.}
#'   \item{days_to_peak}{Day on which `peak_infected` occurs.}
#'   \item{cum_infected}{Cumulative infections over the horizon.}
#'   \item{attack_rate}{`cum_infected` divided by the initial density `N`.}
#' }
#'
#' @examples
#' set.seed(42)
#' sim <- simulate_epidemic(build_community(4), duration = 200)
#' species_metrics(sim)
#'
#' @seealso [simulate_epidemic()]
#' @export
species_metrics <- function(sim) {
  if (!inherits(sim, "multihost_epidemic")) {
    stop("`sim` must be created by simulate_epidemic().", call. = FALSE)
  }
  community <- sim$community
  n <- nrow(community)
  out <- sim$out

  # Columns are addressed by name, so the result is unaffected by any change to
  # the ordering deSolve returns.
  I_mat <- as.matrix(out[, paste0("I", seq_len(n)), drop = FALSE])
  C_final <- as.numeric(out[nrow(out), paste0("C", seq_len(n))])

  peak_infected <- apply(I_mat, 2L, max)
  days_to_peak <- out$time[apply(I_mat, 2L, which.max)]

  data.frame(
    species       = community$species %||% seq_len(n),
    type          = community$type,
    mass          = community$mass,
    peak_infected = as.numeric(peak_infected),
    days_to_peak  = as.numeric(days_to_peak),
    cum_infected  = C_final,
    attack_rate   = C_final / community$N,
    stringsAsFactors = FALSE,
    row.names = NULL)
}

#' Community-level metrics for one epidemic
#'
#' @param sim An object from [simulate_epidemic()].
#' @return A one-row data frame plus the daily total-prevalence series.
#' @noRd
community_metrics <- function(sim) {
  community <- sim$community
  n <- nrow(community)
  out <- sim$out
  duration <- max(out$time)

  I_mat <- as.matrix(out[, paste0("I", seq_len(n)), drop = FALSE])
  total_I_daily <- rowSums(I_mat)
  total_C_final <- sum(as.numeric(out[nrow(out), paste0("C", seq_len(n))]))

  is_pred <- community$type == "Predator"
  total_N <- sum(community$N)
  total_births <- sum(community$mu * community$N) * duration

  peak_row <- which.max(total_I_daily)
  I_at_peak <- I_mat[peak_row, ]

  counts <- table(community$competence_category)
  composition <- paste0(
    "R", n, "_",
    paste(sprintf("%s:%d", names(counts), as.integer(counts)), collapse = "_"))

  n_high <- if ("High" %in% names(counts)) as.integer(counts[["High"]]) else NA_integer_
  n_low  <- if ("Low"  %in% names(counts)) as.integer(counts[["Low"]])  else NA_integer_

  row <- data.frame(
    richness               = n,
    composition            = composition,
    avg_competence         = mean(community$competence),
    prop_high_competence   = n_high / n,
    prop_low_competence    = n_low / n,
    avg_body_mass          = mean(community$mass),
    pop_weighted_body_mass = sum(community$mass * community$N) / total_N,
    avg_trophic_level      = mean(community$trophic_weight),
    n_predators            = sum(is_pred),
    predator_proportion    = sum(is_pred) / n,
    total_density          = total_N,
    prey_density           = sum(community$N[!is_pred]),
    predator_density       = sum(community$N[is_pred]),
    R0                     = sim$R0,
    peak_prevalence        = max(total_I_daily),
    days_to_peak           = out$time[peak_row],
    cum_infected_prop      = total_C_final / (total_N + total_births),
    prey_infected_at_peak  = sum(I_at_peak[!is_pred]),
    pred_infected_at_peak  = sum(I_at_peak[is_pred]),
    seed_species_type      = community$type[sim$seed_species],
    stringsAsFactors       = FALSE)

  list(row = row, prevalence = total_I_daily, time = out$time)
}
