#' Rules for assembling random host communities
#'
#' Bundles every choice that governs how a random community is drawn: the
#' probability that a species is prey, the body-mass pools available to each
#' trophic group, the thinning applied to predator densities, and the ranges
#' from which host competence is sampled.
#'
#' @param prey_prob Probability that any one species is prey. Predators make up
#'   the remainder. Must lie in `[0, 1]`.
#' @param prey_masses,predator_masses Numeric vectors of candidate body masses
#'   (kg) for prey and predators. Masses are sampled from these pools with
#'   replacement.
#' @param predator_density_factor Predator densities are multiplied by this
#'   factor after allometric scaling, representing top-down thinning of upper
#'   trophic levels. `0.1` means predators occur at one tenth of their
#'   allometric baseline density.
#' @param competence_ranges A named list of length-2 numeric vectors giving the
#'   range of the within-species competence parameter \eqn{\beta} for each
#'   competence category. Names become the category labels.
#' @param competence_log_scale If `TRUE` (default), competence is drawn
#'   uniformly on the log scale and back-transformed, which spreads draws evenly
#'   across orders of magnitude. If `FALSE`, competence is drawn uniformly on
#'   the linear scale.
#'
#' @return An object of class `"multihost_community_spec"`.
#'
#' @examples
#' community_spec()
#'
#' # Larger-bodied prey pool, predators at 25% of baseline density
#' community_spec(prey_masses = c(5, 20, 80), predator_density_factor = 0.25)
#'
#' @seealso [build_community()], [allometry()]
#' @export
community_spec <- function(prey_prob = 0.8,
                           prey_masses = c(0.02, 0.1, 2, 10, 50),
                           predator_masses = c(1, 5, 25),
                           predator_density_factor = 0.1,
                           competence_ranges = list(Low  = c(0.01, 0.1),
                                                    Med  = c(0.10, 0.3),
                                                    High = c(0.30, 0.8)),
                           competence_log_scale = TRUE) {

  check_scalar_num(prey_prob, "prey_prob", lo = 0, hi = 1)
  check_scalar_num(predator_density_factor, "predator_density_factor", lo = 0)
  check_positive_vector(prey_masses, "prey_masses")
  check_positive_vector(predator_masses, "predator_masses")
  check_flag(competence_log_scale, "competence_log_scale")

  if (!is.list(competence_ranges) || length(competence_ranges) < 1L ||
      is.null(names(competence_ranges)) || any(names(competence_ranges) == "")) {
    stop("`competence_ranges` must be a named list of length-2 numeric vectors.",
         call. = FALSE)
  }
  for (nm in names(competence_ranges)) {
    rng <- competence_ranges[[nm]]
    if (!is.numeric(rng) || length(rng) != 2L || any(!is.finite(rng)) ||
        any(rng <= 0) || rng[1] > rng[2]) {
      stop(sprintf(paste0("`competence_ranges$%s` must be two finite positive ",
                          "numbers in increasing order."), nm), call. = FALSE)
    }
  }

  structure(
    list(prey_prob = prey_prob,
         prey_masses = prey_masses,
         predator_masses = predator_masses,
         predator_density_factor = predator_density_factor,
         competence_ranges = competence_ranges,
         competence_log_scale = competence_log_scale),
    class = "multihost_community_spec")
}

#' @param x An object of class `"multihost_community_spec"`.
#' @param ... Ignored.
#' @rdname community_spec
#' @export
print.multihost_community_spec <- function(x, ...) {
  cat("<multihost community specification>\n")
  cat(sprintf("  P(prey)                 : %.2f\n", x$prey_prob))
  cat(sprintf("  prey masses (kg)        : %s\n",
              paste(format(x$prey_masses, trim = TRUE), collapse = ", ")))
  cat(sprintf("  predator masses (kg)    : %s\n",
              paste(format(x$predator_masses, trim = TRUE), collapse = ", ")))
  cat(sprintf("  predator density factor : %.3g\n", x$predator_density_factor))
  cat(sprintf("  competence sampling     : %s scale\n",
              if (x$competence_log_scale) "log" else "linear"))
  for (nm in names(x$competence_ranges)) {
    cat(sprintf("    %-5s : [%.3g, %.3g]\n", nm,
                x$competence_ranges[[nm]][1], x$competence_ranges[[nm]][2]))
  }
  invisible(x)
}

#' Draw one random host community
#'
#' Assembles a community of `richness` species. Each species is assigned a
#' trophic type, a body mass drawn from the matching mass pool, allometrically
#' scaled demographic rates, and a competence category with an associated
#' within-species transmission rate.
#'
#' @param richness Number of host species. Must be a whole number `>= 1`.
#' @param spec A specification from [community_spec()].
#' @param allo An allometry object from [allometry()].
#'
#' @return A data frame with one row per species and the columns:
#' \describe{
#'   \item{species}{Integer species index.}
#'   \item{type}{`"Prey"` or `"Predator"`.}
#'   \item{mass}{Body mass (kg).}
#'   \item{trophic_weight}{`1` for prey, `2` for predators.}
#'   \item{mu, gamma, alpha}{Mortality, recovery and virulence rates (day\eqn{^{-1}}).}
#'   \item{N}{Carrying-capacity density (ind km\eqn{^{-2}}), after predator thinning.}
#'   \item{competence_category}{Label drawn from `spec$competence_ranges`.}
#'   \item{competence}{Within-species competence \eqn{\beta}.}
#'   \item{beta_ii}{Within-species transmission rate, `competence / N`.}
#' }
#'
#' @examples
#' set.seed(1)
#' comm <- build_community(4)
#' comm
#'
#' # Densities fall with body mass; predator rows are thinned further
#' comm[order(comm$mass), c("type", "mass", "N")]
#'
#' @seealso [transmission_matrix()], [community_R0()]
#' @export
build_community <- function(richness, spec = community_spec(),
                            allo = allometry()) {

  check_scalar_num(richness, "richness", lo = 1, integer = TRUE)
  if (!inherits(spec, "multihost_community_spec")) {
    stop("`spec` must be created by community_spec().", call. = FALSE)
  }
  if (!inherits(allo, "multihost_allometry")) {
    stop("`allo` must be created by allometry().", call. = FALSE)
  }
  richness <- as.integer(richness)

  type <- sample(c("Prey", "Predator"), richness, replace = TRUE,
                 prob = c(spec$prey_prob, 1 - spec$prey_prob))

  # Draw each species' mass from the pool matching its trophic type. Sampling
  # both pools at full length and selecting with ifelse() keeps a single RNG
  # path regardless of how many predators are drawn.
  #
  # sample_pool() is used rather than sample(): sample(x, k) treats a length-1
  # numeric x as 1:x, so a single-mass pool such as c(50) would silently draw
  # masses from 1:50.
  mass <- ifelse(type == "Prey",
                 sample_pool(spec$prey_masses, richness),
                 sample_pool(spec$predator_masses, richness))

  rates <- apply_allometry(mass, allo)
  is_predator <- type == "Predator"
  rates$N[is_predator] <- rates$N[is_predator] * spec$predator_density_factor

  categories <- names(spec$competence_ranges)
  category <- categories[sample.int(length(categories), richness,
                                    replace = TRUE)]
  competence <- vapply(category, function(ct) {
    rng <- spec$competence_ranges[[ct]]
    if (spec$competence_log_scale) {
      exp(runif(1L, log(rng[1L]), log(rng[2L])))
    } else {
      runif(1L, rng[1L], rng[2L])
    }
  }, numeric(1L), USE.NAMES = FALSE)

  out <- data.frame(
    species             = seq_len(richness),
    type                = type,
    mass                = mass,
    trophic_weight      = ifelse(is_predator, 2, 1),
    mu                  = rates$mu,
    gamma               = rates$gamma,
    alpha               = rates$alpha,
    N                   = rates$N,
    competence_category = factor(category, levels = categories),
    competence          = competence,
    stringsAsFactors    = FALSE)

  out$beta_ii <- out$competence / out$N
  out
}
