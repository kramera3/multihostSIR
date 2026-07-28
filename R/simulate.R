#' Simulate one random community from end to end
#'
#' Draws a community, builds its transmission matrix, computes \eqn{R_0}, solves
#' the SIR system, and returns a one-row summary together with the daily
#' community-wide prevalence series. This is the unit of work replicated by
#' [run_multihost_sir()].
#'
#' @param richness Number of host species.
#' @param spec A specification from [community_spec()].
#' @param allo An allometry object from [allometry()].
#' @param coefficients Trophic coefficients from [trophic_coefficients()].
#' @param duration Simulation horizon in days.
#' @param ... Passed to [simulate_epidemic()].
#'
#' @return A list with elements `row` (one-row data frame), `prevalence`
#'   (numeric vector of daily total infected) and `time`.
#'
#' @examples
#' set.seed(42)
#' one <- simulate_once(5, duration = 200)
#' one$row[, c("richness", "R0", "peak_prevalence", "days_to_peak")]
#'
#' @seealso [run_multihost_sir()]
#' @export
simulate_once <- function(richness,
                          spec = community_spec(),
                          allo = allometry(),
                          coefficients = trophic_coefficients(),
                          duration = 365,
                          ...) {
  community <- build_community(richness, spec = spec, allo = allo)
  beta <- transmission_matrix(community, coefficients = coefficients)
  sim <- simulate_epidemic(community, beta = beta, duration = duration, ...)
  community_metrics(sim)
}

#' Replicate multi-host epidemics across community sizes
#'
#' Runs `n_iterations` independent random communities at each requested richness
#' level, recording epidemiological outcomes, the mean epidemic trajectory, and
#' the wall-clock time taken by each level. The timing table supports the
#' run-time scaling study described in [timing_complexity()].
#'
#' @param richness_levels Integer vector of community sizes to simulate. Each
#'   value must be at least 1; at least 2 species are needed for any
#'   between-species transmission to occur.
#' @param n_iterations Number of stochastic replicate communities per richness
#'   level.
#' @param duration Simulation horizon in days.
#' @param spec A specification from [community_spec()].
#' @param allo An allometry object from [allometry()].
#' @param coefficients Trophic coefficients from [trophic_coefficients()].
#' @param seed Optional integer seed applied before the timed loop, for
#'   reproducibility. `NULL` (default) leaves the RNG untouched.
#' @param parallel If `TRUE`, distribute iterations across a PSOCK cluster.
#'   Defaults to `FALSE`, which keeps the default run single-threaded.
#' @param n_cores Number of worker processes when `parallel = TRUE`. Defaults to
#'   `2`, the maximum permitted in CRAN checks. Raise it for real work.
#' @param warmup If `TRUE` (default), a short untimed batch is run first so that
#'   one-off costs such as cluster start-up and JIT compilation do not inflate
#'   the timing of the first richness level.
#' @param progress If `TRUE`, report progress after each richness level.
#' @param ... Passed to [simulate_epidemic()].
#'
#' @return An object of class `"multihost_sir"`: a list with elements
#' \describe{
#'   \item{results}{Data frame with one row per iteration.}
#'   \item{trajectories}{Data frame of mean daily prevalence per richness level.}
#'   \item{timing}{Data frame of elapsed time per richness level.}
#'   \item{settings}{The inputs used.}
#'   \item{call}{The matched call.}
#' }
#'
#' @examples
#' # Small, fast example. Use n_iterations in the hundreds or thousands for
#' # publishable results.
#' res <- run_multihost_sir(richness_levels = 2:4, n_iterations = 5,
#'                          duration = 120, seed = 1)
#' res
#' summary(res)
#'
#' \donttest{
#' # A larger run, in parallel
#' big <- run_multihost_sir(richness_levels = 3:8, n_iterations = 200,
#'                          seed = 42, parallel = TRUE, n_cores = 2)
#' plot(big, which = "r0")
#' }
#'
#' @seealso [simulate_once()], [timing_complexity()], [export_results()]
#' @export
run_multihost_sir <- function(richness_levels = 3:8,
                              n_iterations = 100,
                              duration = 365,
                              spec = community_spec(),
                              allo = allometry(),
                              coefficients = trophic_coefficients(),
                              seed = NULL,
                              parallel = FALSE,
                              n_cores = 2L,
                              warmup = TRUE,
                              progress = FALSE,
                              ...) {

  cl_call <- match.call()

  if (!is.numeric(richness_levels) || length(richness_levels) < 1L ||
      any(!is.finite(richness_levels)) ||
      any(richness_levels != round(richness_levels)) ||
      any(richness_levels < 1)) {
    stop("`richness_levels` must be a vector of whole numbers >= 1.",
         call. = FALSE)
  }
  richness_levels <- as.integer(richness_levels)
  check_scalar_num(n_iterations, "n_iterations", lo = 1, integer = TRUE)
  check_scalar_num(duration, "duration", lo = 1, integer = TRUE)
  check_scalar_num(n_cores, "n_cores", lo = 1, integer = TRUE)
  check_flag(parallel, "parallel")
  check_flag(warmup, "warmup")
  check_flag(progress, "progress")
  if (!is.null(seed)) check_scalar_num(seed, "seed", integer = TRUE)

  n_iterations <- as.integer(n_iterations)
  n_cores <- as.integer(n_cores)

  # Capture dots once, as a plain list. Closures sent to PSOCK workers must not
  # depend on a `...` binding in an enclosing frame, which does not serialise
  # reliably.
  dots <- list(...)

  cl <- NULL
  if (parallel && n_cores > 1L) {
    cl <- parallel::makePSOCKcluster(n_cores)
    on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
    parallel::clusterEvalQ(cl, {
      suppressPackageStartupMessages(requireNamespace("multihostSIR",
                                                      quietly = TRUE))
      NULL
    })
  }

  run_batch <- function(richness, k) {
    args <- c(list(richness = richness, spec = spec, allo = allo,
                   coefficients = coefficients, duration = duration), dots)
    fn <- function(i) do.call(multihostSIR::simulate_once, args)
    if (!is.null(cl)) {
      parallel::parLapply(cl, seq_len(k), fn)
    } else {
      lapply(seq_len(k), fn)
    }
  }

  # Untimed warm-up so the first level is not penalized by start-up costs.
  if (warmup) invisible(run_batch(min(richness_levels), 2L))

  # Seed after the warm-up so the recorded results are reproducible.
  if (!is.null(seed)) {
    if (!is.null(cl)) {
      parallel::clusterSetRNGStream(cl, iseed = as.integer(seed))
    } else {
      set.seed(as.integer(seed))
    }
  }

  results_list <- vector("list", length(richness_levels))
  traj_list    <- vector("list", length(richness_levels))
  timing_list  <- vector("list", length(richness_levels))

  for (k in seq_along(richness_levels)) {
    richness <- richness_levels[k]
    elapsed <- system.time(batch <- run_batch(richness, n_iterations))

    rows <- do.call(rbind, lapply(batch, `[[`, "row"))
    rows <- cbind(iteration = seq_len(n_iterations), rows)
    prev_mat <- vapply(batch, `[[`, numeric(duration + 1L), "prevalence")

    results_list[[k]] <- rows
    traj_list[[k]] <- data.frame(
      richness = richness,
      time = batch[[1L]]$time,
      mean_prevalence = rowMeans(prev_mat),
      stringsAsFactors = FALSE)
    timing_list[[k]] <- data.frame(
      richness = richness,
      n_iterations = n_iterations,
      state_dim = 4L * richness,
      elapsed_sec = as.numeric(elapsed[["elapsed"]]),
      sec_per_iter = as.numeric(elapsed[["elapsed"]]) / n_iterations,
      stringsAsFactors = FALSE)

    if (progress) {
      message(sprintf("  richness %2d: %8.2f s total (%.5f s/iter)",
                      richness, elapsed[["elapsed"]],
                      elapsed[["elapsed"]] / n_iterations))
    }
  }

  results <- do.call(rbind, results_list)
  rownames(results) <- NULL

  structure(
    list(results = results,
         trajectories = do.call(rbind, traj_list),
         timing = do.call(rbind, timing_list),
         settings = list(richness_levels = richness_levels,
                         n_iterations = n_iterations,
                         duration = duration,
                         spec = spec, allometry = allo,
                         coefficients = coefficients,
                         seed = seed, parallel = parallel,
                         n_cores = if (parallel) n_cores else 1L),
         call = cl_call),
    class = "multihost_sir")
}

#' @param x A `"multihost_sir"` object.
#' @param ... Ignored.
#' @rdname run_multihost_sir
#' @export
print.multihost_sir <- function(x, ...) {
  s <- x$settings
  cat("<multihost_sir>\n")
  cat(sprintf("  richness levels : %s\n",
              paste(range(s$richness_levels), collapse = " to ")))
  cat(sprintf("  iterations/level: %d\n", s$n_iterations))
  cat(sprintf("  horizon         : %d days\n", s$duration))
  cat(sprintf("  total runs      : %d\n", nrow(x$results)))
  cat(sprintf("  total time      : %.1f s\n", sum(x$timing$elapsed_sec)))
  cat(sprintf("  median R0       : %.3f (%.1f%% above 1)\n",
              median(x$results$R0), 100 * mean(x$results$R0 > 1)))
  cat("\n")
  cat(describe_assumptions(allo = s$allometry, spec = s$spec,
                           coefficients = s$coefficients), sep = "\n")
  cat("\nUse summary() for a table by richness, plot() for figures.\n")
  invisible(x)
}
