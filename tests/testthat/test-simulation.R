test_that("simulate_epidemic returns a well-formed solution", {
  set.seed(21)
  comm <- build_community(4)
  sim <- simulate_epidemic(comm, duration = 100)

  expect_s3_class(sim, "multihost_epidemic")
  expect_equal(nrow(sim$out), 101L)
  expect_true(all(c("time", "S1", "I1", "R1", "C1", "C4") %in% names(sim$out)))
  expect_equal(sim$out$time, 0:100)
})

test_that("state variables stay non-negative", {
  set.seed(22)
  sim <- simulate_epidemic(build_community(5), duration = 150)
  state_cols <- setdiff(names(sim$out), "time")
  expect_true(all(as.matrix(sim$out[, state_cols]) >= -1e-6))
})

test_that("cumulative incidence is non-decreasing", {
  set.seed(23)
  comm <- build_community(3)
  sim <- simulate_epidemic(comm, duration = 120)
  for (i in seq_len(nrow(comm))) {
    cvec <- sim$out[[paste0("C", i)]]
    expect_true(all(diff(cvec) >= -1e-6))
  }
})

test_that("the epidemic is seeded in the most competent species", {
  set.seed(24)
  comm <- build_community(5)
  sim <- simulate_epidemic(comm, duration = 50)
  expect_equal(sim$seed_species, which.max(comm$competence))

  sim2 <- simulate_epidemic(comm, duration = 50, seed_species = 2)
  expect_equal(sim2$seed_species, 2L)
})

test_that("with no infection seeded nothing happens", {
  set.seed(25)
  comm <- build_community(4)
  sim <- simulate_epidemic(comm, duration = 60, seed_size = 0)
  metrics <- species_metrics(sim)
  expect_true(all(metrics$peak_infected < 1e-6))
})

test_that("species_metrics has one row per species", {
  set.seed(26)
  comm <- build_community(6)
  sim <- simulate_epidemic(comm, duration = 80)
  m <- species_metrics(sim)
  expect_equal(nrow(m), 6L)
  expect_true(all(m$days_to_peak >= 0 & m$days_to_peak <= 80))
  expect_true(all(m$attack_rate >= 0))
})

test_that("simulate_once produces one summary row and a prevalence series", {
  set.seed(27)
  one <- simulate_once(4, duration = 90)
  expect_equal(nrow(one$row), 1L)
  expect_equal(length(one$prevalence), 91L)
  expect_equal(one$row$richness, 4L)
  expect_true(is.finite(one$row$R0))
})

test_that("run_multihost_sir assembles results, trajectories and timings", {
  res <- run_multihost_sir(richness_levels = 2:4, n_iterations = 3,
                           duration = 60, seed = 99, warmup = FALSE)

  expect_s3_class(res, "multihost_sir")
  expect_equal(nrow(res$results), 9L)
  expect_equal(sort(unique(res$results$richness)), 2:4)
  expect_equal(nrow(res$timing), 3L)
  expect_equal(nrow(res$trajectories), 3L * 61L)
  expect_true(all(res$timing$elapsed_sec >= 0))
})

test_that("a fixed seed reproduces the run exactly", {
  a <- run_multihost_sir(3:4, n_iterations = 4, duration = 60, seed = 7,
                         warmup = FALSE)
  b <- run_multihost_sir(3:4, n_iterations = 4, duration = 60, seed = 7,
                         warmup = FALSE)
  expect_equal(a$results$R0, b$results$R0)
  expect_equal(a$results$peak_prevalence, b$results$peak_prevalence)
})

test_that("summary and helper tables behave", {
  res <- run_multihost_sir(3:5, n_iterations = 4, duration = 60, seed = 5,
                           warmup = FALSE)
  s <- summary(res)
  expect_s3_class(s, "summary.multihost_sir")
  expect_equal(nrow(s$by_richness), 3L)
  expect_true(all(s$by_richness$n_sims == 4L))

  expect_s3_class(predator_summary(res), "data.frame")
  expect_s3_class(competence_summary(res), "data.frame")
  expect_true(is.numeric(timing_complexity(res)$exponent_richness))
})

test_that("invalid driver arguments are caught", {
  expect_error(run_multihost_sir(richness_levels = c(2, 2.5)), "whole numbers")
  expect_error(run_multihost_sir(2:3, n_iterations = 0), "n_iterations")
  expect_error(run_multihost_sir(2:3, n_iterations = 2, duration = 0),
               "duration")
})

test_that("results can be exported to CSV", {
  res <- run_multihost_sir(3:4, n_iterations = 2, duration = 50, seed = 3,
                           warmup = FALSE)
  dir <- file.path(tempdir(), "multihost_export_test")
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  paths <- export_results(res, dir = dir, formats = "csv")
  expect_true(all(file.exists(paths)))
  expect_true(any(grepl("results\\.csv$", paths)))
})
