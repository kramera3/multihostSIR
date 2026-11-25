res <- run_multihost_sir(3:5, n_iterations = 4, duration = 80, seed = 12,
                         warmup = FALSE)

test_that("plotting functions return ggplot objects", {
  expect_s3_class(plot_r0_richness(res), "ggplot")
  expect_s3_class(plot_peak_prevalence(res), "ggplot")
  expect_s3_class(plot_trajectories(res), "ggplot")
  expect_s3_class(plot_timing(res), "ggplot")
  expect_s3_class(plot_timing(res, log_scale = TRUE), "ggplot")
  expect_s3_class(plot_allometry(), "ggplot")
})

test_that("the plot method dispatches on `which`", {
  expect_s3_class(plot(res, which = "r0"), "ggplot")
  expect_s3_class(plot(res, which = "trajectories"), "ggplot")
  expect_error(plot(res, which = "nonsense"))
})

test_that("the transmission heat map builds for a mixed community", {
  set.seed(31)
  comm <- build_community(6)
  expect_s3_class(plot_transmission_matrix(comm), "ggplot")
  expect_s3_class(plot_transmission_matrix(comm, show_values = TRUE), "ggplot")
})

test_that("plots actually render without error", {
  skip_on_cran()
  p <- plot_r0_richness(res)
  tmp <- tempfile(fileext = ".png")
  on.exit(unlink(tmp), add = TRUE)
  ggplot2::ggsave(tmp, p, width = 6, height = 4, dpi = 72)
  expect_true(file.exists(tmp))
})

test_that("richness levels are thinned for legibility", {
  levels_out <- multihostSIR:::display_levels(3:60, max_levels = 8L)
  expect_lte(length(levels_out), 8L)
  expect_equal(min(levels_out), 3)
  expect_equal(max(levels_out), 60)
})

test_that("plot helpers reject non-result objects", {
  expect_error(plot_r0_richness(mtcars), "run_multihost_sir")
  expect_error(plot_allometry(list()), "allometry")
})
