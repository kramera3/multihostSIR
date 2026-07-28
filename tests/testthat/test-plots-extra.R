test_that("extra richness plots return ggplot objects", {
  res <- run_multihost_sir(3:6, n_iterations = 6, duration = 100, seed = 3,
                           warmup = FALSE)
  expect_s3_class(plot_cumulative_richness(res), "ggplot")
  expect_s3_class(plot_days_to_peak_richness(res), "ggplot")
})

test_that("macro-driver scatter plots return ggplot objects", {
  res <- run_multihost_sir(3:6, n_iterations = 6, duration = 100, seed = 3,
                           warmup = FALSE)
  expect_s3_class(plot_r0_trophic(res), "ggplot")
  expect_s3_class(plot_r0_body_mass(res), "ggplot")
  expect_s3_class(plot_r0_density(res), "ggplot")
  expect_s3_class(plot_r0_competence(res), "ggplot")
  expect_s3_class(plot_cumulative_competence(res), "ggplot")
  expect_s3_class(plot_peak_density(res), "ggplot")
  expect_s3_class(plot_days_to_peak_r0(res), "ggplot")
})

test_that("predator-structure plots return ggplot objects", {
  res <- run_multihost_sir(3:6, n_iterations = 6, duration = 100, seed = 3,
                           warmup = FALSE)
  expect_s3_class(plot_peak_predator(res), "ggplot")
  expect_s3_class(plot_r0_predator(res), "ggplot")
})

test_that("plot_gallery returns the full named list of ggplots", {
  res <- run_multihost_sir(3:6, n_iterations = 6, duration = 100, seed = 3,
                           warmup = FALSE)
  g <- plot_gallery(res)
  expect_type(g, "list")
  expect_true(all(vapply(g, inherits, logical(1), what = "ggplot")))
  expect_true(all(c("r0_richness", "r0_trophic", "peak_predator",
                    "timing", "allometry") %in% names(g)))
  expect_length(g, 16L)
})

test_that("the extra plots reject a non-result object", {
  expect_error(plot_r0_trophic(mtcars), "run_multihost_sir")
  expect_error(plot_gallery(mtcars), "run_multihost_sir")
})

test_that("the transmission heatmap shows actual values by default", {
  set.seed(1)
  comm <- build_community(4)
  p <- plot_transmission_matrix(comm)
  # a geom_text layer is present when values are shown
  has_text_layer <- any(vapply(p$layers, function(l)
    inherits(l$geom, "GeomText"), logical(1)))
  expect_true(has_text_layer)

  # and the default is not the relative/percentage scale
  p_nolabel <- plot_transmission_matrix(comm, show_values = FALSE)
  expect_s3_class(p_nolabel, "ggplot")
})
