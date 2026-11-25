complete_species <- function() {
  data.frame(
    species    = c("mouse", "deer", "fox"),
    type       = c("Prey", "Prey", "Predator"),
    mass       = c(0.02, 50, 5),
    mu         = c(0.013, 0.002, 0.003),
    gamma      = c(0.27, 0.04, 0.07),
    alpha      = c(0.027, 0.004, 0.007),
    N          = c(9000, 27, 6),
    competence = c(0.5, 0.1, 0.2),
    stringsAsFactors = FALSE)
}

test_that("a complete species table reads as a complete community", {
  comm <- read_community(complete_species())
  expect_true(attr(comm, "complete"))
  expect_equal(nrow(comm), 3L)
  expect_true(all(c("type", "mass", "mu", "gamma", "alpha", "N",
                    "competence", "beta_ii") %in% names(comm)))
  expect_equal(comm$beta_ii, comm$competence / comm$N)
})

test_that("a partial table is marked incomplete and filled from allometry", {
  partial <- complete_species()[, c("species", "type", "mass")]
  comm <- read_community(partial)
  expect_false(attr(comm, "complete"))
  expect_true(length(attr(comm, "filled")) > 0)
  # filled rates should match the allometric relationships
  allo <- allometry()
  expect_equal(comm$mu, allo$a_mort * comm$mass^allo$b_mort)
})

test_that("bad type values are rejected", {
  bad <- complete_species()
  bad$type[1] <- "Plant"
  expect_error(read_community(bad), "Prey")
})

test_that("a missing type column is rejected", {
  no_type <- complete_species()
  no_type$type <- NULL
  expect_error(read_community(no_type), "type")
})

test_that("read_community reads a CSV file", {
  f <- tempfile(fileext = ".csv")
  on.exit(unlink(f), add = TRUE)
  write.csv(complete_species(), f, row.names = FALSE)
  comm <- read_community(f)
  expect_true(attr(comm, "complete"))
  expect_equal(nrow(comm), 3L)
})

test_that("the transmission matrix reader validates shape and values", {
  m <- matrix(c(0.05, 0.01, 0.00,
                0.01, 0.03, 0.00,
                0.02, 0.02, 0.04), nrow = 3, byrow = TRUE)
  beta <- read_transmission_matrix(m)
  expect_equal(dim(beta), c(3L, 3L))

  expect_error(read_transmission_matrix(matrix(1:6, nrow = 2)), "square")
  neg <- m; neg[1, 1] <- -1
  expect_error(read_transmission_matrix(neg), "non-negative")
})

test_that("matrix dimension is checked against the community", {
  comm <- read_community(complete_species())
  wrong <- matrix(0.1, nrow = 2, ncol = 2)
  expect_error(read_transmission_matrix(wrong, community = comm), "species")
})

test_that("a fully specified community runs deterministically", {
  comm <- read_community(complete_species())
  beta <- read_transmission_matrix(
    matrix(c(0.05, 0.01, 0.00,
             0.01, 0.03, 0.00,
             0.02, 0.02, 0.04), nrow = 3, byrow = TRUE),
    community = comm)

  fit <- run_fixed_community(comm, beta, duration = 100)
  expect_s3_class(fit, "multihost_fixed")
  expect_true(fit$assumptions$deterministic)
  expect_equal(nrow(fit$results), 1L)
  expect_true(is.finite(fit$R0))
  expect_equal(nrow(fit$species), 3L)
})

test_that("the deterministic run is reproducible without a seed", {
  comm <- read_community(complete_species())
  beta <- read_transmission_matrix(
    matrix(c(0.05, 0.01, 0.00,
             0.01, 0.03, 0.00,
             0.02, 0.02, 0.04), nrow = 3, byrow = TRUE),
    community = comm)
  a <- run_fixed_community(comm, beta, duration = 100)
  b <- run_fixed_community(comm, beta, duration = 100)
  expect_equal(a$results$R0, b$results$R0)
  expect_equal(a$results$peak_prevalence, b$results$peak_prevalence)
})

test_that("without a supplied matrix the run is not flagged deterministic", {
  comm <- read_community(complete_species())
  fit <- run_fixed_community(comm, beta = NULL, duration = 80)
  expect_false(fit$assumptions$deterministic)
  expect_s3_class(fit, "multihost_fixed")
})

test_that("print method shows the assumptions block", {
  comm <- read_community(complete_species())
  beta <- read_transmission_matrix(
    matrix(c(0.05, 0.01, 0.00,
             0.01, 0.03, 0.00,
             0.02, 0.02, 0.04), nrow = 3, byrow = TRUE),
    community = comm)
  fit <- run_fixed_community(comm, beta, duration = 60)
  expect_output(print(fit), "Assumptions used")
  expect_output(print(fit), "Fully specified")
})

test_that("fixed-run plots return ggplot objects", {
  comm <- read_community(complete_species())
  beta <- read_transmission_matrix(
    matrix(c(0.05, 0.01, 0.00,
             0.01, 0.03, 0.00,
             0.02, 0.02, 0.04), nrow = 3, byrow = TRUE),
    community = comm)
  fit <- run_fixed_community(comm, beta, duration = 100)

  expect_s3_class(plot_fixed_trajectory(fit), "ggplot")
  expect_s3_class(plot_fixed_trajectory(fit, by_species = TRUE), "ggplot")
  expect_s3_class(plot(fit), "ggplot")
  expect_s3_class(plot(fit, which = "species"), "ggplot")
  expect_s3_class(plot(fit, which = "matrix"), "ggplot")
})

test_that("species_curves has one column per species plus time", {
  comm <- read_community(complete_species())
  beta <- read_transmission_matrix(
    matrix(c(0.05, 0.01, 0.00,
             0.01, 0.03, 0.00,
             0.02, 0.02, 0.04), nrow = 3, byrow = TRUE),
    community = comm)
  fit <- run_fixed_community(comm, beta, duration = 80)

  expect_equal(ncol(fit$species_curves), 4L)   # time + 3 species
  expect_equal(nrow(fit$species_curves), 81L)  # 0..80
  expect_true(all(c("mouse", "deer", "fox") %in% names(fit$species_curves)))
})

test_that("export_results writes a fixed run to CSV", {
  comm <- read_community(complete_species())
  beta <- read_transmission_matrix(
    matrix(c(0.05, 0.01, 0.00,
             0.01, 0.03, 0.00,
             0.02, 0.02, 0.04), nrow = 3, byrow = TRUE),
    community = comm)
  fit <- run_fixed_community(comm, beta, duration = 60)

  dir <- file.path(tempdir(), "fixed_export_test")
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  paths <- export_results(fit, dir = dir, formats = "csv")
  expect_true(all(file.exists(paths)))
  expect_true(any(grepl("results\\.csv$", paths)))
  expect_true(any(grepl("species_curves\\.csv$", paths)))
})

test_that("plot_fixed_trajectory rejects the wrong object", {
  expect_error(plot_fixed_trajectory(mtcars), "run_fixed_community")
})
