make_toy <- function() {
  comm <- data.frame(
    species = 1:4,
    type = c("Prey", "Prey", "Predator", "Predator"),
    mass = c(0.02, 2, 1, 25),
    trophic_weight = c(1, 1, 2, 2),
    stringsAsFactors = FALSE)
  rates <- multihostSIR:::apply_allometry(comm$mass, allometry())
  comm <- cbind(comm, rates)
  comm$competence <- c(0.2, 0.4, 0.3, 0.5)
  comm$beta_ii <- comm$competence / comm$N
  comm
}

test_that("the transmission matrix has the right shape and diagonal", {
  comm <- make_toy()
  beta <- transmission_matrix(comm)
  expect_true(is.matrix(beta))
  expect_equal(dim(beta), c(4L, 4L))
  expect_equal(diag(beta), comm$beta_ii, ignore_attr = TRUE)
})

test_that("predator-to-prey transmission is blocked by default", {
  comm <- make_toy()
  beta <- transmission_matrix(comm)
  prey <- comm$type == "Prey"
  pred <- comm$type == "Predator"
  # Rows are recipients, columns are sources.
  expect_true(all(beta[prey, pred] == 0))
  expect_true(all(beta[pred, prey] > 0))
})

test_that("off-diagonal entries are scaled geometric means", {
  comm <- make_toy()
  co <- trophic_coefficients()
  beta <- transmission_matrix(comm, co)
  expect_equal(beta[1, 2],
               co$prey_prey * sqrt(comm$beta_ii[1] * comm$beta_ii[2]))
  expect_equal(beta[3, 1],
               co$prey_to_pred * sqrt(comm$beta_ii[3] * comm$beta_ii[1]))
  expect_equal(beta[3, 4],
               co$pred_pred * sqrt(comm$beta_ii[3] * comm$beta_ii[4]))
})

test_that("the matrix is symmetric when all coefficients are equal", {
  comm <- make_toy()
  co <- trophic_coefficients(0.4, 0.4, 0.4, 0.4)
  beta <- transmission_matrix(comm, co)
  expect_equal(beta, t(beta))
})

test_that("R0 matches the single-species analytic value", {
  comm <- make_toy()[1, ]
  comm$species <- 1L
  r0 <- community_R0(comm)
  expected <- comm$beta_ii * comm$N / (comm$gamma + comm$mu + comm$alpha)
  expect_equal(r0, expected, ignore_attr = TRUE)
})

test_that("R0 is non-negative and scales with density", {
  set.seed(11)
  comm <- build_community(5)
  r0 <- community_R0(comm)
  expect_true(is.finite(r0))
  expect_gte(r0, 0)

  doubled <- comm
  doubled$N <- doubled$N * 2
  doubled$beta_ii <- doubled$competence / doubled$N
  # beta_ii falls as 1/N while K carries a factor N, so R0 is unchanged under
  # this frequency-dependent parameterisation.
  expect_equal(community_R0(doubled), r0)
})

test_that("the next-generation matrix has the documented structure", {
  comm <- make_toy()
  beta <- transmission_matrix(comm)
  K <- next_generation_matrix(comm, beta)
  duration <- 1 / (comm$gamma + comm$mu + comm$alpha)
  expect_equal(K[2, 3], beta[2, 3] * comm$N[2] * duration[3])
})

test_that("malformed input is rejected", {
  comm <- make_toy()
  expect_error(transmission_matrix(comm[, c("type", "mass")]), "missing")
  expect_error(next_generation_matrix(comm, diag(3)), "4 x 4")
})
