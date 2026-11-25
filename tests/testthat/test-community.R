test_that("build_community returns the requested number of species", {
  set.seed(1)
  comm <- build_community(6)
  expect_s3_class(comm, "data.frame")
  expect_equal(nrow(comm), 6L)
  expect_true(all(c("type", "mass", "mu", "gamma", "alpha", "N", "beta_ii") %in%
                    names(comm)))
})

test_that("allometric rates follow the specified power laws", {
  allo <- allometry()
  set.seed(2)
  comm <- build_community(5, allo = allo)
  expect_equal(comm$mu, allo$a_mort * comm$mass^allo$b_mort)
  expect_equal(comm$gamma, allo$a_recovery * comm$mass^allo$b_recovery)
  expect_equal(comm$alpha, allo$a_virulence * comm$mass^allo$b_virulence)
})

test_that("predator densities are thinned by predator_density_factor", {
  spec <- community_spec(prey_prob = 0, predator_density_factor = 0.1)
  allo <- allometry()
  set.seed(3)
  comm <- build_community(4, spec = spec, allo = allo)
  expect_true(all(comm$type == "Predator"))
  expect_equal(comm$N, allo$a_density * comm$mass^allo$b_density * 0.1)
})

test_that("larger animals are rarer and longer-lived", {
  allo <- allometry()
  small <- apply_allometry(0.02, allo)
  large <- apply_allometry(50, allo)
  expect_gt(small$N, large$N)      # Damuth's law
  expect_gt(small$mu, large$mu)    # quarter-power lifespan scaling
})

test_that("beta_ii is competence divided by density", {
  set.seed(4)
  comm <- build_community(5)
  expect_equal(comm$beta_ii, comm$competence / comm$N)
})

test_that("competence stays inside its category range", {
  spec <- community_spec()
  set.seed(5)
  comm <- build_community(30, spec = spec)
  for (i in seq_len(nrow(comm))) {
    rng <- spec$competence_ranges[[as.character(comm$competence_category[i])]]
    expect_gte(comm$competence[i], rng[1])
    expect_lte(comm$competence[i], rng[2])
  }
})

test_that("invalid inputs are rejected", {
  expect_error(build_community(0), "richness")
  expect_error(build_community(2.5), "whole number")
  expect_error(community_spec(prey_prob = 1.5), "prey_prob")
  expect_error(community_spec(prey_masses = c(-1, 2)), "prey_masses")
  expect_error(allometry(a_mort = -1), "a_mort")
})
