test_that("simulate_hybrid_cox returns expected structure", {
  set.seed(1)
  sim <- simulate_hybrid_cox(nI1 = 50, nI0 = 50, nE = 100,
                             theta0 = log(0.8), delta0 = 0, p = 3)
  expect_named(sim, c("data", "truth", "settings"))
  expect_true(all(c("time", "status", "T", "Z") %in% names(sim$data)))
  expect_equal(nrow(sim$data), 50 + 50 + 100)
  expect_equal(sum(sim$data$T), 50)
  expect_equal(sum(sim$data$Z), 100)
  expect_true(all(sim$data$time > 0))
  expect_true(all(sim$data$status %in% c(0L, 1L)))
})

test_that("simulate_hybrid_cox validates inputs", {
  expect_error(simulate_hybrid_cox(nI1 = 50, nI0 = 50, nE = 100,
                                   p = 3, beta = rep(0.2, 2)),
               "Length of beta")
  expect_error(simulate_hybrid_cox(nI1 = 0, nI0 = 50, nE = 100))
  expect_error(simulate_hybrid_cox(nI1 = 50, nI0 = 50, nE = 100,
                                   target_cens = 1.5))
})

test_that("simulate_hybrid_cox with nonzero cov_shift differs from zero shift", {
  set.seed(2)
  sim_a <- simulate_hybrid_cox(nI1 = 50, nI0 = 50, nE = 100, p = 3,
                               cov_shift = rep(0, 3))
  set.seed(2)
  sim_b <- simulate_hybrid_cox(nI1 = 50, nI0 = 50, nE = 100, p = 3,
                               cov_shift = c(0.5, 0.5, 0))
  # External-control covariate means should differ between the two
  mu_a <- mean(sim_a$data$X1[sim_a$data$Z == 1])
  mu_b <- mean(sim_b$data$X1[sim_b$data$Z == 1])
  expect_false(isTRUE(all.equal(mu_a, mu_b)))
})
