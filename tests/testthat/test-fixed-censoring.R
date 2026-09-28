test_that("population censoring matches an analytical exponential benchmark", {
  # With shape=1, T~Exp(lambda) and P(C<T)=r/(r+lambda).
  expect_equal(calibrate_censor_rate(p = 0, shape = 1, lambda = 0.02,
                                     target_cens = 0.4), 0.02 * 0.4 / 0.6,
               tolerance = 1e-8)
  r <- calibrate_censor_rate(300, 150, log(0.7), p = 0, shape = 1,
                             target_cens = 0.4)
  expect_equal((2/3) * r/(r + 0.014) + (1/3) * r/(r + 0.02), 0.4,
               tolerance = 1e-8)
})

test_that("censoring calibration is deterministic and does not consume RNG", {
  set.seed(551); before <- .Random.seed
  r <- calibrate_censor_rate(300, 150, p = 0, target_cens = 0.4)
  expect_identical(.Random.seed, before)
  set.seed(5)
  a <- simulate_hybrid_cox(300, 150, 450, p = 0, target_cens = 0.4)
  set.seed(99)
  b <- simulate_hybrid_cox(300, 150, 900, p = 0, delta0 = log(1.25), target_cens = 0.4)
  expect_identical(a$settings$rate_c, r)
  expect_identical(b$settings$rate_c, r)
})

test_that("explicit and zero censoring rates are honored", {
  a <- simulate_hybrid_cox(30, 30, 40, p = 0, target_cens = 0)
  expect_equal(a$settings$rate_c, 0)
  expect_true(all(a$data$status == 1L))
  b <- simulate_hybrid_cox(30, 30, 40, p = 0, censor_rate = 0.015)
  expect_equal(b$settings$rate_c, 0.015)
  for (r in list(-1, NA_real_, Inf, c(1, 2))) {
    expect_error(simulate_hybrid_cox(30, 30, 40, censor_rate = r), "censor_rate")
  }
})

test_that("large independent samples reproduce target censoring", {
  for (theta in c(0, log(0.7))) {
    set.seed(971)
    x <- simulate_hybrid_cox(40000, 20000, 100, theta0 = theta,
                             p = 0, target_cens = 0.4)
    observed <- mean(x$data$status[x$data$Z == 0] == 0)
    expect_lt(abs(observed - 0.4), 0.01)
  }
  set.seed(63)
  x <- simulate_hybrid_cox(30000, 30000, 100, theta0 = log(0.7),
                           p = 2, beta = c(0.4, -0.2), rho = 0.3,
                           target_cens = 0.4)
  expect_lt(abs(mean(x$data$status[x$data$Z == 0] == 0) - 0.4), 0.01)
})

test_that("drivers pass the fixed rate into calibration and evaluation replicates", {
  sc <- scenario_S1
  sc$p <- 0; sc$beta <- numeric(); sc$cov_shift <- numeric()
  sc$nI1 <- 40; sc$nI0 <- 40; sc$nE <- 60
  sc$censor_rate <- 0
  out <- run_simulation(2, sc, lambdas_default, seed = 121)
  expect_equal(out$scenario$censor_rate, 0)
  fun <- .evaluate_calibration_drift_cached
  env <- new.env(parent = environment(fun)); environment(fun) <- env
  original <- simulate_hybrid_cox
  passed <- numeric()
  env$simulate_hybrid_cox <- function(..., censor_rate) {
    passed <<- c(passed, censor_rate)
    original(..., censor_rate = censor_rate)
  }
  fun("Li", 0.1, sc, 0, 2, 0.025, 1, FALSE, 1, FALSE,
      0.001, 1, 1.64, 0.25, 3, c(-2, 2), 21)
  expect_identical(passed, c(0, 0))
})

test_that("labels match manuscript and legacy labels can be translated", {
  expected <- c("Internal-only", "Adaptive lasso", "Precision-weighted L1",
                "Integrated-gate", "Information-adaptive MCP", "LR-weighted L1",
                "Naive pooled")
  expect_identical(method_labels(), expected)
  expect_identical(method_labels(c("IO", "Li", "P1", "P2", "P3", "P4", "NP")), expected)
  expect_identical(method_labels(c("P2_GatedL1", "InternalOnly")), expected[c(4, 1)])
  set.seed(5)
  fits <- fit_all_methods(simulate_hybrid_cox(50, 50, 100, p = 0)$data)
  expect_setequal(fits$method, expected)
})
