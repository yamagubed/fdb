test_that("summaries retain missing methods and report valid denominators", {
  raw <- data.frame(method = c("A", "A", "B", "B"),
                    theta_hat = c(-0.2, 0.2, NA, NA),
                    se_theta = c(0.1, NA, NA, NA), se_sand = c(0.1, NA, NA, NA),
                    reject_mod = c(1, NA, NA, NA), reject_sand = c(1, NA, NA, NA))
  s <- .simulation_summary(raw, 0, "sandwich")
  expect_equal(s$n_estimates, c(2, 0))
  expect_equal(s$n_valid, c(1, 0))
  expect_equal(s$n_missing, c(1, 2))
  expect_equal(s$rej_rate[1], 1)
  expect_true(is.na(s$rej_rate[2]))
  expect_equal(s$coverage_95[1], 0)
  expect_true(is.na(s$coverage_95[2]))
})

test_that("drift curves optionally retain paired replicate rows", {
  fun <- run_drift_curve
  env <- new.env(parent = environment(fun)); environment(fun) <- env
  env$run_simulation <- function(scenario, ...) {
    raw <- data.frame(method = rep(c("InternalOnly", "P2_GatedL1"), each = 3),
                      sim = rep(1:3, 2), theta_hat = c(1:3, 2:4),
                      theta0 = scenario$theta0, delta0 = scenario$delta0)
    list(raw = raw, summary = data.frame(method = c("InternalOnly", "P2_GatedL1")))
  }
  sc <- list(nI1 = 10, nI0 = 10)
  a <- fun(0, c(0, log(1.1)), sc, list(), nsim = 3, seed = NULL)
  b <- fun(0, c(0, log(1.1)), sc, list(), nsim = 3, seed = NULL, keep_raw = TRUE)
  expect_identical(a, b$summary)
  expect_equal(nrow(b$raw), 12)
  expect_false(anyDuplicated(b$raw[c("method", "sim", "drift_index")]) > 0)
  expect_equal(sort(unique(b$raw$driftHR)), c(1, 1.1))
})

test_that("study wrapper refuses uncalibrated fallback tuning", {
  fun <- run_fdb_study
  env <- new.env(parent = environment(fun)); environment(fun) <- env
  env$calibrate_all_lambdas <- function(...) list(lambda_star = matrix(
    NA_real_, 5, 2, dimnames = list(c("Li", "P1", "P2", "P3", "P4"),
                                   c("model_based", "sandwich"))))
  expect_error(fun(), "No uncalibrated fallback")
})

test_that("confirmation honors an explicit drift set and supports NULL seeds", {
  fun <- calibrate_lambda_grid_two_stage
  env <- new.env(parent = environment(fun)); environment(fun) <- env
  calls <- list()
  env$calibrate_lambda_grid <- function(drift_set, seed, ...) {
    calls[[length(calls) + 1L]] <<- list(drift = drift_set, seed = seed)
    list(summary = data.frame(lambda = 0.1), calibration_table = data.frame(lambda = 0.1),
         lambda_star = data.frame(inference = c("model_based", "sandwich"),
                                  lambda_star = c(0.1, 0.1)))
  }
  fun("Li", c(0.1, 0.2), list(), 0, c(0, log(1.2)), seed = NULL)
  expect_equal(calls[[3]]$drift, c(0, log(1.2)))
  expect_true(all(vapply(calls, function(x) is.null(x$seed), logical(1))))
  calls <- list()
  fun("Li", c(0.1, 0.2), list(), log(1.1), seed = NULL)
  expect_equal(calls[[3]]$drift, log(1.1))
})

test_that("core defaults and count validation are explicit", {
  expect_equal(.resolve_ncores(NULL), 2L)
  expect_error(.resolve_ncores(0), "integer")
  expect_error(.resolve_ncores(1.5), "integer")
  expect_error(.validate_count(NA_real_, "nsim"), "integer")
  expect_null(.offset_seed(NULL, 100))
})

test_that("simulation metadata reports usable estimates", {
  sc <- scenario_S1
  sc$nI1 <- 30; sc$nI0 <- 30; sc$nE <- 60
  sc$p <- 0; sc$beta <- numeric(); sc$cov_shift <- numeric()
  out <- run_simulation(nsim = 2, scenario = sc, lambdas = lambdas_default, seed = 814)
  expect_true(all(out$summary$n_valid == 2))
  expect_equal(out$summary$mcse_rej_rate,
               sqrt(out$summary$rej_rate * (1 - out$summary$rej_rate) / 2))
})
