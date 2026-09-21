test_that("calibration preserves errors and non-finite statistics", {
  fun <- .calibration_rejects_for_dataset
  env <- new.env(parent = environment(fun))
  environment(fun) <- env
  run <- function() fun(data.frame(), "Li", c(0.1, 0.2), -1.96,
                       1, 1.64, 0.25, 3, c(-2, 2), FALSE, 0.001, 21)
  env$cox_fit_full <- function(...) stop("deliberate Cox failure")
  out <- run()
  expect_true(all(is.na(out$reject_mod)))
  expect_true(all(grepl("deliberate Cox failure", out$fit_error)))
  env$cox_fit_full <- function(...) list()
  env$fit_one_penalized_method_cached <- function(...) stop("deliberate penalty failure")
  expect_true(all(grepl("deliberate penalty failure", run()$fit_error)))
  env$fit_one_penalized_method_cached <- function(...) list(z = NA_real_, z_sand = -3)
  out <- run()
  expect_true(all(is.na(out$reject_mod)))
  expect_equal(out$reject_sand, c(1, 1))
})

test_that("calibration resolves caller settings before replicate evaluation", {
  fun <- .evaluate_calibration_drift_cached
  env <- new.env(parent = environment(fun))
  environment(fun) <- env
  caller <- new.env(parent = env)
  caller$fun <- fun
  caller$CONFIG <- list(eps = 0.001, gamma = 1, rho = 0.1)
  env$simulate_hybrid_cox <- function(...) {
    # Model a worker where the caller's global configuration is absent.
    rm("CONFIG", envir = caller)
    list(data = data.frame())
  }
  env$.calibration_rejects_for_dataset <- function(...) {
    args <- list(...)
    stopifnot(args$eps == 0.001, args$gamma_li == 1, args$rho_mcp == 0.1)
    data.frame(lambda = 0.1, reject_mod = 0, reject_sand = 0,
               fit_error = NA_character_)
  }
  out <- evalq(fun("Li", 0.1, list(), 0, 1, 0.025, 1, FALSE, NULL,
                  FALSE, CONFIG$eps, CONFIG$gamma, 1.64, 0.25, 3,
                  c(-2, 2), 21, rho_mcp = CONFIG$rho), caller)
  expect_equal(out$n_nonmissing_model_based, 1)
})

test_that("all failed calibration stops with the original error", {
  fun <- .evaluate_calibration_drift_cached
  env <- new.env(parent = environment(fun))
  environment(fun) <- env
  env$simulate_hybrid_cox <- function(...) list(data = data.frame())
  env$.calibration_rejects_for_dataset <- function(...) {
    data.frame(lambda = 0.1, reject_mod = NA_real_, reject_sand = NA_real_,
               fit_error = "original fitting error")
  }
  expect_error(fun("Li", 0.1, list(), 0, 2, 0.025, 1, FALSE, NULL,
                   FALSE, 0.001, 1, 1.64, 0.25, 3, c(-2, 2), 21),
               "original fitting error")
})

test_that("missing or incomplete drift results cannot qualify a lambda", {
  details <- data.frame(method = "Li", lambda = 0.1, delta0 = c(0, 0.1),
                        type1_model_based = c(0.01, NA_real_),
                        type1_sandwich = c(0.01, 0.01),
                        n_nonmissing_model_based = c(10, 0),
                        n_nonmissing_sandwich = c(10, 9))
  out <- .summarize_calibration_details(details, 10)
  expect_true(is.na(out$worst_type1_model_based))
  expect_true(is.na(out$worst_type1_sandwich))
  expect_true(is.na(.select_lambda_star(out, 0.04)))
  details$n_nonmissing_sandwich <- c(10, 10)
  expect_equal(.summarize_calibration_details(details, 10)$worst_type1_sandwich, 0.01)
})
