# Deterministic stage results reproduce the failed P4 refinement without
# relying on a particular random simulation sample.
run_search_fixture <- function(fine_all_fail = FALSE, confirm_all_fail = FALSE,
                               confirm = TRUE) {
  fun <- calibrate_lambda_grid_two_stage
  env <- new.env(parent = environment(fun))
  environment(fun) <- env
  calls <- list()
  env$calibrate_lambda_grid <- function(method, lambda_grid, ...) {
    stage <- length(calls) + 1L
    calls[[stage]] <<- lambda_grid
    limit <- c(2.272, 0.5, 0.1)[stage]
    rate <- ifelse(lambda_grid <= limit, 0.03, 0.05)
    if ((stage == 2 && fine_all_fail) || (stage == 3 && confirm_all_fail)) {
      rate[] <- 0.05
    }
    summary <- data.frame(method = method, lambda = lambda_grid,
                          worst_type1_model_based = 0.08,
                          worst_type1_sandwich = rate)
    stars <- data.frame(method = method,
                        inference = c("model_based", "sandwich"),
                        lambda_star = c(NA_real_,
                          .select_lambda_star(summary, 0.04, "sandwich", "point")))
    list(summary = summary, calibration_table = summary, lambda_star = stars)
  }
  coarse <- c(0.02, 0.044, 0.097, 0.213, 0.469, 1.032, 2.272, 5)
  result <- fun(method = "P4", lambda_grid_coarse = coarse,
                scenario_base = list(), drift_set_cal = c(0, log(1.1)),
                primary_inference = "sandwich", confirm_full_drift = confirm)
  list(result = result, calls = calls, coarse = coarse)
}

test_that("smaller coarse candidates survive refinement and confirmation", {
  x <- run_search_fixture()
  expect_true(all(x$coarse %in% x$calls[[2]]))
  expect_true(all(x$coarse[x$coarse <= 0.469] %in% x$calls[[3]]))
  expect_equal(x$result$fine$lambda_star$lambda_star[2], 0.469)
  expect_equal(x$result$lambda_star$lambda_star[2], 0.097)
  expect_identical(x$result$final_stage, x$result$confirm)
})

test_that("all failed fine candidates still receive independent confirmation", {
  x <- run_search_fixture(fine_all_fail = TRUE)
  expect_true(is.na(x$result$fine$lambda_star$lambda_star[2]))
  expect_equal(x$calls[[3]], x$calls[[2]])
  expect_equal(x$result$lambda_star$lambda_star[2], 0.097)
})

test_that("failed confirmation never falls back to a coarse-stage success", {
  x <- run_search_fixture(confirm_all_fail = TRUE)
  expect_true(is.finite(x$result$coarse$lambda_star$lambda_star[2]))
  expect_true(is.finite(x$result$fine$lambda_star$lambda_star[2]))
  expect_true(is.na(x$result$lambda_star$lambda_star[2]))
})

test_that("without confirmation selection uses the expanded fine grid", {
  x <- run_search_fixture(confirm = FALSE)
  expect_length(x$calls, 2)
  expect_null(x$result$confirm)
  expect_equal(x$result$lambda_star$lambda_star[2], 0.469)
})
