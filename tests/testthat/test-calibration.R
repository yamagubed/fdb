test_that("make_drift_set returns expected values", {
  d <- make_drift_set(drift_hr_range = c(0.8, 1.2), by_hr = 0.1,
                      include_zero = TRUE)
  expect_true(0 %in% d)
  expect_true(all(is.finite(d)))
  expect_true(any(abs(d - log(0.8)) < 1e-8))
  expect_true(any(abs(d - log(1.2)) < 1e-8))
})

test_that("make_drift_set_from_values handles HR values correctly", {
  d <- make_drift_set_from_values(c(0.9, 1.0, 1.1), include_zero = TRUE)
  expect_true(0 %in% d)
  expect_true(any(abs(d - log(0.9)) < 1e-8))
})

test_that("calibrate_lambda_grid runs end-to-end on a tiny example", {
  skip_on_cran()
  small_scenario <- list(
    nI1 = 40, nI0 = 40, nE = 80,
    theta0 = 0, delta0 = 0,
    p = 2, beta = c(0.2, 0.2),
    rho = 0, cov_shift = c(0, 0),
    shape = 1.2, lambda = 0.02,
    target_cens = 0.2
  )
  cal <- calibrate_lambda_grid(
    method = "P1",
    lambda_grid = c(0.1, 0.5),
    scenario_base = small_scenario,
    drift_set = c(0, log(1.1)),
    nsim = 10,
    alpha = 0.025,
    seed = 1
  )
  expect_named(cal,
               c("method", "details", "summary", "calibration_table",
                 "lambda_star"))
  expect_s3_class(cal$summary, "data.frame")
  expect_s3_class(cal$details, "data.frame")
  expect_equal(nrow(cal$lambda_star), 2)
  expect_true(all(cal$lambda_star$inference %in%
                  c("model_based", "sandwich")))
})
