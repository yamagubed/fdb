test_that("fit_internal_only and fit_naive_pooled return basic structure", {
  set.seed(11)
  sim <- simulate_hybrid_cox(nI1 = 60, nI0 = 60, nE = 100,
                             theta0 = log(0.8), delta0 = 0, p = 3)
  fit_int <- fit_internal_only(sim$data)
  fit_pool <- fit_naive_pooled(sim$data)
  for (fit in list(fit_int, fit_pool)) {
    expect_true(is.finite(fit$theta_hat))
    expect_true(is.finite(fit$se_theta))
    expect_true(is.finite(fit$z))
    expect_true(fit$se_theta > 0)
  }
})

test_that("All penalized methods run and return required fields", {
  set.seed(12)
  sim <- simulate_hybrid_cox(nI1 = 80, nI0 = 80, nE = 150,
                             theta0 = log(0.8), delta0 = 0, p = 3)
  required <- c("method", "theta_hat", "se_theta", "se_sand",
                "z", "z_sand", "delta_hat", "pen_curv")

  fit_li <- fit_li_adaptive_lasso(sim$data, lambda = 0.2)
  expect_true(all(required %in% names(fit_li)))
  expect_equal(fit_li$method, "LiAdaptiveLasso")

  fit_p1 <- fit_P1_precision_L1(sim$data, lambda = 0.2)
  expect_true(all(required %in% names(fit_p1)))
  expect_equal(fit_p1$method, "P1_SEScaledL1")

  fit_p2 <- fit_P2_gated_L1(sim$data, lambda = 0.2)
  expect_true(all(required %in% names(fit_p2)))
  expect_equal(fit_p2$method, "P2_GatedL1")
  expect_true(fit_p2$g_hat >= 0 && fit_p2$g_hat <= 1)

  fit_p3 <- fit_P3_info_MCP(sim$data, lambda = 0.2)
  expect_true(all(required %in% names(fit_p3)))
  expect_equal(fit_p3$method, "P3_SEScaledMCP")

  fit_p4 <- fit_P4_LRweighted_L1(sim$data, lambda = 0.2)
  expect_true(all(required %in% names(fit_p4)))
  expect_equal(fit_p4$method, "P4_LRWeightedL1")
  expect_true(is.finite(fit_p4$LR0))
})

test_that("fit_one_penalized_method dispatches correctly", {
  set.seed(13)
  sim <- simulate_hybrid_cox(nI1 = 60, nI0 = 60, nE = 120,
                             theta0 = log(0.8), delta0 = 0, p = 3)
  for (m in c("Li", "P1", "P2", "P3", "P4")) {
    f <- fit_one_penalized_method(sim$data, method = m, lambda = 0.2)
    expect_true(is.finite(f$theta_hat))
    expect_true(is.finite(f$se_theta))
  }
})

test_that("fit_all_methods returns a 7-row data frame with required columns", {
  set.seed(14)
  sim <- simulate_hybrid_cox(nI1 = 60, nI0 = 60, nE = 120,
                             theta0 = log(0.8), delta0 = 0, p = 3)
  res <- fit_all_methods(sim$data)
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 7)
  expect_true(all(c("method", "theta_hat", "se_theta", "se_sand",
                    "z", "z_sand", "delta_hat",
                    "pen_curv", "pen_curv_raw") %in% names(res)))
  expect_true(all(is.finite(res$theta_hat)))
})

test_that("As lambda -> infinity, P1 shrinks delta toward zero", {
  set.seed(15)
  sim <- simulate_hybrid_cox(nI1 = 80, nI0 = 80, nE = 150,
                             theta0 = log(0.8), delta0 = log(1.5), p = 3)
  fit_int <- fit_internal_only(sim$data)
  fit_p1_big <- fit_P1_precision_L1(sim$data, lambda = 50)
  # With very strong shrinkage, delta_hat should collapse near 0.
  expect_lt(abs(fit_p1_big$delta_hat), 0.05)
  # The theta_hat then approaches the internal-only estimator, but not
  # exactly: at delta = 0 the external controls are still in the risk
  # set, so finite-sample bias prevents exact equality. A 0.3 tolerance
  # is comfortable for the default sample sizes used here.
  expect_lt(abs(fit_p1_big$theta_hat - fit_int$theta_hat), 0.3)
})
