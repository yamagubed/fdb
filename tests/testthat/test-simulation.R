test_that("run_simulation returns expected structure", {
  skip_on_cran()
  small_scenario <- list(
    nI1 = 40, nI0 = 40, nE = 80,
    theta0 = log(0.8), delta0 = 0,
    p = 2, beta = c(0.2, 0.2),
    rho = 0, cov_shift = c(0, 0),
    shape = 1.2, lambda = 0.02,
    target_cens = 0.2
  )
  sim_out <- run_simulation(nsim = 5, scenario = small_scenario,
                            lambdas = lambdas_default,
                            alpha = 0.025, seed = 1)
  expect_named(sim_out,
               c("raw", "summary", "scenario", "lambdas", "settings"))
  expect_s3_class(sim_out$raw, "data.frame")
  expect_s3_class(sim_out$summary, "data.frame")
  expect_true("inference" %in% names(sim_out$summary))
  expect_true(all(c("model_based", "sandwich") %in%
                  sim_out$summary$inference))
})

test_that("compute_ess_from_raw produces nonnegative ESS for InternalOnly", {
  skip_on_cran()
  small_scenario <- list(
    nI1 = 40, nI0 = 40, nE = 80,
    theta0 = log(0.8), delta0 = 0,
    p = 2, beta = c(0.2, 0.2),
    rho = 0, cov_shift = c(0, 0),
    shape = 1.2, lambda = 0.02,
    target_cens = 0.2
  )
  sim_out <- run_simulation(nsim = 10, scenario = small_scenario,
                            lambdas = lambdas_default,
                            alpha = 0.025, seed = 2)
  ess <- compute_ess_from_raw(sim_out$raw, NS = 80,
                              ref_method = "InternalOnly")
  expect_s3_class(ess, "data.frame")
  expect_true("ESS" %in% names(ess))
  expect_equal(ess$ESS[ess$method == "InternalOnly"], 0)
})
