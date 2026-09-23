test_that("installed PSOCK workers evaluate caller settings reproducibly", {
  skip_if(Sys.getenv("FDB_RUN_PARALLEL_TESTS") != "true",
          "Opt-in native PSOCK integration test")
  cfg <- list(eps = 0.001, gamma = 1, rho = 0.1)
  sc <- scenario_S1
  sc$nI1 <- 40; sc$nI0 <- 40; sc$nE <- 80
  sc$p <- 0; sc$beta <- numeric(); sc$cov_shift <- numeric()
  run <- function() calibrate_lambda_grid("Li", c(0.02, 0.2), sc, 0,
    nsim = 2, parallel = TRUE, ncores = 2, seed = 51,
    eps = cfg$eps, gamma_li = cfg$gamma, rho_mcp = cfg$rho)
  a <- run(); b <- run()
  expect_identical(a$details, b$details)
  expect_true(all(a$details$n_nonmissing_sandwich == 2))
})
