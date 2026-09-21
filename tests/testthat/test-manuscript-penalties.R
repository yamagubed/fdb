test_that("integrated P2 matches quadrature, derivatives and monotonicity", {
  lambda <- 0.7; se <- 0.3; eps <- 0.001; c <- 1.64; tau <- 0.25
  p <- function(d) pen_integrated_gate(d, lambda, se, c, tau, eps)
  dp <- function(d) pen_derivative_gated_L1(d, lambda, se, c, tau, eps)
  for (d in c(0, 1e-5, 0.1, 0.5, 1, 3, -0.5)) {
    ref <- lambda * integrate(function(u) plogis((c-u/se)/tau),
                               eps, sqrt(d*d+eps*eps), rel.tol=1e-10)$value
    expect_equal(p(d), ref, tolerance=1e-9)
    step <- 1e-6
    expect_equal((p(d+step)-p(d-step))/(2*step), dp(d), tolerance=1e-5)
    expect_equal((dp(d+step)-dp(d-step))/(2*step),
                 pen_curvature_gated_L1(d, lambda, se, c, tau, eps),
                 tolerance=1e-5)
    expect_equal(p(-d), p(d))
  }
  expect_equal(p(0), 0)
  expect_true(all(diff(vapply(seq(0, 3, length.out=101), p, numeric(1))) >= -1e-14))
  expect_equal(pen_integrated_gate(0.8, 0, se), 0)
  expect_true(is.finite(pen_integrated_gate(1e6, 1, 0.001)))
  expect_true(is.finite(pen_integrated_gate(1, 1, 1, c=1e6, tau=1e-6)))
})

test_that("MCP matches the integrated slope and is C2 at both joins", {
  a <- 0.4; gamma <- 3; eps <- 0.001; rho <- 0.1
  b <- gamma*a; h <- rho*b
  p <- function(d) pen_MCP(d, a, gamma, eps, rho)
  dp <- function(d) pen_derivative_MCP(d, a, gamma, eps, rho)
  ddp <- function(d) pen_curvature_MCP(d, a, gamma, eps, rho)
  slope <- function(u) ifelse(u <= b-h, a-u/gamma,
                            ifelse(u < b+h, (b+h-u)^2/(4*gamma*h), 0))
  points <- c(0, 0.1, sqrt((b-h+eps)^2-eps^2),
              sqrt((b+eps)^2-eps^2), sqrt((b+h+eps)^2-eps^2), 2)
  for (d in c(points, -points)) {
    t <- sqrt(d*d+eps*eps)-eps
    cuts <- sort(unique(c(0, t, c(b-h,b+h)[c(b-h,b+h) < t])))
    ref <- if(t == 0) 0 else sum(vapply(seq_len(length(cuts)-1L),
      function(i) integrate(slope, cuts[i], cuts[i+1], rel.tol=1e-10)$value, numeric(1)))
    expect_equal(p(d), ref, tolerance=1e-9)
    step <- 1e-6
    expect_equal((p(d+step)-p(d-step))/(2*step), dp(d), tolerance=1e-5)
    expect_equal((dp(d+step)-dp(d-step))/(2*step), ddp(d), tolerance=1e-5)
    expect_equal(p(-d), p(d))
  }
  for (u in c(b-h, b+h)) {
    d <- sqrt((u+eps)^2-eps^2)
    expect_equal(ddp(d-1e-8), ddp(d+1e-8), tolerance=1e-7)
  }
  expect_equal(c(p(0), dp(0), ddp(0)), c(0, 0, a/eps))
  expect_equal(p(3), p(2))
  expect_equal(c(dp(3), ddp(3)), c(0, 0))
  expect_true(all(diff(vapply(seq(0,3,length.out=101), p, numeric(1))) >= -1e-14))
  for (d in c(0, 0.01, 1)) {
    expect_equal(pen_MCP(d, 0), 0)
    expect_equal(pen_curvature_MCP(d, 0, 3), 0)
    original <- if(abs(d)<=b) a*abs(d)-d*d/(2*gamma) else gamma*a*a/2
    small_eps <- 1e-8
    small_rho <- 1e-6
    actual <- pen_MCP(d, a, gamma, eps=small_eps, rho=small_rho)
    # Finite smoothing is not exact equality. Moving |d| by at most eps
    # changes MCP by at most a*eps; rounding the slope at b adds at most
    # h^2/(6*gamma). Test this absolute bound, not a platform-dependent
    # relative tolerance on penalty values near zero.
    bound <- a * small_eps + (small_rho * gamma * a)^2 / (6 * gamma)
    rounding <- 10 * .Machine$double.eps * max(1, abs(original))
    expect_lte(abs(actual - original), bound + rounding)
  }
})

test_that("no-covariate fits support the manuscript and zero penalty", {
  set.seed(208)
  dat <- simulate_hybrid_cox(nI1=80,nI0=60,nE=100,p=0)$data
  expect_named(dat, c("time","status","T","Z"))
  full <- cox_fit_full(dat, character())
  for (method in c("P2", "P3")) {
    f <- fit_one_penalized_method(dat, method, lambda=0)
    expect_equal(f$delta_hat, full$delta_hat, tolerance=1e-3)
    expect_equal(f$theta_hat, full$theta_hat, tolerance=1e-3)
    expect_equal(f$pen_curv_raw, 0)
    expect_true(is.finite(f$se_sand) && f$se_sand > 0)
  }
  all <- fit_all_methods(dat, rho_mcp=0.25)
  expect_true(all(is.finite(all$se_sand)))
  direct <- fit_P3_info_MCP(dat, lambda=0.2, rho_mcp=0.25)
  dispatched <- fit_one_penalized_method(dat, "P3", lambda=0.2, rho_mcp=0.25)
  expect_equal(direct, dispatched)
  expect_equal(direct$rho_mcp, 0.25)
  expect_equal(direct$h, 0.25*3*direct$lambda_eff)
  expect_equal(direct$pen_curv, max(direct$pen_curv_raw,0))
  expect_equal(all$theta_hat[all$method=="P3_SEScaledMCP"], direct$theta_hat)
  expect_error(fit_P3_info_MCP(dat, lambda=0.2, rho_mcp=1), "rho_mcp")
  expect_error(fit_P2_gated_L1(dat, lambda=0.2, tau=0), "tau")
})

test_that("sandwich uses score residuals at the fitted coefficients, including ties", {
  set.seed(209)
  dat <- simulate_hybrid_cox(nI1=80,nI0=60,nE=100,p=0)$data
  dat$time <- pmax(0.1, round(dat$time,1))
  fit <- fit_P3_info_MCP(dat, lambda=0.2, rho_mcp=0.2)
  eta <- c(fit$theta_hat,fit$delta_hat)
  ev <- survival::coxph(survival::Surv(time,status)~T+Z, data=dat,
                       ties="efron", init=eta,
                       control=survival::coxph.control(iter.max=0), x=TRUE)
  expect_equal(unname(coef(ev)), eta)
  B <- crossprod(residuals(ev,type="score"))
  A <- solve(vcov(ev)); A[2,2] <- A[2,2]+fit$pen_curv
  inv <- solve(A)
  expected <- inv %*% B %*% t(inv)
  expect_equal(fit$se_sand, sqrt(expected[1,1]), tolerance=1e-8)
})

test_that("rho_mcp propagates through cached calibration", {
  set.seed(210)
  dat <- simulate_hybrid_cox(nI1=60,nI0=60,nE=100,p=0)$data
  direct <- fit_P3_info_MCP(dat, lambda=0.2, rho_mcp=0.3)
  cached <- fit_one_penalized_method_cached(dat, character(),
    cox_fit_full(dat, character()), method="P3", lambda=0.2, rho_mcp=0.3)
  expect_equal(cached, direct)
  sc <- scenario_S1
  sc$nI1 <- 40; sc$nI0 <- 40; sc$nE <- 60
  sc$p <- 0; sc$beta <- numeric(); sc$cov_shift <- numeric()
  cal <- calibrate_lambda_grid("P3", c(0,0.2), sc, drift_set=0,
                                nsim=2, rho_mcp=0.3)
  expect_true(all(cal$details$n_nonmissing_sandwich == 2))
  tuning <- lambdas_default; tuning$rho_mcp <- 0.3
  sim <- run_simulation(nsim=2, scenario=sc, lambdas=tuning)
  expect_equal(sim$lambdas$rho_mcp, 0.3)
  expect_true(all(is.finite(sim$raw$se_sand)))
})
