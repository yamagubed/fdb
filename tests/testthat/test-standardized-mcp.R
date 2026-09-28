test_that("standardized MCP has the correct scale, joins, and derivatives", {
  lambda <- 0.25; gamma <- 3; rho <- 0.1; eps <- 0.001
  for (se in c(0.05, 0.122, 0.4, 1.5)) {
    b <- gamma * lambda; h <- rho * b
    p <- function(d) pen_MCP(d, lambda, gamma, eps, rho, se)
    dp <- function(d) pen_derivative_MCP(d, lambda, gamma, eps, rho, se)
    ddp <- function(d) pen_curvature_MCP(d, lambda, gamma, eps, rho, se)
    joins <- sqrt((se * c(b-h, b+h) + eps)^2 - eps^2)
    for (d in c(0, eps, 0.1*se, joins, 2*se, -joins)) {
      step <- 1e-9
      expect_equal((p(d+step)-p(d-step))/(2*step), dp(d), tolerance=1e-5)
      expect_equal((dp(d+step)-dp(d-step))/(2*step), ddp(d), tolerance=1e-5)
      expect_equal(p(d), p(-d))
      # Independent quadrature on the standardized scale.
      v <- (sqrt(d*d+eps*eps)-eps)/se
      slope <- function(u) ifelse(u <= b-h, lambda-u/gamma,
        ifelse(u < b+h, (b+h-u)^2/(4*gamma*h), 0))
      cuts <- sort(unique(c(0,v,c(b-h,b+h)[c(b-h,b+h)<v])))
      integral <- if (v==0) 0 else sum(vapply(seq_len(length(cuts)-1L),
        function(i) integrate(slope,cuts[i],cuts[i+1],rel.tol=1e-10)$value,numeric(1)))
      expect_equal(p(d), integral, tolerance=1e-10)
    }
    expect_equal(c(p(0),dp(0),ddp(0)), c(0,0,lambda/(se*eps)))
    expect_equal(c(dp(2*se),ddp(2*se)), c(0,0))
    expect_equal(p(2*se), gamma*lambda^2/2+h^2/(6*gamma))
    for (v in c(0.1, 0.5, 0.75, 1.5)) {
      target <- if (v <= b) lambda*v-v^2/(2*gamma) else gamma*lambda^2/2
      # Absolute approximation bound: origin displacement plus rounded tail.
      bound <- lambda*1e-10/se + (1e-6*b)^2/(6*gamma)
      expect_lte(abs(pen_MCP(se*v,lambda,gamma,1e-10,1e-6,se)-target),
                 bound+10*.Machine$double.eps)
    }
  }
})

test_that("P3 fitting minimizes standardized MCP and uses its raw curvature", {
  set.seed(8112)
  for (drift in log(c(0.8, 1, 1.25))) {
    dat <- simulate_hybrid_cox(nI1=80,nI0=60,nE=100,p=0,delta0=drift)$data
    full <- cox_fit_full(dat,character())
    for (lambda in c(0,0.03,0.25,0.8)) {
      fit <- fit_P3_info_MCP(dat,lambda=lambda,full_fit=full,n_grid_opt=5)
      se <- full$se_delta; eps <- 0.001
      obj <- function(d) -cox_fit_profile(dat,d,character())$loglik +
        pen_MCP(d,lambda,3,eps,0.1,se)
      # Independent dense-grid search, refining every sampled local minimum.
      grid <- seq(-0.8,0.8,length.out=161)
      vals <- vapply(grid,obj,numeric(1))
      idx <- which(vals <= c(Inf,head(vals,-1)) & vals <= c(tail(vals,-1),Inf))
      refs <- vapply(idx,function(i) optimize(obj,
        grid[c(max(1,i-1),min(length(grid),i+1))],tol=1e-9)$objective,numeric(1))
      expect_lte(obj(fit$delta_hat), min(refs)+1e-5)
      expect_equal(fit$pen_curv_raw,
        pen_curvature_MCP(fit$delta_hat,lambda,3,eps,0.1,se))
      expect_equal(fit$pen_curv,max(0,fit$pen_curv_raw))
      expect_equal(fit$h,0.1*3*lambda)
      expect_identical(fit$p3_parameterization,"standardized_drift_v1")
      if (lambda==0) expect_equal(fit$delta_hat,full$delta_hat,tolerance=1e-3)
    }
  }
})

test_that("P3 reaches its flat tail at drift measured in standard-error units", {
  se <- 0.122; lambda <- 0.25
  expect_equal(pen_derivative_MCP(0.2,lambda,se_delta=se),0)
  expect_lt(pen_curvature_MCP(0.05,lambda,3,se_delta=se),0)
  # The old raw-scale parameterization would still shrink this drift.
  expect_gt(pen_derivative_MCP(0.2,lambda/se),0)
})
