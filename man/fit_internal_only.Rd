\name{fit_internal_only}
\alias{fit_internal_only}
\title{Internal-control-only Cox analysis}
\usage{
fit_internal_only(dat, xnames = NULL, robust = FALSE)
}
\arguments{
\item{dat}{A data frame from \code{\link{simulate_hybrid_cox}} or
  conforming to its column conventions (\code{time}, \code{status},
  \code{T}, \code{Z}, and covariates named \code{X1}, ..., \code{Xp}).}
\item{xnames}{Character vector of covariate names. If \code{NULL},
  automatically detected as columns matching \code{"^X\\d+$"}.}
\item{robust}{Use robust (Lin-Wei) Cox standard errors.}
}
\value{
A list of estimates and inference quantities including \code{method},
\code{theta_hat}, \code{se_theta}, \code{se_sand}, \code{z},
\code{z_sand}, \code{delta_hat}, and \code{pen_curv}.
}
\description{
Fits a Cox proportional hazards model on the internal (\code{Z = 0})
subjects only, without any borrowing from external controls.
}
