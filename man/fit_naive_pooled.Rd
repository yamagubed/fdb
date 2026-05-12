\name{fit_naive_pooled}
\alias{fit_naive_pooled}
\title{Naive pooled Cox analysis}
\usage{
fit_naive_pooled(dat, xnames = NULL, robust = FALSE)
}
\arguments{
\item{dat}{A data frame from \code{\link{simulate_hybrid_cox}} or
  conforming to its column conventions.}
\item{xnames}{Character vector of covariate names. If \code{NULL},
  automatically detected.}
\item{robust}{Use robust (Lin-Wei) Cox standard errors.}
}
\value{
A list of estimates and inference quantities.
}
\description{
Fits a Cox proportional hazards model treating concurrent and
external controls as exchangeable (full pooling), without any
indicator for external-control membership.
}
