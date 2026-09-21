\name{lambdas_default}
\alias{lambdas_default}
\docType{data}
\title{Default tuning parameters}
\format{
A list with elements \code{lambda_li}, \code{gamma_li},
\code{lambda_p1}, \code{lambda_p2}, \code{gate_c}, \code{gate_tau},
\code{lambda_p3}, \code{gamma_mcp}, \code{rho_mcp} (default 0.1), \code{lambda_p4}, and
\code{delta_bounds}.
}
\usage{
lambdas_default
}
\description{
A reference tuning list for the borrowing methods, used when no
calibrated tuning is supplied. The default penalty strengths
(\eqn{\lambda = 0.20}) are reasonable starting values for the
default scenario but should generally be replaced by calibrated
values for confirmatory use.
}
\keyword{datasets}
