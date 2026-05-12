\name{add_ess_to_simulation_result}
\alias{add_ess_to_simulation_result}
\title{Append ESS to a simulation summary}
\usage{
add_ess_to_simulation_result(simres, NS, ref_method = "InternalOnly")
}
\arguments{
\item{simres}{Output of \code{\link{run_simulation}}.}
\item{NS}{Internal-control sample size used as the scaling factor.}
\item{ref_method}{Reference method (default \code{"InternalOnly"}).}
}
\value{
The same \code{simres} list with an \code{ESS} column merged into
its \code{summary} component.
}
\description{
Convenience wrapper around \code{\link{compute_ess_from_raw}}.
}
