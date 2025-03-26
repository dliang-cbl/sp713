\name{extract_log_lik-methods}
\docType{methods}
\alias{extract_log_lik}
\alias{extract_log_lik-methods}
\alias{extract_log_lik,inla-methods}

\title{Extract log likelihood for inla models}
\description{
  Extract log likelihood for inla models from \code{inla} samples.
}
\usage{
extract_log_lik(object , n=1000 , refresh=0.1 , pointwise=FALSE ,log_like = "log_lik", ...)
\S4method{extract_log_lik}{inla}(object , n=1000 , refresh=0 , pointwise=FALSE , ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{Object of class \code{inla}}
  \item{n}{Number of samples to use}
  \item{refresh}{TBD}
  \item{pointwise}{TBD}
  \item{...}{Other parameters to pass to someone}
}
\details{
  TBD
}
\value{
a data frame of \code{n} rows of log likelihood samples from \code{inla}.
}
\references{
}
\author{Dong Liang}
\examples{
\dontrun{
}
}
\keyword{}

