\name{extract.samples-methods}
\docType{methods}
\alias{extract.samples}
\alias{extract.samples-methods}
\alias{extract.samples,inla-method}
\title{Predictions for inla models}
\description{
    Collect posterior samples from a \code{inla} model.
}
\usage{
extract.samples(object,n=10000,pars,... )
\S4method{extract.samples}{inla}(object,n=10000,pars,... )
}
\arguments{
    \item{object}{Object of class \code{inla}}
    \item{n}{Number of samples to simulate}
    \item{pars}{Character vector of parameters to return}
    \item{...}{Other parameters to pass}
}
\details{
 This function return lists of samples for posterior 
 distributions, respectively.
}
\value{
a data frame of \code{n} rows of posterior samples from \code{inla}.
}
\references{
}
\author{Dong Liang}
\examples{
\dontrun{
}
}
\keyword{}
