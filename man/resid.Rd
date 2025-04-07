\name{resid-methods}
\docType{methods}
\alias{resid}
\alias{resid-methods}
\alias{resid,inla-method}
\title{Residuals for inla models}
\description{
  Computes residuals for \code{inla} objects.
}
\usage{
resid( fit,type="working",... )
\S4method{resid}{inla}(  fit,type="working",... )
}
\arguments{
  \item{fit}{Object of class \code{inla}}
  \item{type}{type of GLM residuals}
  }
\details{
  This function computes the residuals from inla objects.

}
\value{
a vector of \code{n} residuals from \code{inla}.
}
\references{
}
\author{Dong Liang}
\examples{
\dontrun{
}
}
\keyword{}