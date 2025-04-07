\name{Variog-methods}
\docType{methods}
\alias{Variog}
\alias{Variog-methods}
\alias{Variog,inla-method}
\title{Residuals for inla models}
\description{
  Computes residuals for \code{inla} objects.
}
\usage{
Variog( fit,plot=TRUE,... )
\S4method{Variog}{inla}( fit,plot=TRUE,... )
}
\arguments{
  \item{fit}{Object of class \code{inla}}
  \item{plot}{whether to plot variogram}
  \item{...}{arguments to \code{geoR::variog}}
  }
\details{
  This function computes the variogram from fitted inla objects.

}
\value{
a \code{geoR} fitted variogram object.
}
\references{
}
\author{Dong Liang}
\examples{
\dontrun{
}
}
\keyword{}