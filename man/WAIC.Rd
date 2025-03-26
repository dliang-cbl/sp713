\name{WAIC-methods}
\docType{methods}
\alias{WAIC}
\alias{WAIC-methods}
\alias{WAIC,list-methods}
\alias{PSIS}
\alias{PSIS-methods}
\alias{PSIS,list-methods}
\alias{LOO-methods}
\alias{LOO,list-methods}

\title{Information Criteria and Pareto-Smoothed Importance Sampling Cross-Validation}
\description{
  Computes WAIC, and PSIS cross validation for \code{inla} model fits.
}
\usage{
WAIC( object , n=1000 , refresh=0.1 , pointwise=FALSE , ... )
PSIS( object , n=1000 , refresh=0.1 , pointwise=FALSE , ... )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{Object of class \code{map} or \code{map2stan}}
  \item{n}{Number of samples to use in computing WAIC. Set to \code{n=0} to use all samples in \code{map2stan} fit}
  \item{refresh}{Refresh interval for progress display. Set to \code{refresh=0} to suppress display.}
  \item{pointwise}{If \code{TRUE}, return a vector of WAIC values for each observation. Useful for computing standard errors.}
  \item{...}{Other parameters to pass to specific methods}
}
\details{
  These functions use the samples and model definition to compute the Widely Applicable Information Criterion (WAIC), or Pareto-smoothed importance-sampling cross-validation estimate (PSIS). 
  
  WAIC is an estimate of out-of-sample relative K-L divergence (KLD), and it is defined as:

  \deqn{WAIC = -2(lppd - pWAIC)}

  Components \code{lppd} (log pointwise predictive density) and \code{pWAIC} (the effective number of parameters) are reported as attributes. See Gelman et al 2013 for definitions and formulas. This function uses the variance definition for \code{pWAIC}.

  PSIS is another estimate of out-of-sample relative K-L divergence. It is computed by the \code{loo} package. See Vehtari et al 2015 for definitions and computation.

  In practice, WAIC and PSIS are extremely similar estimates of KLD.

}
\value{
}
\references{
Watanabe, S. 2010. Asymptotic equivalence of Bayes cross validation and Widely Applicable Information Criterion in singular learning theory. Journal of Machine Learning Research 11:3571-3594.

Gelman, A., J. Hwang, and A. Vehtari. 2013. Understanding predictive information criteria for Bayesian models.

Vehtari, A., A. Gelman, and J. Gabry. 2015. Efficient implementation of leave-one-out cross-validation and WAIC for evaluating fitted Bayesian models.
}
\author{Dong Liang}
\seealso{\code{\link{link}}}
\examples{
\dontrun{
}
}
\keyword{ }

