\name{link-methods}
\docType{methods}
\alias{link}
\alias{link-methods}
\alias{link,inla-method}

\title{Predictions for inla models}
\description{
  Computes inverse-link linear model values for \code{inla} samples.
}
\usage{
link( fit,data,n=1000,type="confidence",E=NULL,Ntrials=NULL,... )
\S4method{link}{inla}( fit,data,n=1000,type="confidence",E=NULL,Ntrials=NULL,... )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{Object of class \code{inla}}
  \item{data}{Optional list of data to compute predictions over. When missing, uses data found inside fit object.}
  \item{n}{Number of samples to use}
  \item{type}{the type of prediction required. The default is the confidence intervals for the mean on the response scale. It could also be predictive interval on the response scale.}
  \item{E}{the offset in original scale for Poisson, negative binomial model, must be equal in length to the number of rows in non missing \code{data}.}
  \item{Ntrials}{the number of trials for Binomial model, must be equal in length to the number of rows in non missing \code{data}.}
  \item{...}{Other parameters to pass to someone}
}
\details{
  This function computes the value of each linear model at each sample for each case in the data. Inverse link functions are applied, so that for example a logit link linear model produces probabilities, using the logistic transform. 

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

