\name{spde.gam}
\alias{spde.gam}
\title{Stochastic Partial Differential Equations approach to geostatistical generalized additive model}
\usage{
  spde.gam(
    formula,data,mesh,spde,xcoord='x',ycoord='y',
    ...)
}
\arguments{
\item{formula}{fixed effects type formula mostly used in \code{gam}}

\item{data}{data frame to evaluate the model}

\item{mesh}{mesh for building the Stochastic Partial Differential Equations}

\item{spde}{the Stochastic Partial Differential Equations model object}

\item{xcoord}{field name of x coordinate}

\item{ycoord}{field name of y coordinate}


\item{...}{arguments for the \code{inla} call}
}
\value{
fitted INLA object
}
\description{
geostatistical GAM with Stochastic Partial Differential Equations using R and Integrated Nested Laplace Approximation.

This function replaces the \code{spde.glm} in earlier version of the package.

The formula can not contain \code{poly}, polynomial regression must be developed by hand.
}
\author{Dong Liang}
