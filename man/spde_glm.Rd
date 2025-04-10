\name{spde.glm}
\alias{spde.glm}
\title{Stochastic Partial Differential Equations approach to geostatistical generalized linear model}
\usage{
  spde.glm(
    formula,data,mesh,spde,query=NULL,xcoord='x',ycoord='y',
    ...)
}
\arguments{
\item{formula}{fixed effects type formula mostly used in \code{lm}}

\item{data}{data frame to evaluate the model}

\item{mesh}{mesh for building the Stochastic Partial Differential Equations}

\item{spde}{the Stochastic Partial Differential Equations model object}

\item{xcoord}{field name of x coordinate}

\item{ycoord}{field name of y coordinate}

\item{query}{prediction data (optional)}


\item{...}{arguments for the \code{inla} call}
}
\value{
fitted INLA object
}
\description{
geostatistical GLM with Stochastic Partial Differential Equations using R and Integrated Nested Laplace Approximation.

The formula can not contain \code{poly}, polynomial regression must be developed by hand.
}
