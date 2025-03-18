\name{sglm.inla}
\alias{sglm.inla}
\title{geostatistical generalized linear model}
\usage{
sglm.inla(formula, data, mesh, spde, random = ~s, family = "gaussian",
  xcoord = "x", ycoord = "y", E = NULL, Ntrials = NULL,
  query = NULL, verbose = F, control = list(nsamples = 0))
}
\arguments{
\item{formula}{fixed effects type formula object}

\item{data}{data frame to evaluate the model}

\item{mesh}{mesh for building the Stochastic Partial Differential Equations}

\item{spde}{the Stochastic Partial Differential Equations model}

\item{random}{formula defining the random term}

\item{family}{family for response variable}

\item{xcoord}{field name of x coordinate}

\item{ycoord}{field name of y coordinate}

\item{E}{offset for Poisson term}

\item{Ntrials}{number of trials for binomial}

\item{query}{prediction data}

\item{verbose}{detailed}

\item{control}{a list of control objects}
}
\value{
fitted INLA object
}
\description{
geostatistical GLM with Stochastic Partial Differential Equations using R and Integrated Nested Laplace Approximation.
}
