\name{plot-methods}
\alias{plot.inla}
\alias{plot-methods}
\alias{plot,inla-methods}
\title{Default GAM plotting}
\usage{
\method{plot}{inla}(x, pages, select,...)
}
\arguments{
\item{x}{a fitted inla and foo object produced by gamm.inla.}

\item{pages}{the number of pages over which to spread the output.}

\item{select}{Allows the plot for a single model term to be selected for printing.}
\item{...}{graphics parameters to pass on to matplot commands}

}

\value{
None
}
\details{
Create a default plot of smooth terms in a fitted INLA model
}
\author{Dong Liang}
