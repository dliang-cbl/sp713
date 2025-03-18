\name{plot.foo}
\alias{plot.foo}
\title{Default GAM plotting}
\usage{
\method{plot}{foo}(x, pages, select,...)
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
\description{
Create a default plot of smooth terms in a fitted INLA model
}
