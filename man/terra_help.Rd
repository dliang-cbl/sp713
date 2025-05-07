\name{raster2nb}
\alias{raster2nb}
\alias{xy2cell}
\alias{cell2xy}
\title{Utility functions for rast objects in CAR models}
\usage{
    raster2nb(object,eps=NULL,drop=NULL,...)
    xy2cell(object,xy,eps=NULL,nna=TRUE)
    cell2xy(object,id,eps=NULL)
}
\arguments{
    \item{object}{a \code{rast} object}
    \item{drop}{an integer vector of cells to drop}
    \item{eps}{a scalar or vector of size 2, resolution for the raster}
    \item{id}{an integer vector of cell id in the raster}
    \item{nna}{logical whether to remove missing cells in the ID}
    \item{xy}{a two column matrix of geographic coordinates}
}
\value{
a \code{nb} object from \code{raster2nb}, cell ID or coordinates from the \code{cell2xy} or \code{xy2cell} functions.
}
\details{
    Helper function to construct data structure for Conditional Autoregressive Models in INLA.
}
\author{Dong Liang}