## inla2raster
## extract the xy coordinates to and from node ID
## ver 3) allow change of resolution
# require(raster)
## arguments:
##  object  : raster* object
##  eps     : resolution (default no change)
##  id      : node id (non-missing)
##  xy      : matrix of x and y coordinates
## nna      : whether to exclude missing values
xy2cell <- function(object,xy,eps=NULL,nna=TRUE){
  if(!is.null(eps)){
    res(object) <- eps
  }
  xyz <- as.data.frame(object,xy=T)      ## extract all non-missing cells
  id0 <- cellFromXY(object,xyz[,1:2]) ## the cell number in the raster 
  id <- cellFromXY(object,xy)
  if(nna){
    id <- match(id,id0)
  }
  id
}
cell2xy <- function(object,id,eps=NULL)
{
  if(!is.null(eps)){
    res(object) <- eps
  }
  xyz <- as.data.frame(object,xy=T)      ## extract all non-missing cells
  xyz[id,1:2]
}

test <- function()
{
  rm(list=ls())
  library(raster)
  source("inla2raster_2.R")
  source("raster2nb_2.R")
  r <- raster(nrows=10, ncols=10)
  xy <- xyFromCell(r,1:100)
  mu <- apply((xy/200)^2,1,sum)
  values(r) <- rnorm(100,mu,0.1)
  
  r[sample(100,30)] <- NA
  plot(r)
  
  obs <- cbind(runif(50,-250,250),runif(50,-80,80))
  plot(r)
  points(obs)
  nodes <- raster2nb(r,obs,"tmp.txt")
  
  obs2 <- inla2raster(r,nodes)
  plot(r)
  for(i in 1:nrow(obs2)){
    Sys.sleep(0.3)
    points(obs[i,,drop=FALSE])
    points(obs2[i,,drop=FALSE],pch="X",cex=2)  
  }
}