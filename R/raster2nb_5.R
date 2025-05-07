## raster2nb
## helper script to make a CAR structure matrix for INLA from a raster object
## assuming NA cells are masked and not in the model
## ver 4) allow removing certain cells
## ver 5) written for terra
#require(spdep)
#require(raster)
## arguments:
## object : raster* object
##          usually a base raster defining the study area
##          missing values were interpreted as outside the study area.
## eps    : chosen raster resolution
## drop   : a logical vector of the length of output
##          cells to remove from the raster output.
## ... : additional arguments to adjacent

## value:
## a nb object defining the neighborhood structures
raster2nb <- function(object,eps=NULL,drop=NULL,...)
{
  ## change object resolution if needed
  if(!is.null(eps)){
    res(object) <- eps
  }
  #browser()
  xyz <- as.data.frame(object,xy=T)      ## extract all non-missing cells
  if(!is.null(drop)){
    xyz <- xyz[!drop,]
  }
  id <- cellFromXY(object,xyz[,1:2]) ## the cell number in the raster 
  adj <- adjacent(object,id,pairs=TRUE,...) ## find adjacent cells to each non-missing cell
  
  ## identify the neighbors for each non-missing cells
  to1 <- match(adj[,2],id)
  from1 <- match(adj[,1],id)

  ## initialize the output list to integer(0)
  lst2 <- vector("list",length(unique(id))) ## some neighbor list is empty, need to account for that
  for(i in 1:length(lst2)){
    lst2[[i]] <- as.integer(0)
  }
  lst3 <- split(to1[!is.na(to1)],from1[!is.na(to1)])
  lst2[as.numeric(names(lst3))] <- lst3
  
  class(lst2) <- "nb"
  
  # ## create the INLA graph
  # nb2INLA(file,lst2)
  # 
  # ## identify the observed cell in the raster
  # obs <- cellFromXY(object,xy)
  # node <- match(obs,id)
  # 
  # node
  lst2
}

test <- function()
{
  rm(list=ls())
  library(terra)
  source("raster2nb_5.R")
  r <- rast(nrows=10, ncols=10)
  xy <- xyFromCell(r,1:100)
  mu <- apply((xy/200)^2,1,sum)
  values(r) <- rnorm(100,mu,0.1)
  plot(r)
  set.seed(112300)
  drop <- rep(F,100)
  drop[sample(100,40)] <- T
  nb1 <- raster2nb(r)  
  library(spdep)
  plot(nb1,xy)
  
  debug(raster2nb)
  nb2 <- raster2nb(r,drop=drop)
  length(nb2)
  summary(nb2)
  sapply(nb2,length)
  tmp <- nb2mat(nb2,style="B",zero.policy = T)
  plot(nb2,xy[!drop,])
  sapply(nb2,max)
  x <- xy[!drop,]
  points(x[3,1],x[3,2],pch="X",col=2)
  length(nb2)
}