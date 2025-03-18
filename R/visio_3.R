#library(classInt)
## Helper function to visualize sp object with legends
## Author: D. Liang CBL/UMCES 12/2016
## Comes with ABSOLUTELY NO WARRANTY
## Arguments:
## obj : a simple feature object to visualize
## field : the field name to visualize
## leg.col: color palette
## choropleth: whether showing choropleth or numeric color
## bins  : for choropleth the binary color dividers
## nquantile.color: default division for binary color definitions
## na.color: The color to return for NA values.
## pt.cex : for point data the size of the point
## pt.col : for point data the color of the outline
## pt.weight: for point data the weight/width of the outline
## project: whether to project raster to EPSG:3857
## legend : whether to include legend
## legend.title: title of the legend
## add    : add to an existing layer, not used 
## ...    : additional arguments to graphic layer code
## More details can be found at
## https://rstudio.github.io/leaflet/
visio <- function(
  obj,leg.col="white",field=NULL,
  choropleth=TRUE,bins=NA,nquantile.color=4,na.color="#808080",
  pt.cex=4,pt.col="black",pt.weight=1,
  project=TRUE,Graticule=FALSE,
  legend=TRUE,legend.position="bottomright",legend.title="",
  add=NULL,
  ...)
{
  if(class(obj)[1]=="SpatRaster"){
    geo_ <- "SpatRaster"
  }else if(class(obj)[1]=="sf"){
    tmp_ <- st_geometry_type(obj)[1] ## assuming same type
    if(length(grep("POINT",tmp_))>0){
      geo_ <- "POINT"
    }else if(length(grep("POLYGON",tmp_))>0){
      geo_ <- "POLYGON"
    }else{
      geo_ <- tmp_
    }
  }else{
    stop("input type not supported.\n")
  }
  ## extract response data fields
  if(geo_=="SpatRaster"){
    resp__ <- (as.vector(obj))
  }
  else if(class(obj)[1] %in% c("sf")){
    ## extract response data fields
    if(is.null(field)){
      stop("field can not be null for spatial objects.")
    }
    resp__ <- st_drop_geometry(obj)[,field]
  } 
  else{
    stop("can't handle this type of spatial object")
  }
  
  ## define filled color palette
  pal <- visio.color.helper(
    x=resp__,leg.col = leg.col,choropleth = choropleth,
    bins=bins,nquantile.color = nquantile.color,na.color = na.color)
  
  ## define base map
  if(is.null(add)){
    if(class(obj)[1]=="SpatRaster"){
      base <- leaflet() %>% addTiles()  ## add base tiles
    }
    else if(class(obj)[1] %in% c("sf")){
      base <- leaflet(obj) %>% addTiles()  ## add base tiles
    }
    map1 <- base %>% addScaleBar() ## add geogrpahic details
    if(Graticule){
      map1 <- map1 %>% addGraticule()
    }
  }else{
    if(class(add)[1] == "leaflet" & class(add)[2]=="htmlwidget"){
      map1 <- add
    }else{
      stop("base map is not leaflet.\n")
    }
  }
  
  ## add features
  if(class(obj)[1]=="SpatRaster"){
    map.out <- map1 %>% addRasterImage(
      x=obj, colors=pal, opacity=1,project=project,...
    )
  }
  else if(geo_ == "POINT"){
    ## spatial point data
    radius__ <- pt.cex  ## size of points
    if(is.na(pt.col)){
      ## outline color transparent
      map.out <- map1%>%addCircleMarkers(
        radius=radius__,
        stroke=FALSE,
        fillColor=~pal(resp__),
        fillOpacity=0.9,
        data=obj,...
      )
    }else{
      map.out <- map1%>%addCircleMarkers(
        radius=radius__,
        stroke=TRUE,
        color=pt.col,
        weight=pt.weight,
        opacity=1,
        fillColor=~pal(resp__),
        fillOpacity=0.9,
        data=obj,...
      )
    }
  }
  else if(geo_ == "POLYGON"){
    if(is.na(pt.col)){
      ## outline color transparent
      map.out <- map1 %>% addPolygons(
        stroke=FALSE, 
        fillColor = ~pal(resp__),
        fillOpacity=0.9,...
      )
    }else{
      ## outline color provided
      map.out <- map1 %>% addPolygons(
        stroke=TRUE,
        color=pt.col,
        weight=pt.weight,
        opacity=1,
        fillColor = ~pal(resp__),
        fillOpacity = 0.9,...
      )
    }
  }

  ## add legend
  if(legend){
    map.out <- map.out %>%
      addLegend(position=legend.position,
                pal=pal,values=resp__,
                opacity=0.9,
                title=legend.title)
  }
 
  map.out
}

## helper function to define color in leaflet
## return a color pallette
visio.color.helper <- function(x,leg.col,choropleth,bins,nquantile.color,na.color){
  if(class(x)=="factor"){
    ## categorical data
    pal <- colorFactor(
      palette = leg.col,
      domain = x,
      na.color=na.color
    )
  }else{
    ## quantitative data
    if(choropleth){
      ## discrete color
      if(any(is.na(bins))){
        ## contain missing values
        if(length(bins)>1){
          ## assume correct
          warning("missing value ignored in bins argument.\n")
          bins <- na.omit(bins)
        }else{
          ## assume using quantiles
          bins <- quantile(x,probs=seq(0,1,length.out = nquantile.color + 1),na.rm=TRUE)
          bins <- unique(bins) ## remove redundant quantiles
          if(length(bins)==1){
            ## a single value
            bins <- c(bins-0.01,bins+0.01)
          }
        }
      }
      pal <- colorBin(
        palette = leg.col,
        domain = x,
        bins = bins,
        na.color = na.color
      )
    }else{
      ## continuous color
      pal <- colorNumeric(
        palette=leg.col,
        domain=x,
        na.color = na.color)
    }
  }
  pal
}

show.knitr <- function(m,file=NULL)
{
  if(is.null(file)){
    file <- tempfile(tmpdir=".",fileext=".png")
  }
  mapshot(m,file=file)
  knitr::include_graphics(substr(file,3,nchar(file)))
}

debug_ <- function(){
  rm(list=ls())
  library(terra)
  library(sf)
  library(leaflet)
  source("visio_3.R")
  rca <- rast("rca.tif")
  source("visio_3.R")
  visio(rca,leg.col="Blues",na.color = NA)
  ## 1. Winter Dredge Survey, an example geo-statistical data
  wds0 <- read.csv("wds.csv")
  ## Make a SF for the point data in Long/Lat
  wds <-  st_as_sf(wds0, coords = c("ALONG", "ALAT"), crs = 4326)
  visio(wds,field = "TOT",leg.col="Blues")
  shoreline <- st_read("CBnet_Shoreline_simplify.shp")
  head(shoreline)
  shoreline2 <- st_transform(shoreline,crs = 4326)
  visio(shoreline2,field = "AREA",leg.col="Blues")
  
}
