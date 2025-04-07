#' geostatistical generalized linear model using SPDE 
#' 
#' spatial GLM with Stochastic Partial Differential Equations using R and Integrated Nested Laplace Approximation. 
#'
#' @param formula fixed effects type formula object
#' @param data data frame to evaluate the model
#' @param mesh mesh for building the Stochastic Partial Differential Equations
#' @param spde the Stochastic Partial Differential Equations model
#' @param xcoord field name of x coordinate
#' @param ycoord field name of y coordinate
#' @param query prediction data
#' @... additional argument to INLA
#' 
#' @return fitted INLA object
#' 
spde.glm <- function(
    formula,data,mesh,spde,query=NULL,xcoord='x',ycoord='y',
    ...){
  
  # input argument check
  
  # Evaluate arguments in data
  attach(data)
  args__ <- list(...)
  detach(data)
  
  #browser()
  ## adjust ... argument to allow posterior simulation
  args__ <- sgam.config.compute(args__)
  ## adjust ... to allow SPDE fitting
  args__[["control.predictor"]] <- NULL
  
  if(!is.null(query)){
    stopifnot(class(query)==class(data))
  }
  
  obs <- nmiss(formula,data)
  
  ## process formula
  l.vars <- getVars(formula)
  zcoord <- l.vars[1]
  if(nchar(zcoord)<1){
    stop("please specify a response in formula.\n")
  }
  
  #browser()
  ## extract fixed effects and error terms
  terms_ <- terms(formula)
  terms_labels <- attr(terms_,"term.labels") ## each term
  if(length(terms_labels)>0){
    terms_vars_unique <- l.vars[-1] ## associated variables
    terms_vars <- rep("",length(terms_labels))
    for(i in 1:length(terms_vars_unique)){
      j <- grep(terms_vars_unique[i],terms_labels)
      terms_vars[j] <- terms_vars_unique[i]
    }
  }else{
    terms_vars <- NULL
  }
  #stopifnot(length(terms_vars)==length(terms_labels))
  is_f <- grep("^f(.+)$",terms_labels)
  if(length(is_f)){
    random <- reformulate(terms_labels[is_f],response = NULL)
    ## f term included
    terms_f <- terms_vars[is_f]
    if(length(terms_vars)==length(terms_f)){
      ## only f term
      fixed <- reformulate("1",response=NULL)
    }else{
      ## fixed term included
      fixed <- reformulate(terms_labels[-is_f],response = NULL)
    }    
    ## look for error terms
    terms_nlevel <- rep(0,length(terms_f))
    for(i in 1:length(terms_f)){
      terms_nlevel[i] <- length(unique(obs[,terms_f[i]]))
    }
    is_err <- terms_nlevel == nrow(obs)
    err_term <- terms_f[is_err]
  }else{
    ## no f term
    random <- reformulate("1",response=NULL)
    fixed <- reformulate(as.character(formula)[3],response=NULL)
    ## no error term
    err_term <- ""
  }
  
  
  ## fixed and random effects observed
  z.obs.lst <- list(b0=rep(1,nrow(obs)))
  
  if(length(c(all.vars(fixed),all.vars(random)))>0){
    z.obs.lst <- c(z.obs.lst,as.list(
      obs[,c(all.vars(fixed),all.vars(random)),drop=FALSE]))
  }
  
  ## fixed effects query
  if(!is.null(query)){
    z.query.lst <- list(b0=rep(1,nrow(query)))
    if(length(all.vars(fixed))>0){
      z.query.lst <- c(z.query.lst,as.list(
        query[,c(all.vars(fixed),all.vars(random)),drop=FALSE]))
    }
  }
  
  ## Projection matrices
  coords <- as.matrix(obs[,c(xcoord,ycoord)])
  A <- inla.spde.make.A(mesh, loc=coords)
  
  ## Define SPDE data stacks -- observed
  stk.z.obs <- inla.stack(tag="obs",
                          data=list(y=obs[,zcoord]),
                          A=list(A,1),
                          effects=list(
                            list(field=1:spde$n.spde),
                            z.obs.lst
                          )
  )
  
  ## Define SPDE data stacks -- query if any
  if(!is.null(query)){
    pcoords <- as.matrix(query[,c(xcoord,ycoord)])
    Ap <- inla.spde.make.A(mesh, loc=pcoords)
    stk.z.query <- inla.stack(tag="query",
                              data=list(y=rep(NA,nrow(query))),
                              A=list(Ap,1),
                              effects=list(
                                list(field=1:spde$n.spde),
                                z.query.lst
                              )
    )
  }
  
  ## SPDE data structure combined
  if(!is.null(query)){
    stk.z <- inla.stack(stk.z.obs,stk.z.query)
  }
  else{
    stk.z <- stk.z.obs
  }
  #browser()
  ## Define INLA formula
  terms_fixed <- attr(terms(fixed),"term.labels")
  interceptonly <- length(terms_fixed)==0
  #interceptonly <- length(getVars(fixed))==1
  if(interceptonly){
    fixed.vars <- c("0+b0")
  }else{
    fixed.vars <- c("0+b0",as.character(fixed)[2])
  }
  
  terms_random <- attr(terms(random),"term.labels")
  if(length(terms_random)>0){
    random.vars <- as.character(random)[2]
  }else{
    random.vars <- NULL
  }
  spatial.vars <- c("f(field,model=spde)")
  formula__ <- reformulate(c(fixed.vars,
                             random.vars,
                             spatial.vars),response="y")
  
  
  ## posterior inference
  #browser()
  mod1 <- try(do.call("inla",c(list(
    formula=formula__,
    data=inla.stack.data(stk.z),
    control.predictor=list(A=inla.stack.A(stk.z),compute=TRUE)),
    args__)))

  ## save arguments
  mod1$.args$formula0 <- formula
  mod1$.args$fixed <- fixed
  mod1$.args$random <- random
  mod1$.args$raw <- obs
  mod1$.args$error <- err_term
  mod1$.args$mesh <- mesh
  mod1$.args$spde <- spde
  mod1$.args$coords <- coords
  
  ## return model
  mod1
}
dev <- function(){
  ## estimation
  rm(list=ls())
  source("sgam_2.R")
  source("nmiss_2.R")
  source("getVars.R")
  source("spde_glm_1.R")
  load(file="Scratch/spde_glm_dev1.rdata")
  
  wds.geo$err <- seq(1,nrow(wds.geo))
  wds.geo$X2 <- as.vector(scale(wds.geo$X))
  wds.geo$Y2 <- as.vector(scale(wds.geo$Y))
  ## poly(X,2) does not work because missing values not allowed
  mod <- spde.glm(TOT~1,data=wds.geo,mesh=mesh1,spde=spde1,
                  xcoord="X",ycoord="Y",family="poisson",
                  E=AREA)
  mod <- spde.glm(TOT~f(err),data=wds.geo,mesh=mesh1,spde=spde1,
                  xcoord="X",ycoord="Y",family="poisson",
                  E=AREA)
  mod <- spde.glm(TOT~X2+I(X2^2)+f(err),data=wds.geo,mesh=mesh1,spde=spde1,
                  xcoord="X",ycoord="Y",family="poisson",
                  E=AREA)
  mod <- spde.glm(TOT~X2+Y2+f(err),data=wds.geo,mesh=mesh1,spde=spde1,
                  xcoord="X",ycoord="Y",family="poisson",
                  E=AREA)
  mod <- spde.glm(TOT~X2+Y2,data=wds.geo,mesh=mesh1,spde=spde1,
                  xcoord="X",ycoord="Y",family="poisson",
                  E=AREA)
  
  summary(mod)
  
  if(F){
    library(INLA)
    mod <- spde.glm(TOT~X2,data=wds.geo,mesh=mesh1,spde=spde1,
                    xcoord="X",ycoord="Y")
    mod <- spde.glm(TOT~X2+I(X2^2),data=wds.geo,mesh=mesh1,spde=spde1,
                    xcoord="X",ycoord="Y")
    wds.geo$yi <- cut(wds.geo$Y,breaks = 13,labels = F)
    wds.geo$xi <- cut(wds.geo$X,breaks = 13,labels = F)
    
    mod <- spde.glm(TOT~X2+I(X2^2)+f(xi,model = "rw1")+f(yi),data=wds.geo,mesh=mesh1,spde=spde1,
                    xcoord="X",ycoord="Y")
    
    wds.geo$N <- as.integer(wds.geo$TOT/runif(nrow(wds.geo)))
    mod <- spde.glm(TOT~1,
                    data=wds.geo,mesh=mesh1,spde=spde1,
                    xcoord="X",ycoord="Y",family="binomial",
                    Ntrials=N)
    mod <- spde.glm(TOT~X2+I(X2^2),
                    data=wds.geo,mesh=mesh1,spde=spde1,
                    xcoord="X",ycoord="Y",family="binomial",
                    Ntrials=N)
    mod <- spde.glm(TOT~X2+I(X2^2)+f(xi,model = "rw1")+f(yi),
                    data=wds.geo,mesh=mesh1,spde=spde1,
                    xcoord="X",ycoord="Y",family="binomial",
                    Ntrials=N)
    
  }
}
test_spde_glm <- function(){
  rm(list=ls())
  library(INLA)
  library(terra)
  data(SPDEtoy)
  pl.dom <- cbind(c(0, 1, 1, 0.7, 0), c(0, 0, 0.7, 1, 1))
  mesh5 <- inla.mesh.2d(loc.domain = pl.dom, max.e = c(0.092, 0.2))
  
  spde5 <- inla.spde2.pcmatern(
    # Mesh and smoothness parameter
    mesh = mesh5, alpha = 1.5,
    # P(practic.range < 0.3) = 0.5
    prior.range = c(0.3, 0.5),
    # P(sigma > 1) = 0.01
    prior.sigma = c(10, 0.01))
  
  source("sgam_2.R")
  source("nmiss_2.R")
  source("getVars.R")
  source("spde_glm_1.R")
  mod <- spde.glm(y~1,data=SPDEtoy,mesh=mesh5,spde=spde5,
                  xcoord="s1",ycoord="s2")
  
  
  rca <- rast("rca.tif")
  
  cru<- read.csv("cruise.csv")
  mesh <- inla.mesh.2d(
    cru[,1:2],max.edge=c(0.1,0.3),cut=0.05)
  
  ## 2) approaximate covariance function for large data
  spde <- inla.spde2.matern(mesh)
  
  grid <- as.data.frame(rca,xy=T)
  
  #library(FNN)
  #names(cru)[1:2] <- c("x","y")
  mod <- spde.glm(Value~1,data=cru,mesh=mesh,spde=spde,
                  xcoord="Longitude",ycoord="Latitude")
  
}
