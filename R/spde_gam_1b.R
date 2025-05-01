#' geostatistical generalized additive model using SPDE 
#' 
#' spatial GAM with Stochastic Partial Differential Equations using R and Integrated Nested Laplace Approximation. 
#'
#' @param formula gam type formula object
#' @param data data frame to evaluate the model
#' @param mesh mesh for building the Stochastic Partial Differential Equations
#' @param spde the Stochastic Partial Differential Equations model
#' @param xcoord field name of x coordinate
#' @param ycoord field name of y coordinate
#' @... additional argument to INLA
#' 
#' @return fitted INLA object
#' 
spde.gam <- function(
    formula,data,mesh,spde,xcoord='x',ycoord='y',
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
  
  obs <- nmiss(formula,data)
  coords <- as.matrix(obs[,c(xcoord,ycoord)])
  
  ## process formula
  l.vars <- getVars(formula)
  zcoord <- l.vars[1]
  if(nchar(zcoord)<1){
    stop("please specify a response in formula.\n")
  }
  
  #browser()
  ## extract fixed effects and random error terms
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
  is_f <- grep("^f\\(.+\\)$",terms_labels)
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
  
  ## check for s terms
  is_s <- grep("^s\\(.+\\)$",terms_labels)
  sout <- NULL
  terms_s <- NULL
  for_char <- ~1
  if(length(is_s)){
    terms_s <- terms_vars[is_s]
    
    ## deal with s terms
    fixed_char <- terms_labels[is_s]
    resp_char <- as.character(formula)[2]
    for_char <- reformulate(fixed_char,response=resp_char)
    #paste(c(resp_char,"~",fixed_char),collapse = "")
    sout <- s2inla(formula(for_char),obs)
    sout$X <- as.list(sout$data[terms_s])
    
    ## update fixed
    terms_fixed <- attr(terms(fixed),"term.labels")
    if(length(terms_fixed[-is_s])){
      fixed <- reformulate(terms_fixed[-is_s])
    }else{
      fixed <- reformulate("1")
    }
  }
  
  #browser()
  offset <- ~1
  if(!is.null(attr(terms_,"offset"))){
    offset_terms <- row.names(attr(terms_,"factors"))[attr(terms_,"offset")]
    offset <- reformulate(offset_terms)
  }
  
  #browser()
  inla_spde_stack_helper <- function(
    fixed,random,formula,obs,
    mesh,spde,xcoord,ycoord,sout,offset){
    ## make a INLA stack
    
    ## translate responses
    formula__ <- reformulate("1",response=as.character(formula)[2])
    obs_fixed_frame <- model.frame(formula__,data=obs)
    obs_fixed_y <- model.response(obs_fixed_frame)

    ## create objects for fixed and random effects observed
    z.obs.lst <- list(b0=rep(1,nrow(obs)))
    
    if(length(all.vars(fixed))){
      obs_fixed_X <- model.matrix(fixed,data=obs)[,-1,drop=F]
      obs_fixed_X <- as.data.frame(obs_fixed_X)
      names(obs_fixed_X) <- make.names(names(obs_fixed_X))
      obs_fixed_X_list <- as.list(obs_fixed_X)
      z.obs.lst <- c(z.obs.lst,obs_fixed_X_list)
    }
    if(!is.null(sout)){
      z.obs.lst <- c(z.obs.lst,sout$Z,sout$X)
    }
    if(length(all.vars(random))){
      z.obs.lst <- c(z.obs.lst,as.list(
        obs[,all.vars(random),drop=FALSE]))
    }
    if(length(all.vars(offset))){
      z.obs.lst <- c(z.obs.lst,as.list(
        obs[,all.vars(offset),drop=FALSE]
      ))
    }
    ## Projection matrices
    coords <- as.matrix(obs[,c(xcoord,ycoord)])
    A <- inla.spde.make.A(mesh, loc=coords)
    
    ## Define SPDE data stacks -- observed
    inla.stack(tag="obs",
             data=list(y=obs_fixed_y),
             A=list(A,1),effects=list(
               list(field=1:spde$n.spde),
               z.obs.lst))
  }
  
  inla_spde_formula_helper <- function(fixed,random,obs,sout,offset){
    ## Define INLA formula
    terms_fixed <- attr(terms(fixed),"term.labels")
    interceptonly <- length(terms_fixed)==0
    if(interceptonly){
      fixed.vars <- c("0+b0")
    }else{
      obs_fixed_X <- model.matrix(fixed,data=obs)[,-1,drop=F]
      obs_fixed_X <- as.data.frame(obs_fixed_X)
      names(obs_fixed_X) <- make.names(names(obs_fixed_X))
      obs_fixed_X_list <- as.list(obs_fixed_X)
      fixed.vars <- paste0(c("0","b0",names(obs_fixed_X_list)),collapse = "+")
    }
    
    terms_random <- attr(terms(random),"term.labels")
    if(length(terms_random)>0){
      random.vars <- as.character(random)[2]
    }else{
      random.vars <- NULL
    }
    ## deal with S terms
    if(is.null(sout)){
      sterms.vars <- NULL
    }else{
      sterms.vars <- as.character(sout$formula)[3]
    }
    ## deal with offset term
    if(length(all.vars(offset))){
      offset.vars <- as.character(offset)[2]
    }else{
      offset.vars <- NULL
    }
    
    spatial.vars <- c("f(field,model=spde)")
    reformulate(c(offset.vars,fixed.vars,random.vars,
                  sterms.vars,spatial.vars),
                response="y")
  }
  
  ## adjust stack data for missing Z terms
  inla_spde_adjustZ_helper <- function(stk.z,sout){
    ## NA in Z terms are imputed with some positive number
    data_ <- inla.stack.data(stk.z)
    if(!is.null(sout)){
      for(zterm in names(sout$Z)){
        Z_ <- data_[[zterm]]
        #Z_[is.na(Z_)] <- 1
        data_[[zterm]] <- na.omit(Z_)
      }
    }
    data_  
  }
  ## translate fixed effects (minus intercepts) and responses
  # formula__ <- reformulate("1",response=as.character(formula)[2])
  # obs_fixed_frame <- model.frame(formula__,data=obs)
  # obs_fixed_y <- model.response(obs_fixed_frame)
  
  ## create objects for fixed and random effects observed
  # z.obs.lst <- list(b0=rep(1,nrow(obs)))
  # 
  # if(length(all.vars(fixed))){
  #   obs_fixed_X <- model.matrix(fixed,data=obs)[,-1]
  #   obs_fixed_X <- as.data.frame(obs_fixed_X)
  #   names(obs_fixed_X) <- make.names(names(obs_fixed_X))
  #   obs_fixed_X_list <- as.list(obs_fixed_X)
  #   z.obs.lst <- c(z.obs.lst,obs_fixed_X_list)
  # }
  # if(length(all.vars(random))){
  #   z.obs.lst <- c(z.obs.lst,as.list(
  #     obs[,all.vars(random),drop=FALSE]))
  # }

  ## Projection matrices
  # coords <- as.matrix(obs[,c(xcoord,ycoord)])
  # A <- inla.spde.make.A(mesh, loc=coords)
  # 
  ## Define SPDE data stacks -- observed
  # stk.z.obs <- inla.stack(tag="obs",
  #                         data=list(y=obs_fixed_y),
  #                         A=list(A,1),
  #                         effects=list(
  #                           list(field=1:spde$n.spde),
  #                           z.obs.lst
  #                         )
  # )
  
  stk.z <- inla_spde_stack_helper(
    fixed,random,formula,obs,
    mesh,spde,xcoord,ycoord,sout,offset
  )
  #browser()
  
  ## Define INLA formula
  # terms_fixed <- attr(terms(fixed),"term.labels")
  # interceptonly <- length(terms_fixed)==0
  # if(interceptonly){
  #   fixed.vars <- c("0+b0")
  # }else{
  #   fixed.vars <- paste0(c("0","b0",names(obs_fixed_X_list)),collapse = "+")
  # }
  
  # terms_random <- attr(terms(random),"term.labels")
  # if(length(terms_random)>0){
  #   random.vars <- as.character(random)[2]
  # }else{
  #   random.vars <- NULL
  # }
  # spatial.vars <- c("f(field,model=spde)")
  # formula__ <- reformulate(c(fixed.vars,
  #                            random.vars,
  #                            spatial.vars),response="y")
  formula__ <- inla_spde_formula_helper(
    fixed,random,obs,sout,offset)
  
  if(F){
    ## debugs
    source("zplot_methods_2.R")
    mod0 <- inla(sout$formula,
                 data=c(list(TOT=data$TOT),sout$X,sout$Z),
                 family = "poisson",
                 E=data$AREA)
    summary(mod0)
    data_ <- inla.stack.data(stk.z)
    A_ <- inla.stack.A(stk.z)
    dim(A_)
    field_idx_ <- which(is.na(data_$b0))
    all(is.na(data_$Z1[field_idx_,]))
    ## how to fill in these NA values with
    ## something reasonable if possible?
    apply(A_[,field_idx],2,sum)
    
    data_$TOT <- data_$y
    Z1 <- data_$Z1
    Z1 <- na.omit(Z1)
    #Z1[is.na(Z1)] <- 1
    data_$Z1 <- Z1
    mod1 <- inla(sout$formula,
                 data=data_,
                 control.predictor=list(A=inla.stack.A(stk.z),compute=TRUE),
                 family = "poisson",
                 E=data$AREA)
    
    plot(mod0$summary.random$x1$mean,
         mod1$summary.random$x1$mean)
    abline(0,1)
  }
  ## posterior inference
  #browser()
  mod1 <- try(do.call("inla",c(list(
    formula=formula__,
    data=inla_spde_adjustZ_helper(stk.z,sout),
    control.predictor=list(A=inla.stack.A(stk.z),compute=TRUE)),
    args__)))

  ## save arguments
  mod1$.args$formula0 <- formula
  mod1$.args$offset <- offset
  mod1$.args$fixed <- fixed
  mod1$.args$random <- random
  mod1$.args$raw <- obs
  mod1$.args$sterms <- terms_s
  mod1$.args$sformula <- formula(for_char)
  mod1$.args$error <- err_term
  mod1$.args$mesh <- mesh
  mod1$.args$spde <- spde
  mod1$.args$coords <- coords
  
  ## return model
  mod1
}
spde.glm <- function(...) {
  cat("Warning: 'spde.glm' has been superseded by 'spde.gam'. Please use 'spde.gam' instead.\n")
  
  # Call the new function with the provided arguments
  spde.gam(...)
}
dev <- function(){
  ## categorical predictor as fixed effects
  rm(list=ls())
  library(INLA)
  source("sgam_2.R")
  source("nmiss_2.R")
  source("getVars.R")
  source("s2inla_2.R")
  load(file="Scratch/spde_glm_dev1.rdata")
  source("spde_gam_1.R")
  
  wds.geo$strata <- gl(3,510)[1:nrow(wds.geo)]
  wds.geo$x1 <- rnorm(nrow(wds.geo))
  ## original response, categorical
  mod <- spde.gam(
    TOT~strata,
    data=wds.geo,mesh=mesh1,spde=spde1,
    xcoord="X",ycoord="Y",family="poisson",
    E=AREA)
  ## transformed repsonse, categorical
  mod <- spde.gam(I(TOT/AREA)~strata,data=wds.geo,
                  mesh=mesh1,spde=spde1,
                  xcoord="X",ycoord="Y")
  
  ## transformed response, categorical + linear
  wds.geo$Xc <- scale(wds.geo$X)
  mod <- spde.gam(I(TOT/AREA)~strata+Xc,data=wds.geo,
                  mesh=mesh1,spde=spde1,
                  xcoord="X",ycoord="Y")

  ## transformed response, categorical + transformed
  mod <- spde.gam(I(TOT/AREA)~strata+Xc+I(Xc^2),data=wds.geo,
                  mesh=mesh1,spde=spde1,
                  xcoord="X",ycoord="Y")
  
  ## transformed response, categorical + single spline term
  library(splines)
  mod <- spde.gam(I(TOT/AREA)~strata+ns(X,df=3),data=wds.geo,
                  mesh=mesh1,spde=spde1,
                  xcoord="X",ycoord="Y")

  ## transformed response, categorical + multiple spline terms
  mod <- spde.gam(I(TOT/AREA)~strata+ns(X,df=3)+ns(Y,df=3),
                  data=wds.geo,
                  mesh=mesh1,spde=spde1,
                  xcoord="X",ycoord="Y")
  
  wds.geo$err <- seq(1,nrow(wds.geo))
  wds.geo$X2 <- as.vector(scale(wds.geo$X))
  wds.geo$Y2 <- as.vector(scale(wds.geo$Y))
  ## poly(X,2) does not work because missing values not allowed
  mod <- spde.gam(TOT~1,data=wds.geo,mesh=mesh1,spde=spde1,
                  xcoord="X",ycoord="Y",family="poisson",
                  E=AREA)
  mod <- spde.gam(TOT~f(err),data=wds.geo,mesh=mesh1,spde=spde1,
                  xcoord="X",ycoord="Y",family="poisson",
                  E=AREA)
  mod <- spde.gam(TOT~X2+I(X2^2)+f(err),data=wds.geo,mesh=mesh1,spde=spde1,
                  xcoord="X",ycoord="Y",family="poisson",
                  E=AREA)
  mod <- spde.gam(TOT~X2+Y2+f(err),data=wds.geo,mesh=mesh1,spde=spde1,
                  xcoord="X",ycoord="Y",family="poisson",
                  E=AREA)
  mod <- spde.gam(TOT~ns(X,df=2)+ns(Y,df=4)+f(err),data=wds.geo,mesh=mesh1,spde=spde1,
                  xcoord="X",ycoord="Y",family="poisson",
                  E=AREA)
  
}
test_spde_gam <- function(){
  rm(list=ls())
  library(INLA)
  library(mgcv)
  source("sgam_2.R")
  source("nmiss_2.R")
  source("getVars.R")
  source("s2inla_2.R")
  load(file="Scratch/spde_glm_dev1.rdata")
  source("spde_gam_1b.R")
  
  wds.geo$strata <- gl(3,510)[1:nrow(wds.geo)]
  wds.geo$x1 <- rnorm(nrow(wds.geo))
  wds.geo$x2 <- rnorm(nrow(wds.geo))
  mod <- spde.gam(
    TOT~x1,
    data=wds.geo,mesh=mesh1,spde=spde1,
    xcoord="X",ycoord="Y",family="poisson",
    E=AREA)
  
  library(splines)
  mod <- spde.gam(
    TOT~ns(x1),
    data=wds.geo,mesh=mesh1,spde=spde1,
    xcoord="X",ycoord="Y",family="poisson",
    E=AREA)
  
  ## original response, categorical
  mod <- spde.gam(
    TOT~strata+s(x1),
    data=wds.geo,mesh=mesh1,spde=spde1,
    xcoord="X",ycoord="Y",family="poisson",
    E=AREA)
  ## transformed repsonse, categorical
  source("lme_inla_formula.R")
  source("proc.formula.R")
  wds.geo$den <- wds.geo$TOT/wds.geo$AREA
  source("zplot_methods_2.R")
  mod <- sgam(den~s(x1),data=wds.geo)
  plot(mod)
  mod <- spde.gam(I(TOT/AREA)~s(x1),data=wds.geo,
                  mesh=mesh1,spde=spde1,
                  xcoord="X",ycoord="Y")
  plot(mod)
  ## transformed response, categorical + linear
  wds.geo$Xc <- scale(wds.geo$X)
  mod <- spde.gam(I(TOT/AREA)~strata+Xc+s(x1),
                  data=wds.geo,
                  mesh=mesh1,spde=spde1,
                  xcoord="X",ycoord="Y")
  
  ## transformed response, categorical + transformed
  mod <- spde.gam(I(TOT/AREA)~strata+Xc+I(Xc^2)+s(x1),
                  data=wds.geo,
                  mesh=mesh1,spde=spde1,
                  xcoord="X",ycoord="Y")
  
  ## transformed response, categorical + single spline term
  library(splines)
  mod <- spde.gam(I(TOT/AREA)~strata+ns(X,df=3)+s(x1),
                  data=wds.geo,
                  mesh=mesh1,spde=spde1,
                  xcoord="X",ycoord="Y")
  
  ## transformed response, categorical + multiple spline terms
  mod <- spde.gam(
    I(TOT/AREA)~strata+ns(X,df=3)+ns(Y,df=3)+s(x1),
    data=wds.geo,mesh=mesh1,spde=spde1,
    xcoord="X",ycoord="Y")
  
  wds.geo$err <- seq(1,nrow(wds.geo))
  wds.geo$X2 <- as.vector(scale(wds.geo$X))
  wds.geo$Y2 <- as.vector(scale(wds.geo$Y))
  ## poly(X,2) does not work because missing values not allowed
  mod <- spde.gam(TOT~s(x1),data=wds.geo,mesh=mesh1,spde=spde1,
                  xcoord="X",ycoord="Y",family="poisson",
                  E=AREA)
  mod <- spde.gam(TOT~s(x1)+f(err),data=wds.geo,mesh=mesh1,spde=spde1,
                  xcoord="X",ycoord="Y",family="poisson",
                  E=AREA)
  mod <- spde.gam(TOT~X2+I(X2^2)+s(x1)+f(err),data=wds.geo,mesh=mesh1,spde=spde1,
                  xcoord="X",ycoord="Y",family="poisson",
                  E=AREA)
  mod <- spde.gam(TOT~X2+Y2+s(x1)+f(err),data=wds.geo,mesh=mesh1,spde=spde1,
                  xcoord="X",ycoord="Y",family="poisson",
                  E=AREA)
  mod <- spde.gam(TOT~ns(X,df=2)+ns(Y,df=4)+s(x1)+f(err),data=wds.geo,mesh=mesh1,spde=spde1,
                  xcoord="X",ycoord="Y",family="poisson",
                  E=AREA)
  
}
test_spde_gam_2 <- function(){
  rm(list=ls())
  library(INLA)
  library(mgcv)
  source("sgam_2.R")
  source("nmiss_2.R")
  source("getVars.R")
  source("s2inla_2.R")
  load(file="Scratch/spde_glm_dev1.rdata")
  source("spde_gam_1.R")
  
  wds.geo$strata <- gl(3,510)[1:nrow(wds.geo)]
  wds.geo$x1 <- rnorm(nrow(wds.geo))
  wds.geo$x2 <- rnorm(nrow(wds.geo))
  ## original response, categorical
  mod <- spde.gam(
    TOT~strata+s(x1)+s(x2),
    data=wds.geo,mesh=mesh1,spde=spde1,
    xcoord="X",ycoord="Y",family="poisson",
    E=AREA)
  ## transformed repsonse, categorical
  source("lme_inla_formula.R")
  source("proc.formula.R")
  wds.geo$den <- wds.geo$TOT/wds.geo$AREA
  source("zplot_methods_2.R")
  mod <- sgam(den~s(x1)+s(x2),data=wds.geo)
  plot(mod,pages = 1)
  mod <- spde.gam(I(TOT/AREA)~s(x1)+s(x2),data=wds.geo,
                  mesh=mesh1,spde=spde1,
                  xcoord="X",ycoord="Y")
  plot(mod,pages = 1)
  ## transformed response, categorical + linear
  wds.geo$Xc <- scale(wds.geo$X)
  mod <- spde.gam(I(TOT/AREA)~strata+Xc+s(x1)+s(x2),
                  data=wds.geo,
                  mesh=mesh1,spde=spde1,
                  xcoord="X",ycoord="Y")
  plot(mod,pages = 1)
  
  ## transformed response, categorical + transformed
  mod <- spde.gam(I(TOT/AREA)~strata+Xc+I(Xc^2)+s(x1)+s(x2),
                  data=wds.geo,
                  mesh=mesh1,spde=spde1,
                  xcoord="X",ycoord="Y")
  plot(mod,pages = 1)
  
  ## transformed response, categorical + single spline term
  library(splines)
  mod <- spde.gam(I(TOT/AREA)~strata+ns(X,df=3)+s(x1)+s(x2),
                  data=wds.geo,
                  mesh=mesh1,spde=spde1,
                  xcoord="X",ycoord="Y")
  plot(mod,pages = 1)
  
  ## transformed response, categorical + multiple spline terms
  mod <- spde.gam(
    I(TOT/AREA)~strata+ns(X,df=3)+ns(Y,df=3)+s(x1)+s(x2),
    data=wds.geo,mesh=mesh1,spde=spde1,
    xcoord="X",ycoord="Y")
  plot(mod,pages = 1)
  
  wds.geo$err <- seq(1,nrow(wds.geo))
  wds.geo$X2 <- as.vector(scale(wds.geo$X))
  wds.geo$Y2 <- as.vector(scale(wds.geo$Y))
  ## poly(X,2) does not work because missing values not allowed
  mod <- spde.gam(TOT~s(x1)+s(x2),data=wds.geo,mesh=mesh1,spde=spde1,
                  xcoord="X",ycoord="Y",family="poisson",
                  E=AREA)
  mod <- spde.gam(TOT~s(x1)+s(x2)+f(err),data=wds.geo,mesh=mesh1,spde=spde1,
                  xcoord="X",ycoord="Y",family="poisson",
                  E=AREA)
  mod <- spde.gam(TOT~X2+I(X2^2)+s(x1)+f(err),data=wds.geo,mesh=mesh1,spde=spde1,
                  xcoord="X",ycoord="Y",family="poisson",
                  E=AREA)
  mod <- spde.gam(TOT~X2+Y2+s(x1)+s(x2)+f(err),data=wds.geo,mesh=mesh1,spde=spde1,
                  xcoord="X",ycoord="Y",family="poisson",
                  E=AREA)
  mod <- spde.gam(TOT~ns(X,df=2)+ns(Y,df=4)+s(x1)+s(x2)+f(err),data=wds.geo,mesh=mesh1,spde=spde1,
                  xcoord="X",ycoord="Y",family="poisson",
                  E=AREA)
  
}
test_spde_glm <- function(){
  rm(list=ls())
  library(INLA)
  library(mgcv)
  source("sgam_2.R")
  source("nmiss_2.R")
  source("getVars.R")
  source("s2inla_2.R")
  load(file="Scratch/spde_glm_dev1.rdata")
  source("spde_gam_1b.R")
  
  wds.geo$strata <- gl(3,510)[1:nrow(wds.geo)]
  wds.geo$x1 <- rnorm(nrow(wds.geo))
  wds.geo$x2 <- rnorm(nrow(wds.geo))
  ## original response, categorical
  mod <- spde.glm(
    TOT~strata+x1+x2,
    data=wds.geo,mesh=mesh1,spde=spde1,
    xcoord="X",ycoord="Y",family="poisson",
    E=AREA)
}
