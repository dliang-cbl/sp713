#' generalized additive model 
#' 
#' GLM using R and Integrated Nested Laplace Approximation. 
#' Dong Liang, UMCES/CBL 2024
#' @param formula  lme type formula object with additive terms
#' @param data data frame to evaluate the model
#' @param newdata data frame to form prediction
#' @param family family for response variable
#' @param E offset term
#' @param Ntrials number of trials INLA term
#' @param param hyper-parameter for random effects
#' @param weights weights INLA term
#' @param verbose whether to show detailed results
#' @return fitted INLA object
#' 
#' @export
sgam.translate <- function(formula,data,param){
  ## translate a formula into the internal objects
  l.vars <- getVars(formula)
  stopifnot(nchar(l.vars)[1]>0)
  stopifnot(all(l.vars %in% names(data)))
  r0 <- lme.inla.formula(formula=formula,data=data,param=param)
  no_fixed <- length(all.vars(r0$fixed))==1
  no_random <- length(r0$formulas)==0
  if(no_fixed){
    if(no_random){
      for0 <- r0$fixed
    }else{
      for0 <- reformulate(r0$formulas,response=all.vars(formula)[1])
    }
  }
  else{
    for0 <- reformulate(r0$formulas,response=all.vars(formula)[1])
  }
  df0 <- data
  if(!is.null(r0$index)){
    df0 <- cbind(data,r0$index)
  }
  
  #browser()
  sout <- s2inla(for0,data=df0)
  #for(i in 1:length(sout$Z)){
  #  assign(names(sout$Z)[i],sout$Z[[i]])
  #}
  #assign(sout$zmat,sout$X)
  data_ <- c(as.list(sout$data),sout$Z)
  
  list(value=sout,data=data_)
}
test_translate <- function(){
  rm(list=ls())
  library(mgcv)
  source("s2inla_2.R")
  source("getVars.R")
  source("nmiss_2.R")
  source("proc.formula.R")
  source("lme_inla_formula.R")
  source("sgam_2.R")
  sim <- gamSim()
  sim$y[1] <- NA
  sim2 <- list(y=sim$y,x1=sim$x1,Z=matrix(rnorm(1200),ncol=3),
               fac=rep(1:4,100))
  sim2$Z[cbind(c(2,4,7),1)] <- NA
  dum <- sgam.translate(y~s(x1)+Z+(1|fac),data=sim2,param = NULL)
  (dum$value$formula)
  str(dum$data)
}
sgam.combine <- function(formula,data,newdata){
  ## combine data with new data and according to the formula
  stopifnot(class(data)==class(newdata))
  vars_ <- getVars(formula)
  resp.var <- vars_[1]
  newdata[[resp.var]] <- rep(NA,length(newdata[[1]])) ## set response as missing
  if(class(data)=="data.frame"){
    return(rbind(data[,vars_],newdata[,vars_]))
  }
  if(class(data)=="list"){
    data_ <- data[vars_]
    newdata_ <- newdata[vars_]
    value_ <- vector("list",length(data_))
    names(value_) <- names(data_)
    for(i in 1:length(value_)){
      if(is.null(dim(data_[[i]]))){
        value_[[i]] <- c(data_[[i]],newdata_[[i]])
      }else{
        value_[[i]] <- rbind(data_[[i]],newdata_[[i]])
      }
    }
    return(value_)
  }
  stop("wrong input type.")
}
test_combine <- function(){
  rm(list=ls())
  library(mgcv)
  source("getVars.R")
  source("sgam_2.R")
  sim <- gamSim()
  sim$y[1] <- NA
  sim2 <- list(y=sim$y,x1=sim$x1,Z=matrix(rnorm(1200),ncol=3),
               fac=rep(1:4,100))
  sim2$Z[cbind(c(2,4,7),1)] <- NA
  sim3 <- list(y=NA,x1=2,Z=sim2$Z[1,drop=F])
  r <- sgam.combine(y~s(x1)+Z,data=sim2,newdata = sim3)
  tail(r$Z)
}
sgam.config.compute <- function(args){
  ## adjustment of the argument list 
  ## to allow posterior simulation
  
  args__ <- args
  #ctrl_compute_ <- list(dic=TRUE,config=TRUE)
  if("control.compute" %in% names(args__)){
    if("config" %in% names(args__[["control.compute"]])){
      if(!(args__[["control.compute"]]$config)){
        warning("reset config\n")
        args__[["control.compute"]]$config <- TRUE
      }
    }else{
      args__[["control.compute"]]$config <- TRUE
    }
  }else{
    args__ <- c(args__,list(control.compute=list(waic=TRUE,config=TRUE)))
  }
  
  return(args__)
  
}
sgam <- function(
  formula,data,newdata=NULL,param=NULL,...)
{
  ## combine data
  if(!is.null(newdata)){
    data <- sgam.combine(formula,data,newdata)
  }
  
  ## translate gamm formula into a list
  internal__ <- sgam.translate(formula = formula,
                               data = data,
                               param = param)
  #sout <- internal__$value
  #data_ <- internal__$data
  #browser()
  ## reset control.compute if provided
  args__ <- sgam.config.compute(list(...))
  # #ctrl_compute_ <- list(dic=TRUE,config=TRUE)
  # if("control.compute" %in% names(args__)){
  #   if("config" %in% names(args__[["control.compute"]])){
  #     if(!(args__[["control.compute"]]$config)){
  #       warning("reset config\n")
  #       args__[["control.compute"]]$config <- TRUE
  #     }
  #   }else{
  #     args__[["control.compute"]]$config <- TRUE
  #   }
  # }else{
  #   args__ <- c(args__,list(control.compute=list(waic=TRUE,config=TRUE)))
  # }
  # 
  ## posterior inference
  result0 <- do.call(
    "inla",
    c(list(formula=internal__$value$formula,
         data=internal__$data),
         args__))
  
  # ## deal with E and NTrials
  # if(!is.null(E)){
  #   result0 <- inla(sout$formula,data=data_,
  #                   family=family,E=E,weights=weights,
  #                   control.predictor = list(link=1),
  #                   control.compute=list(dic=TRUE,config=TRUE),
  #                   verbose=verbose)
  # }
  # if(!is.null(Ntrials)){
  #   result0 <- inla(sout$formula,data=data_,
  #                   family=family,Ntrials = Ntrials,weights=weights,
  #                   control.predictor = list(link=1),
  #                   control.compute=list(dic=TRUE,config=TRUE),
  #                   verbose=verbose)
  # }
  # if(is.null(E)&is.null(Ntrials)){
  #   result0 <- inla(sout$formula,data=data_,
  #                   family=family,weights=weights,
  #                   control.predictor = list(link=1),
  #                   control.compute=list(dic=TRUE,config=TRUE),
  #                   verbose=verbose)
  # }
  

  #result0$DIC <- with(result0$dic,c(Dbar=mean.deviance,pD=p.eff,DIC=dic))
  result0$.args$sterms <- internal__$value$sterms
  result0$.args$raw <- internal__$value$raw
  result0$.args$formula0 <- formula
  result0$.args$param <- param
  #result0$.args$weights <- weights
  #result0$.args
  #class(result0) <- c("foo",class(result0))
  result0
}
report_ <- function(){
  rm(list=ls())
  library(INLA)
  library(mgcv)
  load(file="scotland2.rData")
  
  unique(data_$x2)
  dim(X[[1]][[1]])
  
  Z <- X[[1]][[1]]
  formula <- Counts ~ f(
    x2, model = "z",Z = Z,
    hyper = list(prec = list(prior = "logtnormal",
                             param = c(0, 1e-05))))
  result0 <- inla(formula,data=data_,
                  family="poisson",E=data_$E,
                  control.predictor = list(link=1),
                  control.compute=list(dic=TRUE,config=TRUE),
                  verbose=F)  
  
  formula2 <- Counts ~ f(
    x2, model = "z",Z = X[[1]][[1]],
    hyper = list(prec = list(prior = "logtnormal",
                             param = c(0, 1e-05))))
  
  try(result1 <- inla(formula2,data=data_,
                  family="poisson",E=data_$E,
                  control.predictor = list(link=1),
                  control.compute=list(dic=TRUE,config=TRUE),
                  verbose=F) )
  
  oo <- order(data_$x2.1)
  matplot(data_$x2.1[oo],
          result0$summary.random$x2[1:56,4:6][oo,],
          type="l",lty=c(2,1,2),col=1)

  ## GAM fit
  mod0 <- gam(Counts~s(x2.1,k=3),data=data_,family="poisson")
  plot(mod0)
  
}
debug_ <- function(){
  rm(list=ls())
  library(mgcv)
  source("lme_inla_formula.R")
  source("getVars.R")
  source("proc.formula.R")
  source("nmiss_2.R")
  source("s2inla_2.R")
  source("zplot_methods.R")
  source("sgam_2.R")
  library(INLA)
  load(file="scotland.rData")
  Scotland$x2 <- Scotland$X/10
  #Scotland$X/10
  
  #getVars(formula)
  #debug(sglm)
  library(mgcv)
  #row.names(Scotland) <- NULL
  m1 <- gam(Counts~s(x2,k=4),data=Scotland)
  #library(splines)
  glmmPC = sgam(Counts~s(x2,k=4),family="gaussian",
                data=Scotland)
  summary(glmmPC)
  plot(glmmPC)
  formula = Counts~f(Region,model="bym2", 
                     graph = "scotland.graph",
                     scale.model = TRUE,constr = TRUE,
                     hyper=list(phi=list(prior="pc",param=c(0.5,2/3),initial=-3),
                                prec=list(prior="logtnormal",param=c(0,0.0001))
                     ))+ s(x2,k=3)
  glmmPC2 = sgam(formula,family="gaussian",
                data=Scotland)
  summary(glmmPC2)
  plot(glmmPC2)
  plot(m1)
  source("extract_ll_methods_2.R")
  source("zWAIC-methods.R")
  source("zLOO-methods.R")
  source("compare.r")
  l1 <- extract_log_lik(glmmPC2)
  WAIC(l1)
  PSIS(l1)
  l2 <- l1
  compare(l1,l2)
  source("Z_link_inla_2b.R")
  p1 <- link(glmmPC2)
  plot(Scotland$Counts,p1[1,])
  p2 <- link(glmmPC2,data = Scotland)
}
sim_ <- function(){
  ## simulate a geostat poisson data
  rm(list=ls())
  library(mgcv)
  library(geoR)
  dat <- gamSim(6,n=200,scale=.2,dist="poisson")
  hist(log(dat$f1))
  sim <- (grf(nrow(dat),cov.pars = c(1,0.3)))
  points(sim)
  hist(sim$data)
  dat$x <- sim$coord[,1]
  dat$y <- sim$coord[,2]
  dat$fs <- sim$data
  mu <- with(dat,fs+log(f0)+log(f1/2)+log(f2/2))
  dat$z <- rpois(nrow(dat),lambda = exp(mu))
  save(dat,file="Scratch/sim1.rData")
}
dev2_ <- function(){
  ## add spde(x,y) to allow spatial random effects
  rm(list=ls())
  library(mgcv)
  source("lme_inla_formula.R")
  source("getVars.R")
  source("proc.formula.R")
  source("nmiss.R")
  source("s2inla.R")
  source("plot_foo.R")
  load(file="Scratch/sim1.rData")
}
dev_ <- function(){
  rm(list=ls())
  library(mgcv)
  source("lme_inla_formula.R")
  source("getVars.R")
  source("proc.formula.R")
  source("nmiss_2.R")
  source("s2inla_2.R")
  source("zplot_methods.R")
  dat <- gamSim(6,n=200,scale=.2,dist="poisson")
  source("sgam_2.R")
  library(INLA)
  # helper <- function(x,y){
  #   print(summary(x)$fixed)
  #   print(summary(x)$hyper)
  #   print(summary(y))
  #   print(summary(y$lme))
  # }
  m1 <- sgam(y~x1,data=dat,family = "poisson")
  m1d <- sgam(y~x1,data=dat,family = "poisson",
              control.compute=list(dic=FALSE,config=TRUE))
  m1c <- sgam(y~x1,data=dat,family = "poisson",
              control.compute=list(dic=FALSE))
  
  source("extract_ll_methods_2.R")
  source("zWAIC-methods.R")
  source("zLOO-methods.R")
  source("compare.r")
  l1 <- extract_log_lik(m1)
  WAIC(l1)
  PSIS(l1)
  # m1b <- glm(y~1,data=dat,family=poisson)
  # helper(m1,m1b)
  
  m2 <- sgam(y~x1,data=dat,family = "poisson")
  l2 <- extract_log_lik(m2)
  compare(l1,l2,func = PSIS)
  
  # m2b <- glm(y~x1,data=dat,family=poisson)
  # helper(m2,m2b)
  #p0 <- link(m2)
  #p1 <- link(m2,data=data.frame(x1=sample(dat$x1,100)))
  m3 <- sgam(y~x1+x2,data=dat,family = "poisson")
  m3b <- sgam(y~x1+x2,data=dat,newdata = dat[1:10,],family = "poisson")
  dim(m3b$summary.fitted.values)
  #l3 <- link(m3)
  WAIC(extract_log_lik(m3))
  WAIC(extract_log_lik(m3b))
  m3b <- glm(y~x1+x2,data=dat,family=poisson)
  helper(m3,m3b)
  
  m4 <- sgam(y~s(x1),data=dat,newdata = dat[1:10,],family = "poisson")
  #l4 <- link(m4)
  WAIC(extract_log_lik(m4))
  
  
  m8 <- sgam(y~s(x1)+(1|fac),data=dat,newdata = dat[1:10,],family = "poisson")
  
  d <- list(y=dat$y,x1=dat$x1,Z=matrix(rnorm(600),ncol=3))
  nd <- list(y=NA,x1=2,Z=matrix(rnorm(3),ncol=3))
  m8b <- sgam(y~s(x1)+Z,data=d,newdata = nd,family = "poisson")
  str(m8b$summary.fitted.values)
  WAIC(extract_log_lik(m8))

  d <- list(y=dat$y,x1=dat$x1,Z=matrix(rnorm(600),ncol=3),fac=dat$fac)
  nd <- list(y=NA,x1=2,Z=matrix(rnorm(3),ncol=3),fac=dat$fac[1])
  m8c <- sgam(y~s(x1)+Z+(1|fac),data=d,newdata = nd,family = "poisson")
  
  coeftab(m8)
  l8 <- link(m8,data=dat[1:30,])
  ll8 <- extract_log_lik(m8)
  WAIC(ll8)
  m8b <- gamm(y~s(x1),random=list(fac=~1),data=dat,
              family = poisson)
  helper(m8,m8b)
  plot(m8)
  plot(m8b$gam)
  
  m9 <- sgam(y~s(x1)+s(x2)+(1|fac),data=dat,family = "poisson")
  m9b <- gamm(y~s(x1)+s(x2),random=list(fac=~1),data=dat,
              family = poisson)
  plot(m9,pages = 1)
  plot(m9b$gam,pages = 1)
  helper(m9,m9b)
  
  
}
