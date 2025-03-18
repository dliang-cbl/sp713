#' generalized linear model 
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

sglm <- function(
  formula,data,newdata=NULL,family="gaussian",
  E=NULL,Ntrials=NULL,param=NULL,weights=NULL,
  verbose=FALSE)
{
  #browser()
  ## get response (if any) and predictor variables
  l.vars <- getVars(formula)
  stopifnot(nchar(l.vars)[1]>0)
  stopifnot(all(l.vars %in% names(data)))
  if(!is.null(newdata)){
    if(!(l.vars[1] %in% colnames(newdata))){
      newdata <- cbind(newdata,NA)
      colnames(newdata)[ncol(newdata)] <- l.vars[1]
    }else{
      newdata[,l.vars[1]] <- NA
    }
    data <- rbind(data[,l.vars,drop=FALSE],
                  newdata[,l.vars,drop=FALSE])
  }
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
                                                                              
  ## deal with E and NTrials
  if(!is.null(E)){
    result0 <- inla(sout$formula,data=data_,
                    family=family,E=E,weights=weights,
                    control.predictor = list(link=1),
                    control.compute=list(dic=TRUE,config=TRUE),
                    verbose=verbose)
  }
  if(!is.null(Ntrials)){
    result0 <- inla(sout$formula,data=data_,
                    family=family,Ntrials = Ntrials,weights=weights,
                    control.predictor = list(link=1),
                    control.compute=list(dic=TRUE,config=TRUE),
                    verbose=verbose)
  }
  if(is.null(E)&is.null(Ntrials)){
    result0 <- inla(sout$formula,data=data_,
                    family=family,weights=weights,
                    control.predictor = list(link=1),
                    control.compute=list(dic=TRUE,config=TRUE),
                    verbose=verbose)
  }
  
  result0$formula <- formula
  result0$DIC <- with(result0$dic,c(Dbar=mean.deviance,pD=p.eff,DIC=dic))
  result0$.args$sterms <- sout$sterms
  result0$.args$raw <- sout$raw
  class(result0) <- c("foo",class(result0))
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
  library(INLA)
  source("lme_inla_formula.R")
  source("getVars.R")
  source("proc.formula.R")
  source("nmiss.R")
  source("s2inla.R")
  source("plot_foo.R")
  source("sglm.R")
  load(file="scotland.rData")
  Scotland$x2 <- rnorm(nrow(Scotland))
  #Scotland$X/10
  formula = Counts~f(Region,model="bym2", 
                     graph = "scotland.graph",
                     scale.model = TRUE,constr = TRUE,
                     hyper=list(phi=list(prior="pc",param=c(0.5,2/3),initial=-3),
                                prec=list(prior="logtnormal",param=c(0,0.0001))
                     ))+ s(x2,k=3)
  #getVars(formula)
  #debug(sglm)
  library(mgcv)
  #row.names(Scotland) <- NULL
  m1 <- gam(Counts~s(x2),data=Scotland)
  #library(splines)
  glmmPC = sglm(Counts~s(x2,k=3),family="gaussian",
                data=Scotland)
  summary(glmmPC)
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
  source("nmiss.R")
  source("s2inla.R")
  source("plot_foo.R")
  dat <- gamSim(6,n=200,scale=.2,dist="poisson")
  source("sglm.R")
  library(INLA)
  helper <- function(x,y){
    print(summary(x)$fixed)
    print(summary(x)$hyper)
    print(summary(y))
    print(summary(y$lme))
  }
  m1 <- sglm(y~1,data=dat,family = "poisson")
  m1b <- glm(y~1,data=dat,family=poisson)
  helper(m1,m1b)
  
  m2 <- sglm(y~x1,data=dat,family = "poisson")
  m2b <- glm(y~x1,data=dat,family=poisson)
  helper(m2,m2b)
  
  m3 <- sglm(y~x1+x2,data=dat,family = "poisson")
  m3b <- glm(y~x1+x2,data=dat,family=poisson)
  helper(m3,m3b)
  
  m4 <- sglm(y~s(x1),data=dat,family = "poisson")
  plot(m4)
  m4b <- gam(y~s(x1),data=dat,family = poisson)
  plot(m4b)
  
  m5 <- sglm(y~s(x1)+s(x2),data=dat,family = "poisson")
  plot(m5,pages=1)
  m5b <- gam(y~s(x1)+s(x2),data=dat,family = "poisson")
  plot(m5b,pages=1)
  
  m6 <- sglm(y~(1|fac),data=dat,family = "poisson")
  m6b <- gamm(y~1,random=list(fac=~1),data=dat,family = poisson)
  helper(m6,m6b)
  
  
  m7 <- sglm(y~x1+(1|fac),data=dat,family = "poisson")
  m7b <- gamm(y~x1,random=list(fac=~1),data=dat,
              family = poisson)
  helper(m7,m7b)

  m8 <- sglm(y~s(x1)+(1|fac),data=dat,family = "poisson")
  m8b <- gamm(y~s(x1),random=list(fac=~1),data=dat,
              family = poisson)
  helper(m8,m8b)
  plot(m8)
  plot(m8b$gam)
  
  m9 <- sglm(y~s(x1)+s(x2)+(1|fac),data=dat,family = "poisson")
  m9b <- gamm(y~s(x1)+s(x2),random=list(fac=~1),data=dat,
              family = poisson)
  plot(m9)
  plot(m9b$gam)
  helper(m9,m9b)
}
