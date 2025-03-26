setGeneric("link",function(fit,data,n=1000,type="confidence",E=NULL,Ntrials=NULL,...){
  message( paste0("No link method for object of class '",
                  class(fit),"'. Returning instead.") )
  return(0)
})

inla_reset_config <- function(fit){
  compute_ <- fit$.args$control.compute
  if(!compute_$config){
    #browser()
    compute_$config <- TRUE
    char_ <- as.character(fit$.args$formula)
    for_ <- as.formula(paste0(char_[2],"~",char_[3]))
    fit_new <- inla(formula = for_,
                    data=fit$.args$data,
                    family = fit$.args$family,
                    E=fit$.args$E,
                    Ntrials=fit$.args$Ntrials,
                    control.compute = compute_, 
                    control.predictor = fit$.args$control.predictor, 
                    control.family = fit$.args$control.family, 
                    control.inla = fit$.args$control.inla, 
                    control.fixed = fit$.args$control.fixed, 
                    control.mode = fit$.args$control.mode, 
                    control.expert = fit$.args$control.expert, 
                    control.lincomb = fit$.args$control.lincomb, 
                    control.update = fit$.args$control.update, 
                    control.lp.scale = fit$.args$control.lp.scale, 
                    control.pardiso = fit$.args$control.pardiso)
    return(fit_new)
  }
  return(fit)
}
setMethod("link","inla",function(fit,data,n=1000,type="confidence",E=NULL,Ntrials=NULL,...){
  #browser()
  stopifnot(type %in% c("confidence","prediction"))
  
  ## generate posterior samples from fitted model
  if(missing(data)){
    ## if no data, use the training samples
    fit_new <- inla_reset_config(fit) ## reset config=TRUE
    samples <- inla.posterior.sample(n=n,result = fit_new)
    
  }else{
    ## extract model formula and training data
    if(is.null(fit$.args$formula0)){
      ## fit from standard inla
      char_ <- as.character(fit$.args$formula)
      for_ <- as.formula(paste0(char_[2],"~",char_[3]))
      raw_ <- fit$.args$data
    }else{
      ## fit from sgam
      for_ <- fit$.args$formula0
      raw_ <- fit$.args$raw
    }
    ## combine new data with the existing data
    data_ <- sgam.combine(for_,raw_,data)
    
    ## translate the sformula
    internal__ <- sgam.translate(formula = for_,
                                 data = data_,
                                 param = fit$.args$param)

    ## refit the model for prediction
    compute_ <- fit$.args$control.compute
    compute_$config <- TRUE ## allow simulation
    predictor_ <- fit$.args$control.predictor
    predictor_$link <- NULL
    
    ## set the E for Poisson
    E_ <- fit$.args$E
    if(!is.null(fit$.args$E)){
      np <- length(data_[[1]]) - length(E_)
      E_ <- c(fit$.args$E,rep(1,np))
    }
    
    ## set the Ntrials for Binomial
    Ntrials_ <- fit$.args$Ntrials
    if(!is.null(fit$.args$Ntrials)){
      np <- length(data_[[1]]) - length(Ntrials_)
      Ntrials_ <- c(fit$.args$Ntrials, rep(1,np))
    }
    #browser()
    args_ <- list(formula = internal__$value$formula,
                  data=internal__$data,
                  family = fit$.args$family,
                  E=E_,
                  Ntrials=Ntrials_,
                  control.compute = compute_, 
                  control.predictor = predictor_, 
                  control.family = fit$.args$control.family, 
                  control.inla = fit$.args$control.inla, 
                  control.fixed = fit$.args$control.fixed, 
                  control.mode = fit$.args$control.mode, 
                  control.expert = fit$.args$control.expert, 
                  control.lincomb = fit$.args$control.lincomb, 
                  control.update = fit$.args$control.update, 
                  control.lp.scale = fit$.args$control.lp.scale, 
                  control.pardiso = fit$.args$control.pardiso)
    if(!is.null(E_)){
      args_$E <- E_
    }
    if(!is.null(Ntrials_)){
      args_$Ntrials <- Ntrials_
    }
    fit_new <- do.call("inla",args = args_)
    # fit_new <- inla(formula = internal__$value$formula,
    #                 data=internal__$data,
    #                 family = fit$.args$family,
    #                 E=E_,
    #                 Ntrials=Ntrials_,
    #                 control.compute = compute_, 
    #                 control.predictor = predictor_, 
    #                 control.family = fit$.args$control.family, 
    #                 control.inla = fit$.args$control.inla, 
    #                 control.fixed = fit$.args$control.fixed, 
    #                 control.mode = fit$.args$control.mode, 
    #                 control.expert = fit$.args$control.expert, 
    #                 control.lincomb = fit$.args$control.lincomb, 
    #                 control.update = fit$.args$control.update, 
    #                 control.lp.scale = fit$.args$control.lp.scale, 
    #                 control.pardiso = fit$.args$control.pardiso)
    samples <- inla.posterior.sample(n=n,result = fit_new)
  }
  ## precision location
  prec.ii <- grep("Precision",names(samples[[1]]$hyperpar))
  name.ii <- gsub("Precision","SD",names(samples[[1]]$hyperpar))[prec.ii]
  ## Predictor locations
  pred.ii <- grep("Predictor",row.names(samples[[1]]$latent))
  if(!missing(data)){
    ## new data, remove training data
    ntrain <- nrow(fit$summary.fitted.values)
    pred.ii <- pred.ii[-c(1:ntrain)]
  }
  ## define link function
  linkinv_ <- switch(fit$.args$family,
                     poisson=poisson()$linkinv,
                     gaussian=gaussian()$linkinv,
                     binomial=binomial()$linkinv,
                     nbinomial=poisson()$linkinv,
                     gaussian()$linkinv)
  #browser()
  # extract linear predictors
  lp_matrix <- t(sapply(samples,function(elmt){
    elmt$latent[pred.ii,1,drop=FALSE]}))
  
  #browser()
  fam_ <- fit$.args$family
  if(type=="confidence"){
    extract <- function(elmt){
      #sigma <- sqrt(elmt$hyperpar[prec.ii][1])
      mu <- elmt$latent[pred.ii,1,drop=FALSE]
      #rnorm(length(mu),mu,sigma)
      linkinv_(mu)
    }
    post <- sapply(samples,extract)
    
    if(fam_=="binomial"){
      if(!is.null(Ntrials)){
        N_ <- cbind(data,Ntrials)[,ncol(data)+1]
        post <- sweep(post,1,STATS = N_,FUN = "*")
      }
    }
    if(fam_ %in% c("poisson","nbinomial")){
      if(!is.null(E)){
        E_ <- cbind(data,E)[,ncol(data)+1]
        post <- sweep(post,1,STATS = E_,FUN = "*")
      }
    }
  }else{
    #browser()
    if(fam_=="gaussian"){
      sim_work <- function(elmt){
        mu <- elmt$latent[pred.ii]
        sigma <- 1/sqrt(elmt$hyperpar["Precision for the Gaussian observations"])
        rnorm(length(mu),mean=mu,sd=sigma)
      }
      post <- sapply(samples,sim_work)
    }
    if(fam_=="binomial"){
      m_matrix <- linkinv_(lp_matrix)
      ## extract N trials
      if(!is.null(Ntrials)){
        N_ <- cbind(data,Ntrials)[,ncol(data)+1]
      }else{
        N_ <- 1
      }
      
      post <- apply(m_matrix,1,rbinom,
                           n=length(pred.ii),
                           size=N_)
    }
    if(fam_=="poisson"){
      #browser()
      ## extract expected
      m_matrix <- linkinv_(lp_matrix)
      if(!is.null(E)){
        E_ <- cbind(data,E)[,ncol(data)+1]
        m_matrix <- sweep(m_matrix,2,STATS = E_,FUN = "*")
      }
      post <- apply(m_matrix,1,rpois,n=length(pred.ii))
    }
    
    if(fam_=="nbinomial"){
      ## extract expected
      m_matrix <- linkinv_(lp_matrix)
      if(!is.null(E)){
        E_ <- cbind(data,E)[,ncol(data)+1]
        m_matrix <- sweep(m_matrix,2,STATS = E_,FUN = "*")
      }
      
      #inla.likelihood.parser
      hyper_ <- sapply(samples,function(elmt) 
        elmt$hyper["size for the nbinomial observations (1/overdispersion)"])
      
      post <- sapply(1:nrow(m_matrix),function(i){
        size_ <- (hyper_[i])
        p_ <- size_/(size_+m_matrix[i,])
        rnbinom(n=length(pred.ii),
                size=size_,prob = p_)
      })
    }
  }
  
  
  if(is.null(dim(post))){
    post <- matrix(post,ncol=1)
  }else{
    post <- t(post)  
  }
  
  as.data.frame(post)
})

dev_Gaus <- function(){
  rm(list=ls())
  library(INLA)
  library(mgcv)
  sim <- gamSim(n=100)
  fit <- inla(y~x1,data=sim)
  source("s2inla_2.R")
  source("getVars.R")
  source("nmiss_2.R")
  source("proc.formula.R")
  source("lme_inla_formula.R")
  source("sgam_2.R")
  source("Z_link_inla_2c.R")
  summary(fit)
  mu <- link(fit,data=sim[1:10,])
  post <- link(fit,data = sim[1:10,],type="prediction")
  plot(apply(mu,2,mean),fit$summary.fitted.values[1:10,"mean"])
  plot(apply(mu,2,mean),apply(post,2,mean))
  plot(apply(mu,2,sd),apply(post,2,sd))
  1/sqrt(fit$summary.hyperpar)
}

dev_Binom <- function(){
  rm(list=ls())
  DeerEcervi <- read.csv(file="DeerEcervi.csv")
  ## Exploratory Data Analyses
  DeerEcervi$fSex <- factor(DeerEcervi$Sex,labels=c("F","M"))
  library(ggplot2)
  ggplot(aes(x=Length,y=Ecervi.01,color=fSex),data=DeerEcervi)+geom_point()+
    geom_smooth(method='glm',se=FALSE,method.args = list(family = "binomial"))
  
  ## center the predictor
  DeerEcervi$Length.cent <- DeerEcervi$Length-mean(DeerEcervi$Length)
  fit <- inla(Ecervi.01~0+fSex+Length.cent,data=DeerEcervi,
                    family = "binomial",Ntrials = 1)
 
  source("s2inla_2.R")
  source("getVars.R")
  source("nmiss_2.R")
  source("proc.formula.R")
  source("lme_inla_formula.R")
  source("sgam_2.R")
  source("Z_link_inla_2c.R") 
  
  mu <- link(fit,data=DeerEcervi[1:10,],Ntrials = 1:10)
  post <- link(fit,data=DeerEcervi[1:10,],type = "prediction",Ntrials = 1:10)
  plot(apply(mu,2,mean)/1:10,fit$summary.fitted.values[1:10,"mean"])
  abline(0,1)
  plot(apply(mu,2,mean),apply(post,2,mean))
  abline(0,1)
  plot(apply(mu,2,sd),apply(post,2,sd))
  table(post[,1])  
}

dev_Poisson <- function(){
  rm(list=ls())
  #library(rethinking)
  library(INLA)
  counts <- c(18,17,15,20,10,20,25,13,12)
  outcome <- gl(3,1,9)
  treatment <- gl(3,3)
  dat <- data.frame(treatment, outcome, counts) # showing data
  ctrl_ <-  list(dic = TRUE, waic = TRUE, config = TRUE)
  fit <- inla(counts ~ outcome + treatment, data=dat,family = "poisson",
                  control.compute = ctrl_,E=1:9)
  source("s2inla_2.R")
  source("getVars.R")
  source("nmiss_2.R")
  source("proc.formula.R")
  source("lme_inla_formula.R")
  source("sgam_2.R")
  source("Z_link_inla_2c.R") 
  mu <- link(fit,data=dat,E=1:9)
  post <- link(fit,data=dat,type = "prediction",E=1:9)
  plot(apply(mu,2,mean)/1:9,fit$summary.fitted.values$mean)
  abline(0,1)
  plot(apply(mu,2,mean),apply(post,2,mean))
  abline(0,1)
  plot(apply(mu,2,sd),apply(post,2,sd))
  table(post[,1])  
  
}

dev_NB <- function(){
  rm(list=ls())
  library(rethinking)
  library(INLA)
  require(foreign)
  require(ggplot2)
  require(MASS)
  dat <- read.dta("https://stats.idre.ucla.edu/stat/stata/dae/nb_data.dta")
  dat <- within(dat, {
    prog <- factor(prog, levels = 1:3, labels = c("General", "Academic", "Vocational"))
    id <- factor(id)
  })
  ctrl_ <-  list(dic = TRUE, waic = TRUE, config = TRUE)
  fit <- inla(daysabs ~ math + prog, data = dat,family = "nbinomial",
             control.compute = ctrl_)#),E = median(dat$daysabs))
  #m1$waic$waic
  source("s2inla_2.R")
  source("getVars.R")
  source("nmiss_2.R")
  source("proc.formula.R")
  source("lme_inla_formula.R")
  source("sgam_2.R")
  source("Z_link_inla_2c.R") 
  
  
  mu <- link(fit,data=dat[1:10,],E=1:10)
  post <- link(fit,data=dat[1:10,],type = "prediction",E=1:10)
  plot(apply(mu,2,mean)/1:10,fit$summary.fitted.values$mean[1:10])
  abline(0,1)
  plot(apply(mu,2,mean),apply(post,2,mean))
  abline(0,1)
  plot(apply(mu,2,sd),apply(post,2,sd))
  table(post[,1])  
}

