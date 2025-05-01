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
inla_spde_stack <- function(fit,data){
  if(F){
    fit <- mod
    data <- wds.geo
  }
  inla_spde_stack_helper <- function(
    fixed,random,formula,obs,
    mesh,spde,xcoord,ycoord,sout,offset,isquery=F){
    ## make a INLA stack
    
    ## translate fixed effects (minus intercepts) and responses
    if(isquery){
      obs_fixed_y <- rep(NA,nrow(obs))
    }else{
      formula__ <- reformulate("1",response=as.character(formula)[2])
      obs_fixed_frame <- model.frame(formula__,data=obs)
      obs_fixed_y <- model.response(obs_fixed_frame)
    }
    
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
    
    ## Define SPDE data stacks
    if(isquery){
      tag_ <- "query"
    }else{
      tag_ <- "obs"
    }
    inla.stack(tag=tag_,
               data=list(y=obs_fixed_y),
               A=list(A,1),effects=list(
                 list(field=1:spde$n.spde),
                 z.obs.lst))
  }
  
  ## adjust stack data for missing Z terms
  inla_spde_adjustZ_helper <- function(stk.z,sout){
    ## NA in Z terms are imputed with some positive number
    data_ <- inla.stack.data(stk.z)
    if(F){
      A_ <- inla.stack.A(stk.z)
    }
    if(!is.null(sout)){
      for(zterm in names(sout$Z)){
        Z_ <- data_[[zterm]]
        #Z_[is.na(Z_)] <- 1
        data_[[zterm]] <- na.omit(Z_)
      }
    }
    data_  
  }
  #browser()
  ## extract arguments from fits
  fixed <- fit$.args$fixed
  random <- fit$.args$random
  obs <- fit$.args$raw
  formula <- fit$.args$formula0
  coords_name <- colnames(fit$.args$coords)
  mesh <- fit$.args$mesh
  spde <- fit$.args$spde
  terms_s <- fit$.args$sterms
  
  ## prepare sobjects jointly for observed
  ## and predicted
  
  ## make sure all variables are in query
  v_ <- all.vars(fit$.args$sformula)
  resp_na_ <- v_[which(!(v_ %in% names(data)))]
  if(length(resp_na_)){
    for(j in 1:length(resp_na_)){
      data$dummy__ <- 0
      names(data)[which(names(data)=="dummy__")] <- resp_na_[j]
    }
  }
  ## prepare Z terms jointly
  joint_ <- rbind(obs[,v_],data[,v_])
  soutJ <- s2inla(fit$.args$sformula,joint_)
  
  ## get the observed and query z terms
  idx.obs <- 1:nrow(obs)
  idx.query <- nrow(obs)+1:nrow(data)
  
  ## prepare the lists
  sout <- list(
    Z=lapply(soutJ$Z,function(mat) mat[idx.obs,]),
    X=as.list(soutJ$data[idx.obs,][terms_s]))
  sout.query <- list(
    Z=lapply(soutJ$Z,function(mat) mat[idx.query,]),
    X=as.list(soutJ$data[idx.query,][terms_s]-nrow(obs)))
  #sout$X <- as.list(sout$data[terms_s])
  offset <- fit$.args$offset
  ## prepare stack for the observed
  stk.z.obs <- inla_spde_stack_helper(
    fixed,random,formula,obs,mesh,spde,
    coords_name[1],coords_name[2],sout,offset,isquery = F)
  
  ## prepare stack for query
  #browser()
  #data$dummyResponse__ <- 1
  #sformula_ <- paste0("dummyResponse__~",as.character(fit$.args$sformula)[3])
  #sout.query <- s2inla(formula(sformula_),data)
  #sout.query$X <- as.list(sout.query$data[terms_s])
  
  stk.z.query <- inla_spde_stack_helper(
    fixed,random,formula,data,mesh,spde,
    xcoord = coords_name[1],ycoord = coords_name[2],
    sout.query,offset,isquery = T
  )
  #browser()
  stk.z <- inla.stack(stk.z.obs,stk.z.query)
  
  ## reset the inla compute argument for prediction
  compute_ <- fit$.args$control.compute
  compute_$config <- TRUE ## allow simulation
  
  formula_ <- reformulate(as.character(fit$.args$formula)[3],
                         response = getVars(fit$.args$formula)[1])
  ## prepare data structure and arguments as a list
  internal__ <- list(
    formula = formula_,
    data=inla_spde_adjustZ_helper(stk.z,sout),
    control.compute = compute_,
    control.predictor=list(A=inla.stack.A(stk.z),compute=TRUE))
  
  ## return
  internal__
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
    ## prepare model formula and training data
    if(!is.null(fit$.args$spde)){
      ## SPDE model
      internal__ <- inla_spde_stack(fit,data)
      N__ <- length(internal__$data$y)
      #browser()
    }else{
      ## regular INLA model
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
      tmp_ <- sgam.translate(formula = for_,
                                   data = data_,
                                   param = fit$.args$param)
      
      ## refit the model for prediction
      compute_ <- fit$.args$control.compute
      compute_$config <- TRUE ## allow simulation
      predictor_ <- fit$.args$control.predictor
      predictor_$link <- NULL
      
      N__ <- length(data_[[1]])
      internal__ <- list(formula = tmp_$value$formula,
                         data=tmp_$data,
                         control.compute = compute_, 
                         control.predictor = predictor_)
    }

    #browser()
    ## set the E for Poisson
    E_ <- fit$.args$E
    if(!is.null(fit$.args$E)){
      #np <- length(data_[[1]]) - length(E_)
      np <- N__ - length(E_)
      E_ <- c(fit$.args$E,rep(1,np))
    }
    
    ## set the Ntrials for Binomial
    Ntrials_ <- fit$.args$Ntrials
    if(!is.null(fit$.args$Ntrials)){
      #np <- length(data_[[1]]) - length(Ntrials_)
      np <- N__ - length(Ntrials_)
      Ntrials_ <- c(fit$.args$Ntrials, rep(1,np))
    }
    #browser()
    args_ <- c(internal__,list(
      family = fit$.args$family,
      E=E_,Ntrials=Ntrials_,
      control.family = fit$.args$control.family,
      control.inla = fit$.args$control.inla,
      control.fixed = fit$.args$control.fixed,
      control.mode = control.mode(theta=fit$mode$theta,restart=FALSE),
#      fit$.args$control.mode,
      control.expert = fit$.args$control.expert,
      control.lincomb = fit$.args$control.lincomb,
      control.update = fit$.args$control.update,
      control.lp.scale = fit$.args$control.lp.scale,
      control.pardiso = fit$.args$control.pardiso))
    
    if(!is.null(E_)){
      args_$E <- E_
    }
    if(!is.null(Ntrials_)){
      args_$Ntrials <- Ntrials_
    }
    #browser()
    
    spde <- fit$.args$spde
    #args_$data$spde <- fit$.args$spde
    fit_new <- do.call("inla",args = args_)
    samples <- inla.posterior.sample(n=n,result = fit_new)
  }
  ## precision location
  prec.ii <- grep("Precision",names(samples[[1]]$hyperpar))
  name.ii <- gsub("Precision","SD",names(samples[[1]]$hyperpar))[prec.ii]
  
  ## Predictor locations
  if(!is.null(fit$.args$spde)){
    ## spde model
    pred.ii <- grep("APredictor",row.names(samples[[1]]$latent))
    if(!missing(data)){
      ## new data, remove training data
      ntrain <- grep("APredictor",row.names(fit$summary.fitted.values))
      #ntrain <- nrow(fit$summary.fitted.values)
      pred.ii <- pred.ii[-c(1:length(ntrain))]
    }
    
  }else{
    pred.ii <- grep("Predictor",row.names(samples[[1]]$latent))
    if(!missing(data)){
      ## new data, remove training data
      ntrain <- nrow(fit$summary.fitted.values)
      pred.ii <- pred.ii[-c(1:ntrain)]
    }
    
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
  if(missing(data)){
    ## for in-sample prediction, use the training data
    data <- fit$.args$raw
  }
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
  source("Z_link_inla_dev_5.R")
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
  source("Z_link_inla_dev_5.R") 
  
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
  source("Z_link_inla_dev_5.R") 
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
  source("Z_link_inla_dev_5.R") 
  
  
  mu <- link(fit,data=dat[1:10,],E=1:10)
  post <- link(fit,data=dat[1:10,],type = "prediction",E=1:10)
  plot(apply(mu,2,mean)/1:10,fit$summary.fitted.values$mean[1:10])
  abline(0,1)
  plot(apply(mu,2,mean),apply(post,2,mean))
  abline(0,1)
  plot(apply(mu,2,sd),apply(post,2,sd))
  table(post[,1])  
}

dev_SPDE <- function(){
  rm(list=ls())
  source("sgam_2.R")
  source("nmiss_2.R")
  source("getVars.R")
  source("s2inla_2.R")
  load(file="Scratch/spde_glm_dev1.rdata")
  source("spde_gam_1.R")
  
  wds.geo$strata <- gl(3,510)[1:nrow(wds.geo)]
  dat <- wds.geo
  source("Z_link_inla_dev_5b.R")
  wds.geo$X2 <- scale(wds.geo$X)
  wds.geo$Y2 <- scale(wds.geo$Y)
  mod <- spde.gam(TOT~strata+ns(X2,df=3)+ns(Y2,df=3),
                  data=wds.geo,mesh=mesh1,spde=spde1,
                  xcoord="X",ycoord="Y",family="poisson",
                  E=AREA)
  mod <- spde.gam(TOT~strata+s(X2),
                  data=wds.geo,mesh=mesh1,spde=spde1,
                  xcoord="X",ycoord="Y",family="poisson",
                  E=AREA)
  mu <- link(mod,n=99,data=wds.geo[1:10,],E=wds.geo$AREA[1:10])
  post <- link(mod,n=99,data=wds.geo[1:10,],
               type = "prediction",
               E=wds.geo$AREA[1:10])
  plot(apply(mu,2,mean)/wds.geo$AREA[1:10],
       mod$summary.fitted.values$mean[1:ncol(mu)])
  abline(0,1)
  plot(apply(mu,2,mean),apply(post,2,mean))
  abline(0,1)
  plot(apply(mu,2,sd),apply(post,2,sd))
  abline(0,1)
  #table(post[,1])  
  
}
