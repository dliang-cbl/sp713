setGeneric("extract_log_lik",
           function( object , n=1000 , refresh=0.1 , pointwise=FALSE , ... ) {
             message( paste0("No extract loglik method for object of class '",
                             class(object),"'. Returning instead.") )
             return(0)
           }
)

setMethod("extract_log_lik", "inla",
          function(object , n=1000 , refresh=0 , pointwise=FALSE , log_lik="log_lik" , ...){
            #browser()
            object <- inla_reset_config(object)
            
            ## extract response
            vars_ <- all.vars(object$.args$formula)
            if(length(vars_)==1){
              ## only response
              formula_ <- as.formula(paste0(vars_,"~1"))
            }else{
              formula_ <- reformulate(vars_[-1],response = vars_[1])
            }
            frame_ <- model.frame(formula_,object$.args$data)
            y_ <- model.response(frame_)
            

            ## posterior simulation
            sim <- inla.posterior.sample(n=n,result=object)
            
            ## extract log likelihood predictors
            lp_matrix <- t(sapply(sim,function(elmt){elmt$latent[1:length(y_)]}))
            
            ## define link function
            linkinv_ <- switch(object$.args$family,
                               poisson=poisson()$linkinv,
                               gaussian=gaussian()$linkinv,
                               binomial=binomial()$linkinv,
                               nbinomial=poisson()$linkinv,
                               gaussian()$linkinv)
            
            fam_ <- object$.args$family
            
            if(fam_=="binomial"){
              ## extract N trials
              if(is.null(object$.args$Ntrials)){
                Ntrials <- 1
              }else{
                Ntrials <- object$.args$Ntrials
              }
              #browser()
              m_matrix <- linkinv_(lp_matrix)
              
              ll_matrix <- t(apply(m_matrix,1,dbinom,x=y_,size=Ntrials,log=TRUE))
              return(list(log_lik=ll_matrix))
            }
            
            if(fam_=="poisson"){
              #browser()
              ## extract expected
              if(is.null(object$.args$E)){
                E_ <- rep(1,length(y_))
              }else{
                E_ <- object$.args$E
              }
              
              m_matrix <- linkinv_(lp_matrix)
              m_matrix_2 <- sweep(m_matrix,2,STATS = E_,FUN = "*")
              ll_matrix <- t(apply(m_matrix_2,1,dpois,x=y_,log=TRUE))
              return(list(log_lik=ll_matrix))
              
            }
            
            if(fam_=="nbinomial"){
              #browser()
              ## extract expected
              if(is.null(object$.args$E)){
                E_ <- rep(1,length(y_))
              }else{
                E_ <- object$.args$E
              }
              m_matrix <- linkinv_(lp_matrix)
              m_matrix_2 <- sweep(m_matrix,2,STATS = E_,FUN = "*")
              
              #inla.likelihood.parser
              hyper_ <- sapply(sim,function(elmt) 
                elmt$hyper["size for the nbinomial observations (1/overdispersion)"])
              llik_ <- t(sapply(1:nrow(m_matrix_2),function(i){
                size_ <- (hyper_[i])
                p_ <- size_/(size_+m_matrix_2[i,])
                dnbinom(y_,size=size_,prob = p_,log = TRUE)
              }))
              
              return(list(log_lik=llik_))
              
            }
            extractll <- function(elmt){
              mu <- elmt$latent[1:length(y_)]
              sigma <- 1/sqrt(elmt$hyperpar["Precision for the Gaussian observations"])
              dnorm(y_,mean=mu,sd=sigma,log=TRUE)
            }
            extractll_t <- function(elmt){
              mu <- elmt$latent[1:length(y_)]
              sigma <- 1/sqrt(elmt$hyperpar["precision for the student-t observations"])
              dt(y_/sigma,df=3,log = TRUE)
            }

            if(fam_=="gaussian"){
              ll_matrix <- t(sapply(sim,extractll))
            }else if(fam_=="t"){
              ll_matrix <- t(sapply(sim,extractll_t))
            }
            list(log_lik=ll_matrix)
          }
)

setMethod("nobs","list",function(object,...){
  return(ncol(object[["log_lik"]]))
})

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
  m1 <- inla(daysabs ~ math + prog, data = dat,family = "nbinomial",
             control.compute = ctrl_)
  m1$waic$waic
  source("extract_ll_methods_2.R")
  l1 <- extract_log_lik(m1)
  WAIC(l1)
}

dev_Poisson <- function(){
  rm(list=ls())
  library(rethinking)
  library(INLA)
  source("extract_ll_methods_2.R")
  counts <- c(18,17,15,20,10,20,25,13,12)
  outcome <- gl(3,1,9)
  treatment <- gl(3,3)
  dat <- data.frame(treatment, outcome, counts) # showing data
  ctrl_ <-  list(dic = TRUE, waic = TRUE, config = TRUE)
  glm.D93 <- inla(counts ~ outcome + treatment, data=dat,family = "poisson",
                  control.compute = ctrl_,E=1:9)
  ll.D93 <- extract_log_lik(glm.D93)
  WAIC(ll.D93)
  glm.D93$waic$waic
}