setGeneric("Variog",function(fit,plot=TRUE,...){
  message( paste0("No Variog method for object of class '",
                  class(fit),"'. Returning instead.") )
  return(0)
})
setMethod("Variog","inla",function(fit,plot=TRUE,...){
  ## return: estimated variogram parameters
  ## side effect: a variogram of the working residuals
  ##  without the spatial effects.
  ## plus the fitted variogram overlaid.
  if(!("spde" %in% names(fit$.args))){
    message( paste0("No spde object for object of class '",
                    class(fit),"'. Returning instead.") )
    return(0)
  }
  ## nugget parameter (if any)
  sigmasq <- rep(NA,3)
  if(fit$.args$family=="gaussian"){
    sigmasq <- 1/fit$summary.hyperpar[1,c(4,3,5)]
  }
  else if(fit$.args$family=="poisson"){
    if(length(fit$.args$err)>0){
      i.err <- grep(fit$.args$err,row.names(fit$summary.hyperpar))
      sigmasq <- 1/fit$summary.hyperpar[i.err,c(4,3,5)]
    }
  }
  else if(fit$.args$family=="nbinomial"){
    ## need checking
    sigmasq <- 1/fit$summary.hyperpar[1,c(4,3,5)]
  }
  
  #browser()
  ## partial sill and range parameter
  tausq <- theta <- rep(NA,3)
  r.f <- try(inla.spde2.result(fit, "field", fit$.args$spde, 
                               do.transform=TRUE),silent=TRUE)
  
  if(class(r.f)!="try-error"){
    ## partial sill
    tausq <- inla.qmarginal(c(0.5,0.025,0.975),
                              r.f$marginals.variance.nominal[[1]])
    ## range
    theta <- inla.qmarginal(c(0.5,0.025,0.975),
                             r.f$marginals.range.nominal[[1]])
  }else{
    stop("fail to estimate variogram.\n")
  }
  
  parm <- as.matrix(rbind(tausq,theta,sigmasq))
  colnames(parm) <- c("fit","lwr","upr")
  rownames(parm) <- c("psill","range","nugget")
  
  ## get the working residuals
  resid <- resid(fit,type="working")
  
  ## sample variogram of the residuals
  data_ <- data.frame(cbind(fit$.args$coords,resid))
  geodata_ <- as.geodata(data_)
  vg_ <- variog(geodata_,...)
  
  ## estimated variogram from INLA
  fit_ <- variofit(vg_,ini.cov.pars = parm[1:2,"fit"] ,
                   cov.model = "matern",
                   nugget = parm[3,"fit"],
                   fix.nugget = TRUE,kappa = 0.5,
                   messages = F)
  fit_$cov.pars <- c(parm[1,"fit"],parm[2,"fit"]/3)
  
  #browser()
  if(plot){
    ## confidence limits
    lwr_ <- fit_
    upr_ <- fit_
    lwr_$cov.pars <- c(parm[1,"lwr"],parm[2,"lwr"]/3)
    lwr_$nugget <- parm[3,"upr"]
    upr_$cov.pars <- c(parm[1,"upr"],parm[2,"upr"]/3)
    upr_$nugget <- parm[3,"lwr"]
    #plot(vg_,ylim=c(0,1.1*sum(parm[c(1,3),3])))
    plot(vg_)
    lines(fit_)
    lines(lwr_,lty=2)
    lines(upr_,lty=2)
  }
  return(fit_)
})

test_variog <- function(){
  rm(list=ls())
  source("sgam_2.R")
  source("nmiss_2.R")
  source("getVars.R")
  source("spde_glm_1.R")
  load(file="Scratch/spde_glm_dev1.rdata")
  library(INLA)
  library(geoR)
  wds.geo$err <- seq(1,nrow(wds.geo))
  wds.geo$X2 <- as.vector(scale(wds.geo$X))
  wds.geo$Y2 <- as.vector(scale(wds.geo$Y))
  
  mod <- spde.glm(TOT~X2+I(X2^2)+f(err),data=wds.geo,mesh=mesh1,spde=spde1,
                  xcoord="X",ycoord="Y",family="poisson",
                  E=AREA)
  source("Z_resid_inla_1.R")
  source("Z_variog_inla_2.R")
  e <- resid(mod,type="working")
  summary(mod)
  Variog(mod)  
}