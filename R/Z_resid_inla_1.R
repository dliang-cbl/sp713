setGeneric("resid",function(fit,type="working",...){
  message( paste0("No resid method for object of class '",
                  class(fit),"'. Returning instead.") )
  return(0)
})
setMethod("resid","inla",function(fit,type="working",...){
  ## ... additional arguments to residuals.glm.
  ## extract the fixed effects
  #browser()
  model.matrix <- model.matrix(fit$.args$fixed,
                               data=fit$.args$raw)
  fixed.eta <- as.vector(model.matrix %*% fit$summary.fixed$mean)
  
  ## Adjust E
  if(!is.null(fit$.args$E)){
    fixed.eta <- fixed.eta + log(fit$.args$E)
  }
  
  ## refit a model conditioning only on the fixed effects
  ## without the spatial term
  for_ <- reformulate("0+offset(E)",
                      response = getVars(fit$.args$formula0)[1])
  ## binomial model testing
  glmfit <- do.call("glm",list(formula=for_,
                               data=cbind(fit$.args$raw,E=fixed.eta),
                               family=fit$.args$family))
  
  ## get the working residuals
  residuals(glmfit,type=type)
  
})

dev_ <- function(){
  rm(list=ls())
  source("sgam_2.R")
  source("nmiss_2.R")
  source("getVars.R")
  source("spde_glm_1.R")
  load(file="Scratch/spde_glm_dev1.rdata")
  library(INLA)
  wds.geo$err <- seq(1,nrow(wds.geo))
  wds.geo$X2 <- as.vector(scale(wds.geo$X))
  wds.geo$Y2 <- as.vector(scale(wds.geo$Y))
  
  mod <- spde.glm(TOT~X2+I(X2^2)+f(err),data=wds.geo,mesh=mesh1,spde=spde1,
                  xcoord="X",ycoord="Y",family="poisson",
                  E=AREA)
  source("Z_resid_inla_1.R")
  e <- resid(mod)
}