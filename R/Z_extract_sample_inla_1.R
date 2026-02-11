setGeneric("extract.samples",function(object,n=10000,pars,...){
  message( paste0("No extract.samples method for object of class '",
                  class(object),"'. Returning instead.") )
  return(0)
})
setMethod("extract.samples","inla",function(object,n=10000,pars,...){
  ## draw posteriors
  samples <- inla.posterior.sample(n=n,result = object)
  ## precision location
  prec.ii <- grep("Precision",names(samples[[1]]$hyperpar))
  name.ii <- gsub("Precision","SD",names(samples[[1]]$hyperpar))[prec.ii]
  ## Predictor locations
  pred.ii <- grep("Predictor",row.names(samples[[1]]$latent))
  #browser()
  ## fixed effects
  names.fixed <- row.names(object$summary.fixed)
  extract <- function(elmt){
    sigma <- sqrt(elmt$hyperpar[prec.ii])
    names(sigma) <- name.ii
    fixed <- elmt$latent[-pred.ii,1]
    c(fixed,sigma)
  }
  post <- t(sapply(samples,extract))
  colnames(post)[1:length(names.fixed)] <- names.fixed
  as.data.frame(post)
})

setMethod("extract.samples","list",function(object,n=10000,pars,...){
  ## draw posteriors
  samples <- object
  ## precision location
  prec.ii <- grep("Precision",names(samples[[1]]$hyperpar))
  name.ii <- gsub("Precision","SD",names(samples[[1]]$hyperpar))[prec.ii]
  ## Predictor locations
  pred.ii <- grep("Predictor",row.names(samples[[1]]$latent))
  ## fixed effects names
  names.fixed <- row.names(object[[1]]$latent)[-pred.ii]
  loc <- gregexpr("\\:",names.fixed)
  out_ <- rep("",length(names.fixed))
  for(i in 1:length(loc)){
    out_[i] <- substr(names.fixed[i],1,loc[[i]]-1)
  }
  extract <- function(elmt){
    sigma <- sqrt(elmt$hyperpar[prec.ii])
    names(sigma) <- name.ii
    fixed <- elmt$latent[-pred.ii,1]
    c(fixed,sigma)
  }
  post <- t(sapply(samples,extract))
  colnames(post)[1:length(names.fixed)] <- out_
  as.data.frame(post)
  
})