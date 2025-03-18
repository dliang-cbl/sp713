#' spatial generalized linear model 
#' 
#' spatial GLM with Stochastic Partial Differential Equations using R and Integrated Nested Laplace Approximation. 
#'
#' @param formula fixed effects type formula object
#' @param data data frame to evaluate the model
#' @param mesh mesh for building the Stochastic Partial Differential Equations
#' @param spde the Stochastic Partial Differential Equations model
#' @param random formula defining the random term
#' @param family family for response variable
#' @param xcoord field name of x coordinate
#' @param ycoord field name of y coordinate
#' @param E offset for Poisson term
#' @param Ntrials number of trials for binomial
#' @param query prediction data
#' @param verbose detailed
#' @param control a list of control objects
#' 
#' @return fitted INLA object
#' 
#' @export
sglm.inla <-function(
  formula,data,mesh,spde,random=~s,family="gaussian",
  xcoord='x',ycoord='y',
  E=NULL,Ntrials=NULL,
  query=NULL,verbose=F,control=list(nsamples=0))
{
  obs <- data
  ## process formula
  l.vars <- getVars(formula)
  zcoord <- l.vars[1]
  if(nchar(zcoord)<1){
    stop("please specify a response in formula.\n")
  }
  fixed <- reformulate(as.character(formula)[3],response=NULL)

  ## fixed effects
  z.obs.lst <- list(b0=rep(1,nrow(obs)))
  if(length(all.vars(fixed))>0){
    z.obs.lst <- c(z.obs.lst,as.list(obs[,all.vars(fixed),drop=FALSE]))
  }
  if(!is.null(query)){
    z.query.lst <- list(b0=rep(1,nrow(query)))
    if(length(all.vars(fixed))>0){
      z.query.lst <- c(z.query.lst,as.list(query[,all.vars(fixed),drop=FALSE]))
    }
  }
  
  ## random effects
  l.terms <- lme.inla.random(formula=random)
  
  coords <- as.matrix(obs[,c(xcoord,ycoord)])
  ## Define INLA data structure
  A <- inla.spde.make.A(mesh, loc=coords)
  stk.z.obs <- inla.stack(tag="obs",
                      data=list(y=obs[,zcoord]),
                      A=list(A,1),
                      effects=list(
                        list(field=1:spde$n.spde),
                        z.obs.lst
                      )
  )
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
  if(!is.null(query)){
    stk.z <- inla.stack(stk.z.obs,stk.z.query)
  }
  else{
    stk.z <- stk.z.obs
  }
  
  if(!is.null(query) & !is.null(E)){
    E <- c(E,rep(1,nrow(query)))
  }
  
  ## Define INLA model
  fixed.vars <- c("0+b0",as.character(fixed)[2])
  spatial.vars <- NULL
  if(l.terms$is.spde){
    spatial.vars <- c(spatial.vars,"f(field,model=spde)")
  }
  formula <- reformulate(c(fixed.vars,spatial.vars),response="y")
  
  if(family=="gaussian"){
    mod1 <- try(inla(formula,family=family,
                     data=inla.stack.data(stk.z),
                     control.compute=list(config=TRUE,dic=TRUE,return.marginals.predictor=FALSE),
                     control.predictor=list(A=inla.stack.A(stk.z),compute=TRUE),
                     control.family=list(hyper=list(prec=list(prior="logtnormal",param=c(0,0.0001)))),
                     control.fixed=list(expand.factor.strategy="inla"),
                     verbose=verbose
    ))
  }
  else if(family=="poisson"){
    mod1 <- try(inla(formula,family=family,E=E,
                     data=inla.stack.data(stk.z),
                     control.compute=list(config=TRUE,dic=TRUE,return.marginals.predictor=FALSE),
                     control.predictor=list(A=inla.stack.A(stk.z),compute=TRUE,link=1,quantiles=NULL),
                     #control.family=list(hyper=list(prec=list(prior="logtnormal",param=c(0,0.0001)))),
                     control.fixed=list(expand.factor.strategy="inla"),
                     verbose=verbose
    ))
  }
  else if(family=="nbinomial"){
    mod1 <- try(inla(formula,family=family,E=E,
                     data=inla.stack.data(stk.z),
                     control.compute=list(config=TRUE,dic=TRUE,return.marginals.predictor=FALSE),
                     control.predictor=list(A=inla.stack.A(stk.z),compute=TRUE,
                                            link=1,quantiles=NULL),
                     #control.family=list(hyper=list(prec=list(prior="logtnormal",param=c(0,0.0001)))),
                     control.fixed=list(expand.factor.strategy="inla"),
                     verbose=verbose
    ))
  }
  else if(family=="binomial"){
    mod1 <- try(inla(formula,family=family,Ntrials = Ntrials,
                     data=inla.stack.data(stk.z),
                     control.compute = list(config=TRUE,dic=TRUE,return.marginals.predictor=FALSE),
                     control.predictor = list(A=inla.stack.A(stk.z),compute=TRUE,
                                              link=1,quantiles=NULL),
                     control.fixed=list(expand.factor.strategy="inla"),
                     verbose = verbose))
  }
  else{
    stop("unimplemented family ",family,"\n")
  }

  
  if(class(mod1)=="try-error"){
    cat("INLA model error, please re-run with verbose=T.\n")
    return(list(est=NULL,inla=NULL,samples=NULL,
                terms=l.terms,predSample=NULL,DIC=NULL))
  }
  
  ## parameter estimation
  if(l.terms$is.spde){
    r.f <- try(inla.spde2.result(mod1, "field", spde, do.transform=TRUE),silent=TRUE)
    ## partial sill
    sigmasq <- inla.qmarginal(c(0.5,0.025,0.975),r.f$marginals.variance.nominal[[1]])
    ## range
    h <- inla.qmarginal(c(0.5,0.025,0.975),r.f$marginals.range.nominal[[1]])
  }
  else{
    sigmasq <- h <- rep(NA,3)
  }
  
  tausq <- rep(NA,3)
  if(family=="gaussian"){
    tausq <- mod1$summary.hyperpar[1,c(4,3,5)]
  }
  else if(family=="nbinomial"){
    tausq <- mod1$summary.hyperpar[1,c(4,3,5)]
  }
  
  
  alpha <- as.matrix(mod1$summary.fixed[,c(4,3,5)])
  parm <- as.matrix(rbind(alpha,sigmasq,h,tausq))
  colnames(parm) <- c("fit","lwr","upr")
  rownames(parm) <- c(row.names(alpha),"sigmasq","range","tausq")
  
  
  predSample <- NULL
  
  if(!is.null(query) & control$nsamples>0){
    cat("## approximate posterior samples ## \n")
    ## posterior samples
    samples <- inla.posterior.sample(control$nsamples,mod1)
  }
  else{
    samples <- NULL
  }
  
  predSample <- NULL
  if(!is.null(query) & !is.null(samples)){ ## generate predictive samples
    predSample <- sapply(samples,lme.inla.pred,
                           zindex=inla.stack.index(stk.z,"query")$data)
  }
  
  ## spatial random effects
  random.1 <- NULL
  if(l.terms$is.spde){
    proj <- inla.mesh.projector(mesh,xlim=range(obs[,xcoord,]),
                                ylim=range(obs[,ycoord]))
    coords. <- as.data.frame(proj$lattice$loc)
    names(coords.) <- c(xcoord,ycoord)
    random. <- data.frame(
      mean=c(inla.mesh.project(proj,mod1$summary.random$field$mean)),
      sd=c(inla.mesh.project(proj,mod1$summary.random$field$sd))
    )
    random.0 <- cbind(coords.,random.)
    
    ## only visualize data 
    ## find nearest neighbor to each grid
    nn <- get.knnx(
      query=random.0[,c(xcoord,ycoord)],
      data=obs[,c(xcoord,ycoord)],
      k=1
    )
    
    random.1 <- subset(random.0,nn$nn.dist[,1]<0.5)
  }
  
  ## prediction
  predict.1 <- NULL
  if(!is.null(query)){
    ii <- inla.stack.index(stk.z,"query")$data
    predict.0 <- data.frame(
      mean=mod1$summary.fitted.values$mean[ii],
      sd=mod1$summary.fitted.values$sd[ii]
    )
    predict.1 <- cbind(query[,c(xcoord,ycoord)],predict.0)
  }
  
  l.dic <- with(mod1$dic,c(Dbar=mean.deviance,pD=p.eff,DIC=dic))
  #browser()
  #l.spde <- paste(deparse(substitute(spde)),"(",deparse(substitute(mesh)),")",sep="")
  l.spde <- deparse(substitute(spde))
  l.formula <- reformulate(c(as.character(fixed)[2],l.spde),response=zcoord)
  list(est=parm,summary=parm,random.=random.1,predict.=predict.1,
       inla=mod1, terms=l.terms,predSample=predSample,DIC=l.dic,
       formula=l.formula)
}

lme.inla.pred <- function(l,zindex){
  ## helper function to derive posterior predictive distribution
  ## l : list of object from INLA approximate posterior sample
  ## zindex : index of catch process
  ## yindex : index of number of catch given catch 
  ## value: a posterior sample draw from m prediction grid
  
  names <- rownames(l$latent)
  
  eta.z <- l$latent[zindex,1]

  #browser()
  zp <- rnorm(length(zindex),eta.z,1/sqrt(l$hyper[1]))
  ## posterior predictive distribution
  zp
}
lme.inla.random <- function(formula){
  l.vars <- all.vars(formula)
  ret <- list(is.spde =FALSE, is.overdisp =FALSE)
  if(length(l.vars)==0){
    return(ret)
  }
  if(any(!(l.vars %in% c("s","e")))){
    stop("random terms must be either s(SPDE) or e(overdisp)\n")
  }
  if(length(l.vars)==1){
    if(l.vars=="s"){
      ret$is.spde <- TRUE
    }
    else{
      ret$is.overdisp <- TRUE
    }
    return(ret)
  }
  if(length(l.vars)==2){
    ret$is.spde <- TRUE
    ret$is.overdisp <- TRUE
  }
  ret
}

## minor edit to make same name object work

