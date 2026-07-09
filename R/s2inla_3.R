
## ver. 7: allow formula and data to be used directly for INLA
## bug: no smooth interface.
## helper function to define smooth term fitting in INLA
## supporting offset through formula instead

## helper function to fix smooth terms in INLA
## formula: the formula containing only smooth terms
## data: the data to define smooth terms
## U  : standard deviations of T Normal priors on standard deviation of random effects
##    : a vector of two, the first entry is for the smooth with rank > 1
##    :                  the second entry is for the smooth with rank 1
##    : NOTE: these numbers appear work with single smooth, but the linear
##    combination part is tricky with z models
## family : glm family to process
## zmat   : character string for zmatrix list name
## verbose : whether to output details
## value: X: the Z matrices for INLA definition
##        formula: the INLA formula including smooth terms
##        data: the expanded data frame for INLA model
##        sterms: smooth terms
##        zmat  : input for zmat name
s2inla <- function(
  formula,data,U=c(1,1),
  family=gaussian,zmat="Z",verbose=FALSE)
{
  if(FALSE){
    rm(list=ls())
    library(mgcv)
    source("../R/getVars.R")
    source("../R/nmiss_2.R")
    source("../R/proc.formula.R")
    source("../R/lme_inla_formula.R")
    d <- gamSim(n=100)

    formula = y~s(x0)+s(x1)+x2
    data=d
    U=c(1,1)
    family=gaussian;zmat="Z";verbose=FALSE
  }
  ##browser()
  l.tt <- terms(formula)
  l.ss <- grep("^s\\(.+\\)",attr(l.tt,"term.labels"))
  l.ff <- grep("^f\\(.+\\)",attr(l.tt,"term.labels"))
  l.oo <- row.names(attr(l.tt,"factors"))[attr(l.tt,"offset")] ## offset terms
  l.terms.nos <- attr(l.tt,"term.labels")[-l.ss] ## identify non-smooth terms
  l.terms.s <- attr(l.tt,"term.labels")[l.ss] ## identify only smooth terms
  l.response <- getVars(formula)[1]              ## identify response variable
  l.intercept <- grep("^intercept$",attr(l.tt,"term.labels")) ## identify intercept term
  
  ## remove the f portion of the formula to be compatible with jagam
  ## also remove the intercept portion of the formula
  l.formula <- formula
  if(length(c(l.intercept,l.ff))>0){
    l.terms.nof <- c(
      getVars(formula)[-1][-c(l.intercept,l.ff,l.ss)],
      attr(l.tt,"term.labels")[l.ss]) ## extract all no f terms
    if(length(l.terms.nof)>0){
      l.formula <- reformulate(l.terms.nof,response=l.response)  
    }else{
      l.formula <- as.formula(paste(l.response,"~1",sep=""))
    }
  }
  
  if(length(l.ss)==0){
    ## the formula contain no s term
    catch <- try(model.matrix(l.formula,data),silent = TRUE)
    frame <- model.frame(l.formula,data)
    offset <- model.offset(frame)
    if(is.null(offset)){
      offset <- rep(0,length(model.response(frame)))
    }
    return(list(y=model.response(frame),offset=offset,
                Design=catch,data=data,formula=formula,
                sterms=NULL,zmat=zmat))
  }
  
  #browser()
  
  ## remove missing values
  ## allow NA in responses for predictions
  s.formula <- reformulate(l.terms.s,response=l.response)
  data_ <- nmiss(l.formula,data,forceResponse = T) 

  ## derive offset
  if(length(l.oo)==0){
    offset <- rep(0,length(data_[[l.response]]))
  }else{
    l.oo.formula <- reformulate(l.oo,response = l.response)
    l.oo.data <- data[all.vars(l.oo.formula)]
    offset <- model.offset(model.frame(l.oo.formula,data=l.oo.data))
  }
  
  jf <- tempfile()
  jd <- jagam(s.formula,data=data_,file=jf,
              diagonalize = TRUE,na.action=na.pass,family=family)
  
  ## creat fixed design matrix
  nms <- names(jd$pregam$cmX) ## if this name is not null, it is a design matrix
  Design <- jd$pregam$X[,nchar(nms)>0,drop=F] ## now without f terms
  #browser()
  
  ## create list for random design matrix, formula and simply IDs
  X <- vector("list",length(jd$pregam$smooth))
  fout <- vector("list",length(X))
  #varnames <- letters[1:length(X)]
  varnames <- sapply(jd$pregam$smooth,function(elmt) elmt$term)
  ones <- vector("list",length(X))
  
  LC <- vector("list",length(X))
  LC_name <- vector("list",length(X))
  
  for(i in 1:length(jd$pregam$smooth)){ ## loop thru smooth term
    ## create id matrices
    tmp <- matrix(1:nrow(Design),nrow(Design),length(jd$pregam$smooth[[i]]$rank))
    colnames(tmp) <- paste(varnames[i],1:length(jd$pregam$smooth[[i]]$rank),sep="")
    colnames(tmp)[1] <- varnames[i]
    #colnames(tmp) <- c(varnames[i],paste(varnames[i],2:length(jd$pregam$smooth[[i]]$rank),sep=""))
    ones[[i]] <- as.data.frame(tmp)
    
    ## create linear combination list
    tmp_ <- list(a__=diag(nrow(ones[[i]])),
                 b__=diag(nrow(ones[[i]])))
    names(tmp_) <- names(ones[[i]])
    #tmp2 <- inla.make.lincombs(x0=diag(100),x02=diag(100))
    LC[[i]] <- do.call(inla.make.lincombs,tmp_) 
    names(LC[[i]]) <- paste0(varnames[i],"_",1:length(LC[[i]]))
    LC_name[[i]] <- rep(varnames[i],nrow(ones[[i]]))
    
    ## create design matrices
    X[[i]] <- vector("list",length(jd$pregam$smooth[[i]]$rank))

    start <- jd$pregam$smooth[[i]]$first.para
    ff <- rep("",length(jd$pregam$smooth[[i]]$rank))
    for(j in 1:length(jd$pregam$smooth[[i]]$rank)){
      idx <- start+seq(0,jd$pregam$smooth[[i]]$rank[j]-1)
      if(verbose){
        cat(i,j,"idx=",idx,"\n")
      }
      if(j==1){ ## smooth rank 
        l.U <- U[1]
      }                  ## end
      else{              ## smooth rank == 1
        l.U <- U[2]
      }                  ## end
      X[[i]][[j]] <- jd$jags.data$X[,idx,drop=FALSE]
      ## INLA formula
      tmpfor <- paste("f(",colnames(tmp)[j],',model="z",',"Z=",zmat,i,"_",j,",",
                      'hyper=list(prec=list(prior="logtnormal",param=c(0,',
                      l.U,'))))',
                      sep="")
      if(j==1){
        ## range space
        ff[j] <- tmpfor
      }
      else{
        ## null space
        ff[j] <- tmpfor
      }
      start <- max(idx)+1
    }
    fout[[i]] <- paste(ff[nchar(ff)>0],collapse = "+")
  }
  
  #browser()
  ## process data for INLA
  #r.data <- cbind(data,do.call(cbind,ones))
  if(class(data)=="data.frame"){
    r.data <- cbind(do.call(cbind,ones),data)
    names(r.data) <- make.names(names(r.data),T)
  }
  if(class(data)=="list"){
    r.data.0 <- do.call(cbind,ones)
    r.data <- c(as.list(r.data.0),data)
    names(r.data) <- make.names(names(r.data),T)
  }
  ## process formula for INLA
  
  ## append intercepts
  if(attr(l.tt,"intercept")==0){
    r.formula <- reformulate(
      c("0",do.call(c,fout),l.terms.nos,l.oo),
      response=l.response)
    
  }else{
    r.formula <- reformulate(
      c(do.call(c,fout),l.terms.nos,l.oo),
      response=l.response)
  }
  #l.ii <- ifelse(attr(l.tt,"intercept")==0,"0",NULL)
  #r.formula <- reformulate(
  #  c(l.ii,do.call(c,fout),l.terms.nos,l.oo),
  #  response=l.response)
  
  ## return both range and null components to be consistent with jagam
  for(i in 1:length(X)){
    names(X[[i]]) <- paste0(zmat,i,"_",seq(1,length(X[[i]])))
  }
  Xout <- do.call(c,X)
  #X <- lapply(X,function(elmt) elmt[[1]])
  #names(X) <- paste0(zmat,seq(1,length(X)))
  
  LCout <- do.call(c,LC)
  LCnameout <- do.call(c,LC_name)
  
  l.r <- list(y=jd$pregam$y,Design=Design,X=Xout,
              LC=LCout,LCname=LCnameout,
              formula=r.formula,data=r.data,
              raw=data,
              pregam=jd$pregam,offset=offset,
              sterms=varnames,zmat=zmat)
  names(l.r)[3] <- zmat
  
  l.r
}

