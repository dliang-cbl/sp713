## helper script to subset to data where the formula is not missing.
nmiss <- function(formula,data,forceResponse=F)
{
  ## forceResponse: whether to keep response regardless of missing
  ## value: subset of data with missing values in formula removed
  
  ## response variable
  pred_var <- labels(terms(formula))
  all_var <- all.vars(formula)
  if(length(pred_var)==length(all_var)){
    resp_var <- NULL
  }else{
    resp_var <- all_var[1] ## univariate assumed
  }
  
  ## variables to look through missing values
  if(forceResponse){
    vars__ <- all_var[-1]
  }else{
    vars__ <- all_var
  }
  
  if(is.null(data)){
    return(data)
  }
  
  if(class(data)=="list"){
    na.action <- rep(F,length(data[[resp_var]]))
    for(var in vars__){
      tmp <- na.omit(data[[var]])
      if(!is.null(attr(tmp,"na.action"))){
        na.action[attr(tmp,"na.action")] <- T
      }
    }
    value <- lapply(data,function(elmt){
      if(is.null(dim(elmt))){
        return(elmt[!na.action])
      }else{
        ## only matrices are allowed
        return(elmt[!na.action,])
      }
    })
    return(value)
  }
  if(class(data)=="data.frame"){
    tmp <- na.omit(data[,vars__])
    #browser()
    if(!is.null(attr(tmp,"na.action"))){
      data <- data[-attr(tmp,"na.action"),]
    }
    return(data)
  }
  return(data)
}

test_ <- function(){
  rm(list=ls())
  source("nmiss_2.R")
  dat <- data.frame(y=rnorm(100),x=rnorm(100),z=rnorm(100))
  dat$x[1] <- NA
  dat$y[2] <- NA
  r <- nmiss(y~x,dat)
  r <- nmiss(y~x,dat,forceResponse = T)
  
  dat <- list(y=rnorm(100),x=rnorm(100),z=matrix(rnorm(200),100,2))
  dat$y[50] <- NA
  dat$x[99] <- NA
  dat$z[cbind(c(1,2,5),c(2,1,2))] <- NA
  r <- nmiss(y~x,dat)
  r <- nmiss(y~x,dat,forceResponse = T)
  r <- nmiss(y~z,dat)    
  r <- nmiss(y~z,dat,forceResponse = T)    
  
}