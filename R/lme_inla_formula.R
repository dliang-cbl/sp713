lme.inla.formula <- function(formula,data,param=NULL)
{
  ## helper function to formulate INLA formula
  ## formula: mixed effect formula in lme4 syntax
  ## data : data frame to evaluate the data
  ## param: a list of hyper-parameters consistent with each random terms
  ##      for univariate random, default to (0,0.00001) for truncated normal prior
  ##      for multivariate random, default to identity metrix and df equal 2 plus dimension
  ##      multivariate random can't have more than five terms due to INLA constraints
  ## value: a vector of INLA formula terms
  ##        index data frame useful for model fitting
  ## assumption: only single group, no nesting allowed
  ## 
  ## this function is based on the proc. formula function in marked package
  ##
  lst <- proc.form(formula) ## process mixed formula
  
  if(length(lst$re.model)==0){
    return(
      list(fixed=as.formula(lst$fix.model),formulas=labels(terms(formula)),index=NULL)
    )
  }
  r <- vector("list",length(lst$re.model)) 
  for(i in 1:length(r)){  ## loop through random terms
    elmt <- lst$re.model[[i]]
    dm <- model.matrix(formula(elmt$model),data)
    v0 <- all.vars(formula(elmt$model))
    g0 <- all.vars(formula(elmt$sub))
    
    ## make copies of index
    ## ASSUMING ONE GROUP ONLY, No NESTING
    if(ncol(dm)==1){
      ind1 <- NULL
    } 
    else{
      ng <- length(unique(data[,g0]))
      ind1 <- matrix(as.numeric(data[,g0]),nrow(data),ncol(dm))
      ind1 <- ind1+rep((0:(ncol(dm)-1))*ng,each=nrow(data))
      colnames(ind1) <- paste(g0,1:ncol(dm),sep="")
    }
    
    ## make formula
    out <- rep("",ncol(dm))
    if(is.null(param)){
      ## default hyper-parameter
      tmp <- diag(ncol(dm))
      parami <- c(ncol(dm)+2,1,1,0)
    }
    else{
      parami <- param[[i]]
    }
    if(ncol(dm)==1){
      out[1] <- paste(
        "f(",g0,",",
        'model="iid",',
        'hyper=list(prec=list(prior="logtnormal",param=c(0,0.00001))))',
        sep="")
    }
    if(ncol(dm)>5){
      stop(ncol(dm)," dimension not supported by INLA\n")
    }
    if(ncol(dm)>1){
      out[1] <- paste(
        "f(",colnames(ind1)[1],",",
        'model="iid',ncol(dm),'d",',
        "n=",ncol(dm)*ng,
        ",hyper=list(theta1=list(param=c(",paste(parami,collapse =","),"))))",
        sep="")
      for(j in 2:ncol(dm)){
        out[j] <- paste("f(",colnames(ind1)[j],",",
                        colnames(dm)[j],',copy="',
                        colnames(ind1)[1],'")',sep="")
      }
    }
    r[[i]] <- list(formula=out,index=ind1)
  }
  
  tmp <- do.call(c,lapply(r,function(elmt) elmt$formula))
  idx <- do.call(cbind,lapply(r,function(elmt) elmt$index))
  
  fix.model <- labels(terms(as.formula(lst$fix.model)))
  #list(formulas=c(lst$fix.model,tmp),index=idx)
  list(fixed=as.formula(lst$fix.model),
       formulas=c(fix.model,tmp),index=idx)
}
