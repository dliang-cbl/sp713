#' deviance information criterion
#' 
#' Model comparison using deviance information criterion
#' 
#' @param ... a list of model objects
#' @return fitted DIC values
#' 
#' @export

DIC <- function(...){
  ## potentially allow inla output directly
  ## and name the rows correctly, not for now.
  names <- deparse(substitute(...))
  lst <- list(...)
  r <- lapply(lst,function(elmt) {
    if(!is.null(elmt$DIC)){
      r <- elmt$DIC
    }
    else if (class(elmt)=="inla"){
      if(!is.null(elmt$dic)){
        r <- with(elmt$dic,c(Dbar=mean.deviance,pD=p.eff,DIC=dic))
      }
      else{
        r <- rep(NA,3)
      }
    }
    else{
      r <- rep(NA,3)
    }
    r
  })
  r2 <- do.call(rbind,r)
  
  #browser()
  #nms <- setdiff(as.character(match.call(expand.dots=TRUE)), 
  #               as.character(match.call(expand.dots=FALSE)))
  nm1 <- as.character(match.call(expand.dots=TRUE))
  nms <- nm1[which(!(nm1 %in% as.character(match.call(expand.dots=FALSE))))]
  
  for(i in 1:length(lst)){
    elmt <- lst[[i]]
    if(!is.null(elmt$formula)){
      cat("model ",i," ",nms[i],":",deparse(elmt$formula),"\n")
    }
    else if (class(elmt)=="inla"){
      cat("model ",i," ",nms[i],":", deparse(elmt$.args$formula),"\n")
    }
    else{
      cat("model ",i,"",nms[i],": formula not found\n")
    }
  }
  row.names(r2) <- make.unique(nms)
  r2
}