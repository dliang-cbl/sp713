coeftab.getfixed <- function(x,intercept=TRUE){
  ## x: inla object 
  ## return: fixed estimates
  if(!intercept){
    if(row.names(x$summary.fixed)[1]=="(Intercept)"){
      fixed <- x$summary.fixed[-1,,drop=F] 
    }else{
      fixed <- x$summary.fixed
    }
  }else{
    fixed <- x$summary.fixed
  }
  return(fixed)
}
coeftab.work <- function(x,intercept=TRUE){
  ## x: inla object 
  ## value: dot chart of the fixed effects
  fixed <- coeftab.getfixed(x,intercept)
  # if(!intercept){
  #   if(row.names(x$summary.fixed)[1]=="(Intercept)"){
  #     fixed <- x$summary.fixed[-1,,drop=F] 
  #   }else{
  #     fixed <- x$summary.fixed
  #   }
  # }else{
  #   fixed <- x$summary.fixed
  # }
  mu <- fixed$`0.5quant`
  left <- fixed$`0.025quant`
  right <- fixed$`0.975quant`
  xlim <- c(min(left),max(right))
  dotchart(mu,labels = row.names(fixed),xlim=xlim)
  for ( i in 1:length(mu) ) lines( c(left[i],right[i]) , c(i,i) , lwd=2 )
}

coeftab <- function(...,intercept=TRUE){
  ## ... a list of R inla objects
  ## value: dot chart of the fixed effects, grouped by models
  #browser()
  names <- as.character(match.call())[-1]
  lst <- list(...)
  if(length(lst)==1){
    coeftab.work(lst[[1]],intercept = intercept)
    return(0)
  }
  fixed_ <- lapply(lst,function(x){
    fixed <-   coeftab.getfixed(x,intercept)
    mu <- fixed$`0.5quant`
    left <- fixed$`0.025quant`
    right <- fixed$`0.975quant`
    data.frame(name=row.names(fixed),mu,left,right)
  })
  for(i in 1:length(lst)) fixed_[[i]]$model <- names[i]
  #browser()
  fixed2 <- do.call(rbind,fixed_)
  
  mu <- pivot_wider(fixed2[,c(1,2,5)], names_from=name,values_from = mu)
  left <- pivot_wider(fixed2[,c(1,3,5)], names_from=name,values_from = left)
  right <- pivot_wider(fixed2[,c(1,4,5)], names_from=name,values_from = right)
  
  xlim <- c(min(left[,-1],na.rm=T),max(right[,-1],na.rm=T))
  dotchart(as.matrix(mu[,-1]),labels = mu$model,xlim=xlim)
  
  x <- mu[,-1]
  ## estimate the y positions on the dotchart
  kn <- nrow(x) ## group size
  base <- (col(x)-1)*(kn+2)
  ypos <- base + row(x)
  ypos2 <- ypos[,ncol(ypos):1]
  #abline(h=c(ypos))
  #browser()
  lwr <- unlist(left[,-1])
  upr <- unlist(right[,-1])
  
  for(i in 1:length(lwr)) {
    if(!is.na(lwr[i])){
      lines( c(lwr[i],upr[i]) , c(ypos2[i],ypos2[i]) , lwd=2  )
    }
  }
}