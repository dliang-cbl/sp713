plot.inla <-  function(
  x,## a fitted inla object producedby gamm.inla()
  pages=0, ## number of pages over which to spread the output
  select=NULL, ## single term to be selected)
  ... # arguments passed to lower level plotting function
){
  args__ <- as.list(match.call())
  
  stopifnot(select %in% seq(1,length(x$.args$sterms)))
  if(is.null(select)){
    sterms <- x$.args$sterms
  }else{
    sterms <- x$.args$sterms[select]
  }
  
  if(pages==1){
    op <- par(mfrow=c(ceiling(length(sterms)/2),2))
  }
  # }else{
  #   op <- par(mfrow=c(1,1))
  # }
  for(i in 1:length(sterms)){
    l.x <- x$.args$raw[,sterms[i]]
    l.fit <- x$summary.random[[sterms[i]]][
      1:nrow(x$.args$raw),
      4:6]
    if("xlab" %in% names(args__)){
      xlab <- args__$xlab
    }else{
      xlab <- sterms[i]
    }
    if("ylab" %in% names(args__)){
      ylab <- args__$ylab
    }else{
      ylab <- paste("s(",xlab,")",sep="")
    }
    plot.foo.work(l.x,l.fit,
                  xlab=xlab,
                  ylab=ylab,...)
  }
  if(pages==1){
    par(op)  
  }
}
plot.foo.work <- function(x,y,...){
  l.oo <- order(x)
  matplot(x[l.oo],y[l.oo,],
          type="l",lty=c(2,1,2),col=1,...)
  rug(x)
}