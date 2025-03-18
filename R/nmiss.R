## helper script to subset to data where the formula is not missing.
nmiss <- function(formula,data)
{
  if(is.null(data)){
    return(data)
  }
  tmp <- na.omit(data[,all.vars(formula)])
  #browser()
  if(!is.null(attr(tmp,"na.action"))){
    data <- data[-attr(tmp,"na.action"),]
  }
  data
}