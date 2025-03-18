## get response variable (if any) and predictor variables from a formula object
getVars <- function(formula)
{
  if(length(formula)<3){
    ## no response
    r <- c("",all.vars(formula))
  }
  else{
    r <- all.vars(formula)
  }
  r
}