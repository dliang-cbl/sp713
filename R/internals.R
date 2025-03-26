col.alpha <- function (acol, alpha = 0.2) 
{
  acol <- col2rgb(acol)
  acol <- rgb(acol[1]/255, acol[2]/255, acol[3]/255, alpha)
  acol
}
set_nice_margins <- function () 
{
  par_mf <- par("mfrow", "mfcol")
  if (all(unlist(par_mf) == 1)) {
    par(mgp = c(1.5, 0.5, 0), mar = c(2.5, 2.5, 2, 1) + 0.1, 
        tck = -0.02)
  }
}