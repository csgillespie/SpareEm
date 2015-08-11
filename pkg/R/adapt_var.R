#' Adapt var
#' @export
adapt_var = function (b, samples, pl, rate.update, 
                      k = 3, c0 = 10, c1 = 0.8, r.opt = 0.234) 
{
  gamma1 <- c0/(b/rate.update + k)^c1
  gamma2 <- 1/(b/rate.update + k)^c1
  r.hat <- mean(!duplicated(samples))
  sigma.m <- exp(log(pl$sigma.m) + gamma1 * (r.hat - r.opt))
  Sigma <- pl$Sigma + gamma2 * (cov(samples) - pl$Sigma)
  pl <- list(sigma.m = sigma.m, Sigma = Sigma)
  pl$C <- pl$sigma.m * chol(pl$Sigma)
  return(pl)
}
