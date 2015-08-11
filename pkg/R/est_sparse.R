#' @export
est_sparse = function (cutoff, x, inflate = 1.2) 
{
  x <- x[sample(1:nrow(x), min(50 * ncol(x), nrow(x))), ]
  d <- mult_dist(x)
  sm <- d[[1]]
  if (length(d) > 1) 
    for (i in 2:length(d)) sm <- sm + d[[i]]
    return(min(mean(sm < cutoff) * inflate, 1))
}
