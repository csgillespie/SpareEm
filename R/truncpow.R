#' @export
truncpow <-
function (d, cutoff, alpha = 1) 
{
    stopifnot(alpha %in% c(1, 3/2, 5/3))
    nu <- switch(match(alpha, c(1, 3/2, 5/3)), 1, 2, 3)
    r <- d/cutoff
    return((1 - r^alpha)^nu * (r <= 1))
}
