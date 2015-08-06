#' @export
bohman <-
function (d, cutoff) 
{
    r <- d/cutoff
    return(((1 - r) * cos(pi * r) + sin(pi * r)/pi) * (r <= 1))
}
