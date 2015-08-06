#' @export
mult_dist <-
function (x) 
{
    if (is.null(ncol(x))) 
        return(list(dist(x)))
    lapply(1:ncol(x), function(i) {
        dist(x[, i])
    })
}
