#' @export
find_tau <-
function (den, dim) 
{
    f.obj <- function(x, dim, den) {
        estden <- x^dim * (2 - x)^dim
        return((estden - den)^2)
    }
    opt <- optimize(f.obj, lower = 0, upper = 1, dim = dim, den = den)
    return(opt$minimum)
}
