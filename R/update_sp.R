#' @export
update_sp <- function (init, cutoff) 
{
    p <- nrow(init$dmat)
    index <- rep(FALSE, init$numlt)
    datalist <- .Fortran("inrangeupdate", p = as.integer(nrow(init$dmat)), 
        numlt = as.integer(init$numlt), dmat = as.single(init$dmat), 
        index = as.logical(index), cf = as.single(cutoff), PACKAGE = "SparseEm")
    rowinds <- init$rowinds[datalist$index]
    rowpointers <- as.integer(cumsum(c(1, tabulate(rowinds, nbins = init$n))))
    return(list(colinds = init$colinds[datalist$index], rowpointers = rowpointers, 
        dmat = init$dmat[, datalist$index, drop = FALSE], numlt = sum(abs(datalist$index)), 
        n = init$n, p = init$p, cutoff = cutoff))
}
