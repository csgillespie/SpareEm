#' @useDynLib SparseEm
#' @export
init_sp <-
function (x, cutoff, sparsity.estimate = est_sparse(cutoff, x)) 
{
    n = nrow(x)
    p = ncol(x)
    if (length(cutoff) > 1) 
        stop("Single cutoff needed for initialization, applied to sum")
    maxnlt = round(n * (n - 1)/2 * sparsity.estimate)
    rowinds = rep(0, maxnlt)
    colinds = rep(0, maxnlt)
    dmat = matrix(0, nrow = p, ncol = maxnlt)
    datalist <- .Fortran("distintri", nrc = as.integer(n), np = as.integer(p), 
        as.single(x), cf = as.single(cutoff), nlt = as.integer(1), 
        nm = as.integer(maxnlt), rc = as.integer(0), rowinds = as.integer(rowinds), 
        colinds = as.integer(colinds), dmat = as.single(dmat), 
        PACKAGE = "SparseEm")
    while (datalist$rc < 0) {
        warning(paste("Overflow: increasing sparsity.estimate up from", 
            round(sparsity.estimate, 3)))
        sparsity.estimate <- min((sparsity.estimate + 0.001) * 
            2, 1)
        maxnlt = round(n * (n - 1)/2 * sparsity.estimate)
        rowinds = rep(0, maxnlt)
        colinds = rep(0, maxnlt)
        dmat = matrix(0, nrow = p, ncol = maxnlt)
        datalist <- .Fortran("distintri", nrc = as.integer(n), 
            np = as.integer(p), as.single(x), cf = as.single(cutoff), 
            nlt = as.integer(1), nm = as.integer(maxnlt), rc = as.integer(0), 
            rowinds = as.integer(rowinds), colinds = as.integer(colinds), 
            dmat = as.single(dmat), PACKAGE = "SparseEm")
    }
    return(list(colinds = datalist$colinds[1:datalist$nlt], rowinds = datalist$rowinds[1:datalist$nlt], 
        dmat = matrix(datalist$dmat[1:(datalist$nlt * p)], nrow = p), 
        numlt = datalist$nlt, n = n, p = p, cutoff = cutoff))
}
