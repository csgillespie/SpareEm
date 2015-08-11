#' @export
polySet <-
function (p, d, w, m) 
{
    tmp <- c(p, d, w, m)
    if (any(tmp < 0)) 
        return(NULL)
    if (any(tmp == 0)) 
        return(matrix(0, 1, p))
    mdm <- 0:min(d, m)
    robj <- vector("list", length(mdm))
    for (i in seq(along = mdm)) {
        j <- mdm[i]
        robj[[i]] <- cbind(as.vector(j), Recall(p - 1, d - j, 
            w - (j > 0), m))
    }
    do.call("rbind", robj)
}
