#' @export
prod_exp <-
function (d, phi, pow = rep(2, length(phi)), eps = 0) 
{
    if (length(phi) != length(d)) 
        stop("Length of phi must match length of d (total dimensions)")
    Sigma <- exp(-phi[1] * d[[1]]^pow[1])
    if (length(phi) > 1) 
        for (i in 2:length(phi)) Sigma <- Sigma * exp(-phi[i] * 
            d[[i]]^pow[i])
    return(as.matrix(Sigma) + (1 + eps) * diag(attr(d[[1]], "Size")))
}
