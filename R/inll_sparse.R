#' @importFrom spam forwardsolve diag
#' @export
inll_sparse <-
function (tau, y, init, cor, cor.args, X) 
{
    n <- length(y)
    setup <- update_sp(init, cutoff = tau)
    G <- corr_sp(setup, cor = cor, cor.args = cor.args)
    Q <- spam::chol(G)
#    rm(G)
#    gc()
    b1 <- spam::forwardsolve(Q, X)
    b2 <- spam::forwardsolve(Q, y)
    beta.hat <- solve(crossprod(b1), crossprod(b1, b2))
    b3 <- b2 - b1 %*% beta.hat
    s2.phi <- drop(crossprod(b3))
    l1 <- 2 * sum(log(spam::diag(Q)))
    l2 <- determinant(crossprod(b1), logarithm = TRUE)$modulus
    return(drop(0.5 * (l1 + l2) + (n - ncol(X))/2 * log(s2.phi)))
}
