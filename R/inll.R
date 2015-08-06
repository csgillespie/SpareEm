#' @export
inll <-
function (phi, y, d, pow, X) 
{
    n <- length(y)
    corr.matrix <- prod_exp(d, phi = phi, pow = rep(pow, length(phi)))
    Q <- chol(corr.matrix)
    b1 <- backsolve(Q, X, transpose = TRUE)
    b2 <- backsolve(Q, y, transpose = TRUE)
    beta.hat <- solve(crossprod(b1), crossprod(b1, b2))
    b3 <- b2 - b1 %*% beta.hat
    s2.phi <- crossprod(b3)
    l1 <- 2 * sum(log(diag(Q)))
    l2 <- determinant(crossprod(b1), logarithm = TRUE)$modulus
    return(drop(0.5 * (l1 + l2) + (n - ncol(X))/2 * log(s2.phi)))
}
