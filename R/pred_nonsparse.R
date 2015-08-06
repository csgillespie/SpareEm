#' @export
pred_nonsparse <-
function (phi, x, y, xpred, pow, degree = 0, maxint = 0, verbose = TRUE) 
{
    dim <- ncol(x)
    n <- nrow(x)
    npred <- nrow(xpred)
    terms <- polySet(p = dim, d = degree, w = maxint, m = degree)
    if (nrow(terms) > n) 
        stop("Specified degree and maxint give more regression terms than observations")
    leg01 <- legFun(0, 1)
    X <- leg01(x, terms = terms)
    Xpred <- leg01(xpred, terms = terms)
    p <- ncol(X)
    d <- mult_dist(rbind(x, xpred))
    m <- v <- matrix(NA, nrow = nrow(phi), ncol = npred)
    if (verbose) 
        start.time <- proc.time()
    for (i in 1:nrow(phi)) {
        if (verbose && i%%10 == 0) 
            print(i)
        corr.matrix <- prod_exp(d, phi = phi[i, ], pow = rep(pow, 
            dim))
        Q <- chol(corr.matrix[1:n, 1:n])
        b1 <- backsolve(Q, X, transpose = TRUE)
        b2 <- backsolve(Q, y, transpose = TRUE)
        beta.hat <- solve(crossprod(b1), crossprod(b1, b2))
        b3 <- b2 - b1 %*% beta.hat
        s2.phi <- drop(crossprod(b3))
        b4 <- backsolve(Q, corr.matrix[1:n, -(1:n)], transpose = TRUE)
        Qb <- chol(crossprod(b1))
        b5 <- backsolve(Qb, t(Xpred), transpose = TRUE)
        m[i, ] <- drop(Xpred %*% beta.hat + crossprod(b4, b3))
        v[i, ] <- 1/(n - p - 2) * s2.phi * (1 - colSums(b4^2) + 
            colSums(b5^2))
        if (verbose && i == 10) {
            est.time <- (proc.time()[3] - start.time[3])/10 * 
                (nrow(phi) - 10)/3600
            print(paste("Estimated running time:", round(est.time, 
                2), "hours"))
        }
    }
    ypred.postmean <- apply(m, 2, mean)
    ypred.postvar <- apply(v, 2, mean) + apply(m, 2, var)
    return(list(mean = ypred.postmean, var = ypred.postvar))
}
