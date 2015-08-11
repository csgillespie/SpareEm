#' @export
pred_sparse <-
function (tau, x, y, xpred, cor = "truncpow", cor.args = NULL, 
    degree = 4, maxint = 2, verbose = TRUE) 
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
    init <- init_sp(x = rbind(x, xpred), cutoff = max(apply(tau, 
        1, sum)))
    m <- v <- matrix(NA, nrow = nrow(tau), ncol = npred)
    if (verbose) 
        start.time <- proc.time()
    for (i in 1:nrow(tau)) {
        if (verbose) 
            print(i)
        setup <- update_sp(init, cutoff = tau[i, ])
        G <- corr_sp(setup, cor = cor, cor.args = cor.args)
        Q <- spam::chol(G[1:n, 1:n])
        b1 <- forwardsolve(Q, X)
        b2 <- forwardsolve(Q, y)
        beta.hat <- solve(t(b1) %*% b1, t(b1) %*% b2)
        b3 <- b2 - b1 %*% beta.hat
        s2.tau <- drop(t(b3) %*% b3)
        b4 <- forwardsolve(Q, G[1:n, -(1:n)])
        Qb <- chol(crossprod(b1))
        b5 <- backsolve(Qb, t(Xpred), transpose = TRUE)
        m[i, ] <- drop(Xpred %*% beta.hat + t(b4) %*% b3)
        v[i, ] <- 1/(n - p - 2) * s2.tau * (1 - colSums(b4^2) + 
            colSums(b5^2))
        if (verbose && i == 10) {
            est.time <- (proc.time()[3] - start.time[3])/10 * 
                (nrow(tau) - 10)/3600
            print(paste("Estimated running time:", round(est.time, 
                2), "hours"))
        }
    }
    ypred.postmean <- apply(m, 2, mean)
    ypred.postvar <- apply(v, 2, mean) + apply(m, 2, var)
    return(list(mean = ypred.postmean, var = ypred.postvar))
}
