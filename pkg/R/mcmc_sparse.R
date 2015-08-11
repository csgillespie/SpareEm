#' @export
mcmc_sparse <-
function (y, x, mc = find_tau(den = 0.1, dim = ncol(x)) * ncol(x), 
    degree = 4, maxint = 2, cor = "truncpow", cor.args = NULL, 
    B = 2000, tau.start = rep(mc/ncol(x)/1.5, ncol(x)), rate.update = 10 * 
        ncol(x), verbose = TRUE) 
{
    dim <- ncol(x)
    n <- nrow(x)
    terms <- polySet(p = dim, d = degree, w = maxint, m = degree)
    if (nrow(terms) > n) 
        stop("Specified degree and maxint give more regression terms than observations")
    leg01 <- legFun(0, 1)
    X <- leg01(x, terms = terms)
    init <- init_sp(x, cutoff = mc)
    loglik.old <- -inll_sparse(tau.start, y, init, cor = cor, 
        cor.args = cor.args, X = X)
    pl <- list(sigma.m = 2.38/sqrt(dim), Sigma = diag(tau.start^2/9, 
        nrow = dim))
    pl$C <- pl$sigma.m * chol(pl$Sigma)
    tau <- matrix(NA, nrow = B, ncol = dim)
    tau[1, ] <- tau.start
    if (verbose) 
        start.time <- proc.time()
    for (b in 2:B) {
        if (verbose && b%%10 == 0) 
            message(b, " iterations completed.")
        prop <- tau[b - 1, ] + rnorm(dim) %*% pl$C
        if (any(prop < 0) | sum(prop) > mc) {
            tau[b, ] <- tau[b - 1, ]
        }
        else {
            loglik.new <- -inll_sparse(prop, y, init, cor, cor.args, X)
            r <- exp(loglik.new - loglik.old)
            if (r > runif(1)) {
                loglik.old <- loglik.new
                tau[b, ] <- prop
            }
            else {
                tau[b, ] <- tau[b - 1, ]
            }
        }
        if (verbose && b == 11) {
            est.time <- (proc.time()[3] - start.time[3])/10 * 
                (B - 11)/3600
            message("Estimated running time: ", round(est.time, 2), " hours.")
        }
        if (b%%rate.update == 0) {
            pl <- adapt_var(b, samples = tau[(b - rate.update + 1):b, ], pl, rate.update)
        }
    }
    return(tau)
}
