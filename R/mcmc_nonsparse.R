#' @export
mcmc_nonsparse <-
function (y, x, degree = 0, maxint = 0, pow = 2, B = 2000, phi.max = 100, 
    phi.start = rep(1, ncol(x)), rate.update = 10 * ncol(x), 
    verbose = TRUE) 
{
    dim <- ncol(x)
    n <- nrow(x)
    terms <- polySet(p = dim, d = degree, w = maxint, m = degree)
    if (nrow(terms) > n) 
        stop("Specified degree and maxint give more regression terms than observations")
    leg01 <- legFun(0, 1)
    X <- leg01(x, terms = terms)
    d <- mult_dist(x)
    loglik.old <- -inll(phi.start, y, d, pow, X)
    pl <- list(sigma.m = 2.38/sqrt(dim), Sigma = diag(phi.start^2/9, 
        nrow = dim))
    pl$C <- pl$sigma.m * chol(pl$Sigma)
    phi <- matrix(NA, nrow = B, ncol = dim)
    phi[1, ] <- phi.start
    if (verbose) 
        start.time <- proc.time()
    for (b in 2:B) {
        if (verbose && b%%10 == 0) 
            print(b)
        prop <- phi[b - 1, ] + rnorm(dim) %*% pl$C
        if (any(prop < 0 | prop > phi.max)) {
            phi[b, ] <- phi[b - 1, ]
        }
        else {
            loglik.new <- -inll(prop, y, d, pow, X)
            r <- exp(loglik.new - loglik.old)
            if (r > runif(1)) {
                loglik.old <- loglik.new
                phi[b, ] <- prop
            }
            else {
                phi[b, ] <- phi[b - 1, ]
            }
        }
        if (verbose && b == 11) {
            est.time <- (proc.time()[3] - start.time[3])/10 * 
                (B - 11)/3600
            print(paste("Estimated running time:", round(est.time, 
                2), "hours"))
        }
        if (b%%rate.update == 0) {
            pl <- adapt.var(b, samples = phi[(b - rate.update + 
                1):b, ], pl, rate.update)
        }
    }
    return(phi)
}
