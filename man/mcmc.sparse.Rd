\name{mcmc.sparse}
\alias{mcmc.sparse}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
temptitle
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mcmc.sparse(y, x, mc = find.tau(den = 0.1, dim = ncol(x)) * ncol(x), degree = 4, maxint = 2, cor = "truncpow", cor.args = NULL, B = 2000, tau.start = rep(mc/ncol(x)/1.5, ncol(x)), rate.update = 10 * ncol(x), verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{mc}{
%%     ~~Describe \code{mc} here~~
}
  \item{degree}{
%%     ~~Describe \code{degree} here~~
}
  \item{maxint}{
%%     ~~Describe \code{maxint} here~~
}
  \item{cor}{
%%     ~~Describe \code{cor} here~~
}
  \item{cor.args}{
%%     ~~Describe \code{cor.args} here~~
}
  \item{B}{
%%     ~~Describe \code{B} here~~
}
  \item{tau.start}{
%%     ~~Describe \code{tau.start} here~~
}
  \item{rate.update}{
%%     ~~Describe \code{rate.update} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (y, x, mc = find.tau(den = 0.1, dim = ncol(x)) * ncol(x), 
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
    init <- init.sp(x, cutoff = mc)
    loglik.old <- -inll.sparse(tau.start, y, init, cor = cor, 
        cor.args = cor.args, X = X)
    pl <- list(sigma.m = 2.38/sqrt(dim), Sigma = diag(tau.start^2/9, 
        nrow = dim))
    pl$C <- pl$sigma.m * chol(pl$Sigma)
    tau <- matrix(NA, nrow = B, ncol = dim)
    tau[1, ] <- tau.start
    if (verbose) 
        start.time <- proc.time()
    for (b in 2:B) {
        if (verbose && b\%\%10 == 0) 
            print(b)
        prop <- tau[b - 1, ] + rnorm(dim) \%*\% pl$C
        if (any(prop < 0) | sum(prop) > mc) {
            tau[b, ] <- tau[b - 1, ]
        }
        else {
            loglik.new <- -inll.sparse(prop, y, init, cor, cor.args, 
                X)
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
            print(paste("Estimated running time:", round(est.time, 
                2), "hours"))
        }
        if (b\%\%rate.update == 0) {
            pl <- adapt.var(b, samples = tau[(b - rate.update + 
                1):b, ], pl, rate.update)
        }
    }
    return(tau)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
