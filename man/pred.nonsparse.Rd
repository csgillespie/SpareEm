\name{pred.nonsparse}
\alias{pred.nonsparse}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
temptitle
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pred.nonsparse(phi, x, y, xpred, pow, degree = 0, maxint = 0, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phi}{
%%     ~~Describe \code{phi} here~~
}
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{xpred}{
%%     ~~Describe \code{xpred} here~~
}
  \item{pow}{
%%     ~~Describe \code{pow} here~~
}
  \item{degree}{
%%     ~~Describe \code{degree} here~~
}
  \item{maxint}{
%%     ~~Describe \code{maxint} here~~
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
    d <- mult.dist(rbind(x, xpred))
    m <- v <- matrix(NA, nrow = nrow(phi), ncol = npred)
    if (verbose) 
        start.time <- proc.time()
    for (i in 1:nrow(phi)) {
        if (verbose && i\%\%10 == 0) 
            print(i)
        corr.matrix <- prod.exp(d, phi = phi[i, ], pow = rep(pow, 
            dim))
        Q <- chol(corr.matrix[1:n, 1:n])
        b1 <- backsolve(Q, X, transpose = TRUE)
        b2 <- backsolve(Q, y, transpose = TRUE)
        beta.hat <- solve(crossprod(b1), crossprod(b1, b2))
        b3 <- b2 - b1 \%*\% beta.hat
        s2.phi <- drop(crossprod(b3))
        b4 <- backsolve(Q, corr.matrix[1:n, -(1:n)], transpose = TRUE)
        Qb <- chol(crossprod(b1))
        b5 <- backsolve(Qb, t(Xpred), transpose = TRUE)
        m[i, ] <- drop(Xpred \%*\% beta.hat + crossprod(b4, b3))
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
