\name{inll}
\alias{inll}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
temptitle
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
inll(phi, y, d, pow, X)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phi}{
%%     ~~Describe \code{phi} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{d}{
%%     ~~Describe \code{d} here~~
}
  \item{pow}{
%%     ~~Describe \code{pow} here~~
}
  \item{X}{
%%     ~~Describe \code{X} here~~
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
function (phi, y, d, pow, X) 
{
    n <- length(y)
    corr.matrix <- prod.exp(d, phi = phi, pow = rep(pow, length(phi)))
    Q <- chol(corr.matrix)
    b1 <- backsolve(Q, X, transpose = TRUE)
    b2 <- backsolve(Q, y, transpose = TRUE)
    beta.hat <- solve(crossprod(b1), crossprod(b1, b2))
    b3 <- b2 - b1 \%*\% beta.hat
    s2.phi <- crossprod(b3)
    l1 <- 2 * sum(log(diag(Q)))
    l2 <- determinant(crossprod(b1), log = TRUE)$modulus
    return(drop(0.5 * (l1 + l2) + (n - ncol(X))/2 * log(s2.phi)))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
