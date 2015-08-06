\name{adapt.var}
\alias{adapt.var}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
temptitle
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
adapt.var(b, samples, pl, rate.update, k = 3, c0 = 10, c1 = 0.8, r.opt = 0.234)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{b}{
%%     ~~Describe \code{b} here~~
}
  \item{samples}{
%%     ~~Describe \code{samples} here~~
}
  \item{pl}{
%%     ~~Describe \code{pl} here~~
}
  \item{rate.update}{
%%     ~~Describe \code{rate.update} here~~
}
  \item{k}{
%%     ~~Describe \code{k} here~~
}
  \item{c0}{
%%     ~~Describe \code{c0} here~~
}
  \item{c1}{
%%     ~~Describe \code{c1} here~~
}
  \item{r.opt}{
%%     ~~Describe \code{r.opt} here~~
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
function (b, samples, pl, rate.update, k = 3, c0 = 10, c1 = 0.8, 
    r.opt = 0.234) 
{
    gamma1 <- c0/(b/rate.update + k)^c1
    gamma2 <- 1/(b/rate.update + k)^c1
    r.hat <- mean(!duplicated(samples))
    sigma.m <- exp(log(pl$sigma.m) + gamma1 * (r.hat - r.opt))
    Sigma <- pl$Sigma + gamma2 * (cov(samples) - pl$Sigma)
    pl <- list(sigma.m = sigma.m, Sigma = Sigma)
    pl$C <- pl$sigma.m * chol(pl$Sigma)
    return(pl)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
