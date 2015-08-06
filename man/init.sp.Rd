\name{init.sp}
\alias{init.sp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
temptitle
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
init.sp(x, cutoff, sparsity.estimate = est.sparse(cutoff, x))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{cutoff}{
%%     ~~Describe \code{cutoff} here~~
}
  \item{sparsity.estimate}{
%%     ~~Describe \code{sparsity.estimate} here~~
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
function (x, cutoff, sparsity.estimate = est.sparse(cutoff, x)) 
{
    n = nrow(x)
    p = ncol(x)
    if (length(cutoff) > 1) 
        stop("Single cutoff needed for initialization, applied to sum")
    maxnlt = round(n * (n - 1)/2 * sparsity.estimate)
    rowinds = rep(0, maxnlt)
    colinds = rep(0, maxnlt)
    dmat = matrix(0, nrow = p, ncol = maxnlt)
    datalist <- .Fortran("distintri", nrc = as.integer(n), np = as.integer(p), 
        as.single(x), cf = as.single(cutoff), nlt = as.integer(1), 
        nm = as.integer(maxnlt), rc = as.integer(0), rowinds = as.integer(rowinds), 
        colinds = as.integer(colinds), dmat = as.single(dmat), 
        PACKAGE = "SparseEm")
    while (datalist$rc < 0) {
        warning(paste("Overflow: increasing sparsity.estimate up from", 
            round(sparsity.estimate, 3)))
        sparsity.estimate <- min((sparsity.estimate + 0.001) * 
            2, 1)
        maxnlt = round(n * (n - 1)/2 * sparsity.estimate)
        rowinds = rep(0, maxnlt)
        colinds = rep(0, maxnlt)
        dmat = matrix(0, nrow = p, ncol = maxnlt)
        datalist <- .Fortran("distintri", nrc = as.integer(n), 
            np = as.integer(p), as.single(x), cf = as.single(cutoff), 
            nlt = as.integer(1), nm = as.integer(maxnlt), rc = as.integer(0), 
            rowinds = as.integer(rowinds), colinds = as.integer(colinds), 
            dmat = as.single(dmat), PACKAGE = "SparseEm")
    }
    return(list(colinds = datalist$colinds[1:datalist$nlt], rowinds = datalist$rowinds[1:datalist$nlt], 
        dmat = matrix(datalist$dmat[1:(datalist$nlt * p)], nrow = p), 
        numlt = datalist$nlt, n = n, p = p, cutoff = cutoff))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
