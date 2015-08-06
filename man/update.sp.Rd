\name{update.sp}
\alias{update.sp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
temptitle
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
update.sp(init, cutoff)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{init}{
%%     ~~Describe \code{init} here~~
}
  \item{cutoff}{
%%     ~~Describe \code{cutoff} here~~
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
function (init, cutoff) 
{
    p <- nrow(init$dmat)
    index <- rep(FALSE, init$numlt)
    datalist <- .Fortran("inrangeupdate", p = as.integer(nrow(init$dmat)), 
        numlt = as.integer(init$numlt), dmat = as.single(init$dmat), 
        index = as.logical(index), cf = as.single(cutoff), PACKAGE = "SparseEm")
    rowinds <- init$rowinds[datalist$index]
    rowpointers <- as.integer(cumsum(c(1, tabulate(rowinds, nbins = init$n))))
    return(list(colinds = init$colinds[datalist$index], rowpointers = rowpointers, 
        dmat = init$dmat[, datalist$index, drop = FALSE], numlt = sum(abs(datalist$index)), 
        n = init$n, p = init$p, cutoff = cutoff))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
