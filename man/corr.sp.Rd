\name{corr.sp}
\alias{corr.sp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
temptitle
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
corr.sp(setup, cor = "truncpow", cor.args = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{setup}{
%%     ~~Describe \code{setup} here~~
}
  \item{cor}{
%%     ~~Describe \code{cor} here~~
}
  \item{cor.args}{
%%     ~~Describe \code{cor.args} here~~
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
function (setup, cor = "truncpow", cor.args = NULL) 
{
    if (setup$numlt == 0) 
        return(diag.spam(setup$n))
    newentries <- do.call(cor, c(list(d = setup$dmat[1, ], cutoff = setup$cutoff[1]), 
        cor.args))
    if (setup$p > 1) 
        for (i in 2:setup$p) newentries <- newentries * do.call(cor, 
            c(list(d = setup$dmat[i, ], cutoff = setup$cutoff[i]), 
                cor.args))
    sp <- new("spam", entries = newentries, colindices = setup$colinds, 
        rowpointers = setup$rowpointers, dimension = rep(setup$n, 
            2))
    sp <- sp + t(sp) + diag.spam(setup$n)
    return(sp)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
