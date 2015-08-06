\name{polySet}
\alias{polySet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
temptitle
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
polySet(p, d, w, m)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{
%%     ~~Describe \code{p} here~~
}
  \item{d}{
%%     ~~Describe \code{d} here~~
}
  \item{w}{
%%     ~~Describe \code{w} here~~
}
  \item{m}{
%%     ~~Describe \code{m} here~~
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
function (p, d, w, m) 
{
    tmp <- c(p, d, w, m)
    if (any(tmp < 0)) 
        return(NULL)
    if (any(tmp == 0)) 
        return(matrix(0, 1, p))
    mdm <- 0:min(d, m)
    robj <- vector("list", length(mdm))
    for (i in seq(along = mdm)) {
        j <- mdm[i]
        robj[[i]] <- cbind(as.vector(j), Recall(p - 1, d - j, 
            w - (j > 0), m))
    }
    do.call("rbind", robj)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
