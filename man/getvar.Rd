\name{getvar}
\alias{getvar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
temptitle
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getvar(name, gs, burnin = 0, drop = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{name}{
%%     ~~Describe \code{name} here~~
}
  \item{gs}{
%%     ~~Describe \code{gs} here~~
}
  \item{burnin}{
%%     ~~Describe \code{burnin} here~~
}
  \item{drop}{
%%     ~~Describe \code{drop} here~~
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
function (name, gs, burnin = 0, drop = 1) 
{
    ok <- !sapply(gs, is.null)
    if (!(length(dim(gs[[1]][[name]])) >= 2)) {
        return(sapply(gs[ok][seq(burnin + 1, sum(ok), by = drop)], 
            function(x) {
                x[[name]]
            }))
    }
    else {
        temp <- lapply(gs[ok][seq(burnin + 1, sum(ok), by = drop)], 
            function(x) {
                x[[name]]
            })
        return(structure(unlist(temp), dim = c(dim(temp[[1]]), 
            length(temp))))
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
