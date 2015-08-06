\name{find.tau}
\alias{find.tau}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
temptitle
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
find.tau(den, dim)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{den}{
%%     ~~Describe \code{den} here~~
}
  \item{dim}{
%%     ~~Describe \code{dim} here~~
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
function (den, dim) 
{
    f.obj <- function(x, dim, den) {
        estden <- x^dim * (2 - x)^dim
        return((estden - den)^2)
    }
    opt <- optimize(f.obj, lower = 0, upper = 1, dim = dim, den = den)
    return(opt$minimum)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
