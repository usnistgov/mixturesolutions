\name{doubleconsensus}
\alias{doubleconsensus}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
doubleconsensus(dataf, preps)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataf}{
%%     ~~Describe \code{dataf} here~~
}
  \item{preps}{
%%     ~~Describe \code{preps} here~~
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
function (dataf, preps) 
{
    fitmodel <- function(indf, preps) {
        fit <- switch(preps$modeltype, twomix = minmodelsolve(indf, 
            preps), threemix = magicmodelsolve(indf, preps))
        return(fit)
    }
    f2opt <- function(par, indata, preps) {
        prepsn <- c(preps, list(mfrac = c(par, 1)))
        (sum(abs(fitmodel(indata, preps = prepsn)[1, 1:length(preps$trueproportions[, 
            1])] - preps$trueproportions[, 1]), abs(fitmodel(indata, 
            preps = prepsn)[2, 1:length(preps$trueproportions[, 
            2])] - preps$trueproportions[, 2])))
    }
    opted <- optim(par = c(1, 1), f2opt, indata = dataf, preps = preps)$par
    return(c(opted, 1)/sum(c(opted, 1)))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
