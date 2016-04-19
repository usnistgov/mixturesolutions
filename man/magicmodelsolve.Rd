\name{magicmodelsolve}
\alias{magicmodelsolve}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
magicmodelsolve(ssx = indf, preps, idvars)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ssx}{
%%     ~~Describe \code{ssx} here~~
}
  \item{preps}{
%%     ~~Describe \code{preps} here~~
}
  \item{idvars}{
%%     ~~Describe \code{idvars} here~~
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
function (ssx = indf, preps, idvars) 
{
    labround <- NULL
    for (idx in 1:length(preps$mixnames)) {
        labround <- rbind(labround, c(coefficients(lm(data = ssx, 
            preps$mixnames[idx] ~ I(preps$componentnames[1] * 
                preps$mfrac[1]) + I(preps$componentnames[2] * 
                preps$mfrac[2]) + I(preps$componentnames[3] * 
                preps$mfrac[3]) + 0))/sum(coefficients(lm(data = ssx, 
            preps$mixnames[idx] ~ I(preps$componentnames[1] * 
                preps$mfrac[1]) + I(preps$componentnames[2] * 
                preps$mfrac[2]) + I(preps$componentnames[3] * 
                preps$mfrac[3]) + 0))), preps$mixnames[idx], 
            idlist))
    }
    return(labround)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
