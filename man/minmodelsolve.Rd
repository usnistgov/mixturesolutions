\name{minmodelsolve}
\alias{minmodelsolve}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
minmodelsolve(ssx = indf, preps)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ssx}{
%%     ~~Describe \code{ssx} here~~
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
function (ssx = indf, preps) 
{
    modelproportions <- NULL
    for (idx in 1:length(preps$mixnames)) {
        form <- paste0(preps$mixnames[idx], "~")
        for (cindex in (1:(length(preps$componentnames)))) {
            form <- paste0(form, "I(", preps$componentnames[cindex], 
                "*", preps$mfrac[cindex], ")+")
        }
        form <- paste0(form, 0)
        modelproportions <- rbind(modelproportions, c(coefficients(lm(data = ssx, 
            form))/sum(coefficients(lm(data = ssx, form))), preps$mixnames[idx]))
    }
    colnames(modelproportions) <- c(preps$componentnames, "mix")
    modelproportions <- as.data.frame(modelproportions)
    for (I in 1:length(preps$componentnames)) {
        modelproportions[, I] <- as.numeric(as.character(modelproportions[, 
            I]))
    }
    return(modelproportions)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
