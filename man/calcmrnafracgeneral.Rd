\name{calcmrnafracgeneral}
\alias{calcmrnafracgeneral}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
calcmrnafracgeneral(dat, spikeID = "ERCC-", spikemassfraction = 0.1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
%%     ~~Describe \code{dat} here~~
}
  \item{spikeID}{
%%     ~~Describe \code{spikeID} here~~
}
  \item{spikemassfraction}{
%%     ~~Describe \code{spikemassfraction} here~~
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
function (dat, spikeID = "ERCC-", spikemassfraction = 0.1) 
{
    countcolumns <- which(unname(unlist(lapply(dat, class)) == 
        "numeric"))
    annotcolumn <- which(unname(unlist(lapply(dat, class)) != 
        "numeric"))
    ercc <- rownames(dat)[which(substr(rownames(dat), 1, 5) == 
        spikeID)]
    if (length(ercc) == 0) {
        ercc <- grep(spikeID, dat[, annotcolumn[1]])
    }
    if (length(ercc) == 0) {
        stop("I can't identify the spike-ins within your count table.   The spikeID variable should be set to something which uniquely identifies spike-ins.   Rownames are first checked for names, then if there are non-numeric columns, only the FIRST is checked for gene names. ")
    }
    nonercc <- !(1:length(dat[, countcolumns[1]])) \%in\% ercc
    count <- rbind(colSums(dat[ercc, countcolumns]), colSums(dat[nonercc, 
        countcolumns]))
    ercc.targ <- spikemassfraction
    mRNA.frac <- ercc.targ * count[2, ]/count[1, ]
    return(mRNA.frac)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
