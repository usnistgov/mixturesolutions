\name{unittests}
\alias{unittests}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
unittests()
}
%- maybe also 'usage' for other objects documented here.
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
function () 
{
    Blmdata = read.csv("Example_2mix_3compBLM.txt")
    preptargetplot(mrnatype = "ercc", prenormalized = TRUE, modeltype = "twomix", 
        indf = Blmdata, trueproportions = data.frame(mix1 = c(0.25, 
            0.25, 0.5), mix2 = c(0.25, 0.5, 0.25)), componentnames = c("bep", 
            "lep", "mep"), mixnames = c("a1", "a2"))
    preptargetplot(mrnatype = "internalconsensus", prenormalized = TRUE, 
        modeltype = "twomix", indf = Blmdata, trueproportions = data.frame(mix1 = c(0.25, 
            0.25, 0.5), mix2 = c(0.25, 0.5, 0.25)), componentnames = c("bep", 
            "lep", "mep"), mixnames = c("a1", "a2"))
    mirnadata = read.csv("Example_2mix_3comp.txt", sep = "\t")
    preptargetplot(indf = mirnadata, mrnatype = "internalconsensus", 
        prenormalized = FALSE, modeltype = "twomix", componentnames = c("Brain", 
            "Liver", "Placenta"), mixnames = c("Mix1", "Mix2"), 
        trueproportions = data.frame(mix1 = c(0.25, 0.25, 0.5), 
            mix2 = c(0.5, 0.25, 0.25)))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
