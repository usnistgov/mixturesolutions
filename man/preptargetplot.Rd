\name{preptargetplot}
\alias{preptargetplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
preptargetplot(mrnatype = "internalconsensus", splittype = "none", prenormalized = TRUE, modeltype = "twomix", indf, retdf = FALSE, idnames = "id", idvars = 1, trueproportions, componentnames = c("Brain", "Liver", "Placenta"), mixnames = c("Mix1", "Mix2"), annot = c("Replicate"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mrnatype}{
%%     ~~Describe \code{mrnatype} here~~
}
  \item{splittype}{
%%     ~~Describe \code{splittype} here~~
}
  \item{prenormalized}{
%%     ~~Describe \code{prenormalized} here~~
}
  \item{modeltype}{
%%     ~~Describe \code{modeltype} here~~
}
  \item{indf}{
%%     ~~Describe \code{indf} here~~
}
  \item{retdf}{
%%     ~~Describe \code{retdf} here~~
}
  \item{idnames}{
%%     ~~Describe \code{idnames} here~~
}
  \item{idvars}{
%%     ~~Describe \code{idvars} here~~
}
  \item{trueproportions}{
%%     ~~Describe \code{trueproportions} here~~
}
  \item{componentnames}{
%%     ~~Describe \code{componentnames} here~~
}
  \item{mixnames}{
%%     ~~Describe \code{mixnames} here~~
}
  \item{annot}{
%%     ~~Describe \code{annot} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
function (mrnatype = "internalconsensus", splittype = "none", 
    prenormalized = TRUE, modeltype = "twomix", indf, retdf = FALSE, 
    idnames = "id", idvars = 1, trueproportions, componentnames = c("Brain", 
        "Liver", "Placenta"), mixnames = c("Mix1", "Mix2"), annot = c("Replicate"), 
    ...) 
{
    normalizedf <- function(indf, idvars) {
        require(edgeR)
        indf <- as.data.frame(indf)
        indf[is.na(indf)] <- 0
        indf <- indf[rowSums(indf[-idvars]) > 0, ]
        indf[-idvars] <- round(indf[-idvars])
        normfac <- calcNormFactors(indf[, -idvars], method = "upperquartile")
        normdf <- indf[-idvars]
        IT <- 0
        for (I in normfac) {
            IT <- IT + 1
            normdf[, IT] <- normdf[, IT] * I
        }
        indf[-idvars] <- normdf
        return(indf)
    }
    preps <- list(modeltype = modeltype, mrnatype = mrnatype, 
        trueproportions = trueproportions, mixnames = mixnames, 
        componentnames = componentnames, idvars = idvars, idnames = idnames)
    getmfrac <- function(indf, type = preps$mrnatype) {
        mfrac <- switch(type, internalconsensus = doubleconsensus(indf, 
            preps), externalagreement = lookupmfrac(mixtureid), 
            none = c(1, 1, 1), ercc = calcmrnafracgeneral(indf[c(1, 
                selectcomponents(indf, preps))]))
        mfrac <- mfrac/sum(mfrac)
        return(mfrac)
    }
    if (prenormalized == FALSE) {
        indf <- normalizedf(indf, idvars)
    }
    collapsereps <- function(indf, preps, annot) {
        attach(preps)
        reslist <- NULL
        for (I in 1:length(componentnames)) {
            if (length(grep(componentnames[I], colnames(indf))) == 
                0) {
                warning(call. = TRUE, ... = paste0("One of your componentnames, ", 
                  componentnames[I], " couldn't be found in the data:  Do colnames match componentnames?"))
            }
            reslist <- c(reslist, grep(componentnames[I], colnames(indf)))
        }
        for (I in 1:length(mixnames)) {
            if (length(grep(mixnames[I], colnames(indf))) == 
                0) {
                warning(call. = TRUE, ... = paste0("One of your mixnames, ", 
                  mixnames[I], " couldn't be found in the data:  Do colnames match mixnames?"))
            }
            reslist <- c(reslist, grep(mixnames[I], colnames(indf)))
        }
        if (length(unique(reslist)) != length(reslist)) {
            warning("One or more entities appear to have non-unique names.  Can you make sure all of your components and mixes have unique names?")
        }
        C1 <- indf[, c(grep(componentnames[1], colnames(indf)))]
        C1l <- NULL
        for (I in 1:ncol(C1)) {
            C1l <- rbind(C1l, data.frame(id = indf[, idvars], 
                data = C1[, I], annot = I))
        }
        colnames(C1l) <- c(idnames, "data", annot)
        C1l$Source <- componentnames[1]
        if (length(componentnames) > 1) {
            C2 <- indf[, c(grep(componentnames[2], colnames(indf)))]
            C2l <- NULL
            for (I in 1:ncol(C2)) {
                C2l <- rbind(C2l, data.frame(id = indf[, idvars], 
                  data = C2[, I], annot = I))
            }
            colnames(C2l) <- c(idnames, "data", annot)
            C2l$Source <- componentnames[2]
        }
        if (length(componentnames) > 2) {
            C3 <- indf[, c(grep(componentnames[3], colnames(indf)))]
            C3l <- NULL
            for (I in 1:ncol(C3)) {
                C3l <- rbind(C3l, data.frame(id = indf[, idvars], 
                  data = C3[, I], annot = I))
            }
            colnames(C3l) <- c(idnames, "data", annot)
            C3l$Source <- componentnames[3]
        }
        if (length(componentnames) > 3) {
            warning("You seem to have input more than 3 component names.  This behavior is not currently supported.")
        }
        if (length(componentnames) > 3) {
            C4 <- indf[, c(grep(componentnames[4], colnames(indf)))]
            C4l <- NULL
            for (I in 1:ncol(C4)) {
                C4l <- rbind(C4l, data.frame(id = indf[, idvars], 
                  data = C4[, I], annot = I))
            }
            colnames(C4l) <- c(idnames, "data", annot)
            C4l$Source <- componentnames[4]
        }
        if (length(componentnames) > 4) {
            C5 <- indf[, c(grep(componentnames[5], colnames(indf)))]
            C5l <- NULL
            for (I in 1:ncol(C5)) {
                C5l <- rbind(C5l, data.frame(id = indf[, idvars], 
                  data = C5[, I], annot = I))
            }
            colnames(C5l) <- c(idnames, "data", annot)
            C5l$Source <- componentnames[5]
        }
        M1 <- indf[, c(grep(mixnames[1], colnames(indf)))]
        M1l <- NULL
        for (I in 1:ncol(M1)) {
            M1l <- rbind(M1l, data.frame(id = indf[, idvars], 
                data = M1[, I], annot = I))
        }
        colnames(M1l) <- c(idnames, "data", annot)
        M1l$Source <- mixnames[1]
        if (length(mixnames) > 1) {
            M2 <- indf[, c(grep(mixnames[2], colnames(indf)))]
            M2l <- NULL
            for (I in 1:ncol(M2)) {
                M2l <- rbind(M2l, data.frame(id = indf[, idvars], 
                  data = M2[, I], annot = I))
            }
            colnames(M2l) <- c(idnames, "data", annot)
            M2l$Source <- mixnames[2]
        }
        if (length(mixnames) > 2) {
            warning("You have input more than 2 mix names.  This behavior is not currently supported")
        }
        if (length(mixnames) > 2) {
            M3 <- indf[, c(grep(mixnames[3], colnames(indf)))]
            M3l <- NULL
            for (I in 1:ncol(M3)) {
                M3l <- rbind(M3l, data.frame(id = indf[, idvars], 
                  data = M3[, I], annot = I))
            }
            colnames(M3l) <- c(idnames, "data", annot)
            M3l$Source <- mixnames[3]
        }
        if (length(mixnames) > 3) {
            M4 <- indf[, c(grep(mixnames[4], colnames(indf)))]
            M4l <- NULL
            for (I in 1:ncol(M4)) {
                M4l <- rbind(M4l, data.frame(id = indf[, idvars], 
                  data = M4[, I], annot = I))
            }
            colnames(M4l) <- c(idnames, "data", annot)
            M4l$Source <- mixnames[4]
        }
        listout <- NULL
        outdf <- NULL
        IT <- 0
        allnames <- c(componentnames, mixnames)
        J <- 0
        require(reshape2)
        detach(preps)
        listout <- c(grep("^C[0-9]l$", ls()), grep("^M[0-9]l$", 
            ls()))
        if (length(listout) != length(allnames)) {
            error("Something went horribly wrong in replicate collapsing.")
        }
        for (J in listout) {
            IT <- IT + 1
            outdf <- rbind(outdf, get(ls()[J]))
        }
        outdf <- dcast(outdf, get(idnames) ~ Source, fun.aggregate = mean, 
            na.rm = TRUE, value.var = "data")
        return(outdf)
    }
    if (length(grep(componentnames[1], colnames(indf))) > 1) {
        indf = collapsereps(indf, preps, annot)
    }
    mfrac <- getmfrac(indf, preps$mrnatype)
    preps <- c(preps, list(mfrac = mfrac))
    splitmultidf <- function(indf, splittype, splitcolumn, idvars) {
        switch(splittype, none = return(indf), globalcon = makelabrounda())
    }
    splitm1 <- splitmultidf(indf, splittype, splitcolumn)
    getm1 <- function(indf, preps) {
        switch(preps$modeltype, twomix = minmodelsolve(ssx = indf, 
            preps = preps), threemix = magicmodelsolve())
    }
    m1 <- getm1(indf, preps = preps)
    return(m1)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
