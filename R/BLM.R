#' Brain/Liver/Muscle (BLM) RNA-sequencing mixture data
#'
#' A dataset containing a set of mixes andraw components.
#'
#' @format Data frame of 16 variables and 28030 rows:
#' \describe{
#'	\item{gene_id}{UCSC gene identifier}
#'	\item{a1}{counts from mix 1, ambion ercc spike pool 'a' 4x - designed proportions 0.25,0.25,0.5}
#'	\item{a1d}{mix 1, ambion ercc spike pool 'a' 1x - designed proportions 0.25,0.25,0.5}
#'	\item{a1u}{mix 1, ambion ercc spike pool 'a' 16x - designed proportions 0.25,0.25,0.5}
#'	\item{b1}{mix 1, ambion ercc spike pool 'b' 4x - designed proportions 0.25,0.25,0.5}
#'	\item{b1d}{mix 1, ambion ercc spike pool 'b' 1x - designed proportions 0.25,0.25,0.5}
#'	\item{b1u}{mix 1, ambion ercc spike pool 'b' 16x - designed proportions 0.25,0.25,0.5}
#'	\item{a2}{mix 2, ambion ercc spike pool 'a' 4x - designed proportions 0.25,0.5,0.25}
#'	\item{b2}{mix 2, ambion ercc spike pool 'b' 4x - designed proportions 0.25,0.5,0.25}
#'	\item{b2d}{mix 2, ambion ercc spike pool 'b' 1x - designed proportions 0.25,0.5,0.25}
#'	\item{b2u}{mix 2, ambion ercc spike pool 'b' 16x - designed proportions 0.25,0.5,0.25}
#'	\item{bep}{Pure brain counts}
#'	\item{lep}{Pure liver counts}
#'	\item{mep}{Pure muscle counts}
#'	\item{norm}{string identifying normalization type:  UQN = upper quartile normalization}
#'	\item{uid}{string identifying the data source:  Illumina sequencing, tophat mapping,
#'	 HTSeq-counts quantification, USCS annotation}
#' }
#' @source \url{http://dx.doi.org/10.1186\%2Fs12864-015-1912-7}
"BLM"