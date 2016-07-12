#' SEQC dataset (Australian Genome Research Facility)
#'
#' A dataset containing a set of mixes andraw components.
#'
#' @format Data frame with 257 columns and 43919 rows:
#' \describe{
#'	\item{Transcript_id}{RefSeq transcript identifier}
#'	\item{Other} {256 columns of data with the following format:  SEQC_ILM_AGR_sample_replicate_lane_barcode_RunID}
#'	\item{Sample} {Sample information:  SEQC sample A is ambion UHRR, B is ambion HBRR, C is a 75/25 mix of A+B, D is a 25/75 mix of A+B}
#'	}
#' @source GEO.GSE47774}
"SEQCagr"