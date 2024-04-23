#' Mock GenomicRanges object in GTF format
#'
#' An object cotaining for 4 hypothetical transcripts from "geneA" with different
#' exon structures
#'
#' @format A GRanges object with 35 ranges and 4 metadata columns:
#' \describe{
#'   \item{ranges}{Chromosome, start, end, and strand info of 4 transcripts
#'   and its exons}
#'   \item{type}{Entry type; transcript or exon}
#'   \item{transcript_id}{Name or ID given to transcripts}
#'   \item{gene_id}{Name or ID given to gene origin of transcripts}
#'   ...
#' }
#' @usage data(test.gtf)
"test.gtf"
