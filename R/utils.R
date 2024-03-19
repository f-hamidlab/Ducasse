.get_coord <- function(x){
  paste0(GenomicRanges::seqnames(x), ":",
         GenomicRanges::start(x), "-",
         GenomicRanges::end(x))
}

.get_coord_strand <- function(x, sep = ";"){
  paste0(GenomicRanges::seqnames(x), ":",
         GenomicRanges::start(x), "-",
         GenomicRanges::end(x),
         sep, GenomicRanges::strand(x))
}

.get_coord_geneid_strand <- function(x, sep = ";"){
  paste0(GenomicRanges::seqnames(x), ":",
         GenomicRanges::start(x), "-",
         GenomicRanges::end(x),
         sep, x$gene_id,
         sep, GenomicRanges::strand(x))
}

is_gtf <- function(x) {
  if (is(x, "GRanges")) {
    x <- as.data.frame(x)
    if (all(c("type", "gene_id", "transcript_id") %in% names(x))) {
      x <- x %>%
        dplyr::select(type, gene_id, transcript_id) %>%
        dplyr::distinct()
      if (nrow(x) > 1) {
        return(TRUE)
      }
    }
  }
  return(FALSE)
}
