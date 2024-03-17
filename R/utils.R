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