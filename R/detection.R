
# TESTDATA:
## TODO: Remove these later
gtf.file <-  system.file("extdata/pb_custom.gtf.gz", package = "factR2")
gtf <- rtracklayer::import(gtf.file)


#' @importFrom dplyr %>%
#TODO: Need a better name for this function
detection <- function(gtf){
    
    # Check inputs
    ## Can be path to gtf file or a GenomicRanges object
    
    ## prefilter for genes with at least 2 multi-exonic transcripts
    
    # get only exon entries from GTF
    exons <- gtf[gtf$type == "exon"]
    
    # Create a GenomicRanges object of all non-redundant introns
    exonsbytx <- S4Vectors::split(exons, ~transcript_id)
    intronsbytx <- GenomicRanges::psetdiff(BiocGenerics::unlist(range(exonsbytx)), exonsbytx)
    introns.nr <- unique(unlist(intronsbytx))
    names(introns.nr) <- NULL
    
    # Create a disjointed version of all exons in each gene family
    disjoint.exons <- .disjoin_by_gene(exons)

    
    # Pair up all exons with flanking introns
    
    # Pair up all exons with "skipping" introns
    
    # Output 1) metadata of all exons, 2) adjacency matrix of exons and flanking introns
    # 3) adjacency matrix of exons and skipped introns
    return(disjoint.exons)
    
}




.disjoin_by_gene <- function(x){
    y <- GenomicRanges::disjoin(S4Vectors::split(x, ~gene_id))
    y <- y %>% 
        as.data.frame() %>% 
        dplyr::mutate(gene_id = group_name) %>% 
        dplyr::select(seqnames:gene_id) %>% 
        GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)
    
    # TODO: add other metadata such as gene names, transcript_ids
    
    return(y)
    
    
}