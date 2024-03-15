
# TESTDATA:
## We have test.gtf (GRanges object) that we can use to test the funcs


#' @importFrom dplyr %>%
#TODO: Need a better name for this function
detection <- function(gtf){
    
    # TODO: Check inputs
    ## Can be path to gtf file or a GenomicRanges object
    
    
    
    # TODO: Prefilter for genes with at least 2 multi-exonic transcripts
    
    
    
    
    # get only exon entries from GTF
    exons <- gtf[gtf$type == "exon"]
    
    # Create a GenomicRanges object of all non-redundant introns
    exonsbytx <- S4Vectors::split(exons, ~transcript_id)
    intronsbytx <- GenomicRanges::psetdiff(BiocGenerics::unlist(range(exonsbytx)), exonsbytx)
    introns.nr <- unique(unlist(intronsbytx))
    names(introns.nr) <- NULL
    
    # Create a disjointed version of all exons in each gene family
    disjoint.exons <- .disjoin_by_gene(exons)

    
    # TODO:  Pair up all exons with flanking introns
    
    # TODO:  Pair up all exons with "skipping" introns
    
    # TODO:  Output
    ## 1) metadata of all exons, 2) adjacency matrix of exons and flanking introns
    ## 3) adjacency matrix of exons and skipped introns
    return()
    
}




.disjoin_by_gene <- function(x){
    y <- GenomicRanges::disjoin(S4Vectors::split(x, ~gene_id))
    y <- y %>% 
        as.data.frame() %>% 
        dplyr::mutate(gene_id = group_name) %>% 
        dplyr::select(seqnames:gene_id) %>% 
        GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)
    
    y$order <- 1:length(y)
    out <- IRanges::findOverlapPairs(y, x, type="within")
    out <- out[S4Vectors::first(out)$gene_id == S4Vectors::second(out)$gene_id ]
    
    out <- out %>% 
        as.data.frame() %>% 
        dplyr::select(order = first.order, gene_name = second.gene_name, 
                      transcript_id = second.transcript_id) %>% 
        dplyr::group_by(order) %>% 
        dplyr::summarise(gene_name = gene_name[1], 
                         transcript_id = paste0(transcript_id, collapse = ";"))
    y$gene_name <- out$gene_name
    y$transcript_ids <- out$transcript_id
    y$order <- NULL
    
    return(y)
    
    
}