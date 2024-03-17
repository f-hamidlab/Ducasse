
# TESTDATA:
## We have test.gtf (GRanges object) that we can use to test the funcs


#' @importFrom dplyr %>%
#TODO: Need a better name for this function
detection <- function(gtf){
    
    # TODO: Check inputs
    ## Can be path to gtf file or a GenomicRanges object
    
    # get only exon entries from GTF
    exons <- gtf[gtf$type == "exon"]
    
    # TODO: Prefilter for genes with at least 2 multi-exonic transcripts
    transcript_counts <- table(GenomicRanges::mcols(exons)$transcript_id)
    transcripts <- gtf[gtf$type == "transcript"]
    multi_transcripts <- transcript_counts[transcript_counts > 2]
    transcripts_filtered <- transcripts[transcripts$transcript_id %in% names(multi_transcripts)]
    gene_ids <- table(GenomicRanges::mcols(transcripts_filtered)$gene_id)
    multi_exonic_genes <- names(gene_ids[gene_ids > 1])
    gtf <- gtf[gtf$gene_id %in% multi_exonic_genes]
    
    # trim ends of transcripts to avoid TS and TE
    
    
    # get only exon entries from prefiltered GTF
    exons <- gtf[gtf$type == "exon"]
    
    # Create a GenomicRanges object of all non-redundant introns
    exonsbytx <- S4Vectors::split(exons, ~transcript_id)
    intronsbytx <- GenomicRanges::psetdiff(BiocGenerics::unlist(range(exonsbytx)), exonsbytx)
    introns.nr <- unique(unlist(intronsbytx))
    names(introns.nr) <- NULL
    
    # Create a disjointed version of all exons in each gene family
    disjoint.exons <- .disjoin_by_gene(exons)

    
    # TODO:  Pair up all exons from disjoint.exons with flanking introns
    ## GenomicRanges::findOverlaps
    ## Output a df with these columns:
    #   1. Exon coordinate
    #   2. Intron coordinate for each hit
    #   3. Position (Upstream or downstream)
    exon.intron.pairs <- .pair_by_exon_intron(disjoint.exons, introns.nr)
    
    # TODO:  Pair up all exons with "skipping" introns
    
    # TODO:  Output
    ## 1) metadata of all exons, 2) adjacency matrix of exons and flanking introns
    ## 3) adjacency matrix of exons and skipped introns
    return()
    
}


.trim_ends_by_gene(x){
  # trim-off TS and TE
  x <- x %>%
    as.data.frame() %>%
    dplyr::group_by(transcript_id) %>%
    dplyr::arrange(start) %>%
    dplyr::mutate(pos = dplyr::row_number()) %>%
    dplyr::mutate(pos = dplyr::case_when(pos == 1 ~ "First",
                                         pos == dplyr::n() ~ "Last",
                                         .default = "a.internal")) %%
    dplyr::group_by(seqnames, end, gene_id) %>%
    dplyr::arrange(pos, dplyr::desc(start)) %>%
    dplyr::mutate(start = ifelse(pos == "First", start[1], start)) %>%
    dplyr::group_by(seqnames, start, gene_id) %>%
    dplyr::arrange(dplyr::desc(pos), dplyr::desc(end)) %>%
    dplyr::mutate(end = ifelse(pos == "Last", end[dplyr::n()], end)) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T)
  
  return(x)
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

.pair_by_exon_intron <- function(x, y){
    x$index <- 1:length(x)
  overlap <- IRanges::findOverlapPairs(x, y, maxgap = 0L)
  adjacent <- subset(overlap, IRanges::start(first) == IRanges::end(second) + 1L | 
                       IRanges::end(first) == IRanges::start(second) - 1L)
  
  pair_df <- as.data.frame(adjacent) %>% 
    dplyr::mutate(position = ifelse(first.X.start == second.end + 1, "Upstream", "Downstream")) %>% 
    dplyr::select(original.index=first.index,  position) %>% 
    dplyr::mutate(overlapped.index = dplyr::row_number())
  GenomicRanges::mcols(adjacent) <- pair_df
  
  # disjoint_df <- as.data.frame(x) %>% 
  #   dplyr::mutate(exonCoord = paste0(seqnames, "_", start, ":", end)) %>% 
  #   dplyr::select(exonCoord)
  # 
  # merge <- dplyr::full_join(x = disjoint_df, y = pair_df, "exonCoord")
  
  return(adjacent)
}        
    
}