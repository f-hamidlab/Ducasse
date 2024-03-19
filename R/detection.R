
# TESTDATA:
## We have test.gtf (GRanges object) that we can use to test the funcs


#' @importFrom dplyr %>%
#TODO: Need a better name for this function
detection <- function(gtf){
    
    # TODO: Check inputs
    ## Can be path to gtf file or a GenomicRanges object
    if(class(gtf) %in% "character"){
      if(file.exists(gtf)){
        gtf <- rtracklayer::import(gtf)
      } else {
        rlang::abort("GTF file does not exist")
      }
    }
  
    # check GTF structure
    if(!is_gtf(gtf)){
      rlang::abort("Input is not a GTF file structure")
    } else{
      # check if gene_name column is present and if not, use "gene_id" column
      if(!"gene_name" %in% colnames(S4Vectors::mcols(gtf))){
        gtf$gene_name <- gtf$gene_id
      }
    }
  
    # TODO: Prefilter for genes with at least 2 multi-exonic transcripts
    transcript_counts <- table(GenomicRanges::mcols(gtf)$transcript_id)
    transcripts <- gtf[gtf$type == "transcript"]
    multi_transcripts <- transcript_counts[transcript_counts > 2]
    transcripts_filtered <- transcripts[transcripts$transcript_id %in% names(multi_transcripts)]
    gene_ids <- table(GenomicRanges::mcols(transcripts_filtered)$gene_id)
    multi_exonic_genes <- names(gene_ids[gene_ids > 1])
    gtf <- gtf[gtf$gene_id %in% multi_exonic_genes]
    
    
    # get only exon entries from prefiltered GTF
    exons <- gtf[gtf$type == "exon"]
    
    # trim ends of transcripts to avoid TS and TE
    exons <- .trim_ends_by_gene(exons)
    
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
    #Find exons that are within intron coordinates
    exon.junction <- .group_by_junction(exon.intron.pairs, introns.nr)
    
    # TODO: Get junctions for Retained introns
    #append the retained introns into previous dataframe
    exon.events <- .append_retained_intron(disjoint.exons, introns.nr, exon.junction)
    
    
    # TODO: Classify events
    exon.classified <- .event_classification(exon.events)
    # TODO:  Output
    ## 1) metadata of all exons, 2) adjacency matrix of exons and flanking introns
    ## 3) adjacency matrix of exons and skipped introns
    
    return()
    
}


.trim_ends_by_gene <- function(x){
  # trim-off TS and TE
  x %>%
    as.data.frame() %>%
    dplyr::group_by(transcript_id) %>%
    dplyr::arrange(start) %>%
    dplyr::mutate(pos = dplyr::row_number()) %>%
    dplyr::mutate(pos = dplyr::case_when(pos == 1 ~ "First",
                                         pos == dplyr::n() ~ "Last",
                                         .default = "a.internal")) %>%
    dplyr::group_by(seqnames, end, gene_id) %>%
    dplyr::arrange(pos, dplyr::desc(start)) %>%
    dplyr::mutate(start = ifelse(pos == "First", start[1], start)) %>%
    dplyr::group_by(seqnames, start, gene_id) %>%
    dplyr::arrange(dplyr::desc(pos), dplyr::desc(end)) %>%
    dplyr::mutate(end = ifelse(pos == "Last", end[dplyr::n()], end)) %>%
    dplyr::ungroup() %>%
    dplyr::select(-pos) %>%
    dplyr::arrange(seqnames, gene_id, transcript_id, start) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T)
  
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
   
.group_by_junction <- function(x, y){
  overlap <- IRanges::findOverlapPairs(x, y, type = "within")
  overlap_df <- as.data.frame(overlap) %>% 
    dplyr::mutate(exon = paste0(first.first.X.seqnames, "_", first.first.X.start, ":", first.first.X.end),
           flankingIntron = paste0(first.second.seqnames,"_",first.second.start,":",first.second.end),
           exonInIntron = paste0(second.seqnames,"_",second.start,":",second.end)) %>% 
    dplyr::select(exon, exonInIntron, flankingIntron, position = first.position, strand = first.first.X.strand)
  
  skipped <- overlap_df %>% 
    dplyr::group_by(exon, exonInIntron) %>% 
    dplyr::filter(dplyr::n()>1) %>% 
    dplyr::mutate(position = "Skipped") %>% 
    dplyr::select(exon, spliceJunction = exonInIntron, position, strand) %>% 
    dplyr::distinct()
  
  upstream_downstream <- overlap_df %>% 
    dplyr::group_by(exon, exonInIntron) %>% 
    dplyr::filter(dplyr::n() == 1) %>% 
    dplyr::ungroup() %>% 
    dplyr::select(exon, spliceJunction = flankingIntron, position, strand)
  
  junction_grp <- rbind(skipped, upstream_downstream)
  
  return(junction_grp)
} 

.append_retained_intron <- function(x, y, z){
  retained_intron <- IRanges::findOverlapPairs(x, y, type = "equal")
  
  ri_skipped <- as.data.frame(retained_intron) %>% 
    dplyr::mutate(exon = paste0(first.X.seqnames, "_", first.X.start, ":", first.X.end),
                  spliceJunction = paste0(second.seqnames,"_",second.start,":",second.end),
                  position = "Skipped",
                  type = "RI") %>% 
    dplyr::select(exon, spliceJunction, position, strand = first.X.strand, type)
  
  ri_upstream <- as.data.frame(retained_intron) %>% 
    dplyr::mutate(exon = paste0(first.X.seqnames, "_", first.X.start, ":", first.X.end),
                  spliceJunction = paste0(second.seqnames,"_",second.start-2,":",second.start-1),
                  position = "Upstream",
                  type = "RI") %>% 
    dplyr::select(exon, spliceJunction, position, strand = first.X.strand, type)
  
  ri_downstream <- as.data.frame(retained_intron) %>% 
    dplyr::mutate(exon = paste0(first.X.seqnames, "_", first.X.start, ":", first.X.end),
                  spliceJunction = paste0(second.seqnames,"_",second.end+1,":",second.end+2),
                  position = "Downstream",
                  type = "RI") %>% 
    dplyr::select(exon, spliceJunction, position, strand = first.X.strand, type)
  
  ri_df <- rbind(ri_skipped,ri_upstream,ri_downstream)
  
  appended <- rbind(z, ri_df)
  
  return(appended)
}

.event_classification <- function(x){
  x$type <- with(x, ifelse(is.na(type) & position == "Upstream", "ALT3",
                           ifelse(is.na(type) & position == "Downstream", "ALT5",
                                  ifelse(is.na(type) & position == "Skipped", "SE", "RI"))))
  
  return(x)
}

}