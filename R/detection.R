
# TESTDATA:
## We have test.gtf (GRanges object) that we can use to test the funcs


#' @importFrom dplyr %>%
#TODO: Need a better name for this function
detection <- function(gtf){
    
    # TODO: Check inputs
    ## Can be path to gtf file or a GenomicRanges object
  if(is_valid_file(gtf)){
    gtf <- rtracklayer::import(gtf)
  } else {
    rlang::abort("GTF file does not exist")
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
    filtered_genes <- gtf %>%
      as.data.frame() %>%
      dplyr::filter(type=="exon") %>%
      dplyr::group_by(gene_id, transcript_id) %>%
      dplyr::tally() %>%
      dplyr::filter(n > 2) %>%
      dplyr::group_by(gene_id) %>%
      dplyr::tally() %>%
      dplyr::filter(n > 1)
    gtf <- gtf[gtf$gene_id %in% filtered_genes$gene_id]
    
    
    # get only exon entries from prefiltered GTF
    exons <- gtf[gtf$type == "exon"]
    
    # trim ends of transcripts to avoid TS and TE
    exons <- .trim_ends_by_gene(exons) 
    
    # classify exons by position (First, Internal, Last)
    exons <- .label_exon_class(exons)
    
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
    exon.juncs <- .get_juncs(disjoint.exons, introns.nr)
    
    
    # TODO: Get junctions for Retained introns
    #generate dataframe of retained introns
    retained.introns <- .find_retained_intron(disjoint.exons, introns.nr)
    
    
    # TODO: Classify events
    exon.classified <- .event_classification(exon.junction)
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

#group by transcript id and label first and last exons
.label_exon_class <- function(x){
  y <- x %>% 
    as.data.frame() %>% 
    dplyr::group_by(transcript_id) %>% 
    dplyr::mutate(exon_pos = dplyr::case_when(
      start==min(start) & strand == "+" ~ "first",
      end==max(end) & strand == "+" ~ "last",
      start==min(start) & strand == "-" ~ "last",
      end==max(end) & strand == "-" ~ "first",
      .default = "internal"
    ))
  x$exon_pos <- y$exon_pos
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
        dplyr::select(order = first.order, 
                      gene_name = second.gene_name, 
                      transcript_id = second.transcript_id,
                      exon_pos = second.exon_pos) %>% 
        dplyr::mutate(exon_pos = factor(exon_pos, c("internal", "first","last"))) %>%
        dplyr::group_by(order) %>% 
        dplyr::arrange(exon_pos) %>%
        dplyr::summarise(gene_name = gene_name[1], 
                         transcript_id = paste0(transcript_id, collapse = ";"),
                         exon_pos = exon_pos[1])
    y$gene_name <- out$gene_name
    y$transcript_ids <- out$transcript_id
    y$exon_pos <- as.character(out$exon_pos)
    y$order <- NULL
    
    return(y)
}

.get_juncs <- function(x, y){
  exon.spljunc <- .get_spljunc(x,y)
  out <- .add_skipjunc(exon.spljunc, y)
  
  return(out)
}


.get_spljunc <- function(x, y){
    x$index <- 1:length(x)
  overlap <- IRanges::findOverlapPairs(x, y, maxgap = 0L)
  adjacent <- subset(overlap, 
                     IRanges::start(first) == IRanges::end(second) + 1L | 
                     IRanges::end(first) == IRanges::start(second) - 1L)
  
  pair_df <- as.data.frame(adjacent) %>% 
    dplyr::mutate(position = ifelse(first.X.start == second.end + 1, "Upstream", "Downstream")) 
  GenomicRanges::mcols(adjacent)$position <- pair_df$position
  
  return(adjacent)
}        
   
.add_skipjunc <- function(x, y){
  overlap <- IRanges::findOverlapPairs(x, y, type = "within")
  
  x.overlap <- S4Vectors::first(overlap)
  y.overlap <- S4Vectors::second(overlap)
  
  exon.w.skipjunc <- as.data.frame(overlap) %>% 
    dplyr::mutate(exon_coord = .get_coord(S4Vectors::first(x.overlap)),
                  junc_coord = .get_coord(y.overlap),
                  junc_type = "Skipped") %>% 
    dplyr::distinct(exon_coord,junc_coord, .keep_all = TRUE) %>%  
    dplyr::select(exon_coord, junc_coord, 
                  gene_id = first.first.gene_id,
                  gene_name = first.first.gene_name,
                  transcript_ids = first.first.transcript_ids,
                  strand = first.first.X.strand, 
                  junc_type,
                  exon_pos=first.first.exon_pos)


  
  
  exon.w.spljunc <- as.data.frame(overlap) %>% 
    dplyr::mutate(exon_coord = .get_coord(S4Vectors::first(x.overlap)),
                  junc_coord = .get_coord(S4Vectors::second(x.overlap))) %>% 
    dplyr::distinct(exon_coord,junc_coord, .keep_all = TRUE) %>% 
    dplyr::select(exon_coord,junc_coord, 
                  gene_id = first.first.gene_id,
                  gene_name = first.first.gene_name,
                  transcript_ids = first.first.transcript_ids,
                  strand=first.first.X.strand,
                  junc_type=first.position,
                  exon_pos=first.first.exon_pos)
  
  return(rbind(exon.w.spljunc, exon.w.skipjunc))
 
} 

.find_retained_intron <- function(x, y){
  retained_intron <- IRanges::findOverlapPairs(x, y, type = "equal")
  
  ri_skipped <- as.data.frame(retained_intron) %>% 
    dplyr::mutate(exon_coord = .get_coord(S4Vectors::first(retained_intron)), 
                  junc_coord = exon_coord,
                  junc_type = "Skipped") %>% 
    dplyr::select(exon_coord,junc_coord, 
                  gene_id = first.gene_id,
                  gene_name = first.gene_name,
                  transcript_ids = first.transcript_ids,
                  strand=first.X.strand,
                  junc_type,
                  exon_pos=first.exon_pos)
  
  ri_resized_up <- GenomicRanges::resize(S4Vectors::first(retained_intron), 1)
  ri_resized_dn <- GenomicRanges::resize(S4Vectors::first(retained_intron), 1, 
                                         fix = "end")
  
  
  ri_spliced <- as.data.frame(retained_intron) %>% 
    dplyr::mutate(exon_coord = .get_coord(S4Vectors::first(retained_intron)), 
                  Upstream = .get_coord(ri_resized_up),
                  Downstream = .get_coord(ri_resized_dn)) %>% 
    dplyr::select(exon_coord,Upstream, Downstream,
                  gene_id = first.gene_id,
                  gene_name = first.gene_name,
                  transcript_ids = first.transcript_ids,
                  strand=first.X.strand,
                  exon_pos=first.exon_pos) %>%
    tidyr::pivot_longer(Upstream:Downstream, names_to = "junc_type", 
                        values_to = "junc_coord")
  
  
  ri_df <- rbind(ri_skipped,ri_spliced)
  ri_df$type <- "RI"
  
  return(ri_df)
}

.event_classification <- function(x,y){
  pivot <- tidyr::pivot_wider(x, names_from = position, values_from = score, values_fill = 0)
  
  classified <- pivot %>% 
    dplyr::mutate(splice_type = ifelse(exon_class=="first", "AF", 
                         ifelse(exon_class=="last", "AL",
                         ifelse(Skipped==1 & Downstream==1 & Upstream==1, "CE",
                         ifelse(strand=="+" & Skipped==1 & Downstream==1, "AD", "AA")))))
  
  return()
}


