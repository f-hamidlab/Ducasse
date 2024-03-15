gr1 <- GenomicRanges::GRanges(
    seqnames = "1", strand = rep("+", 8),
    ranges = IRanges::IRanges(
        start = c(10,80,140,200,300,400,500,620),
        end = c(30,100,160,240,340,440,540,650)
    ),
    transcript_id = "transcript1",
    gene_id = "geneA",
    gene_name = "gene001"
)

gr2 <- GenomicRanges::GRanges(
    seqnames = "1", strand = rep("+", 7),
    ranges = IRanges::IRanges(
        start = c(1,80,140,200,400,500,620),
        end = c(30,100,160,240,440,540,660)
    ),
    transcript_id = "transcript2",
    gene_id = "geneA",
    gene_name = "gene001"
)

gr3 <- GenomicRanges::GRanges(
    seqnames = "1", strand = rep("+", 8),
    ranges = IRanges::IRanges(
        start = c(10,80,140,180,300,400,500,620),
        end = c(30,100,160,240,340,460,540,650)
    ),
    transcript_id = "transcript3",
    gene_id = "geneA",
    gene_name = "gene001"
)

gr4 <- GenomicRanges::GRanges(
    seqnames = "1", strand = rep("+", 7),
    ranges = IRanges::IRanges(
        start = c(40,80,200,280,400,500,560),
        end = c(60,160,240,340,460,540,580)
    ),
    transcript_id = "transcript4",
    gene_id = "geneA",
    gene_name = "gene001"
)

test.gtf <- c(gr1,gr2,gr3,gr4)
test.gtf$type <- "exon"

txs <- range(S4Vectors::split(test.gtf, ~transcript_id)) %>% 
    as.data.frame() %>% 
    dplyr::mutate(transcript_id = group_name, gene_id = "geneA", 
                  gene_name = "gene001", type = "transcript") %>% 
    dplyr::select(seqnames:type) %>% 
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T)

gn <- range(test.gtf) %>% 
    as.data.frame() %>% 
    dplyr::mutate( gene_id = "geneA", 
                  gene_name = "gene001", type = "gene") %>% 
    dplyr::select(seqnames:type) %>% 
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T)


test.gtf <- c(gn, txs, test.gtf)

usethis::use_data(test.gtf, overwrite = T)