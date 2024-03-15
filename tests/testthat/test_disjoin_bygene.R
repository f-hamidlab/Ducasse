# setup GRanges 
gr2 <- test.gtf
gr2$gene_id <- "geneB"
gr2$transcript_id <- paste0(gr2$transcript_id, "_1")

comb.gr <- c(gr2, test.gtf)


test_that("Test disjoin function", {
    out <- .disjoin_by_gene(test.gtf[test.gtf$type == "exon"])
    expect_equal(length(out), 16)
    expect_equal(GenomicRanges::start(out),
                 c(1,10,40,80,101,140,180,200,280,300,400,441,500,560,620,651))
    
})

test_that("Test disjoin by gene", {
    out <- .disjoin_by_gene(comb.gr[comb.gr$type == "exon"])
    expect_equal(length(out), 32)
    expect_equal(GenomicRanges::start(out)[1:16],
                 c(1,10,40,80,101,140,180,200,280,300,400,441,500,560,620,651))
    expect_equal(GenomicRanges::start(out)[17:32],
                 c(1,10,40,80,101,140,180,200,280,300,400,441,500,560,620,651))
    
})