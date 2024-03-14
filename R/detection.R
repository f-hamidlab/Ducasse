
# TESTDATA:
gtf.file <-  system.file("extdata/pb_custom.gtf.gz", package = "factR2")
gtf <- rtracklayer::import(gtf.file)



#TODO: Need a better name for this function
detection <- function(gtf){
    
    # Check inputs
    ## Can be path to gtf file or a GenomicRanges object
    
    # get only exon entries from GTF
    exons <- gtf[gtf$type == "exon"]
    
    # Create a GenomicRanges object of all non-redundant introns
    
    # Create a disjointed version of all exons
    
    # Pair up all exons with flanking introns
    
    # Pair up all exons with "skipping" introns
    
    # Output 1) metadata of all exons, 2) adjacency matrix of exons and flanking introns
    # 3) adjacency matrix of exons and skipped introns
    
}