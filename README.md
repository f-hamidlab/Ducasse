# **Ducasse**

## Detection and classification of alternative splicing events in R

**Ducasse** is an R package that detects alternatively-
spliced events from a transcriptome annotation file (GTF) and classifies them into
the major splicing categories (Casette Exons [CE], Alternative Donor [AD],
Alternative Acceptor [AA], Alternative First Exon [AF], Alternative Last Exon [AL]
and Retained Intron [RI]). This tool outputs exon coordinates of all alternatively-
spliced events as well as the intron coordinates that splices or skips the
spliced events. The latter data facilitates the quantification of exon
inclusion levels using the Percent-Spliced In (PSI) method.
  
## How to install
The development version can be installed using devtools:
```r
# install.packages("devtools")
devtools::install_github("f-hamidlab/Ducasse")
```

## How to use

The easiest way to use Ducasse is to supply the path to a GTF transcriptome
file to the `findASevents` function:

```r
output <- Ducasse::findASevents("path/to/gtf")
```
