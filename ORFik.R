# https://bioconductor.org/packages/release/bioc/vignettes/ORFik/inst/doc/ORFikOverview.html

library(ORFik)                        # This package
library(GenomicFeatures)              # For basic transcript operations
library(data.table)                   # For fast table operations
library(BSgenome.Hsapiens.UCSC.hg19)  # Human genome

txdbFile <- system.file("extdata", "hg19_knownGene_sample.sqlite", 
                        package = "GenomicFeatures")

txdb <- loadTxdb(txdbFile)
fiveUTRs <- loadRegion(txdb, "leaders")
fiveUTRs[1]

# Extract sequences of fiveUTRs.
tx_seqs <- extractTranscriptSeqs(BSgenome.Hsapiens.UCSC.hg19::Hsapiens, 
                                 fiveUTRs) 
tx_seqs[1]

fiveUTR_ORFs <- findMapORFs(fiveUTRs, tx_seqs, groupByTx = FALSE)
fiveUTR_ORFs[1:2]

txNames(fiveUTR_ORFs[1:2]) # <- Which transcript

orf_seqs <- extractTranscriptSeqs(BSgenome.Hsapiens.UCSC.hg19::Hsapiens,
                                  fiveUTR_ORFs[1:3])
orf_seqs

orf_aa_seq <- Biostrings::translate(orf_seqs)
orf_aa_seq
