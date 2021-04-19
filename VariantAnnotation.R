# https://bioconductor.org/packages/release/bioc/vignettes/VariantAnnotation/inst/doc/VariantAnnotation.pdf

library(VariantAnnotation)
library(tidyverse)
library(ape)
library(Biostrings)
library(ggseqlogo)

fl <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
vcf <- readVcf(fl, "hg19")

head(rowRanges(vcf), 3)
ref(vcf)[1:5]
alt(vcf)[1:5]
geno(vcf)$DS
