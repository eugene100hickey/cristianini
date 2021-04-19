library(rentrez)
library(ape)

m_genitalium <- "NC_000908"

z <- read.GenBank(m_genitalium, as.character = T)
z1 <- Biostrings::DNAString(z[[1]] %>% str_c(collapse = ""))

ntd <- "H2N3"
h2n3_search <- entrez_search(db="nuccore", term=ntd, retmax=40) 
h2n3_seqs <- entrez_fetch(db="nuccore", 
                                 id=h2n3_search$ids[1], 
                                 rettype="fasta")
