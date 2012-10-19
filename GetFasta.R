# Get fasta files from GenBank based on genebank ID
# e.g. viral genomes for RNAseq
#-----------------------------------------------------------------------------
library(stringr)
library(plyr)
library(Biostrings)
#-----------------------------------------------------------------------------
# viral sequences found here:
# http://viralzone.expasy.org/
# http://www.ncbi.nlm.nih.gov/genomes/GenomesHome.cgi?taxid=10239
# Let's get fasta for HPV16 and HIV1
ids <- data.frame(V1=c("NC_001526", "NC_001802"))
ids$urls <- str_c("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=", ids$V1, "&rettype=fasta&retmode=text")
#-----------------------------------------------------------------------------
saveMyFasta <- function(url) {
  x <- readDNAStringSet(filepath = url, format = "fasta", nrec = 1L, skip = 0L, use.names = TRUE)
  writeXStringSet(x, str_c(ids[ids$urls == url, ]$V1,".fa"))
}
#-----------------------------------------------------------------------------
l_ply(ids$urls, saveMyFasta, .progress = "text")
#-----------------------------------------------------------------------------