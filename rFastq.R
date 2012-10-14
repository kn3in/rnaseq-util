# R CMD BATCH rFastq.R

# Sample 1e6 reads and save into tmp folder
# should be used for QA, aligner run etc

library(ShortRead)
library(parallel)
library(stringr)

fastqFiles <- dir(pattern = "*.fastq.gz$")
sampl_name <- unique(str_split_fixed(fastqFiles,"_[1,2].fastq.gz$" , n=2)[ ,1])

# Be nice
ncores <- min(length(fastqFiles), 40) 

mclapply(sampl_name, function(fl) {
  base <- paste("tmp", fl, sep="/")
  left <- paste(fl, "_1.fastq.gz", sep="")
  right <- paste(fl, "_2.fastq.gz", sep="")
  dir.create("tmp", showWarnings=TRUE, recursive=TRUE, mode="0755")
  f1 <- FastqSampler(left, n=1e6)
  f2 <- FastqSampler(right, n=1e6)
  set.seed(123L)
  p1 <- yield(f1)
  set.seed(123L)
  p2 <- yield(f2)
  writeFastq(p1, paste("tmp", str_c(fl, "_1.fastq"), sep="/"))
  writeFastq(p2, paste("tmp", str_c(fl, "_2.fastq"), sep="/"))
}, mc.cores=ncores)