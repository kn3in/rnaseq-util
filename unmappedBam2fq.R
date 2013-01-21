#!/usr/bin/env Rscript

# Take tophat unmapped reads, keep only paired reads(i.e. both reads are unmapped)
# add back \1 \2 to names of reads
# filter them based on quality and complexity of a sequence
# save as fastq files

file <- commandArgs(trailingOnly=TRUE)
if(length(file) != 1) stop("Usage: unmappedBam2fq.R file")

library(Biostrings)
library(ShortRead)
library(stringr)

# reads in read name, sequence, quality and bam flag
read_my <- function(path) {
  bum <- scanBam(path, param=ScanBamParam(what=c("qname", "seq", "qual", "flag"), flag=scanBamFlag(isUnmappedQuery=TRUE)))
  bum <- bum[[1]]
}

subByind <- function(my_list, ind) lapply(my_list, function(x) x[ind])

# given we operating on unmapped reads we are not
# expecting to see multiple appearance of a read (contrast with multi-hits reads in accepted_hits.bam)
# tophat clips off \1 and \2 from the end of reads names (the only difference in names of a pair)
# hence we want reads whose names present exactly twice in whole unmapped.bam file
# as control you can always check table(table(bam$qname))
# bam flag isPaired is set to TRUE even if only one read is present in the unmapped.bam file

my_paired <- function(bum) {
  my_cnt <- table(bum$qname)
  ind <- my_cnt == 2
  my_names <- names(my_cnt[ind])
  ind2 <- bum$qname %in% my_names
  bum <- subByind(bum, ind2)
}

left_ind  <- function(bum) {
  bamFlagAsBitMatrix(bum$flag, bitnames=FLAG_BITNAMES)[ ,"isFirstMateRead"] == 1
}

right_ind <- function(bum) {
  bamFlagAsBitMatrix(bum$flag, bitnames=FLAG_BITNAMES)[ ,"isSecondMateRead"] == 1
}

bam2fq <- function(bam) {
  my_fq <- ShortReadQ()
  my_fq@id <- BStringSet(bam$qname)
  my_fq@sread <- bam$seq
  my_fq@quality <- FastqQuality(BStringSet(bam$qual))
  my_fq
}

reads2paired_list <- function(paired_reads) {
  all_fq   <- bam2fq(paired_reads)
  left_fq  <- all_fq[left_ind(paired_reads)]
  right_fq <- all_fq[right_ind(paired_reads)]
  
  left_fq  <- left_fq[order(as.character(id(left_fq)))]
  right_fq <- right_fq[order(as.character(id(right_fq)))]
  stopifnot(
    all.equal(as.character(id(right_fq)), as.character(id(left_fq)))
            )
  list(left=left_fq, right=right_fq)
}

add_pair_name <- function(pair_list) {
  pair_list$left@id  <- BStringSet(str_c(as.character(id(pair_list$left)), "/1"))
  pair_list$right@id <- BStringSet(str_c(as.character(id(pair_list$right)), "/2"))
  pair_list
}

filter_quality <- function(pair_list, cut_off=25) {
  left_q <- rowMeans(as(quality(pair_list$left), "matrix"))
  right_q <- rowMeans(as(quality(pair_list$right), "matrix"))
  index <- left_q > cut_off & right_q > cut_off
  pair_list <- subByind(pair_list, index)
}

sweep_floor <- function(pair_list, cut_off=70) {
  dust_left  <- dustyScore(pair_list$left)
  dust_right <- dustyScore(pair_list$right)
  index <- dust_left > cut_off & dust_right > cut_off
  pair_list <- subByind(pair_list, index)
}

rdw <- function(paired_list) {
  writeFastq(paired_list$left,  "unmapped_1.fq")
  writeFastq(paired_list$right, "unmapped_2.fq")
}

rdw(sweep_floor(filter_quality(add_pair_name(reads2paired_list(my_paired(read_my(file)))))))