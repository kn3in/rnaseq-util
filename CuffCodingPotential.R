# Given gtf file from cufflinks/cuffmerge
# predict protein coding potential of every transcript
# using CPAT http://dldcc-web.brc.bcm.edu/lilab/liguow/CGI/cpat/_build/html/index.html
# save as bed file
# Map color to coding potential via itemRgb
# Darkgreen->Coding; Indianred->Non-coding

# set itemRgb="On" when feeding bed to UCSC
# consider using bedToBigBed from 
# http://genomewiki.ucsc.edu/index.php/Kent_source_utilities

#---------------------------------------------------------------
library(stringr)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg19) # change accordingly
library(GenomicFeatures)
#---------------------------------------------------------------
cuffMerged <- import.gff2("merged.gtf", asRangedData=FALSE)

#---------------------------------------------------------------
# cufflinks going to produce a few single-exon transcripts,
# if RNA-seq protocol is not strand specific there is no way
# to tell strand of a single-exon transcript.
# (you decide whether you belive in those single-exon models)
# one way is to drop all single exon "new genes"
# (flag 'u' in cufflinks terminology) or
# we have to have sense/asense pair for each single-exon transcript
#---------------------------------------------------------------

# Here we're going to add _Sense _ASense pair per each un-stranded
# transcript regardless of exon composition

# track strand information as metadata
values(cuffMerged)$strnd <- strand(cuffMerged)
ind <- values(cuffMerged)$strnd == "*"

# tmp copy is going to have "+" strand
tmp <- cuffMerged[ind]
strand(tmp) <- "+"
values(tmp)$transcript_id <- paste(values(tmp)$transcript_id, "_Sense", sep="")

# fix strand to "-" for the rest
strand(cuffMerged[ind]) <- "-"
values(cuffMerged[ind])$transcript_id <- paste(values(cuffMerged[ind])$transcript_id, "_ASense", sep="")

# Merge together
cuffMerged <- append(cuffMerged, tmp)

#---------------------------------------------------------------
# Extract fasta sequence given genome and transcript models
#---------------------------------------------------------------
trx <- split(cuffMerged, values(cuffMerged)$transcript_id)
mySeq <- extractTranscriptsFromGenome(Hsapiens, trx, decreasing.rank.on.minus.strand=TRUE)
writeXStringSet(mySeq, "Merged.fa")

#---------------------------------------------------------------
# Run CPAT prediction
#---------------------------------------------------------------
# Given python and cpat installed
system("cpat.py -g Merged.fa -o Mergedcpat -d ~/src/CPAT-1.2.1/dat/Human_train.RData -x ~/src/CPAT-1.2.1/dat/Human_Hexamer.tab")

#---------------------------------------------------------------
# Get CPAT prediction and set itemRgb accordingly
#---------------------------------------------------------------
cpat <- read.table("Mergedcpat", sep="\t", header=TRUE)
# !!!Human cut-off, change for your species accordingly!!!
# Darkgreen->Coding; Indianred->Non-coding
cpat$col <- ifelse(cpat$coding_prob >=0.363, "darkgreen", "indianred")
cpat$trx_id <- rownames(cpat)
rownames(cpat) <- NULL
cpat<-cpat[ ,c("trx_id", "col")]

# cpat messed up names, fixing it!
cpat$trx_id <- str_replace(cpat$trx_id, "ENSE$", "ense")

# set itemRgb
values(trx)$itemRgb <- cpat[match(names(trx), cpat$trx_id), ]$col
my_bed <- asBED(trx)

# set thick range to zero; default is whole transcript range
# would be nice to actually pull CDS from cpat and set thick range accordingly
my_bed$thick <- IRanges(end(ranges(my_bed)), width=0)

export(my_bed, "merged_coding_pot_paint.bed", format = "bed")

#---------------------------------------------------------------
# Extra
# You may want do this as well
# From command line:
# sort bed
bedSort merged_coding_pot_paint.bed merged_coding_pot_paint.bed
# convert bed to bigbed
bedToBigBed merged_coding_pot_paint.bed ~/seq/hg19/hg19Kent.genome merged_coding_pot_paint.bb
# help: http://genomewiki.ucsc.edu/index.php/Kent_source_utilities
#---------------------------------------------------------------
