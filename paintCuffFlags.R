# Paint cuffmerged/cufflinks gtf according to cufflinks transcript flags
# http://cufflinks.cbcb.umd.edu/manual.html
# output merged.bed file

# set itemRgb="On" when feeding bed to UCSC
# consider using bedToBigBed from 
# http://genomewiki.ucsc.edu/index.php/Kent_source_utilities

library(rtracklayer)

cuffMerged <- import.gff2("merged.gtf", asRangedData=FALSE)

my_flags <- c("=", "c", "j", "e", "i", "o", "p", "r", "u", "x", "s", ".")
my_colors <- c("black", "gray", "steelblue", "darkmagenta", "darkolivegreen", "pink4", "darkorange1", "cadetblue4", "darkred", "bisque4", "chartreuse", "darkorchid4")

# flag to color mapping
lupt <- as.data.frame(cbind(flags=my_flags, col=my_colors))

#transcript_id to class_code(flag) mapping
tr2code <- as.data.frame(unique(values(cuffMerged)[ ,c("transcript_id", "class_code")]))

tr2col <- merge(tr2code, lupt, by.x="class_code", by.y="flags", all.x=TRUE)
grl_cm <- split(cuffMerged, values(cuffMerged)$transcript_id)
values(grl_cm)$itemRgb <- tr2col[match(names(grl_cm), tr2col$transcript_id), "col"]

export(grl_cm, "merged.bed", "bed")