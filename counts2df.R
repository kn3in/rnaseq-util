#!/usr/bin/env Rscript

# aggregate htseq-counts output into a data frame
# discard genes with 0 counts across all samples
# save resulting dataframe as 'my_counts.csv' file

files <- dir(pattern = "_count.txt")
names <- as.vector(sapply(files, function(x) strsplit(x, "_count.txt")[[1]]))
data <- lapply(files, read.table, sep = "\t", header = FALSE)
# drop last 5 summary rows
data <- lapply(data, function(x) head(x, -5))
names(data) <- names

data <- lapply(names, function(x) {
  colnames(data[[x]])[2] <- x
  data[[x]] 
})

data <- Reduce(function(x, y) merge(x, y, by.x = "V1", by.y = "V1"), data, accumulate = FALSE)
names(data)[1] <- "gene_id"
data$filter <- rowSums(data[ ,-1])
data <- subset(data, filter != 0)
data$filter <- NULL
write.csv(data, file = "my_counts.csv", quote = FALSE, row.names = FALSE)