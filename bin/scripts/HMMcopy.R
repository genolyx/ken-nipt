args <- commandArgs(TRUE)

library(HMMcopy)
library(MASS)

rfile <- args[1]
gfile <- args[2]
mfile <- args[3]
resolution <- args[4]  # resolution argument ("50kb" or "10mb")

normal <- correctReadcount(wigsToRangedData(rfile, gfile, mfile))
na.filter <- normal[-which(is.na(normal$copy))] # R>=3.6
all <- as.data.frame(na.filter) # R>=3.6

for (i in c(1:22,"X","Y")){
  assign(paste0("chr",i), subset(all, all[,1]==paste0("chr",i)))
}

pattern <- ifelse(resolution == "50kb", ".50kb.wig", ".10mb.wig")
output_file <- paste0(gsub(pattern, "", rfile), pattern, ".Normalization.txt")

write.table(normal, output_file,
            quote=FALSE, sep='\t', col.names=TRUE, row.names=FALSE)
