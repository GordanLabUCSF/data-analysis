# MIBS
# extract Gene names and log2FC
rm(list = ls())
setwd("C:/Users/Gordan Lab/Box/Gordan Lab FGFR/Network Propagation/ht-kam/");
data <- read.table("13-7MIBsFinalComparisons.txt", header = TRUE, sep = "\t", quote = "", stringsAsFactors=FALSE, fill=TRUE)

# extract first gene name and log2FC, normalize to highest value
Genes <- sub("(\\w+).*", "\\1", data$Gene.names)
log2FC <- data$log2FC
highest <- max(data$log2FC)
norm <- log2FC/highest

newfile <- cbind(Genes,log2FC,norm)
#write.table(newfile, "mibs.txt", sep="\t", row.names = FALSE)

# PHOSFATE
phosfate <- read.table("13-7_24h_kinase_activity.txt", header = TRUE, sep = "\t", quote = "", stringsAsFactors=FALSE, check.names=FALSE, fill=TRUE)

# extract first gene name and log2FC, normalize to highest value
Genes <- phosfate$kinase
value <- phosfate$activity
highest <- max(phosfate$activity)
norm <- value/highest

newphosfate <- cbind(Genes,value,norm)
#write.table(newphosfate, "phosfate.txt", sep="\t", row.names = FALSE)

# combine normalized mibs and phosfate data
both <- as.data.frame(rbind(newfile,newphosfate), stringsAsFactors = FALSE)

# get duplicated rows, store them, but delete from the combined table
duprows <- both[duplicated(both$Genes)|duplicated(both$Genes, fromLast=TRUE),]
toAdd <- both[duplicated(both$Genes),]
toAdd <- toAdd[order(toAdd$Genes),]

both <- both[!(duplicated(both$Genes) | duplicated(both$Genes, fromLast = TRUE)), ]

duprows <-duprows[order(duprows$Genes),]
#both <- aggregate(both$norm,by=list(Genes=both$Genes),data=both,FUN=mean)
# 
i <- 1
for (p in 1:nrow(duprows)) {
  if((p %% 2) == 0) {
    next
  }
  average <- (as.numeric(duprows$norm[p])+as.numeric(duprows$norm[p+1]))/2
  toAdd$average[i] <- average
  i <- i+1
}

final <- as.data.frame(both$Genes)
final$norm <- both$norm
finalToAdd <- as.data.frame(toAdd$Genes)
finalToAdd$norm <- toAdd$average
names(final)[1] <- "Genes"
names(finalToAdd)[1] <- "Genes"

propfile <- rbind(final,finalToAdd)
mean <- mean(as.numeric(propfile$norm))
sd <- sd(as.numeric(propfile$norm))
propfile$original_z_score <- (as.numeric(propfile$norm)-mean)/sd
propfile$z_score <- abs(propfile$original_z_score)

# write table containing z-scores of union of MIBs and phosfate, with any duplicates being averaged
write.table(propfile, "13-7phosmibs.txt", sep="\t", row.names = FALSE)
