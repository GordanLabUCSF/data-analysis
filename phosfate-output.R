# create phosfate output file
# no header
# uniprot id,modification site,log2fc
# on excel, replace all "#NAME?"s with "NA" first
rm(list = ls())
data <- read.table("annotated-YAPCsu360.txt", header = TRUE, sep = "\t", quote = "", stringsAsFactors=FALSE, fill=TRUE)

uniprot_id <- data[,c("Protein")]
mod_sites <- data[,c("mod_sites")]
log2fc <- data[,c("LESS.MORE_log2FC")]

phosfate_output <- cbind(uniprot_id,mod_sites,log2fc)
phosfate_output <- as.data.frame(phosfate_output)
phosfate_output[] <- lapply(phosfate_output, gsub, pattern='"', replacement='')

phosfate_output$uniprot_id <- substr(phosfate_output$uniprot_id, 0, 6)
phosfate_output$mod_sites <- substr(phosfate_output$mod_sites, 9, length(phosfate_output$mod_sites))
for (i in 1:nrow(phosfate_output)){
  comma <- regexpr(pattern=",", phosfate_output$mod_sites[i])[1]
  if (comma != -1){
    replace <- substr(phosfate_output$mod_sites[i], 1, comma-1)
    phosfate_output$mod_sites[i] <- replace
  }
}
write.table(phosfate_output, file="annotated-YAPCsu360-LESS.MORE-phosfate.txt", row.names=FALSE, col.names=FALSE, sep=",")
