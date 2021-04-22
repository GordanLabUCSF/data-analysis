library(MSstats)

#read in file
SKY.data <- read.table("MSstats Input.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)

#convert to MSstats format
Skyfiltered.data <- SkylinetoMSstatsFormat(SKY.data, annotation = NULL, removeiRT = TRUE, useUniquePeptide = TRUE, removeOxidationMpeptides = FALSE, filter_with_Qvalue = FALSE)

#generate quant pdfs before running normalization
SKYnone.TMP <- dataProcess(raw = Skyfiltered.data, fillIncompleteRows = TRUE, normalization = FALSE,  summaryMethod = "TMP", censoredInt = "NA", cutoffCensored = "minFeature", MBimpute = TRUE)
dataProcessPlots(data = SKYnone.TMP, type= "ProfilePlot", ylimUp=35, text.angle = 90, address = "SKYfiltered_before")
dataProcessPlots(data = SKYnone.TMP, type= "QCplot", ylimUp=35, text.angle = 90, address = "SKYfiltered_before")

#generate quant pdfs after running normalization
SKYquant.TMP <- dataProcess(raw = Skyfiltered.data, fillIncompleteRows = TRUE, normalization = "equalizeMedians", summaryMethod = "TMP", censoredInt = "NA", cutoffCensored = "minFeature", MBimpute = TRUE)
dataProcessPlots(data = SKYquant.TMP, type= "ProfilePlot", ylimUp=35, text.angle = 90, address = "SKYfiltered_after_quant_")
dataProcessPlots(data = SKYquant.TMP, type= "QCplot", ylimUp=35, text.angle = 90, address = "SKYfiltered_after_quant_")

#print to console the conditions in the dataset
levels(SKYquant.TMP$ProcessedData$GROUP_ORIGINAL)
#[1] "DMSO" "DOX" 

#creating a comparison matrix
comparison1<-matrix(c(-1,1),nrow=1)

#merge all comparisions into single matrix
comparison<-rbind(comparison1)

#give comparison a name
row.names(comparison) <- c("DMSO-DOX")

#run comparision of different conditions
SKYquant.comparisons <- groupComparison(contrast.matrix = comparison, data = SKYquant.TMP)
write.table(SKYquant.comparisons$ComparisonResult, "SKYfiltered_ComparisonResults.txt", sep="\t", row.names=FALSE)

#print out sample log2 quantification values per condition and replicate
sampleQuant <- quantification(SKYquant.TMP, type = "Sample")
write.table(sampleQuant, "SKYfiltered_sampleQuant.txt", sep="\t", row.names=FALSE)

#print out sample log2 quantification values per condition
sampleQuant <- quantification(SKYquant.TMP, type = "Group")
write.table(sampleQuant, "SKYfiltered_GroupQuant.txt", sep="\t", row.names=FALSE)

#append gene names
data <- read.table("SKYfiltered_ComparisonResults.txt", header = TRUE, sep = "\t", quote = "", check.names=FALSE, stringsAsFactors=FALSE)
additional <- read.table("uniprot_protein_descriptions_HUMAN.txt", header = TRUE, sep = "\t", quote = "", check.names=FALSE, stringsAsFactors=FALSE)

results <- data[-c(1:nrow(data)),, drop = FALSE] # create empty set
toAppend <- additional[-c(1:nrow(additional)),, drop = FALSE]
for (k in 1:nrow(data)) { # loop through each row of dataset
  prot_id <- data[k,1] # get protein ID, used to find unique peptides of that protein
  for (u in 1:nrow(additional)) { # loop through all proteins
    found <- grepl(additional$Entry[u], prot_id) # check if protein k of dataset matches protein u
    if (found == TRUE) {
      break # break out of inner forloop and save u = location (row #)
    }
  }
  toAppend <- rbind(toAppend, additional[u,]) # put protein info to append in correct info
}
results <- cbind(data, toAppend$`Gene names`) # append additional info to filtered dataset

write.table(results, "639VPRKACA-MIBsFinalComparisons.txt", sep="\t", row.names = FALSE)
