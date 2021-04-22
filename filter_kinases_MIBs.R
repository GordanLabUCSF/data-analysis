
filter_kinases_MIBs <- function(path_to_data, data_file_string, path_to_database) {
  # read in files
  data <- read.table(paste0(path_to_data, data_file_string), header = TRUE, sep = "\t", check.names=FALSE, stringsAsFactors=FALSE)
  data_file <- gsub(".txt", "", data_file_string)
  kinases <- read.table(paste0(path_to_database, "uniprot_pkinfam.txt"), header = FALSE, sep="\t", stringsAsFactors=FALSE)
  kinases <- append(kinases[,1], kinases[,3])
  kinases <- kinases[kinases != ""]
  #proteins <- read.table(paste0(path_to_data, "proteinGroups.txt"), header = TRUE, sep="\t", stringsAsFactors=FALSE)

  data1 <- data[!grepl(";", data$Proteins), ] # check for unique peptides that only match to single protein

  # filter for only kinases
  data2 <- data1[data1$Proteins %in% kinases, ]

  write.table(data2, paste0(path_to_data, data_file, "_filtered.txt"), sep="\t", row.names = FALSE)
  
  # filter out one-hit wonders - we don't use this as of yet
  # finaldata <- data2[-c(1:nrow(data2)),, drop = FALSE] # create empty set to add final filtered data to
  # for (k in 1:nrow(data2)) { # loop through each row of dataset
  #   prot_id <- data2[k,13] # get protein ID, used to find unique peptides of that protein from proteinGroups.txt file
  #   for (u in 1:nrow(proteins)) { # loop through all proteins
  #     prot_list <- proteins$Protein.IDs[u]
  #     proteinfound <- grepl(prot_id, proteins$Protein.IDs[u]) # check if protein k of dataset matches protein u of proteinGroups.txt
  #     if (proteinfound == TRUE) {
  #       break # break out of inner forloop and save u = location (row #) in proteinGroups.txt of corresponding protein from dataset
  #     }
  #   }
  #   if (proteins[u,12] > 1){ # if there are more than 1 unique peptides for the protein
  #     finaldata <- rbind(finaldata, data2[k,]) # add corresponding row from dataset into created dataset
  #   }
  # }
  
  # filter out modified proteins - we don't use this as of yet either
  #finaldata <- data[grepl("Unmodified", data$Modifications), ]

  # sort data by id to compare files if necessary
  #sortdata2 <- finaldata[-c(1:nrow(finaldata)),, drop = FALSE]
  #sortdata2 <- with(finaldata, finaldata[order(finaldata$id), ])
  #write.table(sortdata2, "data_sorted.txt", sep="\t", row.names = FALSE)
}