
## pre flight: break up reference genome into 200 bp windows using bedtools & use bedtools map to get average per read level PDN values for each 200 bp window 

##handling the per genomic window levels of PDN - each of the 200 bp windows of the genome has the average PDN value of the reads covering it (reads have to be 51% or more within the region) and
##the number of reads which cover that region

library(dplyr)
library(stringr)
library(readr)

long_data<- readRDS("/scratch/ckelsey4/Cayo_meth/long_data_adjusted")
long_data<- long_data %>%
  dplyr::rename(sex = individual_sex) %>%
  filter(n > 1) %>%
  arrange(lid_pid)

lid_pid<- long_data$lid_pid
lid_pid<- paste(lid_pid, ".meanPDN.bed", sep = "")

##set unique file ID's
unique_fileID <- list.files(path = "/scratch/nsnyderm/cayo_rrbs/discordance/",
                            pattern = "*.meanPDN.bed")

overlap <- read.table(filePaths[1], header = FALSE, fill = TRUE)

unique_fileID<- unique_fileID[unique_fileID %in% lid_pid]


##set file paths
filePaths <- paste0("/scratch/nsnyderm/cayo_rrbs/discordance/", unique_fileID)

##loop to calculate PDN
for(i in 1:length(unique_fileID)){
  ##read in overlap file written as a plain text file with headers removed
  overlap <- read.table(filePaths[i], header = FALSE, fill = TRUE)
  ##replace any regions with less than 5 reads covering it (number of reads is V5)
  overlap2 <- overlap %>%  mutate_at(vars(V4, V5), ~replace(., V5 <= 4, NA)) ##less than or equal to 4
  ##rename the avg PDN with sample name
  assign(paste(unique_fileID[i], "PDN", sep = "_"), overlap2$V4)
  assign(paste(unique_fileID[i], "reads", sep = "_"), overlap2$V5)
  rm(overlap) ##remove large files from memory to reduce load
  rm(overlap2)
}

##now bind all sample PDN values together with original coordinate values
coord <- read.table("/scratch/emb19132/GENOME_PDR/RAT/PDN_output/SRR13012807_PDN_overlap.txt",  header = FALSE)
##give each set of coordinates a name (chrom - V1, start - V2, end - V3)
coord$region.name <- paste(coord$V1,coord$V2,coord$V3, sep=":")
assign(paste("genomic_region"), coord$region.name)

compiled <- cbind(genomic_region, )

##remove rows with more than 20% NAs
compiled <- as.data.frame(compiled)
row.names(compiled) <- compiled$genomic_region
compiled <- compiled[,c(-1)]
compiled <- compiled[!rowSums(is.na(compiled)) > ncol(compiled)*.2,] ##remove regions not covered by at least 20% of samples

write.csv(compiled, file="/scratch/emb19132/GENOME_PDR/RAT/PDN_output/rat_merged_region_PDN.csv")