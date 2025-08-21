library(tidyverse)
library(GENESIS)

path<- "/scratch/ckelsey4/Cayo_meth/"

#FILTER METADTA FOR REPEATED SAMPLES--------------------------------------------
#import metadata
blood_metadata<- read.table(paste0(path, "metadata_temp_clean_241106.txt", sep = ""), header = T, fill = T)

#Filter for whole blood and add pid col
blood_metadata<- blood_metadata %>%
  filter(grantparent_tissueType == "whole_blood") %>%
  mutate(pid = paste0("PID_", str_split_i(lid_pid, "_", 4))) %>%
  relocate(pid, .after = lid_pid)

#Add university prepped
blood_metadata<- blood_metadata %>%
  mutate(prep_year = year(prep_date))

blood_metadata$university<- "uw"
blood_metadata$university[blood_metadata$prep_year > 2019]<- "asu"

#filter for samples with n>2 per id
long_metadata<- blood_metadata %>%
  group_by(monkey_id) %>%
  mutate(n = n()) %>%
  filter(n >= 2)

#Add mean age and within age col
long_metadata<- long_metadata %>%
  group_by(monkey_id) %>%
  mutate(mean.age = mean(age_at_sampling),
         within.age = age_at_sampling - mean.age) %>%
  ungroup() %>%
  arrange(monkey_id)

#Select important cols
long_metadata<- long_metadata %>%
  dplyr::select(monkey_id, lid_pid, pid, age_at_sampling, mean.age, within.age, individual_sex, n, 
                processing_timestamp, prep_date, university, unique)

write.table(long_metadata_short, 'long_data_adjusted.txt',
            quote=F)
write.table(blood_metadata, "blood_metadata_full.txt", 
            quote=F)

#Organize BSseq data in chrs
#Autosomes/ChrX
cayo<- readRDS("/scratch/nsnyderm/cayo_rrbs/bismarkBSseq.rds")
colnames(cayo)=gsub(".CpG_report.merged_CpG_evidence.cov.gz","",str_split_i(colnames(cayo),"\\/",6))

#ChrY
cayo_Y<- readRDS("bismarkBSseq_Y.rds")
colnames(cayo_Y)=gsub(".CpG_report.merged_CpG_evidence.cov.gz","",colnames(cayo_Y))

selected_chrs<- c(1:20, "X")

#Arrange metadata by lid
blood_metadata<- blood_metadata %>%
  arrange(lid_pid)

#Subset cayo bsseq with lids from longitudinal metadata
cayo_blood<- cayo[, cayo@colData@rownames %in% blood_metadata$lid_pid]
cayo_Y<- cayo_Y[, cayo_Y@colData@rownames %in% blood_metadata$lid_pid]
all.equal(cayo_blood@colData@rownames, blood_metadata$lid_pid)

#Split bsseq by chromosome
cayo_blood_list<- parallel::mclapply(selected_chrs,function(x){
  chrSelectBSseq(cayo_blood, seqnames = x, order = TRUE)
}, mc.cores=21)

names(cayo_blood_list)<- 1:21

saveRDS(cayo_blood_list, paste0(path, "cayo_blood_full.rds", sep=""))

########################### GENERATE KINSHIP MATRIX ############################
#Generate vector of animal ids
monkey_vector<- blood_metadata$monkey_id

#Import king output file
file.king <- c("king.kin0")

#Generate kin matrix
kin.matrix<- kingToMatrix(file.king, estimator = "Kinship", sample.include=monkey_vector)
kinmat<- as.matrix(kin.matrix)

#Arrange kinmat colnames by metadata to match r_matrix
kinmat<- kinmat[unique(blood_metadata$monkey_id),unique(blood_metadata$monkey_id)]

#Generate id vectors for Z-matrix of samples(lids) x individuals(monkey_id)
monkey_ids<- blood_metadata %>%
  dplyr::select(monkey_id) %>%
  unique()
lids<- blood_metadata %>%
  ungroup() %>%
  dplyr::select(lid_pid)

#generate empty z-matrix
r_matrix<- data.frame(matrix(ncol = nrow(monkey_ids), nrow=nrow(lids)))
colnames(r_matrix)<- monkey_ids$monkey_id
rownames(r_matrix)<- lids$lid_pid

#Add monkey_id column to match colnames to
all.equal(rownames(r_matrix), blood_metadata$lid_pid)
r_matrix$ids<- blood_metadata$monkey_id

#Assign 1's to colnames that match ids in the id column (i.e. 1's for the same id)
r_matrix[sapply(colnames(r_matrix), `==`, r_matrix$ids)] <- 1

#Remove id column and assign class matrix
r_matrix<- r_matrix %>%
  dplyr::select(-length(r_matrix))
r_matrix<- as.matrix(r_matrix)

#Replace NAs with 0s
r_matrix[is.na(r_matrix)]<- 0

#Multiply matrices together to get full kinship matrix
full_kin<- r_matrix %*% kinmat %*% t(r_matrix)

#Save output file as rds
saveRDS(full_kin, paste0(path, "full_kin_matrix.rds", sep=""))
