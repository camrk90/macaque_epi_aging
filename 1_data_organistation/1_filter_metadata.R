library(tidyverse)
library(GENESIS)
setwd("/scratch/ckelsey4/Cayo_meth")
setwd("/Users/cameronkelsey/Documents/smack_lab/cayo_data/cayo_meth copy")

#FILTER METADTA FOR REPEATED SAMPLES--------------------------------------------
#import metadata
blood_metadata<- read.table("metadata_temp_clean_241106.txt", header = T, fill = T)

#Filter for whole blood and add pid col
blood_metadata<- blood_metadata %>%
  filter(grantparent_tissueType == "whole_blood") %>%
<<<<<<< HEAD
  mutate(pid = str_split_i(lid_pid, "_", 4)) %>%
=======
  mutate(pid = paste0("PID_", str_split_i(lid_pid, "_", 4))) %>%
>>>>>>> 4f8b7c952b676820ba3dd924efca5fe1003eb690
  relocate(pid, .after = lid_pid)

#Add university prepped
blood_metadata<- blood_metadata %>%
  mutate(prep_year = year(prep_date))

blood_metadata$university<- "uw"
blood_metadata$university[blood_metadata$prep_year > 2019]<- "asu"

#filter for samples with n>2 per id
<<<<<<< HEAD
long_metadata<- blood_metadata %>%
=======
long_metadata<- long_metadata %>%
>>>>>>> 4f8b7c952b676820ba3dd924efca5fe1003eb690
  group_by(monkey_id) %>%
  mutate(n = n()) %>%
  filter(n >= 2)

#Add mean age col
long_metadata<- long_metadata %>%
  group_by(monkey_id) %>%
  mutate(mean.age = mean(age_at_sampling))

#Add mean-centred age col
long_metadata<- long_metadata %>%
  mutate(within.age = age_at_sampling - mean.age)

#Arrange by id
long_metadata<- long_metadata %>%
  ungroup() %>%
  arrange(monkey_id)
<<<<<<< HEAD

#Filter out ids that have two entries at the same age
#long_metadata_short<- long_metadata_short %>%
  #group_by(monkey_id) %>%
  #distinct(age_at_sampling, .keep_all = T) %>%
  #mutate(n = n())
=======

#Select important cols
long_metadata<- long_metadata %>%
  dplyr::select(monkey_id, lid_pid, pid, age_at_sampling, mean.age, within.age, individual_sex, n, 
                processing_timestamp, prep_date, university)

#Filter out ids that have two entries at the same age
#MAYBE NEED TO ADJUST THIS TO FILTER OUT THE LID WITH THE WORST COVERAGE
long_metadata<- long_metadata %>%
  group_by(monkey_id) %>%
  distinct(age_at_sampling, .keep_all = T) %>%
  mutate(n = n())
>>>>>>> 4f8b7c952b676820ba3dd924efca5fe1003eb690

#Drop NA rows
#long_metadata_short<- long_metadata_short %>%
  #drop_na()

write.table(long_metadata_short, 'long_data_adjusted.txt',
            quote=F)
write.table(blood_metadata, "blood_metadata_full.txt", 
            quote=F)

########################### GENERATE KINSHIP MATRIX ############################
#Generate vector of animal ids
monkey_vector<- blood_metadata$monkey_id

#Import king output file
file.king <- c("king.kin0")

#Generate kin matrix
kin.matrix<-kingToMatrix(file.king, estimator = "Kinship", sample.include=monkey_vector)
kinmat<-as.matrix(kin.matrix)

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
saveRDS(full_kin, "full_kin_matrix")
