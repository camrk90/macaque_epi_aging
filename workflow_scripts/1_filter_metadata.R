library(tidyverse)
library(GENESIS)
setwd("/scratch/ckelsey4/Cayo_meth")

#FILTER METADTA FOR REPEATED SAMPLES--------------------------------------------
#import metadata
blood_metadata<- read.table("metadata_temp_clean_241106.txt", header = T, fill = T)

#Filter for whole blood and add pid col
blood_metadata<- blood_metadata %>%
  filter(grantparent_tissueType == "whole_blood") %>%
  mutate(pid = lid_pid) %>%
  relocate(pid, .after = lid_pid)

#separate pid col by '_' to isolate pid number
blood_metadata<- blood_metadata %>%
  separate_wider_delim(pid, delim = "_",
                       names = c("a", "b", "c", "pid")) %>%
  dplyr::select(-c(a, b, c))

#filter for samples with n>2 per id
blood_metadata<- blood_metadata %>%
  group_by(monkey_id) %>%
  mutate(n = n()) %>%
  filter(n >= 2)

#Add mean age col
blood_metadata<- blood_metadata %>%
  group_by(monkey_id) %>%
  mutate(mean.age = mean(age_at_sampling))

#Add mean-centred age col
blood_metadata<- blood_metadata %>%
  mutate(within.age = age_at_sampling - mean.age)

#Arrange by id
blood_metadata<- blood_metadata %>%
  ungroup() %>%
  arrange(monkey_id)

#Select important cols
blood_metadata_short<- blood_metadata %>%
  dplyr::select(monkey_id, lid_pid, pid, age_at_sampling, mean.age, within.age, individual_sex, n)

#Filter out ids that have two entries at the same age
blood_metadata_short<- blood_metadata_short %>%
  group_by(monkey_id) %>%
  distinct(age_at_sampling, .keep_all = T) %>%
  mutate(n = n())

#Drop NA rows
blood_metadata_short<- blood_metadata_short %>%
  drop_na()

saveRDS(blood_metadata_short, 'long_data_adjusted')

########################### GENERATE KINSHIP MATRIX ############################
#Generate vector of animal ids
monkey_vector<- blood_metadata_short$monkey_id

#Import king output file
file.king <- c("king.kin0")

#Generate kin matrix
kin.matrix<-kingToMatrix(file.king, estimator = "Kinship", sample.include=monkey_vector)
kinmat<-as.matrix(kin.matrix)

#Arrange kinmat colnames by metadata to match r_matrix
kinmat<- kinmat[unique(blood_metadata_short$monkey_id),unique(blood_metadata_short$monkey_id)]

#Generate id vectors for Z-matrix of samples(lids) x individuals(monkey_id)
monkey_ids<- blood_metadata_short %>%
  dplyr::select(monkey_id) %>%
  unique()
lids<- blood_metadata_short %>%
  ungroup() %>%
  dplyr::select(lid_pid)

#generate empty z-matrix
r_matrix<- data.frame(matrix(ncol = nrow(monkey_ids), nrow=nrow(lids)))
colnames(r_matrix)<- monkey_ids$monkey_id
rownames(r_matrix)<- lids$lid_pid

#Add monkey_id column to match colnames to
all.equal(rownames(r_matrix), blood_metadata_short$lid_pid)
r_matrix$ids<- blood_metadata_short$monkey_id

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

#IT'S NOT CLEAR THIS IS NEEDED AS MOST OF THE DUPLICATED IDS THAT GET SUBSETTED BY
#MEAN COV END UP BEING ONE SAMPLE AND GET FILTERED OUT BEFORE GLMER ANYWAY
#KEEPING THE CODE HERE JUST IN CASE
#Import m/cov rds------------------------------------------------------------
# load region lists that have been filtered for 5x coverage in 90% of samples
regions_cov_list<- readRDS("regions_cov_list")
regions_m_list<- readRDS("regions_m_list")

#Subset regions list for lids with repeated measures
regions_cov<- lapply(names(regions_cov_list), function(x){
  regions_cov<- subset(regions_cov_list[[x]], select=duplicate_lids$lid_pid)
  return(regions_cov)
})

regions_m<- lapply(names(regions_m_list), function(x){
  regions_m<- subset(regions_m_list[[x]], select=duplicate_lids$lid_pid)
  return(regions_m)
})

regions_cov<- do.call(rbind, regions_cov)
regions_m<- do.call(rbind, regions_m)
ratio<- as.matrix(regions_m/regions_cov)
ratio[is.nan(ratio)]<- 0

m_sum<- as.data.frame(colMeans(regions_m))
all.equal(rownames(m_sum), duplicate_lids$lid_pid)
duplicate_lids<- cbind(duplicate_lids, m_sum)

cov_sum<- as.data.frame(colMeans(regions_cov))
all.equal(rownames(cov_sum), duplicate_lids$lid_pid)
duplicate_lids<- cbind(duplicate_lids, cov_sum)

perc_meth<- as.data.frame(colMeans(ratio))
all.equal(rownames(perc_meth), duplicate_lids$lid_pid)
duplicate_lids<- cbind(duplicate_lids, perc_meth)

duplicate_lids_min<- duplicate_lids %>%
  group_by(monkey_id) %>%
  slice_min(`colMeans(regions_cov)`, n = 1)

#Subset out lids with multiple entries at the same age
blood_metadata_short<- blood_metadata_short[!blood_metadata_short$lid_pid %in% duplicate_lids_min$lid_pid,]






