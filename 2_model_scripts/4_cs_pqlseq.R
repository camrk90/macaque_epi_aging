#!/usr/bin/env /packages/apps/spack/18/opt/spack/gcc-11.2.0/r-4.2.2-kpl/bin/Rscript

#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ckelsey4@asu.edu

SAMP <- Sys.getenv("SLURM_ARRAY_TASK_ID")
SAMP <- as.integer(SAMP)

#### THIS SCRIPT RUNS PQLseq FOR THE CROSS SECTIONAL DATA ####

library(tidyverse)
library(PQLseq2)
setwd("/scratch/ckelsey4/Cayo_meth/cross_models")

#Import metadata----------------------------------------------------------------
blood_metadata<- read.table("/scratch/ckelsey4/Cayo_meth/blood_metadata_full.txt")
long_data<- read.table("/scratch/ckelsey4/Cayo_meth/long_data_adjusted.txt")

blood_metadata<- blood_metadata[!blood_metadata$lid_pid %in% long_data$lid_pid,]

blood_metadata<- blood_metadata %>%
  drop_na(age_at_sampling)

#Import kinship matrix----------------------------------------------------------
kinship<- readRDS("/scratch/ckelsey4/Cayo_meth/full_kin_matrix")

#Subset and rearrange kinship rows and cols to match metadata
kinship<- kinship[blood_metadata$lid_pid, blood_metadata$lid_pid]

#Check kinship matrix dims match metadata
all.equal(nrow(blood_metadata), ncol(kinship))
all.equal(nrow(blood_metadata), nrow(kinship))

#Import m/cov rds------------------------------------------------------------
regions_cov<- readRDS("/scratch/ckelsey4/Cayo_meth/regions_cov_filtered.rds")

#Filter metadata to lids in regions list
blood_metadata<- blood_metadata[blood_metadata$lid_pid %in% colnames(regions_cov[[1]]),]

regions_cov<- lapply(names(regions_cov), function(x){
  regions_cov<- subset(regions_cov[[x]], select=blood_metadata$lid_pid)
  return(regions_cov)
})

regions_m<- readRDS("/scratch/ckelsey4/Cayo_meth/regions_m_filtered.rds")

regions_m<- lapply(names(regions_m), function(x){
  regions_m<- subset(regions_m[[x]], select=blood_metadata$lid_pid)
  return(regions_m)
})

names(regions_cov)<- 1:21 #turn all chroms into integers (X = 21)
names(regions_m)<- 1:21 #turn all chroms into integers (X = 21)

rm(long_data)

#Check metadata lids match the lids (cols) of a random chromosome
if (all.equal(blood_metadata$lid_pid, colnames(regions_cov[[runif(1, 1, 21)]])) == T) {
  
  #Model Vectors for lme4---------------------------------------------------------
  cov<- regions_cov[[SAMP]]
  meth<- regions_m[[SAMP]]
  
  ###################################
  #####        Run PQLseq       #####
  ###################################
  #Run PQLseq for within_age-------------------------------------------------------------
  #Generate model matrix
  predictor_matrix<- model.matrix(~ age_at_sampling + individual_sex + unique, data = blood_metadata)
  age_pheno<- predictor_matrix[, 2]
  age_cov<- predictor_matrix[, 3:4]
  
  #Run pqlseq model
  age_model<- pqlseq2(Y = meth, x = age_pheno, 
                                K = kinship, W = age_cov, 
                                lib_size = cov, model="BMM")
  
  #Save pqlseq models
  saveRDS(age_model, paste("wb", "cs", SAMP, ".rds", sep = "_"))
  
} else {
  
  print("metadata lids did not match cov matrix lids")
  
}



