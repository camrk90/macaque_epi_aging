#!/usr/bin/env /packages/apps/spack/18/opt/spack/gcc-11.2.0/r-4.2.2-kpl/bin/Rscript

SAMP <- Sys.getenv("SLURM_ARRAY_TASK_ID")
SAMP <- as.integer(SAMP)

library(tidyverse)
library(PQLseq2)
setwd("/scratch/ckelsey4/Cayo_meth/glmer_model_compare")

#Import metadata----------------------------------------------------------------
long_data<- read.csv("/scratch/ckelsey4/Cayo_meth/long_data_adjusted.txt", sep="\t")

long_data<- long_data %>%
  arrange(monkey_id) %>%
  filter(n > 1)

#Import kinship matrix----------------------------------------------------------
kinship<- readRDS("/scratch/ckelsey4/Cayo_meth/full_kin_matrix")

#Subset and rearrange kinship rows and cols to match metadata
kinship<- kinship[long_data$lid_pid, long_data$lid_pid]

#Import m/cov rds------------------------------------------------------------
# load region lists that have been filtered for 5x coverage in 90% of samples
regions_cov<- readRDS("/scratch/ckelsey4/Cayo_meth/regions_cov_filtered")
regions_cov<- regions_cov[1:21]
regions_m<- readRDS("/scratch/ckelsey4/Cayo_meth/regions_m_filtered")
regions_m<- regions_m[1:21]

#Filter metadata to lids in regions list
long_data<- long_data[long_data$lid_pid %in% colnames(regions_cov[[1]]),]

regions_cov<- lapply(names(regions_cov), function(x){
  regions_cov<- subset(regions_cov[[x]], select=long_data$lid_pid)
  return(regions_cov)
})

regions_m<- lapply(names(regions_m), function(x){
  regions_m<- subset(regions_m[[x]], select=long_data$lid_pid)
  return(regions_m)
})

names(regions_cov)<- 1:21 #turn all chroms into integers (X = 21)
names(regions_m)<- 1:21 #turn all chroms into integers (X = 21)

if (all.equal(long_data$lid_pid, colnames(regions_cov[[runif(1, 1, 21)]])) == T) {
  
  #Check metadata lids match the lids (cols) of a random chromosome
  all.equal(long_data$lid_pid, colnames(regions_cov[[runif(1, 1, 21)]]))
  
  #Model Vectors for lme4---------------------------------------------------------
  cov<- regions_cov[[SAMP]]
  meth<- regions_m[[SAMP]]
  
  ###################################
  #####        Run PQLseq       #####
  ###################################
  #Run PQLseq for males:withinage ------------------------------------------------
  #Generate model matrix
  predictor_matrix<- model.matrix(~ individual_sex:within.age + mean.age + individual_sex, data = long_data)
  
  #Generate predictor and covariate matrices
  m_phenotype<- predictor_matrix[, 5]
  m_covariates<- predictor_matrix[, c(2:4)]
  
      #Run pqlseq model
      male_age_pqlseq2<- pqlseq2(Y = meth, x = m_phenotype, 
                                 K = kinship, W = m_covariates, 
                                 lib_size = cov, model="BMM", verbose=F)
  
  #Run PQLseq for females:withinage-----------------------------------------------
  #Generate model matrix
  f_phenotype<- predictor_matrix[, 4]
  f_covariates<- predictor_matrix[, c(2:3, 5)]
  
  #Run pqlseq model
  female_age_pqlseq2<- pqlseq2(Y = meth, x = f_phenotype, 
                               K = kinship, W = f_covariates, 
                               lib_size = cov, model="BMM")
  
  #Run PQLseq for mean_age-------------------------------------------------------------
  #Generate model matrix
  m.age_phenotype<- predictor_matrix[, 2]
  m.age_covariates<- predictor_matrix[, c(3:5)]
  
  #Run pqlseq model
  mean_age_pqlseq2<- pqlseq2(Y = meth, x = m.age_phenotype, 
                             K = kinship, W = m.age_covariates, 
                             lib_size = cov, model="BMM")
  
  #Save pqlseq models
  saveRDS(male_age_pqlseq2, paste("nested", "model", "m", SAMP, sep = "_"))
  saveRDS(female_age_pqlseq2, paste("nested", "model", "f",  SAMP, sep = "_"))
  saveRDS(mean_age_pqlseq2, paste("nested", "model", "mean", SAMP, sep = "_"))
  
} else {
  
  print("long_data lids did not match cov matrix lids")
  
}



