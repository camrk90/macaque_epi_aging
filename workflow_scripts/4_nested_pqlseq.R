#!/usr/bin/env /packages/apps/spack/18/opt/spack/gcc-11.2.0/r-4.2.2-kpl/bin/Rscript

SAMP <- Sys.getenv("SLURM_ARRAY_TASK_ID")
SAMP <- as.integer(SAMP)

library(tidyverse)
library(PQLseq2)
setwd("/scratch/ckelsey4/Cayo_meth/glmer_model_compare")

#Import metadata----------------------------------------------------------------
long_data<- readRDS("/scratch/ckelsey4/Cayo_meth/long_data_adjusted")
long_data<- long_data %>%
  arrange(lid_pid) %>%
  filter(n > 1)
kinship<- readRDS("/scratch/ckelsey4/Cayo_meth/full_kin_matrix")

#Import m/cov rds------------------------------------------------------------
regions_cov<- readRDS("/scratch/ckelsey4/Cayo_meth/regions_cov_filtered")
regions_m<- readRDS("/scratch/ckelsey4/Cayo_meth/regions_m_filtered")

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

#Check metadata lids match the lids (cols) of a random chromosome
all.equal(long_data$lid_pid, colnames(regions_cov[[runif(1, 1, 21)]]))

names(regions_cov)<- 1:21 #turn all chroms into integers (X = 21)
names(regions_m)<- 1:21 #turn all chroms into integers (X = 21)

#Model Vectors for lme4---------------------------------------------------------
cov<- regions_cov[[SAMP]]
meth<- regions_m[[SAMP]]

###################################
#####        Run PQLseq       #####
###################################
#Run PQLseq for males:withinage ------------------------------------------------
#Generate model matrix
predictor_matrix<- model.matrix(~ individual_sex:within.age + mean.age +  pid, data = long_data)
m_phenotype<- predictor_matrix[, 12]
m_covariates<- predictor_matrix[, c(2,11)]

#Run pqlseq model
male_age_pqlseq2<- pqlseq2(Y = meth, x = m_phenotype, 
                              K = kinship, W = m_covariates, 
                              lib_size = cov, model="BMM")

#Run PQLseq for females:withinage-----------------------------------------------
#Generate model matrix
predictor_matrix<- model.matrix(~ individual_sex:within.age + mean.age +  pid, data = long_data)
f_phenotype<- predictor_matrix[, 11]
f_covariates<- predictor_matrix[, c(2,12)]

#Run pqlseq model
female_age_pqlseq2<- pqlseq2(Y = meth, x = f_phenotype, 
                            K = kinship, W = f_covariates, 
                            lib_size = cov, model="BMM")

#Run PQLseq for mean_age-------------------------------------------------------------
#Generate model matrix
predictor_matrix<- model.matrix(~ individual_sex:within.age + mean.age + pid, data = long_data)
m.age_phenotype<- predictor_matrix[, 2]
m.age_covariates<- predictor_matrix[, c(11,12)]

#Run pqlseq model
mean_age_pqlseq2<- pqlseq2(Y = meth, x = m.age_phenotype, 
                              K = kinship, W = m.age_covariates, 
                              lib_size = cov, model="BMM")

#Save pqlseq models
saveRDS(male_age_pqlseq2, paste("male", "age", "pqlseq2", SAMP, sep = "_"))
saveRDS(female_age_pqlseq2, paste("female", "age", "pqlseq2",  SAMP, sep = "_"))
saveRDS(mean_age_pqlseq2, paste("mean", "age", "pqlseq2", SAMP, sep = "_"))

