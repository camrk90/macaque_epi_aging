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
#Run PQLseq for within_age-------------------------------------------------------------
#Generate model matrix
predictor_matrix<- model.matrix(~ within.age + mean.age + individual_sex + pid, data = long_data)
w.age_phenotype<- predictor_matrix[, 2]
w.age_covariates<- predictor_matrix[, 3:12]

#Run pqlseq model
w.age_pqlseq2_model<- pqlseq2(Y = meth, x = w.age_phenotype, 
                              K = kinship, W = w.age_covariates, 
                              lib_size = cov, model="BMM")

#Run PQLseq for mean_age-------------------------------------------------------------
#Generate model matrix
predictor_matrix<- model.matrix(~ within.age + mean.age + individual_sex + pid, data = long_data)
m.age_phenotype<- predictor_matrix[, 3]
m.age_covariates<- predictor_matrix[, c(2,4)]

#Run pqlseq model
m.age_pqlseq2_model<- pqlseq2(Y = meth, x = m.age_phenotype, 
                              K = kinship, W = m.age_covariates, 
                              lib_size = cov, model="BMM")

#Run PQLseq for sex-------------------------------------------------------------
#Generate model matrix
predictor_matrix<- model.matrix(~ within.age + mean.age + individual_sex + pid, data = long_data)
sex_phenotype<- predictor_matrix[, 4]
sex_covariates<- predictor_matrix[, 2:3]

#Run pqlseq model
sex_pqlseq2_model<- pqlseq2(Y = meth, x = sex_phenotype, 
                            K = kinship, W = sex_covariates, 
                            lib_size = cov, model="BMM")

#Save pqlseq models
saveRDS(w.age_pqlseq2_model, paste("wb", "pqlseq2", "within", "age", SAMP, sep = "_"))
saveRDS(m.age_pqlseq2_model, paste("wb", "pqlseq2", "mean", "age", SAMP, sep = "_"))
saveRDS(sex_pqlseq2_model, paste("wb", "pqlseq2", "sex", SAMP, sep = "_"))

