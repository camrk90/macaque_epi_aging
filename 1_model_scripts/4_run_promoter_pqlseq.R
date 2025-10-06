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

#Load promoters-----------------------------------------------------------------
prom_cov<- readRDS("/scratch/ckelsey4/Cayo_meth/prom_cov_filtered")
prom_cov<- prom_cov[c(1:21)]
prom_m<- readRDS("/scratch/ckelsey4/Cayo_meth/prom_m_filtered")
prom_m<- prom_m[c(1:21)]

#Filter metadata to lids in regions list
long_data<- long_data[long_data$lid_pid %in% colnames(prom_cov[[1]]),]

prom_cov<- lapply(names(prom_cov), function(x){
  prom_cov<- subset(prom_cov[[x]], select=long_data$lid_pid)
  return(prom_cov)
})

prom_m<- lapply(names(prom_m), function(x){
  prom_m<- subset(prom_m[[x]], select=long_data$lid_pid)
  return(prom_m)
})

names(prom_cov)<- 1:21 #turn all chroms into integers (X = 21)
names(prom_m)<- 1:21 #turn all chroms into integers (X = 21)

#Check metadata lids match the lids (cols) of a random chromosome
if (all.equal(long_data$lid_pid, colnames(prom_cov[[runif(1, 1, 21)]]))) {
  
  ###################################
  #####        Prom PQLseq      #####
  ###################################
  #Model Vectors for genes pqlseq-------------------------------------------------
  cov<- prom_cov[[SAMP]]
  meth<- prom_m[[SAMP]]
  
  #Generate model matrix
  predictor_matrix<- model.matrix(~ within.age + mean.age + individual_sex + university, data = long_data)
  w.age_phenotype<- predictor_matrix[, 2]
  w.age_covariates<- predictor_matrix[, 3:4]
  
  #Run PQLseq for within_age ---------------------------------------------------
  prom_w_age<- pqlseq2(Y = meth, x = w.age_phenotype, 
                           K = kinship, W = w.age_covariates, 
                           lib_size = cov, model="BMM")
  
  #Run PQLseq for mean_age------------------------------------------------------
  m.age_phenotype<- predictor_matrix[, 3]
  m.age_covariates<- predictor_matrix[, c(2,4)]
  
  prom_m_age<- pqlseq2(Y = meth, x = m.age_phenotype, 
                           K = kinship, W = m.age_covariates, 
                           lib_size = cov, model="BMM")
  
  #Run PQLseq for sex-----------------------------------------------------------
  sex_phenotype<- predictor_matrix[, 4]
  sex_covariates<- predictor_matrix[, c(2:3)]
  
  prom_sex<- pqlseq2(Y = meth, x = sex_phenotype, 
                              K = kinship, W = sex_covariates, 
                              lib_size = cov, model="BMM")
  
  #Save pqlseq models
  saveRDS(prom_w_age, paste("prom", "pqlseq2", "agew", SAMP, sep = "_"))
  saveRDS(prom_m_age, paste("prom", "pqlseq2", "agem", SAMP, sep = "_"))
  saveRDS(prom_sex, paste("prom", "pqlseq2", "sex", SAMP, sep = "_"))
  
  ###################################
  #####      Nested PQLseq      #####
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
  saveRDS(male_age_pqlseq2, paste("prom", "nested", "m", SAMP, sep = "_"))
  saveRDS(female_age_pqlseq2, paste("prom", "nested", "f",  SAMP, sep = "_"))
  saveRDS(mean_age_pqlseq2, paste("prom", "nested", "mean", SAMP, sep = "_"))
  
} else {
  
  print("long_data lids did not match cov matrix lids")
  
}






