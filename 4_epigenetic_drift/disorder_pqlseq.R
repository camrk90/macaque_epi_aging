#!/usr/bin/env /packages/apps/spack/18/opt/spack/gcc-11.2.0/r-4.2.2-kpl/bin/Rscript

SAMP <- Sys.getenv("SLURM_ARRAY_TASK_ID")
SAMP <- as.integer(SAMP)

library(tidyverse)
library(PQLseq2)

#Import metadata----------------------------------------------------------------
long_data<- readRDS("/scratch/ckelsey4/Cayo_meth/long_data_adjusted")
long_data<- long_data %>%
  dplyr::rename(sex = individual_sex) %>%
  filter(n > 1) %>%
  arrange(lid_pid)

kinship<- readRDS("/scratch/ckelsey4/Cayo_meth/full_kin_matrix")

#Subset and rearrange kinship rows and cols to match metadata
kinship<- kinship[long_data$lid_pid, long_data$lid_pid]

#Import PDN and COV lists
avg_pdn<- readRDS("/scratch/ckelsey4/Cayo_meth/epigenetic_drift/avg_pdn/avg_pdn_list")

clean_list<- function(x){
  x<- x %>%
    mutate_if(is.character, as.numeric) %>%
    select(-chrom)
  x<- x[!rowSums(is.na(x)) >= 0.25*ncol(x),]
}

avg_pdn<- lapply(avg_pdn, clean_list)

read_cov<- readRDS("/scratch/ckelsey4/Cayo_meth/epigenetic_drift/avg_pdn/read_cov_list")
subset_cov<- function(x){
  df<- read_cov[[x]]
  df<- df %>%
    mutate_if(is.character, as.numeric) %>%
    select(-chrom)
  df<- df[rownames(df) %in% rownames(avg_pdn[[x]]),]
}

read_cov<- lapply(names(read_cov), subset_cov)

names(read_cov)<- names(avg_pdn)

#Multiply PDN values by cov as PQLseq weights the input by coverage
avg_pdn<- as.data.frame(do.call(rbind, avg_pdn))
read_cov<- as.data.frame(do.call(rbind, read_cov))
avg_pdn<- avg_pdn*read_cov

###################################
#####        Run PQLseq       #####
###################################

if (SAMP == 1) {
  
  #Run PQLseq for males:withinage ----------------------------------------------
  #Generate model matrix
  predictor_matrix<- model.matrix(~ sex:within.age + mean.age + sex, data = long_data)
  
  #Generate predictor and covariate matrices
  m_phenotype<- predictor_matrix[, 5]
  m_covariates<- predictor_matrix[, c(2:4)]
  
  #Run pqlseq model
  male_age_pqlseq2<- pqlseq2(Y = avg_pdn, x = m_phenotype, 
                             K = kinship, W = m_covariates, 
                             lib_size = read_cov, model="BMM", verbose=F)
  
  saveRDS(male_age_pqlseq2, 
          paste("/scratch/ckelsey4/Cayo_meth/epigenetic_drift/avg_pdn/", "pdn", "m", sep = "_"))
  
} else if (SAMP == 2) {
  
  #Run PQLseq for females:withinage---------------------------------------------
  #Generate model matrix
  predictor_matrix<- model.matrix(~ sex:within.age + mean.age + sex, data = long_data)
  
  #Generate model matrix
  f_phenotype<- predictor_matrix[, 4]
  f_covariates<- predictor_matrix[, c(2:3, 5)]
  
  #Run pqlseq model
  female_age_pqlseq2<- pqlseq2(Y = avg_pdn, x = f_phenotype, 
                               K = kinship, W = f_covariates, 
                               lib_size = read_cov, model="BMM")
  
  saveRDS(female_age_pqlseq2, 
          paste("/scratch/ckelsey4/Cayo_meth/epigenetic_drift/avg_pdn/", "pdn", "f", sep = "_"))
} else {
  
  #Run PQLseq for mean_age------------------------------------------------------
  #Generate model matrix
  predictor_matrix<- model.matrix(~ sex:within.age + mean.age + sex, data = long_data)
  
  #Generate model matrix
  m.age_phenotype<- predictor_matrix[, 2]
  m.age_covariates<- predictor_matrix[, c(3:5)]
  
  #Run pqlseq model
  mean_age_pqlseq2<- pqlseq2(Y = avg_pdn, x = m.age_phenotype, 
                             K = kinship, W = m.age_covariates, 
                             lib_size = read_cov, model="BMM")
  
  saveRDS(mean_age_pqlseq2, 
          paste("/scratch/ckelsey4/Cayo_meth/epigenetic_drift/avg_pdn/", "pdn", "mean", "age",  sep = "_"))
}



