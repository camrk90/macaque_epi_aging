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
#load genes
genes_cov<- readRDS("/scratch/ckelsey4/Cayo_meth/genes_cov_filtered")
genes_cov<- genes_cov[c(1:21)]
genes_m<- readRDS("/scratch/ckelsey4/Cayo_meth/genes_m_filtered")
genes_m<- genes_m[c(1:21)]

#Filter metadata to lids in regions list
long_data<- long_data[long_data$lid_pid %in% colnames(genes_cov[[1]]),]

genes_cov<- lapply(names(genes_cov), function(x){
  genes_cov<- subset(genes_cov[[x]], select=long_data$lid_pid)
  return(genes_cov)
})

genes_m<- lapply(names(genes_m), function(x){
  genes_m<- subset(genes_m[[x]], select=long_data$lid_pid)
  return(genes_m)
})

#Check metadata lids match the lids (cols) of a random chromosome
all.equal(long_data$lid_pid, colnames(genes_cov[[runif(1, 1, 21)]]))

names(genes_cov)<- 1:21 #turn all chroms into integers (X = 21)
names(genes_m)<- 1:21 #turn all chroms into integers (X = 21)

#Load promoters-----------------------------------------------------------------
prom_cov<- readRDS("/scratch/ckelsey4/Cayo_meth/prom_cov_filtered")
prom_cov<- prom_cov[c(1:21)]
prom_m<- readRDS("/scratch/ckelsey4/Cayo_meth/prom_m_filtered")
prom_m<- prom_m[c(1:21)]

prom_cov<- lapply(names(prom_cov), function(x){
  prom_cov<- subset(prom_cov[[x]], select=long_data$lid_pid)
  return(prom_cov)
})

prom_m<- lapply(names(prom_m), function(x){
  prom_m<- subset(prom_m[[x]], select=long_data$lid_pid)
  return(prom_m)
})

#Check metadata lids match the lids (cols) of a random chromosome
all.equal(long_data$lid_pid, colnames(prom_cov[[runif(1, 1, 21)]]))

names(prom_cov)<- 1:21 #turn all chroms into integers (X = 21)
names(prom_m)<- 1:21 #turn all chroms into integers (X = 21)

###################################
#####      Gene PQLseq        #####
###################################
#Model Vectors for genes pqlseq-------------------------------------------------
cov<- genes_cov[[SAMP]]
meth<- genes_m[[SAMP]]

#Run PQLseq for males:withinage ------------------------------------------------
#Generate model matrix
predictor_matrix<- model.matrix(~ individual_sex:within.age + mean.age +  pid, data = long_data)
m_phenotype<- predictor_matrix[, 12]
m_covariates<- predictor_matrix[, c(2,11)]

#Run pqlseq model
genes_m_pqlseq2<- pqlseq2(Y = meth, x = m_phenotype, 
                           K = kinship, W = m_covariates, 
                           lib_size = cov, model="BMM")

#Run PQLseq for females:withinage-----------------------------------------------
f_phenotype<- predictor_matrix[, 11]
f_covariates<- predictor_matrix[, c(2,12)]

#Run pqlseq model
genes_f_pqlseq2<- pqlseq2(Y = meth, x = f_phenotype, 
                             K = kinship, W = f_covariates, 
                             lib_size = cov, model="BMM")

#Run PQLseq for mean_age-------------------------------------------------------------
#Generate model matrix
mean_phenotype<- predictor_matrix[, 2]
mean_covariates<- predictor_matrix[, c(11,12)]

#Run pqlseq model
genes_mean_pqlseq2<- pqlseq2(Y = meth, x = mean_phenotype, 
                           K = kinship, W = mean_covariates, 
                           lib_size = cov, model="BMM")

#Save pqlseq models
saveRDS(genes_m_pqlseq2, paste("gene", "pqlseq2", "m", SAMP, sep = "_"))
saveRDS(genes_f_pqlseq2, paste("gene", "pqlseq2", "f", SAMP, sep = "_"))
saveRDS(genes_mean_pqlseq2, paste("gene", "pqlseq2", "mean", SAMP, sep = "_"))

###################################
#####        Prom PQLseq      #####
###################################
#Model Vectors for genes pqlseq-------------------------------------------------
cov<- prom_cov[[SAMP]]
meth<- prom_m[[SAMP]]

#Run PQLseq for males:withinage ------------------------------------------------
prom_m_pqlseq2<- pqlseq2(Y = meth, x = m_phenotype, 
                          K = kinship, W = m_covariates, 
                          lib_size = cov, model="BMM")

#Run PQLseq for females:withinage-----------------------------------------------
prom_f_pqlseq2<- pqlseq2(Y = meth, x = f_phenotype, 
                          K = kinship, W = f_covariates, 
                          lib_size = cov, model="BMM")

#Run PQLseq for mean_age--------------------------------------------------------
prom_mean_pqlseq2<- pqlseq2(Y = meth, x = mean_phenotype, 
                             K = kinship, W = mean_covariates, 
                             lib_size = cov, model="BMM")

#Save pqlseq models
saveRDS(prom_m_pqlseq2, paste("prom", "pqlseq2", "m", SAMP, sep = "_"))
saveRDS(prom_f_pqlseq2, paste("prom", "pqlseq2", "f", SAMP, sep = "_"))
saveRDS(prom_mean_pqlseq2, paste("prom", "pqlseq2", "mean", SAMP, sep = "_"))

