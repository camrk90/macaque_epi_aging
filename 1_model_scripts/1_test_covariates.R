#!/usr/bin/env /packages/apps/spack/18/opt/spack/gcc-11.2.0/r-4.2.2-kpl/bin/Rscript

#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ckelsey4@asu.edu
#SBATCH --mem=50G 
#SBATCH --array=1-21

SAMP <- Sys.getenv("SLURM_ARRAY_TASK_ID")
SAMP <- as.integer(SAMP)

library(tidyverse)
library(PQLseq2)
setwd("/home/ckelsey4/macaque_epi_aging/models_out/")

#Generate function--------------------------------------------------------------
run_pqlseq<- function(pheno, covariates){
  
  mod_df<- pqlseq2(Y = meth, x = pheno, 
                   K = kinship, W = covariates, 
                   lib_size = cov, model="BMM")
  
  mod_df<- mod_df %>%
    filter(converged == TRUE) %>%
    mutate(fdr = p.adjust(pvalue, method = "fdr")) %>%
    relocate(fdr, .after = pvalue) %>%
    dplyr::select(-c(converged, elapsed_time))
  
  return(mod_df)
  
}

#Import metadata----------------------------------------------------------------
long_data<- read.table("/scratch/ckelsey4/Cayo_meth/long_data_adjusted.txt")

long_data$pid_cov<- 0
long_data$pid_cov[long_data$pid == 'PID_10123']<- 1

#Import kinship matrix----------------------------------------------------------
kinship<- readRDS("/scratch/ckelsey4/Cayo_meth/full_kin_matrix")

#Subset and rearrange kinship rows and cols to match metadata
kinship<- kinship[long_data$lid_pid, long_data$lid_pid]

#Import m/cov rds------------------------------------------------------------
# load region lists that have been filtered for 5x coverage in 90% of samples
regions_cov<- readRDS("/scratch/ckelsey4/Cayo_meth/regions_cov_filtered.rds")
regions_m<- readRDS("/scratch/ckelsey4/Cayo_meth/regions_m_filtered.rds")

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
  
  #Model Vectors for lme4---------------------------------------------------------
  cov<- regions_cov[[SAMP]]
  meth<- regions_m[[SAMP]]
  
  #Select vars to model
  vars <- c("age_at_sampling:individual_sexF", "age_at_sampling:individual_sexM")
  
  #No batch effects model-------------------------------------------------------
  #Generate model matrix
  no_matrix<- model.matrix(~ age_at_sampling:individual_sex + mean.age, data = long_data)
  
  no_model <- lapply(setNames(vars, vars), function(i) {
    
    no_phenotype <- no_matrix[, i]
    no_covariates <- no_matrix[, setdiff(colnames(no_matrix), i)]
    
    run_pqlseq(no_phenotype, no_covariates)
    
  })
  
  #Save pqlseq model
  saveRDS(no_model, paste("dnam_no_model", SAMP, sep = "_"))
  
  #PID model--------------------------------------------------------------------
  #Generate model matrix
  pid_matrix<- model.matrix(~ age_at_sampling:individual_sex + mean.age + pid_cov, data = long_data)
  
  pid_model <- lapply(setNames(vars, vars), function(i) {
    
    pid_phenotype <- pid_matrix[, i]
    pid_covariates <- pid_matrix[, setdiff(colnames(pid_matrix), vars[1])]
    
    run_pqlseq(pid_phenotype, pid_covariates)
    
  })
  
  #Save pqlseq model
  saveRDS(pid_model, paste("dnam_pid_model", SAMP, sep = "_"))
  
  #Perc unique model------------------------------------------------------------
  #Generate model matrix
  perc_matrix<- model.matrix(~ age_at_sampling:individual_sex + mean.age + perc_unique, data = long_data)
  
  perc_model <- lapply(setNames(vars, vars), function(i) {
    
    perc_phenotype <- perc_matrix[, i]
    perc_covariates <- perc_matrix[, setdiff(colnames(perc_matrix), i)]
    
    run_pqlseq(perc_phenotype, perc_covariates)
    
  })
  
  #Save pqlseq model
  saveRDS(perc_model, paste("dnam_perc_model", SAMP, sep = "_"))
  
  #University model------------------------------------------------------------
  #Generate model matrix
  uni_matrix<- model.matrix(~ age_at_sampling:individual_sex + mean.age + university, data = long_data)
  
  uni_model <- lapply(setNames(vars, vars), function(i) {
    
    uni_phenotype <- uni_matrix[, i]
    uni_covariates <- uni_matrix[, setdiff(colnames(uni_matrix), i)]
    
    run_pqlseq(uni_phenotype, uni_covariates)
    
  })
  
  #Save pqlseq model
  saveRDS(uni_model, paste("dnam_uni_model", SAMP, sep = "_"))
  
} else {
  
  print("long_data lids did not match cov matrix lids")
  
}

