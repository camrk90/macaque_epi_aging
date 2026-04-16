#!/usr/bin/env /packages/apps/spack/18/opt/spack/gcc-11.2.0/r-4.2.2-kpl/bin/Rscript

#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ckelsey4@asu.edu
#SBATCH --mem=50G 
#SBATCH --array=1-21

SAMP <- Sys.getenv("SLURM_ARRAY_TASK_ID")
SAMP <- as.integer(SAMP)

library(tidyverse)
library(PQLseq2)
setwd("/scratch/ckelsey4/Cayo_meth/glmer_model_compare")

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

#Import kinship matrix----------------------------------------------------------
kinship<- readRDS("/scratch/ckelsey4/Cayo_meth/full_kin_matrix")

#Subset and rearrange kinship rows and cols to match metadata
kinship<- kinship[long_data$lid_pid, long_data$lid_pid]

#Import m/cov rds------------------------------------------------------------
# load region lists that have been filtered for 5x coverage in 90% of samples
regions_cov<- readRDS("/scratch/ckelsey4/Cayo_meth/regions_cov_filtered2.rds")
regions_m<- readRDS("/scratch/ckelsey4/Cayo_meth/regions_m_filtered2.rds")

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
  
  #Model Vectors for lme4-------------------------------------------------------
  cov<- regions_cov[[SAMP]]
  meth<- regions_m[[SAMP]]
  
  ###################################
  #####        Run PQLseq       #####
  ###################################
  #Additive model---------------------------------------------------------------
  #Generate model matrix
  additive_matrix<- model.matrix(~ age_at_sampling + mean.age + individual_sex + university, data = long_data)
  
  vars <- c("age_at_sampling", "individual_sexM")
  
  additive_model <- lapply(setNames(vars, vars), function(i) {
    
    additive_phenotype <- additive_matrix[, i]
    additive_covariates <- additive_matrix[, setdiff(colnames(additive_matrix), i)]
    
    run_pqlseq(additive_phenotype, additive_covariates)
    
  })
  
  #Save pqlseq model
  saveRDS(additive_model, paste("dnam_additive_model", SAMP, sep = "_"))
  
  #Eq.2 Nested model-----------------------------------------------------------------
  #Generate model matrix
  nested_matrix<- model.matrix(~ age_at_sampling:individual_sex + mean.age + university, data = long_data)
  
  vars <- c("age_at_sampling:individual_sexF", "age_at_sampling:individual_sexM")
  
  nested_model <- lapply(setNames(vars, vars), function(i) {
    
    nested_phenotype <- nested_matrix[, i]
    nested_covariates <- nested_matrix[, setdiff(colnames(nested_matrix), i)]
    
    run_pqlseq(nested_phenotype, nested_covariates)
    
  })
  
  #Save pqlseq model
  saveRDS(nested_model, paste("dnam_nested_model", SAMP, sep = "_"))
  
  #Eq.1 Nested model------------------------------------------------------------
  #Generate model matrix
  eq1_nested_matrix<- model.matrix(~ age_at_sampling:individual_sex + university, data = long_data)
  
  vars <- c("age_at_sampling:individual_sexF", "age_at_sampling:individual_sexM")
  
  eq1_nested_model <- lapply(setNames(vars, vars), function(i) {
    
  eq1_nested_phenotype <- eq1_nested_matrix[, i]
  eq1_nested_covariates <- eq1_nested_matrix[, setdiff(colnames(eq1_nested_matrix), i)]
    
    run_pqlseq(eq1_nested_phenotype, eq1_nested_covariates)
    
  })
  
  #Save pqlseq model
  saveRDS(eq1_nested_model, paste("dnam_nestedEq1_model", SAMP, sep = "_"))
  
} else {
  
  print("long_data lids did not match cov matrix lids")
  
}



