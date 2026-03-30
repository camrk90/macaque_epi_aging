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

long_data<- long_data %>%
  group_by(monkey_id) %>%
  mutate(n = n()) %>%
  ungroup()

long_data<- long_data %>%
  filter(age_at_sampling > 1) %>%
  filter(n > 1) %>%
  dplyr::rename(perc_unique = unique) %>%
  drop_na() %>%
  arrange(lid_pid)

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
  
  ###################################
  #####        Run PQLseq       #####
  ###################################
  #Age*sex interaction----------------------------------------------------------
  #Generate model matrix
  interaction_matrix<- model.matrix(~ age_at_sampling*mean.age + individual_sex + perc_unique, data = long_data)
  
  vars <- c("age_at_sampling", "individual_sexM", "age_at_sampling:mean.age")
  
  interaction_model <- lapply(setNames(vars, vars), function(i) {
    
    interaction_phenotype <- interaction_matrix[, i]
    interaction_covariates <- interaction_matrix[, setdiff(colnames(interaction_matrix), i)]
    
    run_pqlseq(interaction_phenotype, interaction_covariates)
    
  })
  
  #Save pqlseq model
  saveRDS(interaction_model, paste("dnam_interaction_model", SAMP, sep = "_"))
  
  #Nested model-----------------------------------------------------------------
  #Generate model matrix
  nested_matrix<- model.matrix(~ age_at_sampling:individual_sex + mean.age + perc_unique, data = long_data)
  
  vars <- c("age_at_sampling:individual_sexF", "age_at_sampling:individual_sexM")
  
  nested_model <- lapply(setNames(vars, vars), function(i) {
    
    nested_phenotype <- nested_matrix[, i]
    nested_covariates <- nested_matrix[, setdiff(colnames(nested_matrix), i)]
    
    run_pqlseq(nested_phenotype, nested_covariates)
    
  })
  
  #Save pqlseq model
  saveRDS(nested_model, paste("dnam_nested_model", SAMP, sep = "_"))
  
} else {
  
  print("long_data lids did not match cov matrix lids")
  
}



