library(tidyverse)
library(PQLseq2)
setwd("/scratch/ckelsey4/Cayo_meth/glmer_model_compare")

#Import metadata----------------------------------------------------------------
long_data<- readRDS("/scratch/ckelsey4/Cayo_meth/long_data_adjusted")

long_data<- long_data %>%
  arrange(lid_pid) %>%
  filter(age_at_sampling > 1)

#Import kinship matrix----------------------------------------------------------
kinship<- readRDS("/scratch/ckelsey4/Cayo_meth/full_kin_matrix")

#Import m/cov rds------------------------------------------------------------
# load region lists that have been filtered for 5x coverage in 90% of samples
cov<- readRDS("/scratch/ckelsey4/Cayo_meth/regions_cov_filtered")
cov<- cov[["Y"]]
meth<- readRDS("/scratch/ckelsey4/Cayo_meth/regions_m_filtered")
meth<- meth[["Y"]]

#Filter metadata to lids in regions list
long_data<- long_data[long_data$lid_pid %in% colnames(cov),]

cov<- subset(cov, select=long_data$lid_pid)
meth<- subset(meth, select=long_data$lid_pid)

#Check metadata lids match the lids (cols) of a random chromosome
all.equal(long_data$lid_pid, colnames(cov))

#Subset and rearrange kinship rows and cols to match metadata
kinship<- kinship[long_data$lid_pid, long_data$lid_pid]

###################################
#####   Run Additive PQLseq   #####
###################################
#Run PQLseq for within_age-------------------------------------------------------------
#Generate model matrix
predictor_matrix<- model.matrix(~ within.age + mean.age + university, data = long_data)
w.age_phenotype<- predictor_matrix[, 2]
w.age_covariates<- predictor_matrix[, 3:4]

#Run pqlseq model
w.age_pqlseq2_model<- pqlseq2(Y = meth, x = w.age_phenotype, 
                              K = kinship, W = w.age_covariates, 
                              lib_size = cov, model="BMM")

#Run PQLseq for mean_age-------------------------------------------------------------
#Generate model matrix
m.age_phenotype<- predictor_matrix[, 3]
m.age_covariates<- predictor_matrix[, c(2,4)]

#Run pqlseq model
m.age_pqlseq2_model<- pqlseq2(Y = meth, x = m.age_phenotype, 
                              K = kinship, W = m.age_covariates, 
                              lib_size = cov, model="BMM")

#Save pqlseq models
saveRDS(w.age_pqlseq2_model, paste("wb", "pqlseq2", "within", "age", "Y", sep = "_"))
saveRDS(m.age_pqlseq2_model, paste("wb", "pqlseq2", "mean", "age", "Y", sep = "_"))

