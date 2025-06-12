#!/usr/bin/env /packages/apps/spack/18/opt/spack/gcc-11.2.0/r-4.2.2-kpl/bin/Rscript

SAMP <- Sys.getenv("SLURM_ARRAY_TASK_ID")
SAMP <- as.integer(SAMP)

library(tidyverse)
library(PQLseq2)
library(lme4)
library(broom)
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
regions_m<- readRDS("/scratch/ckelsey4/Cayo_meth/regions_m_filtered")

#Filter metadata to lids in regions list
long_data<- long_data[long_data$lid_pid %in% colnames(regions_cov[[1]]),]

#Check metadata lids match the lids (cols) of a random chromosome
test<- all.equal(long_data$lid_pid, colnames(regions_cov[[runif(1, 1, 21)]]))

#Subset regions list for lids n > 1
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

#Model Vectors for lme4---------------------------------------------------------
cov<- t(regions_cov[[SAMP]])
meth<- t(regions_m[[SAMP]])
age.w<- long_data$within.age
age.m<- long_data$mean.age
ids<- long_data$monkey_id
sex<- long_data$individual_sex
#batch<- long_data$pid

###################################
#####        Run GLMM         #####
###################################
#Run GLMMM----------------------------------------------------------------------
#Generate blank objects for glms
df<- data.frame(matrix(nrow=0, ncol=21))
colnames(df)<- c("estimate_(Intercept)", "estimate_age.w", "estimate_age.m", "estimate_sexM",
                 "se_(Intercept)", "se_age.w", "se_age.m", "se_sexM", 
                 "z_(Intercept)", "z_age.w", "z_age.m", "z_sexM",
                 "pval_(Intercept)", "pval_age.w", "pval_age.m", "pval_sexM",
                 "aic", "bic", "loglik", "converged")

model_list_both<- list()
sum_list<- list()
ranef_list<- list()

#Run slope/intercept glm
for(i in 1:ncol(meth)){
  
  rfx=glmer(cbind(meth[, i], cov[, i]) ~ age.w + age.m + sex + (1 + age.w|ids), 
            family = binomial(link = "logit"))
  
  slopes<- as.data.frame(ranef(rfx))
  rfx_sum<- summary(rfx)
  cfs<- rfx_sum[["coefficients"]]
  cfs<- as.data.frame(cfs)
  cfs$term<- rownames(cfs)
  #cfs<- cfs[1:4, 1:5] #this is for when the batch variable gets incorporated
  colnames(cfs)<- c("estimate", "se", "z", "pval", "term")
  rfx_row<- pivot_wider(cfs, names_from = term, values_from = c(estimate, se, z, pval))
  rfx_row$aic<- rfx_sum[["AICtab"]][["AIC"]]
  rfx_row$bic<- rfx_sum[["AICtab"]][["BIC"]]
  rfx_row$loglik<- rfx_sum[["logLik"]]
  
  if(length(rfx@optinfo[["conv"]][["lme4"]]) == 2){
    rfx_row$converged<- "FALSE"
  }else{
    rfx_row$converged<- "TRUE"
  }
  
  df<- rbind(df, rfx_row)
  model_list_both[[length(model_list_both)+1]] = rfx
  sum_list[[length(sum_list)+1]] = rfx_sum
  ranef_list[[length(ranef_list)+1]] = slopes
}

df$region<- colnames(meth)
names(model_list_both)<- colnames(meth)
names(sum_list)<- colnames(meth)
names(ranef_list)<- colnames(meth)

#Save outputs to RDS
saveRDS(df, paste("glmer", "output", SAMP, sep = "_"))
saveRDS(model_list_both, paste("glmer",  "models", SAMP, sep = "_"))
saveRDS(sum_list, paste("glmer", "summaries", SAMP, sep = "_"))
saveRDS(ranef_list, paste("glmer", "ranefx", SAMP, sep = "_"))

###################################
#####        Run PQLseq       #####
###################################
#Run PQLseq for within_age-------------------------------------------------------------
#Generate model matrix
predictor_matrix<- model.matrix(~ within.age + mean.age + individual_sex + pid, data = long_data)
w.age_phenotype<- predictor_matrix[, 2]
w.age_covariates<- predictor_matrix[, 3:4]

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

