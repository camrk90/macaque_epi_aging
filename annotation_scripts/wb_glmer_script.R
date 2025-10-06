#!/usr/bin/env /packages/apps/spack/18/opt/spack/gcc-11.2.0/r-4.2.2-kpl/bin/Rscript

SAMP <- Sys.getenv("SLURM_ARRAY_TASK_ID")
SAMP <- as.integer(SAMP)

library(tidyverse)
library(lme4)
library(broom)
setwd("/scratch/ckelsey4/Cayo_meth/glmer_model_compare")

#Import metadata----------------------------------------------------------------
long_data<- readRDS("/scratch/ckelsey4/Cayo_meth/long_data_adjusted")

#Import m/cov rds------------------------------------------------------------
# load region lists that have been filtered for 5x coverage in 90% of samples
regions_cov_list<- readRDS("/scratch/ckelsey4/Cayo_meth/regions_cov_list")
regions_m_list<- readRDS("/scratch/ckelsey4/Cayo_meth/regions_m_list")

#Subset regions list for lids with repeated measures
regions_cov<- lapply(names(regions_cov_list), function(x){
  regions_cov<- subset(regions_cov_list[[x]], select=long_data$lid_pid)
  return(regions_cov)
})

regions_m<- lapply(names(regions_m_list), function(x){
  regions_m<- subset(regions_m_list[[x]], select=long_data$lid_pid)
  return(regions_m)
})

names(regions_cov)<- 1:21 #turn all chroms into integers (X = 21)
names(regions_m)<- 1:21 #turn all chroms into integers (X = 21)


#Vectors for lme4---------------------------------------------------------------
cov<- t(regions_cov[[SAMP]])
meth<- t(regions_m[[SAMP]])
ratio<- meth/cov
ratio[is.nan(ratio)]<- 0
age.w<- long_data$within.age
age.m<- long_data$mean.age
age<- long_data$age_at_sampling
ids<- long_data$monkey_id
sex<- long_data$sex
#batch<- long_data$university

#Run glm with random intercept--------------------------------------------------
#Generate blank objects for glms
df<- data.frame(matrix(nrow=0, ncol=20))
colnames(df)<- c("estimate_(Intercept)", "estimate_age.w", "estimate_age.m", "estimate_sexM",
                 "se_(Intercept)", "se_age.w", "se_age.m", "se_sexM", 
                 "z_(Intercept)", "z_age.w", "z_age.m", "z_sexM",
                 "pval_(Intercept)", "pval_age.w", "pval_age.m", "pval_sexM",
                 "aic", "bic", "loglik", "converged")

model_list_intercept<- list()
#sum_list_intercept<- list()
#ranef_list_intercept<- list()

#Run glm for regions within chromosomes
for(i in 1:ncol(ratio)){
  
  rfx=glmer(ratio[, i] ~ age.w + age.m + sex + (1|ids), 
            family = binomial(link = "logit"), 
            weights = cov[, i])
  
  slopes<- as.data.frame(ranef(rfx))
  rfx_sum<- summary(rfx)
  cfs<- rfx_sum[["coefficients"]]
  cfs<- as.data.frame(cfs)
  cfs$term<- rownames(cfs)
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
  model_list_intercept[[length(model_list_intercept)+1]] = rfx
  #sum_list[[length(sum_list)+1]] = rfx_sum
  #ranef_list[[length(ranef_list)+1]] = slopes
}

rownames(df)<- colnames(ratio)
names(model_list_intercept)<- colnames(ratio)
#names(sum_list)<- colnames(ratio)
#names(ranef_list)<- colnames(ratio)

#Save outputs to RDS
saveRDS(df, paste("glmer", "intercept", "output", SAMP, sep = "_"))
#saveRDS(model_list_intercept, paste("glmer", "intercept", "models", SAMP, sep = "_"))
#saveRDS(sum_list, paste("glmer", "intercept", "summaries", SAMP, sep = "_"))
#saveRDS(ranef_list, paste("glmer", "intercept", "ranefx", SAMP, sep = "_"))

#Run glm with random intercept and slope----------------------------------------
#Generate blank objects for glms
df<- data.frame(matrix(nrow=0, ncol=20))
colnames(df)<- c("estimate_(Intercept)", "estimate_age.w", "estimate_age.m", "estimate_sexM",
                 "se_(Intercept)", "se_age.w", "se_age.m", "se_sexM", 
                 "z_(Intercept)", "z_age.w", "z_age.m", "z_sexM",
                 "pval_(Intercept)", "pval_age.w", "pval_age.m", "pval_sexM",
                 "aic", "bic", "loglik", "converged")

model_list_both<- list()
#sum_list<- list()
#ranef_list<- list()

#Run glm for regions within chromosomes
for(i in 1:ncol(ratio)){
  
  rfx=glmer(ratio[, i] ~ age.w + age.m + sex + (1 + age.w|ids), 
            family = binomial(link = "logit"), 
            weights = cov[, i])
  
  slopes<- as.data.frame(ranef(rfx))
  rfx_sum<- summary(rfx)
  cfs<- rfx_sum[["coefficients"]]
  cfs<- as.data.frame(cfs)
  cfs$term<- rownames(cfs)
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
  #sum_list[[length(sum_list)+1]] = rfx_sum
  #ranef_list[[length(ranef_list)+1]] = slopes
}

rownames(df)<- colnames(ratio)
names(model_list_both)<- colnames(ratio)
#names(sum_list)<- colnames(ratio)
#names(ranef_list)<- colnames(ratio)

#Save outputs to RDS
saveRDS(df, paste("glmer", "both", "output", SAMP, sep = "_"))
#saveRDS(model_list_both, paste("glmer", "both",  "models", SAMP, sep = "_"))
#saveRDS(sum_list, paste("glmer", "both", "summaries", SAMP, sep = "_"))
#saveRDS(ranef_list, paste("glmer", "both", "ranefx", SAMP, sep = "_"))

#Run anova to determine which model fits best-----------------------------------
#Generate blank list to save output
model_compare<- list()
pvals<- data.frame(matrix(nrow=0, ncol=1))
colnames(pvals)<- c("pval")

#Run for loop comparing models
for(i in 1:length(model_list_intercept)){
  test<- anova(model_list_intercept[[i]], model_list_both[[i]], 
                        test = "LTR")
  
  test_tidy<- broom::tidy(test)
  
  p<- test_tidy[2, 9]
  
  pvals<- rbind(pvals, p)
  model_compare[[length(model_compare)+1]] = test_tidy
}

colnames(pvals)<- "pval"
rownames(pvals)<- names(model_list_intercept)
names(model_compare)<- names(model_list_intercept)

saveRDS(model_compare, paste("model", "compare", SAMP, sep = "_"))
saveRDS(pvals, paste("model", "compare", "pvals", SAMP, sep = "_"))

