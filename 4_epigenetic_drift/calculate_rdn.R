library(tidyverse)
library(stringr)
library(readr)

#Import kinship matrix----------------------------------------------------------
kinship<- readRDS("/scratch/ckelsey4/Cayo_meth/full_kin_matrix")

#Subset and rearrange kinship rows and cols to match metadata
kinship<- kinship[long_data$lid_pid, long_data$lid_pid]

#Import metadata----------------------------------------------------------------
long_data<- readRDS("/scratch/ckelsey4/Cayo_meth/long_data_adjusted")
long_data<- long_data %>%
  dplyr::rename(sex = individual_sex) %>%
  filter(n > 1) %>%
  arrange(lid_pid)

lid_pid<- long_data$lid_pid
lid_pid<- paste(lid_pid, "_meanPDN.bed", sep = "")

##set unique file ID's
unique_fileID <- list.files(path = "/scratch/ckelsey4/Cayo_meth/epigenetic_drift/avg_pdn/",
                            pattern = "*_meanPDN.bed")

unique_fileID<- unique_fileID[unique_fileID %in% lid_pid]

##set file paths
filePaths <- paste0("/scratch/ckelsey4/Cayo_meth/epigenetic_drift/avg_pdn/", unique_fileID)

##loop to calculate PDN
avg_pdn<- data.frame(matrix(nrow=190910, ncol=0))
read_cov<- data.frame(matrix(nrow=190910, ncol=0))
overlap<- read.table(filePaths[1], header = FALSE, fill = TRUE)

for(i in 1:5){
  if (i == 1){
    print(i)
    #Read in file
    overlap<- read.table(filePaths[i], header = FALSE, fill = TRUE)
    overlap<- overlap %>%
      dplyr::rename(chrom = V1,
                    start = V2,
                    end = V3)
    overlap$region<- paste(overlap$chrom, overlap$start, overlap$end, sep = "_")
    overlap<- overlap %>%
      relocate(region, .before = chrom)
    
    pdn<- overlap[c(1, 6)]
    cov<- overlap[c(1,5)]
    
    #Generate unique id
    id<- paste(unique_fileID[i])
    id<- gsub("_meanPDN.bed", "", id)
    
    #Rename PDN col with unique id
    pdn<- dplyr::rename_with(pdn, ~ paste0(id), V5)
    cov<- dplyr::rename_with(cov, ~ paste0(id), V4)
    
    #Bind to df
    avg_pdn<- cbind(avg_pdn, pdn)
    read_cov<- cbind(read_cov, cov)
    
  } else {
    #Read in file
    overlap<- read.table(filePaths[i], header = FALSE, fill = TRUE)
    pdn<- overlap[5]
    cov<- overlap[4]
    
    #Generate unique id
    id<- paste(unique_fileID[i])
    id<- gsub("_meanPDN.bed", "", id)
    
    #Rename PDN col with unique id
    pdn<- dplyr::rename_with(pdn, ~ paste0(id), V5)
    cov<- dplyr::rename_with(cov, ~ paste0(id), V4)
    
    #Bind to df
    avg_pdn<- cbind(avg_pdn, pdn)
    read_cov<- cbind(read_cov, cov)
    print(paste(i, "done!", sep = " "))
  }
}

#Convert periods to NA
avg_pdn[avg_pdn == "."]<- NA

##Filter NAs
#avg_pdn<- avg_pdn[!rowSums(is.na(avg_pdn)) > ncol(avg_pdn)*.2,] #this filters out rows with more than 20% NAs
avg_pdn<- drop_na(avg_pdn) #lme4 drops rows with NAs anyway so this code just removes them ahead of time
read_cov<- read_cov[read_cov$region %in% avg_pdn$region,]

rownames(avg_pdn)<- avg_pdn$region
avg_pdn<- avg_pdn %>%
  select(-region) %>%
  t()

rownames(read_cov)<- read_cov$region
read_cov<- read_cov %>%
  select(-region) %>%
  t()

write.csv(compiled, file="/scratch/ckelsey4/Cayo_meth/epigenetic_drift/mean_pdn_all.csv")

#Model Vectors for lme4---------------------------------------------------------
cov<- regions_cov[[SAMP]]
meth<- regions_m[[SAMP]]

###################################
#####        Run PQLseq       #####
###################################
#Run PQLseq for within_age-------------------------------------------------------------
#Generate model matrix
predictor_matrix<- model.matrix(~ within.age + mean.age + individual_sex + university, data = long_data)
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

#Run PQLseq for sex-------------------------------------------------------------
#Generate model matrix
predictor_matrix<- model.matrix(~ within.age + mean.age + individual_sex + university, data = long_data)
sex_phenotype<- predictor_matrix[, 4]
sex_covariates<- predictor_matrix[, c(2:3)]

#Run pqlseq model
sex_pqlseq2_model<- pqlseq2(Y = meth, x = sex_phenotype, 
                            K = kinship, W = sex_covariates, 
                            lib_size = cov, model="BMM")



