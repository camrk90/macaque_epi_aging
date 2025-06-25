#!/usr/bin/env /packages/apps/spack/18/opt/spack/gcc-11.2.0/r-4.2.2-kpl/bin/Rscript

library(tidyverse)
library(stringr)
library(readr)

##set unique file ID's
unique_fileID <- list.files(path = "/scratch/ckelsey4/Cayo_meth/epigenetic_drift/avg_pdn/",
                            pattern = "*_meanPDN.bed")

##set file paths
filePaths <- paste0("/scratch/ckelsey4/Cayo_meth/epigenetic_drift/avg_pdn/", unique_fileID)

##loop to calculate PDN
avg_pdn<- data.frame(matrix(nrow=190910, ncol=0))
read_cov<- data.frame(matrix(nrow=190910, ncol=0))

for(i in 1:length(unique_fileID)){
  if (i == 1) {
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
    
    pdn<- overlap[c(1:2, 6)]
    cov<- overlap[c(1:2,5)]
    
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
  rm(overlap)
  rm(cov)
  rm(pdn)
}

#Convert periods to NA
avg_pdn[avg_pdn == "."]<- NA

##Filter NAs
#avg_pdn<- avg_pdn[!rowSums(is.na(avg_pdn)) > ncol(avg_pdn)*.2,] #this filters out rows with more than 20% NAs
avg_pdn<- drop_na(avg_pdn) #lme4 drops rows with NAs anyway so this code just removes them ahead of time
read_cov<- read_cov[read_cov$region %in% avg_pdn$region,]

chrs<- str_sort(c(1:20, "X"), numeric = TRUE)

rownames(avg_pdn)<- avg_pdn$region
avg_pdn<- avg_pdn %>%
  mutate(chrom = factor(chrom, levels = chrs)) %>%
  arrange(chrom) %>%
  select(-region)
avg_pdn_list<- split(avg_pdn, avg_pdn$chrom)

rownames(read_cov)<- read_cov$region
read_cov<- read_cov %>%
  mutate(chrom = factor(chrom, levels = chrs)) %>%
  arrange(chrom) %>%
  select(-region)
read_cov_list<- split(read_cov, read_cov$chrom)

saveRDS(avg_pdn_list, "/scratch/ckelsey4/Cayo_meth/epigenetic_drift/avg_pdn/avg_pdn_list")
saveRDS(read_cov_list, "/scratch/ckelsey4/Cayo_meth/epigenetic_drift/avg_pdn/read_cov_list")




