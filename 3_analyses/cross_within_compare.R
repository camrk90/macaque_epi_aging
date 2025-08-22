library(tidyverse)
library(ggplot2)

#Define import function
import_pqlseq<- function(x, y){

    #Generate list of file names
    file_list<- list.files(pattern = x)
    file_order<- str_split_i(file_list, "_", y)
    
    #Import glm models as list
    model_list<- lapply(file_list, readRDS)
    
    #Rename list elements
    names(model_list)<- file_order
    model_list<- model_list[1:21]
    
    #Bind model list to df and add rownames
    model<- do.call(rbind, model_list)
    model$outcome<- str_split_i(model$outcome, "\\.", 3)
    model$outcome2<- model$outcome
    
    #Separate region coordinates into start and end, delete the chr col, and move region col to front
    model<- model %>% 
      separate_wider_delim(outcome2, names=c("chr", "chromStart", "chromEnd"), delim = "_") %>%
      relocate(c(chr, chromStart, chromEnd), .after = outcome)
    
    #Add length col and filter by length
    model<- model %>%
      mutate(length = 1+(as.numeric(chromEnd) - as.numeric(chromStart))) %>%
      relocate(length, .after=outcome)
    
    #Filter for true convergences
    model<- model %>%
      filter(converged == "TRUE")
    
    #Generate df of adjusted pvalues
    model_fdr<- p.adjust(model$pvalue, method = "fdr")
    
    #Bind padj cols to model df and relocate
    model<- cbind(model, model_fdr)
    model<- model %>%
      dplyr::rename(fdr = model_fdr) %>%
      relocate(fdr, .after = pvalue) %>%
      dplyr::select(-c(elapsed_time, converged, h2, sigma2))
}

#Import metadata
blood_metadata<- read.table("/scratch/ckelsey4/Cayo_meth/blood_metadata_full.txt")
long_data<- read.table("/scratch/ckelsey4/Cayo_meth/long_data_adjusted.txt")

blood_metadata<- blood_metadata[!blood_metadata$lid_pid %in% long_data$lid_pid,]

blood_metadata<- blood_metadata %>%
  drop_na(age_at_sampling)

long_data<- long_data %>%
  dplyr::rename(sex = individual_sex) %>%
  group_by(monkey_id) %>%
  mutate(n = n()) %>%
  filter(n > 1) %>%
  ungroup()

######################################
###      Import Within Models      ###
######################################
setwd('/scratch/ckelsey4/Cayo_meth/glmer_model_compare')

#Import longitudinal pqlseq files-----------------------------------------------
#Chronological Age
chron_age_files<- 'wb_pqlseq2_agechron'
chron_age_pqlseq<- import_pqlseq(chron_age_files, y = 4)

#Age Within
age_w_files<- 'wb_pqlseq2_within_age'
age_w_pqlseq<- import_pqlseq(age_w_files, y = 5)

#Mean Age
age_m_files<- 'wb_pqlseq2_mean_age'
age_m_pqlseq<- import_pqlseq(age_m_files, y = 5)

#Rename cols for each df to indicate variable
colnames(age_w_pqlseq)<- c("outcome", "length", "chr", "chromStart", "chromEnd", "n", 
                           paste(names(age_w_pqlseq[,7:12]), "age", sep = "_"))
colnames(age_m_pqlseq)<- c(paste(names(age_m_pqlseq), "mean_age", sep = "_"))
colnames(chron_age_pqlseq)<- c(paste(names(chron_age_pqlseq), "chron_age", sep = "_"))

#Cbind cols for age and sex dfs
pqlseq_model<- cbind(age_w_pqlseq, age_m_pqlseq[,7:12], chron_age_pqlseq[,7:12])

# Import Cross_Se Models--------------------------------------------------------
setwd('/scratch/ckelsey4/Cayo_meth/cross_models')

cross_files<- 'wb_cs'
cross_pqlseq<- import_pqlseq(cross_files, y = 3)

pqlseq_model<- pqlseq_model[pqlseq_model$outcome %in% cross_pqlseq$outcome,]

#Sort chromosome factors
sorted_labels<- str_sort(unique(pqlseq_model$chr), numeric=T)

pqlseq_model<- pqlseq_model %>%
  mutate(chr = factor(chr, levels = sorted_labels)) %>%
  arrange(chr)

rm(age_w_pqlseq);rm(age_m_pqlseq);rm(chron_age_pqlseq)

######################################
###     Import Cross_Se Models     ###
######################################
setwd('/scratch/ckelsey4/Cayo_meth/cross_models')

cross_files<- 'wb_cs'
cross_pqlseq<- import_pqlseq(cross_files, y = 3)


