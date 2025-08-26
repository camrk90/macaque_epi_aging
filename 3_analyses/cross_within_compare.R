library(tidyverse)
library(ggplot2)
library(ggcorrplot)

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
  arrange(lid_pid) %>%
  filter(age_at_sampling > 1) %>%
  dplyr::rename(perc_unique = unique) %>%
  drop_na()

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
                           paste(names(age_w_pqlseq[,7:12]), "within", "age", sep = "_"))
colnames(age_m_pqlseq)<- c(paste(names(age_m_pqlseq), "mean_age", sep = "_"))
colnames(chron_age_pqlseq)<- c(paste(names(chron_age_pqlseq), "chron_age", sep = "_"))

#Cbind cols for age and sex dfs
age_w_pqlseq<- age_w_pqlseq[age_w_pqlseq$outcome %in% chron_age_pqlseq$outcome_chron_age,]
age_m_pqlseq<- age_m_pqlseq[age_m_pqlseq$outcome %in% chron_age_pqlseq$outcome_chron_age,]
all.equal(age_w_pqlseq$outcome, chron_age_pqlseq$outcome_chron_age)
pqlseq_model<- cbind(age_w_pqlseq, age_m_pqlseq[,7:12], chron_age_pqlseq[,7:12])

# Import Cross_Se Models--------------------------------------------------------
setwd('/scratch/ckelsey4/Cayo_meth/cross_models')

long_cross_files<- '_long'
long_cross_pqlseq<- import_pqlseq(long_cross_files, y = 3)

short_cross_files<- '_short'
short_cross_pqlseq<- import_pqlseq(short_cross_files, y = 3)

pqlseq_model<- pqlseq_model[pqlseq_model$outcome %in% long_cross_pqlseq$outcome,]
long_cross_pqlseq<- long_cross_pqlseq[long_cross_pqlseq$outcome %in% pqlseq_model$outcome,]
short_cross_pqlseq<- short_cross_pqlseq[short_cross_pqlseq$outcome %in% long_cross_pqlseq$outcome,]

all.equal(short_cross_pqlseq$outcome, long_cross_pqlseq$outcome)

long_cross_pqlseq<- long_cross_pqlseq %>%
  select(-c(outcome, length, chr, chromStart, chromEnd, n)) %>%
  mutate_if(is.character, as.numeric)
colnames(long_cross_pqlseq)<- paste(colnames(long_cross_pqlseq), "long", "cross", sep = "_")

short_cross_pqlseq<- short_cross_pqlseq %>%
  select(-c(outcome, length, chr, chromStart, chromEnd, n))
colnames(short_cross_pqlseq)<- paste(colnames(short_cross_pqlseq), "short", "cross", sep = "_")

age_full<- cbind(pqlseq_model, long_cross_pqlseq, short_cross_pqlseq)

#Sort chromosome factors
sorted_labels<- str_sort(unique(age_full$chr), numeric=T)

age_full<- age_full %>% 
  mutate(chr = factor(chr, levels = sorted_labels)) %>%
  arrange(chr)

#rm(age_w_pqlseq);rm(age_m_pqlseq);rm(chron_age_pqlseq)

age_full<- age_full %>%
  mutate(beta_cross_age = as.numeric(beta_cross_age),
         cross_diff = abs(beta_cross_age) - abs(beta_chron_age),
         long_diff = abs(beta_within_age) - abs(beta_chron_age),
         aging_diff = abs(beta_cross_age) - abs(beta_within_age))

######################################
###       Descriptive Stats        ###
######################################

age_full %>%
  select(c(beta_within_age, beta_chron_age, beta_long_cross, beta_short_cross)) %>%
  cor(use="pairwise.complete.obs") %>%
  ggcorrplot(show.diag=FALSE, type="lower", lab=TRUE, lab_size=5, sig.level = 0.05, insig = "blank")

######################################
###       Plot Distributions       ###
######################################

compare_plot<- function(df, fdr1, fdr2, var1, var2, c1, c2, lab1, lab2) {
  
  df<- df %>%
    filter(fdr1 < .05 | fdr2 < .05) %>%
    mutate(diff = abs(var1) - abs(var2))
  
  df_lm<- lm(var1 ~ var2, data=df)
    
  df %>%
    ggplot(aes(var1, var2, colour = diff)) +
    geom_point() +
    geom_abline() +
    geom_abline(slope = test_lm[["coefficients"]][[2]], 
                intercept = test_lm[["coefficients"]][[1]],
                colour = "red") +
    geom_vline(xintercept=0, linetype="dashed") +
    geom_hline(yintercept=0, linetype="dashed") +
    scale_color_gradient2(low = c1, mid = "white", high = c2, midpoint = 0, name = "") +
    theme_classic(base_size=32) +
    theme(legend.key.height= unit(2, 'cm')) +
    xlab(lab1) +
    ylab(lab2) 
}

#NEED TO FIX THIS FUNCTION
compare_plot(age_full, 'fdr_short_cross', "fdr_chron_age",
             beta_short_cross, beta_chron_age,
             "darkgoldenrod2", "hotpink3",
             "Short Cross Age", "Chronological Age")

age_full %>%
  filter(fdr_short_cross < .05 | fdr_chron_age < .05) %>%
  mutate(diff = abs(beta_short_cross) - abs(beta_chron_age)) %>%
  ggplot(aes(beta_short_cross, beta_chron_age, colour = diff)) +
  geom_point() +
  geom_abline() +
  geom_smooth(method="lm") +
  geom_vline(xintercept=0, linetype="dashed") +
  geom_hline(yintercept=0, linetype="dashed") +
  scale_color_gradient2(low = "darkgoldenrod2", mid = "white", high = "hotpink3", midpoint = 0, name = "") +
  theme_classic(base_size=32) +
  theme(legend.key.height= unit(2, 'cm')) +
  ylab("Chronological Age") +
  xlab("Short Cross Age")

age_full %>%
  filter(fdr_long_cross < .05 | fdr_chron_age < .05) %>%
  mutate(diff = abs(beta_long_cross) - abs(beta_chron_age)) %>%
  ggplot(aes(beta_long_cross, beta_chron_age, colour = diff)) +
  geom_point() +
  geom_abline() +
  geom_smooth(method="lm") +
  geom_vline(xintercept=0, linetype="dashed") +
  geom_hline(yintercept=0, linetype="dashed") +
  scale_color_gradient2(low = "darkgoldenrod2", mid = "white", high = "steelblue2", midpoint = 0, name = "") +
  theme_classic(base_size=32) +
  theme(legend.key.height= unit(2, 'cm')) +
  ylab("Chronological Age") +
  xlab("Long Cross Age")

test<- age_full %>% filter(fdr_short_cross < .05 | fdr_within_age < .05)
test_lm<- lm(beta_short_cross ~ beta_within_age, data = test)
test_lm[["coefficients"]][1]

age_full %>%
  filter(fdr_short_cross < .05 | fdr_within_age < .05) %>%
  mutate(diff = abs(beta_short_cross) - abs(beta_within_age)) %>%
  ggplot(aes(beta_short_cross, beta_within_age, colour = diff)) +
  geom_point() +
  geom_abline() +
  geom_abline(slope = test_lm[["coefficients"]][[2]], 
              intercept = test_lm[["coefficients"]][[1]],
              colour = "red",
              linetype="dashed") +
  #geom_smooth(method="lm") +
  geom_vline(xintercept=0, linetype="dashed") +
  geom_hline(yintercept=0, linetype="dashed") +
  scale_color_gradient2(low = "purple", mid = "white", high = "hotpink3", midpoint = 0, name = "") +
  theme_classic(base_size=32) +
  theme(legend.key.height= unit(2, 'cm')) +
  ylab("Within Age") +
  xlab("Short Cross Age")

age_full %>%
  filter(fdr_long_cross < .05 | fdr_within_age < .05) %>%
  mutate(diff = abs(beta_long_cross) - abs(beta_within_age)) %>%
  ggplot(aes(beta_long_cross, beta_within_age, colour = diff)) +
  geom_point() +
  geom_abline() +
  geom_smooth(method="lm") +
  geom_vline(xintercept=0, linetype="dashed") +
  geom_hline(yintercept=0, linetype="dashed") +
  scale_color_gradient2(low = "purple", mid = "white", high = "steelblue2", midpoint = 0, name = "") +
  theme_classic(base_size=32) +
  theme(legend.key.height= unit(2, 'cm')) +
  ylab("Within Age") +
  xlab("Long Cross Age")

age_full %>%
  filter(fdr_chron_age < .05 | fdr_within_age < .05) %>%
  mutate(diff = abs(beta_chron_age) - abs(beta_within_age)) %>%
  ggplot(aes(beta_chron_age, beta_within_age, colour = diff)) +
  geom_point() +
  geom_abline() +
  geom_smooth(method="lm") +
  geom_vline(xintercept=0, linetype="dashed") +
  geom_hline(yintercept=0, linetype="dashed") +
  scale_color_gradient2(low = "purple", mid = "white", high = "darkgoldenrod2", midpoint = 0, name = "") +
  theme_classic(base_size=32) +
  theme(legend.key.height= unit(2, 'cm')) +
  ylab("Within Age") +
  xlab("Chronological Age")


