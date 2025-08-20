library(tidyverse)
library(GenomicFeatures)
library(GenomicRanges)
library(rtracklayer)
library(biomaRt)
library(bsseq)
library(fgsea)
library(broom)
library(ggpubr)
library(ggcorrplot)
library(ggborderline)
library(ggrepel)
library(ggridges)
library(ggvenn)
library(msigdbr)
library(effectsize)
options(scipen = 999)
setwd('/scratch/ckelsey4/Cayo_meth')
load("wb_age_analysis.RData")

######################################
####        IMPORT FILES           ###
######################################
#Import metadata----------------------------------------------------------------
long_data<- readRDS("long_data_adjusted")
long_data<- long_data %>%
  dplyr::rename(sex = individual_sex) %>%
  filter(n > 1)

#Generate %methylation matrix
regions_cov<- readRDS("regions_cov_filtered")
regions_cov<- regions_cov[1:21]
regions_m<- readRDS("regions_m_filtered")
regions_m<- regions_m[1:21]

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


######################################
###  Plot metadata distributions   ###
######################################
#No. of individuals by sex------------------------------------------------------
#Total unique individuals
length(unique(long_data$monkey_id))

long_data %>%
  distinct(monkey_id, .keep_all = T) %>%
  ggplot(aes(sex, fill = sex)) +
  geom_bar(colour="black", alpha = 0.8) + 
  scale_fill_manual(values = c("royalblue1", "orangered1"), name = "Sex") +
  geom_text(stat = "count", aes(label = after_stat(count)),
            position = position_dodge(),
            color = "black",
            size = 8,
            vjust = -0.2) +
  ylim(0, 100) +
  theme_classic(base_size = 24)

#Mean age by sex----------------------------------------------------------------
long_data %>% 
  distinct(monkey_id, .keep_all = T) %>%
  ggplot(aes(x=age_at_sampling, fill = sex)) +
  geom_histogram(alpha=0.7, position = position_dodge(width = 1.5), bins=10, colour="black") +
  scale_fill_manual(values = c("darkolivegreen", "darkmagenta"), name = "Sex") +
  scale_x_continuous(breaks = seq(0, 30, by=5)) +
  theme_classic(base_size = 24) +
  labs(y = "Count", x = "Mean Age")

length(unique(long_data$monkey_id))

#N samples by sex---------------------------------------------------------------
long_data %>% 
  distinct(monkey_id, .keep_all = T) %>%
  ggplot(aes(x=n, fill = sex)) +
  geom_bar(alpha = 0.7, position = position_dodge(width = 0.5), colour = "black") +
  scale_fill_manual(values = c("darkolivegreen", "darkmagenta"), name = "Sex") +
  theme_classic(base_size = 32) +
  labs(y = "Count", x = "N")

#N samples per individual-------------------------------------------------------
long_data<- long_data %>%
  group_by(monkey_id) %>%
  mutate(min_age = min(age_at_sampling))
long_data$age_at_sampling<- round(long_data$age_at_sampling, 0)

long_data %>%
  ggplot(aes(x=age_at_sampling, y=reorder(monkey_id, min_age), colour=sex)) +
  geom_path(linewidth = 1.5, alpha = 0.8) +
  geom_point(colour="black") +
  scale_x_continuous(breaks = seq(0, 30, by=2)) +
  scale_colour_manual(values = c("darkolivegreen", "darkmagenta"), name = "Sex") +
  ylab("Individual") +
  xlab("Age") +
  theme_classic(base_size = 24) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

#N pids-------------------------------------------------------------------------
long_data %>%
  ggplot(aes(pid, fill=pid)) +
  geom_bar() +
  geom_text(stat='count', aes(label = after_stat(count)), vjust=-1) +
  ylim(0, 120) +
  theme_classic(base_size=24) +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

######################################
### Generate PCA for %Meth matrix  ###
######################################
#Generate PCA-------------------------------------------------------------------
pca_wb<- prcomp(cor(ratio, use="pairwise.complete.obs"))
pcs<- as.data.frame(pca_wb$x) 
pcs$lid_pid<- rownames(pcs)
pcs<- pcs %>%
  dplyr::relocate(lid_pid, .before = PC1)
pcs<- pcs %>%
  arrange(lid_pid)
long_data<- long_data %>%
  arrange(lid_pid)
all.equal(pcs$lid_pid, long_data$lid_pid)

#Check which pca's explain the most variance
summary(pca_wb)$importance[2, ]

pcs<- cbind(pcs[1:10], dplyr::select(long_data, c("age_at_sampling", "within.age", "mean.age", "sex", "pid")))

pcs %>% 
  ggplot(aes(PC1, within.age)) + 
  geom_point() +
  theme_classic(base_size=24)

pc_lm<- lm(PC1 ~ within.age + mean.age + sex + pid, data = pcs)
predicted_pcs<- predict.lm(pc_lm)
summary(pc_lm)

pcs<- cbind(pcs, predicted_pcs)

pcs %>%
  ggplot(aes(within.age, PC1)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_classic(base_size = 24)

pc.matrix<- model.matrix(~ PC1 + PC2 + PC3 + age_at_sampling + within.age + sex + mean.age, data = pcs)
pc.matrix %>% 
  cor(use="pairwise.complete.obs") %>%
  ggcorrplot(show.diag=FALSE, type="lower", lab=TRUE, lab_size=5, sig.level = 0.05, insig = "blank")

######################################
###      Import GLMER Models       ###
######################################
#Import glmer model data--------------------------------------------------------
setwd('/scratch/ckelsey4/Cayo_meth/glmer_model_compare')
#Generate list of file names
glm_list<- list.files(pattern = 'glmer_output')
glm_order<- str_split_i(glm_list, "_", 3)

#Import glm models as list
glm_model_list<- lapply(glm_list, readRDS)

#Rename list elements
names(glm_model_list)<- glm_order

#Bind model list to df and add rownames
glm_model<- do.call(rbind, glm_model_list)
glm_model$region<- str_split_i(glm_model$region, "\\.", 3)
rownames(glm_model)<- glm_model$region

#Separate region coordinates into start and end, delete the chr col, and move region col to front
glm_model<- glm_model %>%
  relocate(region, .before = `estimate_(Intercept)`)
glm_model$region2<- glm_model$region
glm_model<- glm_model %>%
  separate_wider_delim(region2, "_", names = c("chrom", "chromStart", "chromEnd"))
glm_model<- glm_model %>%
  relocate(c(chrom, chromStart, chromEnd), .before = `estimate_(Intercept)`)
glm_model<- glm_model %>%
  mutate(chromStart = as.numeric(chromStart))

#Filter for true convergences
glm_model.t<- glm_model %>%
  filter(converged == "TRUE")

#Generate df of adjusted pvalues
glm_fdr<- lapply(glm_model.t[, c("pval_(Intercept)", "pval_age.w", "pval_age.m", "pval_sexM")], function(x){
  fdr<- p.adjust(x, method = "fdr")
})

glm_fdr<- do.call(cbind, glm_fdr)

#Rename padj cols
colnames(glm_fdr)<- c("intercept.fdr", "age.w.fdr", "age.m.fdr", "sex.fdr")

#Bind padj cols to glm_model df and relocate
glm_model.t<- cbind(glm_model.t, glm_fdr)
glm_model.t<- glm_model.t %>%
  relocate(c("intercept.fdr", "age.w.fdr", "age.m.fdr", "sex.fdr"), .after = pval_sexM)

age<- glm_model.t %>%
  filter(age.w.fdr < 0.10)



######################################
###      Import PQLseq Models      ###
######################################
setwd('/scratch/ckelsey4/Cayo_meth/glmer_model_compare')
#Import additive pqlseq files---------------------------------------------------
import_pqlseq<- function(x, v, y){
  
  if (v == "age") {
    
    #Generate list of file names
    file_list<- list.files(pattern = x)
    
    if (y == 5){
      file_order<- str_split_i(file_list, "_", 5)
    } else {
      file_order<- str_split_i(file_list, "_", 4)
    }
    
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
      mutate(region_range = paste(as.character(chromStart), "-", as.character(chromEnd))) %>% 
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
      dplyr::select(-c(elapsed_time, converged))
    
  } else if (v == "sex") {
    
    #Generate list of file names
    file_list<- list.files(pattern = x)
    
    if (y == 5){
      file_order<- str_split_i(file_list, "_", 5)
    } else {
      file_order<- str_split_i(file_list, "_", 4)
    }
    
    #Import glm models as list
    model_list<- lapply(file_list, readRDS)
    
    #Rename list elements
    names(model_list)<- file_order
    
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
      mutate(region_range = paste(as.character(chromStart), "-", as.character(chromEnd))) %>%
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
      dplyr::select(-c(elapsed_time, converged))
  }
}

#Chronological Age
chron_age_files<- 'wb_pqlseq2_agechron'
chron_age_pqlseq<- import_pqlseq(chron_age_files, v = "age", y = 4)

#Age Within
age_w_files<- 'wb_pqlseq2_within_age'
age_w_pqlseq<- import_pqlseq(age_w_files, v = "age", y = 5)

#Mean Age
age_m_files<- 'wb_pqlseq2_mean_age'
age_m_pqlseq<- import_pqlseq(age_m_files, v = "age", y = 5)

#Sex
sex_files<- 'wb_pqlseq2_sex'
sex_pqlseq<- import_pqlseq(sex_files, v = "sex", y = 4)

#Rename cols for each df to indicate variable
colnames(age_w_pqlseq)<- c("outcome", "length", "chr", "chromStart", "chromEnd", "n", 
                           paste(names(age_w_pqlseq[,7:14]), "age", sep = "_"), "region_range")
colnames(age_m_pqlseq)<- c(paste(names(age_m_pqlseq), "mean_age", sep = "_"))
colnames(chron_age_pqlseq)<- c(paste(names(chron_age_pqlseq), "chron_age", sep = "_"))

#Cbind cols for age and sex dfs
pqlseq_model<- cbind(age_w_pqlseq, age_m_pqlseq[,7:14], chron_age_pqlseq[,7:14])

#Sort chromosome factors
sorted_labels<- str_sort(unique(pqlseq_model$chr), numeric=T)

pqlseq_model<- pqlseq_model %>%
  mutate(chr = factor(chr, levels = sorted_labels)) %>%
  arrange(chr)

pqlseq_model<- pqlseq_model %>%
  mutate(age_ratio = abs(beta_age)/abs(beta_chron_age),
         age_diff = abs(beta_age) - abs(beta_chron_age))

sex_pqlseq<- sex_pqlseq %>%
  mutate(chr = factor(chr, levels = sorted_labels[1:21])) %>%
  arrange(chr)
  
rm(age_w_pqlseq);rm(age_m_pqlseq);rm(chron_age_pqlseq)

#Compare age-within and chron-age-----------------------------------------------
#Scatterplot of chron-age vs age-within
pqlseq_model %>%
  ggplot(aes(beta_chron_age, beta_age, colour = age_diff)) +
  geom_point() +
  geom_abline() +
  geom_smooth(method = "lm") +
  geom_vline(xintercept=0, linetype="dashed") +
  geom_hline(yintercept=0, linetype="dashed") +
  scale_color_gradient2(low = , mid = "white", high = , midpoint = 0, name = "") +
  theme_classic(base_size=32) +
  theme(legend.key.height= unit(2, 'cm')) +
  ylab("Age Within") +
  xlab("Chronological Age")

pqlseq_model %>%
  ggplot(aes(log(age_ratio), fill=after_stat(x))) +
  geom_histogram(bins = 50, colour = "black", position = 'identity') +
  geom_vline(xintercept = 0, linetype = "dashed", colour = 'red') +
  scale_fill_gradient2(low = "steelblue2", mid = "white", high = "hotpink3", midpoint = 0, name = "Beta Difference") +
  xlab("Log(|Beta (Age Within)|/|Beta (Chron Age)|)") +
  ylab("Hypomethylated Regions") +
  theme_classic(base_size = 32) +
  theme(legend.position = "none")

#Comparison of distributions
age_betas<- pqlseq_model %>%
  dplyr::select(c(outcome, beta_age,, beta_mean_age, beta_chron_age, )) %>%
  pivot_longer(!outcome, names_to = "var", values_to = "beta")

age_fdr<- pqlseq_model %>%
  dplyr::select(c(outcome, fdr_age, fdr_mean_age, fdr_chron_age)) %>%
  pivot_longer(!outcome, names_to = "var", values_to = "fdr")

age_effects<- cbind(age_betas, age_fdr[,3])

age_effects$var<- gsub("beta_", "", age_effects$var)

rm(age_betas);rm(age_fdr)

age_effects %>%
  filter(var != "mean_age") %>%
  ggplot(aes(beta, fill = var)) +
  geom_density(colour="black", position = "identity", alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "red") +
  scale_fill_manual(values = c("hotpink3", "steelblue2")) +
  theme_classic(base_size = 32) +
  xlab("Beta") +
  ylab("Density")

age_effects %>%
  filter(var != "mean_age") %>%
  ggplot(aes(var, beta, fill = var)) +
  geom_boxplot(width=0.3) +
  geom_hline(yintercept = 0, colour = "black", linetype = "dashed") +
  scale_fill_manual(values = c("hotpink3", "steelblue2")) +
  theme_classic(base_size = 32) +
  xlab("Density") +
  ylab("Beta")

#Plot pqlseq features-----------------------------------------------------------
age.w.count<- nrow(pqlseq_model[pqlseq_model$fdr_age < 0.05,])
age.m.count<- nrow(pqlseq_model[pqlseq_model$fdr_mean_age < 0.05,])
sex.count<- nrow(sex_pqlseq[sex_pqlseq$fdr < 0.05,])
counts<- data.frame(count = c(age.w.count, age.m.count, sex.count),
                    predictor = as.factor(c('Age Within', 'Age Between', 'Sex')))
counts<- counts %>%
  mutate(predictor = as.factor(predictor)) %>%
  mutate(predictor=fct_reorder(predictor, count, .desc=T))

rm(age.w.count);rm(age.m.count);rm(sex.count)

counts %>%
  #filter(predictor != "Sex") %>%
  ggplot(aes(predictor, count, fill = predictor)) +
  geom_bar(stat = 'identity', colour="black") +
  geom_text(label=counts$count, vjust=-1, size =5) +
  theme_classic(base_size = 32) +
  xlab("Predictor") +
  ylab("Count") +
  scale_fill_manual(values = c("hotpink3", "darkgoldenrod2", "seagreen4"))

## Significant regions Venn diagram
age.w<- pqlseq_model$outcome[pqlseq_model$fdr_age < 0.05 & pqlseq_model$chrom != "Y"]
age.m<- pqlseq_model$outcome[pqlseq_model$fdr_mean_age < 0.05 & pqlseq_model$chrom != "Y"]
sex<- sex_pqlseq$outcome[sex_pqlseq$fdr < 0.05]

venn_list<- list(age.w, age.m, sex)

rm(age.w);rm(age.m);rm(sex)

ggvenn(venn_list,
       fill_color = c("hotpink3", "seagreen4", "darkgoldenrod2"),
       text_size = 8,
       show_percentage = F)

#Distribution of effect sizes and across the chromosomes
#Age
pqlseq_model %>%
  filter(fdr_age < 0.05) %>%
  ggplot(aes(beta_age, fill = beta_age>0)) +
  geom_histogram(bins=100, colour="black") +
  geom_vline(xintercept=0, linetype="dashed") +
  scale_fill_manual(values = c("hotpink3", "hotpink")) +
  xlab("Age Estimate") +
  ylab("Count") +
  theme_classic(base_size=32)

pqlseq_model %>%
  filter(fdr_age < .05 | fdr_mean_age < .05) %>%
  ggplot(aes(beta_age, beta_sex, colour = abs_diff)) +
  geom_point() +
  geom_abline() +
  geom_vline(xintercept=0, linetype="dashed") +
  geom_hline(yintercept=0, linetype="dashed") +
  scale_color_gradient2(low = "darkgoldenrod2", mid = "white", high = "hotpink3", midpoint = 0, name = "") +
  theme_classic(base_size=32) +
  theme(legend.key.height= unit(2, 'cm')) +
  ylab("Sex") +
  xlab("Age Within")

pqlseq_full %>%
  filter(fdr_age < .05 & fdr_sex < .05) %>%
  ggplot(aes(beta_age, beta_sex)) +
  geom_point() +
  geom_abline() +
  geom_smooth(method = "lm") +
  geom_vline(xintercept=0, linetype="dashed") +
  geom_hline(yintercept=0, linetype="dashed") +
  #scale_color_gradient2(low = "darkgoldenrod2", mid = "white", high = "hotpink3", midpoint = 0, name = "") +
  theme_classic(base_size=32) +
  theme(legend.key.height= unit(2, 'cm')) +
  ylab("Age Within") +
  xlab("Sex")

log_age<- pqlseq_model %>%
  filter(fdr_age < .05 | fdr_mean_age < .05) %>%
  mutate(log_age = log(abs(beta_age)/abs(beta_mean_age))) %>%
  dplyr::select(log_age)

mean(log_age$log_age)

log_age %>%
  ggplot(aes(log_age, fill = log_age<0)) +
  geom_histogram(colour="black", position = "identity", bins = 100) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "red") +
  geom_vline(xintercept = median(log_age$log_age), linetype = "dashed", colour = "green") +
  scale_fill_manual(values = c("hotpink3", "seagreen4")) +
  theme_classic(base_size = 32) +
  theme(legend.position = "none") +
  ylab("Count") +
  xlab("Log(Age Within/Age Between)")

df<- pqlseq_model %>%
  dplyr::select(c(chr, beta_age, fdr_age))
df$type<- "autosomes"
df$type[df$chr == "X"]<- "X"
df$type[df$chr == "Y"]<- "Y"

df %>%
  filter(fdr_age < 0.05) %>%
  ggplot(aes(type, beta_age, colour=beta_age>0)) +
  geom_jitter(alpha=0.5) +
  geom_hline(yintercept = 0, linetype="dashed") +
  scale_colour_manual(values = c("slateblue3", "orangered2")) +
  theme_classic(base_size = 32) +
  theme(legend.position = "none") +
  ylab("Age Estimate (FDR < .05)") +
  xlab("Chromosome Type")

#Sex
sex_pqlseq %>%
  #filter(fdr < 0.05) %>%
  ggplot(aes(beta, fill=beta>0)) +
  geom_histogram(bins=100, colour="black") +
  geom_vline(xintercept=0, linetype="dashed") +
  scale_fill_manual(values = c("darkmagenta", "darkolivegreen")) +
  xlab("Sex Estimate") +
  ylab("Count") +
  theme_classic(base_size=32) +
  theme(legend.position = "none")

#Without X-chrom
sex_pqlseq %>%
  filter(chr != "X") %>%
  ggplot(aes(beta, fill=beta>0)) +
  geom_histogram(bins=100, colour="black") +
  geom_vline(xintercept=0, linetype="dashed") +
  scale_fill_manual(values = c("darkmagenta", "darkolivegreen")) +
  xlab("Sex Estimate") +
  ylab("Count") +
  theme_classic(base_size=32) +
  theme(legend.position = "none")

nrow(sex_pqlseq %>%
  filter(fdr < 0.05 & chr != "X"))

nrow(sex_pqlseq %>%
       filter(fdr < 0.05))

df<- sex_pqlseq %>%
  dplyr::select(c(chr, beta, fdr))
df$type<- "autosomes"
df$type[df$chr == "X"]<- "X"

df %>%
  filter(fdr < 0.05) %>%
  ggplot(aes(type, beta, colour=beta>0)) +
  geom_jitter(alpha=0.5) +
  geom_hline(yintercept = 0, linetype="dashed") +
  scale_colour_manual(values = c("darkolivegreen", "darkmagenta")) +
  theme_classic(base_size = 32) +
  theme(legend.position = "none") +
  ylab("Sex Estimate (FDR < .05)") +
  xlab("Chromosome Type")

rm(df)
  
#Distribution of coefficients for age and sex combined
age_sex<- pqlseq_model %>% 
  filter(chr != "Y") %>%
  dplyr::select(c(chr, fdr_age, beta_age))
age_sex<- cbind(age_sex, sex_pqlseq[,c(9,12)])

age_sex$direction[age_sex$beta_age > 0 & age_sex$beta > 0]<- "Both Positive"
age_sex$direction[age_sex$beta_age > 0 & age_sex$beta < 0]<- "Age Positive, Sex Negative"
age_sex$direction[age_sex$beta_age < 0 & age_sex$beta > 0]<- "Age Negative, Sex Positive"
age_sex$direction[age_sex$beta_age < 0 & age_sex$beta < 0]<- "Both Negative"

age_sex<- age_sex %>%
  mutate(direction = as.factor(direction)) %>%
  mutate(direction=reorder(direction, direction, FUN=length))

age_sex %>%
  filter(fdr_age < 0.05 & fdr < 0.05 & chr != "X") %>%
  ggplot(aes(beta_age, beta)) +
  geom_point(aes(colour=direction)) +
  geom_smooth(method = "lm") +
  theme_classic(base_size=32) +
  theme(legend.position = "none") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_colour_brewer(palette = "Set3") +
  ylab("Sex Estimate") +
  xlab("Age Estimate")

summary(lm(beta ~ beta_age, data = age_sex))
cor.test(age_sex$beta, age_sex$beta_age)

age_sex %>%
  filter(fdr_age < 0.05 & fdr < 0.05) %>%
  ggplot(aes(direction, fill=direction)) +
  geom_bar(colour="black") +
  geom_text(aes(label=after_stat(count)), stat="count", vjust=-1) +
  theme_classic(base_size = 32) +
  theme(axis.text.x=element_blank()) +
  theme(legend.position = "none") +
  scale_fill_brewer(palette = "Set3") +
  ylab("Count") +
  xlab("Direction")

#Import glm random effects summaries--------------------------------------------
#Generate list of file names
glm_rfx_list<- list.files('/scratch/ckelsey4/Cayo_meth/glmer_model_compare', 
                          pattern = 'glmer_ranefx', 
                          full.names = T)

#Import glm models as list
glm_rfx<- lapply(glm_rfx_list, readRDS)
rfx_order<- gsub("/scratch/ckelsey4/Cayo_meth/glmer_model_compare/glmer_ranefx_", "", glm_rfx_list)
names(glm_rfx)<- rfx_order

#This code block filters out unnecessary rows and cols from each region model within a single chromosome, have to scale up to all chrs
df1<- glm_rfx[[1]]
df1<- sapply(names(df1), function(x){
  df<- df1[[x]]
  df<- df %>%
    filter(term == "age.w") %>%
    dplyr::select(-c(grpvar, term, condsd))
}, USE.NAMES = T, simplify = F)

df2<- do.call(rbind, df1)
df2$region<- str_split_i(rownames(df2), "\\.", 3)

df3<- df2 %>%
  pivot_wider(names_from = region, values_from = condval)

sex<- long_data %>%
  distinct(monkey_id, .keep_all = T) %>%
  dplyr::select(monkey_id, sex) %>%
  arrange(monkey_id)

df3<- cbind(df3, sex[,2])
df3<- df3 %>%
  relocate(sex, .after = grp)

#Remove intercept random effects

for (i in names(glm_rfx)){
  
  rfx_chr<- glm_rfx[[i]]
  
  rfx_chr<- lapply(rfx_chr, function(x){
    x %>%
      filter(term == "age.w") %>%
      dplyr::select(-c(grpvar, term, condsd))
  })
  
}
rm(df1);rm(df2)


######################################
###      JOIN INTERSECT FILES      ###   
######################################
##Import annotation files-------------------------------------------------------
re_anno<- read_csv("/scratch/ckelsey4/Cayo_meth/re_annotations.csv")
re_anno<- re_anno %>%
  filter(chr != "Y")
chmm_intersect<- read_csv("/scratch/ckelsey4/Cayo_meth/chmm_annotations.csv")
chmm_intersect<- chmm_intersect %>%
  filter(chr != "Y")
promoters<- read_csv("/scratch/ckelsey4/Cayo_meth/promoters.csv")
promoters<- promoters %>%
  filter(chr != "Y")

#Promoters----------------------------------------------------------------------
pqlseq_prom<- left_join(promoters, pqlseq_model, by = c("region_range", "chr"))
pqlseq_prom<- pqlseq_prom %>%
  drop_na() %>%
  distinct(anno, .keep_all = T)
pqlseq_prom$anno_class<- "Promoter"
pqlseq_prom<- pqlseq_prom %>%
  dplyr::select(-c(chromStart, chromEnd, outcome))
  
#CHMM---------------------------------------------------------------------------
#Join pqlseq model and chmm 
pqlseq_chmm<- left_join(chmm_intersect, pqlseq_model, by = c("region_range", "chr"))

#Filter out regions with models that didn't converge resulting in NAs in the annotation join
pqlseq_chmm<- pqlseq_chmm %>%
  drop_na()

#Select out unnecessary cols and rename
pqlseq_chmm<- pqlseq_chmm %>%
  dplyr::select(-c(outcome, chromStart, chromEnd))

#Set annotations as factor and reorder
annotations_ordered<- str_sort(unique(pqlseq_chmm$anno), numeric = TRUE)
pqlseq_chmm$anno<- factor(pqlseq_chmm$anno, levels = rev(annotations_ordered))

#Create column of broad categories for annotations
#pqlseq_chmm$anno_class<- "A"
pqlseq_chmm$anno_class<- "Transcription Start Sites"
pqlseq_chmm$anno_class[pqlseq_chmm$anno %in% annotations_ordered[3:5]]<- "Active Transcription"
pqlseq_chmm$anno_class[pqlseq_chmm$anno %in% annotations_ordered[6:8]]<- "Enhancer Regions"
pqlseq_chmm$anno_class[pqlseq_chmm$anno %in% annotations_ordered[9:15]]<- "Quiescent States"

#Set classes as factors
class_factors<- c("Transcription Start Sites", "Active Transcription", "Enhancer Regions", "Quiescent States")
pqlseq_chmm$anno_class<- factor(pqlseq_chmm$anno_class, levels = rev(class_factors))


sex_pqlseq_chmm<- left_join(chmm_intersect, sex_pqlseq, by = c("region_range", "chr"))

sex_pqlseq_chmm<- sex_pqlseq_chmm %>%
  drop_na()

sex_pqlseq_chmm<- sex_pqlseq_chmm %>%
  dplyr::select(-c(outcome, chromStart, chromEnd))

sex_pqlseq_chmm$anno<- factor(sex_pqlseq_chmm$anno, levels = rev(annotations_ordered))

#Create column of broad categories for annotations
sex_pqlseq_chmm$anno_class<- "A"
sex_pqlseq_chmm$anno_class[sex_pqlseq_chmm$anno %in% annotations_ordered[1:2]]<- "Transcription Start Sites"
sex_pqlseq_chmm$anno_class[sex_pqlseq_chmm$anno %in% annotations_ordered[3:5]]<- "Active Transcription"
sex_pqlseq_chmm$anno_class[sex_pqlseq_chmm$anno %in% annotations_ordered[6:8]]<- "Enhancer Regions"
sex_pqlseq_chmm$anno_class[sex_pqlseq_chmm$anno %in% annotations_ordered[9:15]]<- "Quiescent States"

sex_pqlseq_chmm$anno_class<- factor(sex_pqlseq_chmm$anno_class, levels = rev(class_factors))

#REPEAT ELEMENTS----------------------------------------------------------------
#Join repeats annotations and glm_models df
pqlseq_re<- left_join(re_anno, pqlseq_model, by = c("region_range", "chr"))

pqlseq_re<- pqlseq_re %>%
  drop_na()

#Remove non-sensical annotations (NA, Unknown etc)
pqlseq_re<- pqlseq_re[!pqlseq_re$repClass == "Unknown",]
pqlseq_re<- pqlseq_re[!pqlseq_re$repClass == "DNA?",]
pqlseq_re<- pqlseq_re[!pqlseq_re$repClass == "LTR?",]
pqlseq_re<- pqlseq_re[!pqlseq_re$repClass == "RC?",]
pqlseq_re<- pqlseq_re[!pqlseq_re$repClass == "Unspecified",]

pqlseq_re$anno_class<- "Simple Repeats"
pqlseq_re$anno_class[pqlseq_re$repClass %in% c("SINE", "LINE", "LTR", "Retroposon")]<- "Transposable Elements Class I"
pqlseq_re$anno_class[pqlseq_re$repClass %in% "DNA"]<- "Transposable Elements Class II"
pqlseq_re$anno_class[pqlseq_re$repClass %in% c("rRNA", "snRNA", "tRNA", "srpRNA", "scRNA")]<- "Structural RNAs"

pqlseq_re<- pqlseq_re %>%
  dplyr::rename(anno = repClass) %>%
  dplyr::select(-c(repName, chromStart, chromEnd)) %>%
  dplyr::relocate(anno, .after = anno_end)

#Select out unnecessary cols and rename
pqlseq_re<- pqlseq_re %>%
  dplyr::select(-c(range, outcome))

sex_pqlseq_re<- left_join(re_anno, sex_pqlseq, by = c("region_range", "chr"))

sex_pqlseq_re<- sex_pqlseq_re %>%
  drop_na()

sex_pqlseq_re<- sex_pqlseq_re %>%
  dplyr::select(-c(range, outcome))

sex_pqlseq_re<- sex_pqlseq_re[!sex_pqlseq_re$repClass == "Unknown",]
sex_pqlseq_re<- sex_pqlseq_re[!sex_pqlseq_re$repClass == "DNA?",]
sex_pqlseq_re<- sex_pqlseq_re[!sex_pqlseq_re$repClass == "LTR?",]
sex_pqlseq_re<- sex_pqlseq_re[!sex_pqlseq_re$repClass == "RC?",]
sex_pqlseq_re<- sex_pqlseq_re[!sex_pqlseq_re$repClass == "Unspecified",]

#Bind annotation dfs together---------------------------------------------------
pqlseq_anno<- rbind(pqlseq_chmm, pqlseq_re, pqlseq_prom)

annos<- c("Promoter", "Transcription Start Sites", "Active Transcription", "Enhancer Regions", "Quiescent States", "Simple Repeats",
          "Transposable Elements Class I", "Transposable Elements Class II", "Structural RNAs")

pqlseq_anno$anno_source<- "Repeat Elements"
pqlseq_anno$anno_source[pqlseq_anno$anno_class == "Transcription Start Sites" | pqlseq_anno$anno_class == "Active Transcription" |
                        pqlseq_anno$anno_class == "Enhancer Regions" | pqlseq_anno$anno_class == "Quiescent States" | 
                        pqlseq_anno$anno_class == "Promoter"]<- "Transcription"

pqlseq_anno<- pqlseq_anno %>%
  arrange(anno_source, anno_class)

#Rearrange factors to sort by type then log_or
pqlseq_anno$anno_class<- factor(pqlseq_anno$anno_class, levels = rev(annos))

pqlseq_anno<- pqlseq_anno %>%
  dplyr::relocate(c(anno_class, anno_source), .after=anno)

pqlseq_anno$agew_signif<- "Non-Significant"
pqlseq_anno$agew_signif[pqlseq_anno$fdr_age < 0.05 & pqlseq_anno$beta_age < 0]<- "Age-Hypomethylated"
pqlseq_anno$agew_signif[pqlseq_anno$fdr_age < 0.05 & pqlseq_anno$beta_age > 0]<- "Age-Hypermethylated"

pqlseq_anno$agew_signif<- factor(pqlseq_anno$agew_signif, 
                                 levels = c("Age-Hypermethylated", "Non-Significant", "Age-Hypomethylated"))

pqlseq_anno$agechron_signif<- "Non-Significant"
pqlseq_anno$agechron_signif[pqlseq_anno$fdr_chron_age < 0.05 & pqlseq_anno$beta_chron_age < 0]<- "Age-Hypomethylated"
pqlseq_anno$agechron_signif[pqlseq_anno$fdr_chron_age < 0.05 & pqlseq_anno$beta_chron_age > 0]<- "Age-Hypermethylated"

pqlseq_anno$agechron_signif<- factor(pqlseq_anno$agechron_signif, 
                                     levels = c("Age-Hypermethylated", "Non-Significant", "Age-Hypomethylated"))


pqlseq_anno$unique_cpg<- paste(pqlseq_anno$chr, pqlseq_anno$cpg_loc, sep="_")

#Plot annotation proportions----------------------------------------------------
agew<- pqlseq_anno %>% 
  distinct(unique_cpg, .keep_all = T) %>%
  group_by(anno_class, agew_signif) %>% 
  dplyr::rename(signif = agew_signif) %>%
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count))

agechron<- pqlseq_anno %>% 
  distinct(unique_cpg, .keep_all = T) %>%
  group_by(anno_class, agechron_signif) %>% 
  dplyr::rename(signif = agechron_signif) %>%
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count))

d3<- pqlseq_anno %>%
  mutate(unique_cpg = paste(chr, cpg_loc, sep="_")) %>%
  distinct(unique_cpg, .keep_all = T) %>%
  group_by(signif) %>%
  summarise(count = n()) %>%
  mutate(perc = count/sum(count))

d3$anno_class<- "All"

agew<- rbind(agew, d3)
annos2<- unique(agew$anno_class)
agew$anno_class<- factor(agew$anno_class, levels = annos2)

agechron<- rbind(agechron, d3)
agechron$anno_class<- factor(agechron$anno_class, levels = annos2)

age_total %>%
  ggplot(aes(x = perc*100, y=anno_class, fill = factor(signif))) +
  geom_bar(stat="identity", width = 0.7, colour="black") +
  #geom_text(label=df$count, hjust=-5) +
  geom_vline(xintercept = 50, linetype = 'dashed') +
  theme_classic(base_size=32) +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("hotpink", "gray90", "hotpink3")) +
  ylab("Annotation") +
  xlab("Percentage") +
  facet_wrap(vars(type))

agechron %>%
  ggplot(aes(x = perc*100, y=anno_class, fill = factor(signif))) +
  geom_bar(stat="identity", width = 0.7, colour="black") +
  #geom_text(label=df$count, hjust=-5) +
  geom_vline(xintercept = 50, linetype = 'dashed') +
  theme_classic(base_size=32) +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("steelblue2", "gray90", "steelblue4")) +
  ylab("Annotation") +
  xlab("Percentage")

######################################
###           ENRICHMENT           ###   
######################################
#Class enrichment---------------------------------------------
enrichment<- function(model_df, opt){
  
  df_list<- list()
  
  for(i in unique(model_df$anno_class)) {
    
    df2<- model_df %>%
      filter(anno_class == i) %>%
      distinct(unique_cpg, .keep_all = T)
    
    df3<- model_df %>%
      distinct(unique_cpg, .keep_all = T) %>%
      filter(!unique_cpg %in% df2$unique_cpg)
    
    if (opt == "hypo"){
      
      #Counts for fdr < 0.05 & beta < 0
      a<- nrow(df2[df2$fdr_age < 0.05 & df2$beta_age < 0,])
      b<- nrow(df3[df3$fdr_age < 0.05 & df3$beta_age < 0,])
      
      #Counts for NOT fdr < 0.05 & beta < 0
      c<- nrow(df2[!(df2$fdr_age < 0.05 & df2$beta_age < 0),])
      d<- nrow(df3[!(df3$fdr_age < 0.05 & df3$beta_age < 0),])
      
      #Generate contingency table
      c_table<- data.frame("sig_hypo_Y" = c(a, b),
                           "sig_hypo_N" = c(c, d),
                           row.names = c(paste(i, "Y", sep=""), paste(i, "N", sep="")))
      
      if (all.equal(sum(c_table), length(unique(model_df$unique_cpg)))){
        print(paste("Contingency table sum for", i, "matches unique cpg_loc length"))
      }
      
      df_list[[length(df_list)+1]] = c_table
      
    } else if (opt == "hyper") {
      
      #Counts for fdr < 0.05 & beta < 0
      a<- nrow(df2[df2$fdr_age < 0.05 & df2$beta_age > 0,])
      b<- nrow(df3[df3$fdr_age < 0.05 & df3$beta_age > 0,])
      
      #Counts for NOT fdr < 0.05 & beta < 0
      c<- nrow(df2[!(df2$fdr_age < 0.05 & df2$beta_age > 0),])
      d<- nrow(df3[!(df3$fdr_age < 0.05 & df3$beta_age > 0),])
      
      #Generate contingency table
      c_table<- data.frame("sig_hyper_Y" = c(a, b),
                           "sig_hyper_N" = c(c, d),
                           row.names = c(paste(i, "Y", sep=""), paste(i, "N", sep="")))
      
      if (all.equal(sum(c_table), length(unique(model_df$unique_cpg)))){
        print(paste("Contingency table sum for", i, "matches unique cpg_loc length"))
      }
      
      df_list[[length(df_list)+1]] = c_table
    }
  }
  #name table list
  names(df_list)<- unique(model_df$anno_class)
  
  #Fisher test for each table and tidy with broom
  ft<- lapply(df_list, fisher.test)
  ft<- lapply(ft, broom::tidy)
  
  ft<- do.call(rbind, ft)
  ft<- ft %>%
    mutate(annotation = rownames(ft))
  
  #FDR p-val adjustment
  ft<- ft %>%
    mutate(padj = p.adjust(p.value)) %>%
    mutate_at(vars(annotation), as.factor)
  
  #Log estimates and CIs
  ft<- ft %>%
    mutate(log_or = log(estimate),
           log_ci.lo = log(conf.low),
           log_ci.hi = log(conf.high))
  
  ft$anno_source<- "Repeat Elements"
  ft$anno_source[ft$annotation == "Transcription Start Sites" | ft$annotation == "Active Transcription" |
                          ft$annotation == "Enhancer Regions" | ft$annotation == "Quiescent States" | 
                          ft$annotation == "Promoter"]<- "Transcription"
  ft<- ft %>%
    arrange(anno_source, log_or)
  hyper_levels<- as.character(ft$annotation)
  
  #Rearrange factors to sort by type then log_or
  ft$annotation<- factor(ft$annotation, levels = rev(annos))
  
  if (opt == "hyper") {
    ft$type<- "hyper"
  } else if (opt == "hypo"){
    ft$type<- "hypo"
  }
  
  return(ft)
}

age_enrichment<- rbind(age_hypo, age_hyper)

age_enrichment %>%
  ggplot(aes(x=annotation, y=estimate, fill=type, alpha=padj<0.05)) +
  geom_col(position = position_dodge(0.5)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_errorbar(ymin = age_enrichment$conf.low, ymax = age_enrichment$conf.high, width = 0.3, position = position_dodge(0.7)) +
  scale_fill_manual(values = c("hotpink", "hotpink3")) +
  theme_classic(base_size = 32) +
  theme(legend.position = "none") +
  #ylim(c(-1, 7)) +
  ylab("Odds Ratio") +
  xlab("Annotation") +
  coord_flip()

#Save workspace image-----------------------------------------------------------
save.image("wb_age_analysis.RData")

