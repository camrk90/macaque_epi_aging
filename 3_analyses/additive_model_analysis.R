library(tidyverse)
library(GenomicFeatures)
library(GenomicRanges)
library(rtracklayer)
library(biomaRt)
library(bsseq)
library(fgsea)
library(broom)
library(lme4)
library(ggpubr)
library(ggcorrplot)
library(ggborderline)
library(ggrepel)
library(PQLseq)
library(ggVennDiagram)
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
####        AVERAGE %METH          ###
######################################
#Generate %methylation matrix (m/cov)-------------------------------------------
ratio<- mapply('/', regions_m, regions_cov, SIMPLIFY = F)
ratio<- do.call(rbind, ratio)
ratio<- ratio[colnames(ratio) %in% long_data$lid_pid]
ratio[is.na(ratio)]<- 0
rownames(ratio)<- str_split_i(rownames(ratio), "\\.", 3)
ratio$chrom<- rownames(ratio)
ratio$chrom<- str_split_i(ratio$chrom, "\\_", 1)

#Generate average %methylation across all regions for each individual
av_meth<- colMeans(ratio)
av_meth<- as.vector(av_meth)
av_meth<- cbind(av_meth, colnames(ratio))
av_meth<- av_meth %>%
  as.data.frame() %>%
  arrange(V2)
av_meth<- cbind(av_meth, long_data$age_at_sampling, long_data$sex, long_data$monkey_id)
colnames(av_meth)<- c('perc_meth', 'lid_pid', 'age', 'sex', "monkey_id")
av_meth$perc_meth<- as.numeric(av_meth$perc_meth)

av_meth<- av_meth %>% 
  group_by(monkey_id) %>%
  mutate(diff = max(perc_meth) - min(perc_meth))

av_meth %>%
  ggplot(aes(age, perc_meth)) +
  geom_point(aes(colour=sex)) +
  geom_path(aes(group = monkey_id), alpha = 0.5) +
  geom_smooth(aes(colour=sex),method = "glm") +
  geom_smooth(method = "glm", colour='black', linetype = "dashed", se = F) +
  scale_colour_manual(values = c("royalblue2", "orangered1"), name = "Sex") +
  scale_x_continuous(breaks = seq(0, 30, by=5)) +
  ylab("Average % Methylation") +
  xlab("Age") +
  theme_classic(base_size = 24)

lm_sum<-av_meth %>%
  distinct(monkey_id, .keep_all = T) %>%
  t.test(diff ~ sex, data = .)

av_meth_diff %>%
  distinct(monkey_id, .keep_all = T) %>%
  ggplot(aes(sex, diff)) +
  geom_dotplot(aes(fill = sex), binaxis = "y", stackdir = "center", dotsize = 0.7) +
  geom_boxplot(width=0.05) +
  scale_fill_manual(values = c("royalblue2", "orangered1"), name = "Sex") +
  ylab("%Methylation Change") +
  xlab("Sex") +
  theme_classic(base_size = 24)

#Import regions info
regions<- readRDS("regions_full")
rownames(regions)<- str_split_i(rownames(regions), "\\.", 2)
regions<- regions_df[rownames(regions) %in% rownames(ratio),]


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
  labs(y = "Count", x = "Age")

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
  ggplot(aes(PC1, PC2, colour=pid)) + 
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

pc.matrix<- model.matrix(~ predicted_pcs + PC1 + PC2 + PC3 + within.age + sex + mean.age + pid, data = pcs)
pc.matrix[, 3:ncol(pc.matrix)] %>% 
  cor(use="pairwise.complete.obs") %>%
  ggcorrplot(show.diag=FALSE, type="lower", lab=TRUE, lab_size=2)

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
import_pqlseq<- function(x, v){
  
  if (v == "age") {
    
    #Generate list of file names
    file_list<- list.files(pattern = x)
    
    if (grepl("age", x,  fixed = T) == T){
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
    model$outcome<- ifelse(startsWith(model$outcome, "Y"), str_split_i(model$outcome, "\\.", 2), str_split_i(model$outcome, "\\.", 3))
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
    
    if (grepl("age", x,  fixed = T) == T){
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

#Age Within
age_w_files<- 'wb_pqlseq2_within_age'
age_w_pqlseq<- import_pqlseq(age_w_files, v = "age")

#Mean Age
age_m_files<- 'wb_pqlseq2_mean_age'
age_m_pqlseq<- import_pqlseq(age_m_files, v = "age")

#Sex
sex_files<- 'wb_pqlseq2_sex'
sex_pqlseq<- import_pqlseq(sex_files, v = "sex")

#Rename cols for each df to indicate variable
colnames(age_w_pqlseq)<- c("outcome", "length", "chrom", "chromStart", "chromEnd", "n", paste(names(age_w_pqlseq[,7:14]), "age", sep = "_"))
colnames(age_m_pqlseq)<- c(paste(names(age_m_pqlseq), "mean_age", sep = "_"))

#Cbind cols for age and sex dfs
pqlseq_model<- cbind(age_w_pqlseq, age_m_pqlseq[,7:14])

#Sort chromosome factors
sorted_labels<- str_sort(unique(pqlseq_model$chr), numeric=T)
pqlseq_model<- pqlseq_model %>%
  mutate(chr = factor(chr, levels = sorted_labels)) %>%
  arrange(chrom)
sex_pqlseq<- sex_pqlseq %>%
  mutate(chr = factor(chr, levels = sorted_labels[1:21])) %>%
  arrange(chr)
  
rm(age_w_pqlseq);rm(age_m_pqlseq)

pqlseq_model<- pqlseq_model %>%
  mutate(abs_diff = abs(beta_age) - abs(beta_mean_age))

pqlseq_model %>% 
  distinct(outcome, .keep_all = T) %>% 
  ggplot(aes(length)) + 
  geom_histogram(bins=100) +
  theme_classic() +
  scale_x_continuous(breaks = seq(0, max(pqlseq_full$length), by=1000))

#Plot pqlseq features------------------------------------------------------------
age.w.count<- nrow(pqlseq_model[pqlseq_model$fdr_age < 0.05 & pqlseq_model$chrom != "Y",])
age.m.count<- nrow(pqlseq_model[pqlseq_model$fdr_mean_age < 0.05 & pqlseq_model$chrom != "Y",])
sex.count<- nrow(sex_pqlseq[sex_pqlseq$fdr < 0.05,])
counts<- data.frame(count = c(age.w.count, age.m.count, sex.count),
                    predictor = as.factor(c('Age Within', 'Age Between', 'Sex')))
counts<- counts %>%
  mutate(predictor = as.factor(predictor)) %>%
  mutate(predictor=fct_reorder(predictor, count, .desc=T))

rm(age.w.count);rm(age.m.count);rm(sex.count)

counts %>%
  filter(predictor != "Sex") %>%
  ggplot(aes(predictor, count, fill = predictor)) +
  geom_bar(stat = 'identity', colour="black") +
  geom_text(label=counts$count[counts$predictor != "Sex"], vjust=-1, size =5) +
  theme_classic(base_size = 32) +
  xlab("Predictor") +
  ylab("Count") +
  scale_fill_manual(values = c("hotpink3", "darkgoldenrod2"))

## Significant regions Venn diagram
age.w<- pqlseq_model$outcome[pqlseq_model$fdr_age < 0.05 & pqlseq_model$chrom != "Y"]
age.m<- pqlseq_model$outcome[pqlseq_model$fdr_mean_age < 0.05 & pqlseq_model$chrom != "Y"]
sex<- sex_pqlseq$outcome[sex_pqlseq$fdr < 0.05]

venn_list<- list(age.w, age.m, sex)

rm(age.w);rm(age.m);rm(sex)

ggVennDiagram(venn_list, 
              category.names = c("Age Within", "Age Between", "Sex"),
              label = 'count',
              label_alpha = 0,
              label_size = 10,
              label_color = "black") +
  scale_x_continuous(expand = expansion(mult = .2)) +
  scale_fill_gradient(low = "white", high = "white")
ggsave('venn.svg', plot = last_plot(), width = 10, height = 10)

#Distribution of effect sizes and across the chromosomes
#Age
pqlseq_model %>%
  filter(fdr_age < 0.05) %>%
  ggplot(aes(beta_age, fill = beta_age>0)) +
  geom_histogram(bins=100, colour="black") +
  geom_vline(xintercept=0, linetype="dashed") +
  scale_fill_manual(values = c("slateblue3", "orangered2")) +
  xlab("Age Estimate") +
  ylab("Count") +
  theme_classic(base_size=32)

pqlseq_model %>%
  filter(fdr_age < .05 | fdr_mean_age < .05) %>%
  ggplot(aes(beta_mean_age, beta_age, colour = abs_diff)) +
  geom_point() +
  geom_abline() +
  geom_vline(xintercept=0, linetype="dashed") +
  geom_hline(yintercept=0, linetype="dashed") +
  scale_color_gradient2(low = "darkgoldenrod2", mid = "white", high = "hotpink3", midpoint = 0, name = "") +
  theme_classic(base_size=32) +
  theme(legend.key.height= unit(2, 'cm')) +
  ylab("Age Within") +
  xlab("Age Between")

log_age<- pqlseq_model %>%
  filter(fdr_age < .05 | fdr_mean_age < .05) %>%
  mutate(log_age = log(abs(beta_age)/abs(beta_mean_age))) %>%
  dplyr::select(log_age)
mean(log_age$log_age)

log_age %>%
  ggplot(aes(log_age, fill = log_age<0)) +
  geom_histogram(colour="black", position = "identity", bins = 100) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "red") +
  scale_fill_manual(values = c("hotpink3", "darkgoldenrod2")) +
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
  filter(fdr < 0.05) %>%
  ggplot(aes(beta, fill=beta>0)) +
  geom_histogram(bins=100, colour="black") +
  geom_vline(xintercept=0, linetype="dashed") +
  scale_fill_manual(values = c("darkolivegreen", "darkmagenta")) +
  xlab("Sex Estimate") +
  ylab("Count") +
  theme_classic(base_size=32)

df<- sex_pqlseq %>%
  dplyr::select(c(chrom, beta, fdr))
df$type<- "autosomes"
df$type[df$chrom == "X"]<- "X"

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
  filter(chrom != "Y") %>%
  dplyr::select(c(fdr_age, beta_age))
age_sex<- cbind(age_sex, sex_pqlseq[,c(9,12)])

age_sex$direction[age_sex$beta_age > 0 & age_sex$beta > 0]<- "Both Positive"
age_sex$direction[age_sex$beta_age > 0 & age_sex$beta < 0]<- "Age Positive, Sex Negative"
age_sex$direction[age_sex$beta_age < 0 & age_sex$beta > 0]<- "Age Negative, Sex Positive"
age_sex$direction[age_sex$beta_age < 0 & age_sex$beta < 0]<- "Both Negative"

age_sex<- age_sex %>%
  mutate(direction = as.factor(direction)) %>%
  mutate(direction=reorder(direction, direction, FUN=length))

age_sex %>%
  filter(fdr_age < 0.05 & fdr < 0.05) %>%
  ggplot(aes(beta_age, beta)) +
  geom_point(aes(colour=direction)) +
  geom_smooth(method = "lm") +
  theme_classic(base_size=32) +
  theme(legend.position = "none") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_colour_brewer(palette = "Set2") +
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
  scale_fill_brewer(palette = "Set2") +
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
setwd("/scratch/ckelsey4/Cayo_meth")
##Import annotation files-------------------------------------------------------
re_anno<- read_csv("re_annotations.csv")
chmm_intersect<- read_csv("chmm_annotations.csv")

#CHMM---------------------------------------------------------------------------
#Join pqlseq model and chmm 
pqlseq_chmm<- left_join(chmm_intersect, pqlseq_model, by = c("region_range", "chr"))
sex_pqlseq_chmm<- left_join(chmm_intersect, sex_pqlseq, by = c("region_range", "chr"))

#Filter out regions with models that didn't converge resulting in NAs in the annotation join
pqlseq_chmm<- pqlseq_chmm %>%
  drop_na()
sex_pqlseq_chmm<- sex_pqlseq_chmm %>%
  drop_na()

#Select out unnecessary cols and rename
pqlseq_chmm<- pqlseq_chmm %>%
  dplyr::select(-c(outcome, chromStart, chromEnd))
sex_pqlseq_chmm<- sex_pqlseq_chmm %>%
  dplyr::select(-c(outcome, chromStart, chromEnd))

#Set annotations as factor and reorder
annotations_ordered<- str_sort(unique(pqlseq_chmm$anno), numeric = TRUE)
pqlseq_chmm$anno<- factor(pqlseq_chmm$anno, levels = rev(annotations_ordered))
sex_pqlseq_chmm$anno<- factor(sex_pqlseq_chmm$anno, levels = rev(annotations_ordered))

#Create column of broad categories for annotations
pqlseq_chmm$chmm_class<- "A"
pqlseq_chmm$chmm_class[pqlseq_chmm$anno %in% annotations_ordered[1:2]]<- "Transcription Start Sites"
pqlseq_chmm$chmm_class[pqlseq_chmm$anno %in% annotations_ordered[3:5]]<- "Active Transcription"
pqlseq_chmm$chmm_class[pqlseq_chmm$anno %in% annotations_ordered[6:8]]<- "Enhancer Regions"
pqlseq_chmm$chmm_class[pqlseq_chmm$anno %in% annotations_ordered[9:15]]<- "Quiescent States"

#Create column of broad categories for annotations
sex_pqlseq_chmm$chmm_class<- "A"
sex_pqlseq_chmm$chmm_class[sex_pqlseq_chmm$anno %in% annotations_ordered[1:2]]<- "Transcription Start Sites"
sex_pqlseq_chmm$chmm_class[sex_pqlseq_chmm$anno %in% annotations_ordered[3:5]]<- "Active Transcription"
sex_pqlseq_chmm$chmm_class[sex_pqlseq_chmm$anno %in% annotations_ordered[6:8]]<- "Enhancer Regions"
sex_pqlseq_chmm$chmm_class[sex_pqlseq_chmm$anno %in% annotations_ordered[9:15]]<- "Quiescent States"

#Set classes as factors
class_factors<- c("Transcription Start Sites", "Active Transcription", "Enhancer Regions", "Quiescent States")
pqlseq_chmm$chmm_class<- factor(pqlseq_chmm$chmm_class, levels = rev(class_factors))
sex_pqlseq_chmm$chmm_class<- factor(sex_pqlseq_chmm$chmm_class, levels = rev(class_factors))

#REPEAT ELEMENTS----------------------------------------------------------------
#Join repeats annotations and glm_models df
pqlseq_re<- left_join(re_anno, pqlseq_model, by = c("region_range", "chr"))
sex_pqlseq_re<- left_join(re_anno, sex_pqlseq, by = c("region_range", "chr"))

pqlseq_re<- pqlseq_re %>%
  drop_na()
sex_pqlseq_re<- sex_pqlseq_re %>%
  drop_na()

#Select out unnecessary cols and rename
pqlseq_re<- pqlseq_re %>%
  dplyr::select(-c(range, outcome))
sex_pqlseq_re<- sex_pqlseq_re %>%
  dplyr::select(-c(range, outcome))

#Remove non-sensical annotations (NA, Unknown etc)
pqlseq_re<- pqlseq_re[!pqlseq_re$repClass == "Unknown",]
pqlseq_re<- pqlseq_re[!pqlseq_re$repClass == "DNA?",]
pqlseq_re<- pqlseq_re[!pqlseq_re$repClass == "LTR?",]
pqlseq_re<- pqlseq_re[!pqlseq_re$repClass == "RC?",]
pqlseq_re<- pqlseq_re[!pqlseq_re$repClass == "Unspecified",]

sex_pqlseq_re<- sex_pqlseq_re[!sex_pqlseq_re$repClass == "Unknown",]
sex_pqlseq_re<- sex_pqlseq_re[!sex_pqlseq_re$repClass == "DNA?",]
sex_pqlseq_re<- sex_pqlseq_re[!sex_pqlseq_re$repClass == "LTR?",]
sex_pqlseq_re<- sex_pqlseq_re[!sex_pqlseq_re$repClass == "RC?",]
sex_pqlseq_re<- sex_pqlseq_re[!sex_pqlseq_re$repClass == "Unspecified",]

#Plot annotation basics---------------------------------------------------------
#Age
pqlseq_re %>%
  filter(fdr_age < .05) %>%
  ggplot(aes(repClass, beta_age, colour=repClass)) +
  geom_jitter(alpha=0.5) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  theme_classic(base_size=32) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  xlab("Repeat Class") +
  ylab("Age Estimate")

pqlseq_chmm %>%
  filter(fdr_age < .05) %>%
  ggplot(aes(anno, beta_age, colour=anno)) +
  geom_jitter(alpha=0.5) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  theme_classic(base_size=32) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  xlab("Annotation") +
  ylab("Age Estimate") +
  coord_flip()

pqlseq_chmm %>%
  filter(fdr_age < .05) %>%
  ggplot(aes(chmm_class, beta_age, colour=chmm_class)) +
  geom_jitter(alpha=0.5) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  theme_classic(base_size=32) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  ylab("Age Estimate") +
  xlab("CHMM Class") +
  coord_flip()

#Sex
sex_pqlseq_re %>%
  filter(fdr < .05) %>%
  ggplot(aes(repClass, beta, colour=repClass)) +
  geom_jitter(alpha=0.5) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  theme_classic(base_size=32) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  xlab("Repeat Class") +
  ylab("Sex Estimate")

sex_pqlseq_re %>%
  filter(fdr < .05 & chr != "X") %>%
  ggplot(aes(repClass, beta, colour=repClass)) +
  geom_jitter(alpha=0.5) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  theme_classic(base_size=32) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  xlab("Repeat Class") +
  ylab("Sex Estimate")

sex_pqlseq_chmm %>%
  filter(fdr < .05) %>%
  ggplot(aes(anno, beta, colour=anno)) +
  geom_jitter(alpha=0.5) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  theme_classic(base_size=32) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  xlab("Annotation") +
  ylab("Sex Estimate") +
  coord_flip()

sex_pqlseq_chmm %>%
  filter(fdr < .05 & chr != "X") %>%
  ggplot(aes(anno, beta, colour=anno)) +
  geom_jitter(alpha=0.5) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  theme_classic(base_size=32) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  xlab("Annotation") +
  ylab("Sex Estimate") +
  coord_flip()

sex_pqlseq_chmm %>%
  filter(fdr < .05) %>%
  ggplot(aes(chmm_class, beta, colour=chmm_class)) +
  geom_jitter(alpha=0.5) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  theme_classic(base_size=32) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  ylab("Sex Estimate") +
  xlab("CHMM Class") +
  coord_flip()

#No x-chrom
sex_pqlseq_chmm %>%
  filter(fdr < .05 & chr != "X") %>%
  ggplot(aes(chmm_class, beta, colour=chmm_class)) +
  geom_jitter(alpha=0.5) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  theme_classic(base_size=32) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  ylab("Sex Estimate") +
  xlab("CHMM Class") +
  coord_flip()



######################################
###           ENRICHMENT           ###   
######################################
#Generate wide dataframe for enrichment (0s,1s for annotations)-----------------
#Select cols
enrich_df<- pqlseq_full %>% 
  dplyr::select(anno, outcome, beta_age, fdr_age, beta_sex, fdr_sex)
#Add col with 1s for pivot_wider values arg and pivot
enrich_df$a<- 1
enrich_df<- enrich_df %>%
  pivot_wider(names_from = anno, values_from = a)
#Assign 0 to NAs (instances where a region does not overlap an annotation)
enrich_df[is.na(enrich_df)]<- 0

chmm<- as.factor(annotations_ordered[1:15])

#Class enrichment---------------------------------------------
class_enrichment<- function(model_df, var_opt){
  
  df_list<- list()
  
  if (var_opt=="sex") {
    
    #filter model for negative sex estimates
    df<- model_df %>%
      filter(fdr_sex < 0.05)
    
    #create list of contingency tables for sex
    for(i in unique(df$class)){
      #Filter for specific class
      cl<- df[df$class == i,]
      
      #Counts for a) negative estimates and b) positive estimates with the focal class
      a<- nrow(cl[cl$beta_sex < 0,])
      b<- nrow(cl[cl$beta_sex > 0,])
      
      #Filter out regions that overlap the focal class
      non_cl<- df[!df$outcome %in% cl$outcome,]
      
      #Counts for a) negative estimates and b) positive estimates for all other classes
      c<- nrow(non_cl[non_cl$beta_sex < 0,])
      d<- nrow(non_cl[non_cl$beta_sex > 0,])
      
      #Generate contingency table - cols are class, rows are direction
      c_table<- data.frame(classY = c(a, b),
                           classN = c(c, d),
                           row.names = c("negY", "negN"))
      
      df_list[[length(df_list)+1]] = c_table
    } 
  } else if (var_opt=="age") {
    
    #filter model for negative sex estimates
    df<- model_df %>%
      filter(fdr_age < 0.05)
    
    #create list of contingency tables for sex
    for(i in unique(df$class)){
      #Filter for specific class
      cl<- df[df$class == i,]
      
      #Counts for a) negative estimates and b) positive estimates with the focal class
      a<- nrow(cl[cl$beta_age < 0,])
      b<- nrow(cl[cl$beta_age > 0,])
      
      #Filter out regions that overlap the focal class
      non_cl<- df[!df$outcome %in% cl$outcome,]
      
      #Counts for a) negative estimates and b) positive estimates for all other classes
      c<- nrow(non_cl[non_cl$beta_age < 0,])
      d<- nrow(non_cl[non_cl$beta_age > 0,])
      
      #Generate contingency table - cols are class, rows are direction
      c_table<- data.frame(classY = c(a, b),
                           classN = c(c, d),
                           row.names = c("negY", "negN"))
      
      df_list[[length(df_list)+1]] = c_table
    }
  } else if (var_opt=="both") {
    
    #filter model for negative sex estimates
    df<- df %>%
      filter(fdr_sex < 0.05 & fdr_sex < 0.05)
    
    #create list of contingency tables for sex
    for(i in unique(df$class)){
      #Filter for specific class
      cl<- df[df$class == i,]
      
      #Counts for a) negative estimates and b) all other estimates with the focal class
      a<- cl[cl$beta_age < 0 & cl$beta_sex < 0,]
      a_count<- nrow(a)
      b<- cl[!cl$outcome %in% a$outcome,]
      b_count<- nrow(b)
      
      #Filter out regions that overlap the focal class
      non_cl<- df[!df$outcome %in% cl$outcome,]
      
      #Counts for a) negative estimates and b) positive estimates for all other classes
      c<- non_cl[non_cl$beta_age < 0 & non_cl$beta_sex < 0,]
      c_count<- nrow(non_cl)
      d<- non_cl[!non_cl$outcome %in% c$outcome,]
      d_count<- nrow(d)
      
      #Generate contingency table - cols are class, rows are direction
      c_table<- data.frame(classY = c(a_count, b_count),
                           classN = c(c_count, d_count),
                           row.names = c("negY", "negN"))
      
      df_list[[length(df_list)+1]] = c_table
    }
  }
  
  #name table list
  names(df_list)<- unique(df$class)
  
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
}

#Age and Sex
age_sex_class<- class_enrichment(autosomes, "both")
age_sex_chmm$annotation<- factor(age_sex_chmm$annotation, levels = rev(annotations_ordered[1:15]))

age_sex_class %>%
  ggplot(aes(x=reorder(annotation, estimate), y=log_or, fill=annotation)) +
  geom_col(colour = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(ymin = age_sex_class$log_ci.lo, ymax = age_sex_class$log_ci.hi, width = 0.3) +
  scale_fill_brewer(palette = "Set3") +
  theme_classic(base_size = 32) +
  theme(legend.position = "none") +
  #ylim(c(-6, 2)) +
  ylab("Log Odds") +
  xlab("Annotation") +
  coord_flip()

sex_class<- class_enrichment(pqlseq_full, "sex")

sex_class %>%
  ggplot(aes(x=reorder(annotation, estimate), y=log_or, fill=annotation)) +
  geom_col(colour = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(ymin = sex_class$log_ci.lo, ymax = sex_class$log_ci.hi, width = 0.3) +
  scale_fill_brewer(palette = "Set3") +
  theme_classic(base_size = 32) +
  theme(legend.position = "none") +
  ylim(c(-2, 3)) +
  ylab("Log Odds") +
  xlab("Annotation") +
  coord_flip()

#chmm enrichment----------------------------------------------------------------
dml_enrichment<- function(model_df, var_opt){
  
  df_list<- list()
  
  df<- model_df
  
  chmm<- unique(df$anno)
  
  if (var_opt=="sex") {

    for (i in chmm) {
      
      #Counts for negative estimates
      a<- nrow(df[df$fdr < 0.05 & df$anno == i,])
      b<- nrow(df[df$fdr < 0.05 & df$anno != i,])
      
      #Counts positive estimates
      c<- nrow(df[df$fdr > 0.05 & df$anno == i,])
      d<- nrow(df[df$fdr > 0.05 & df$anno != i,])
      
      #Generate contingency table - cols are class, rows are direction
      c_table<- data.frame(negY = c(a, b),
                           negN = c(c, d),
                           row.names = c(paste(i, "Y"), paste(i,"N")))
      
      df_list[[length(df_list)+1]] = c_table
    } 
  } else if (var_opt=="age") {
    
    for(i in chmm){
      
      #Counts for negative estimates
      a<- nrow(df[df$fdr_age < 0.05 & df$anno == i,])
      b<- nrow(df[df$fdr_age < 0.05 & df$anno != i,])
      
      #Counts positive estimates
      c<- nrow(df[df$fdr_age > 0.05 & df$anno == i,])
      d<- nrow(df[df$fdr_age > 0.05 & df$anno != i,])
      
      #Generate contingency table - cols are class, rows are direction
      c_table<- data.frame(negY = c(a, b),
                           negN = c(c, d),
                           row.names = c(paste(i, "Y"), paste(i,"N")))
      print(c_table)
      
      df_list[[length(df_list)+1]] = c_table
      
    }
  } 
  
  #name table list
  names(df_list)<- chmm
  
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
}

age_dml<- dml_enrichment(pqlseq_chmm, "age")
age_dml$annotation<- factor(age_dml$annotation, levels = rev(annotations_ordered[1:15]))

age_dml %>%
  ggplot(aes(x=annotation, y=log_or, fill=padj<0.05)) +
  geom_col() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(ymin = age_dml$log_ci.lo, ymax = age_dml$log_ci.hi, width = 0.3) +
  theme_classic(base_size = 32) +
  theme(legend.position = "none") +
  ylim(c(-1, 6)) +
  ylab("Log Odds") +
  xlab("Annotation") +
  coord_flip()

directional_enrichment<- function(model_df, var_opt){
  
  df_list<- list()
  
  df<- model_df
  
  if (var_opt=="sex") {
    
    #filter model for negative sex estimates
    df<- df %>%
      filter(fdr < 0.05)
    
    #create list of contingency tables for sex
    for (i in unique(df$anno)) {
      
      #Counts for negative estimates
      a<- nrow(df[df$beta < 0 & df$anno == i,])
      b<- nrow(df[df$beta < 0 & df$anno != i,])
      
      #Counts positive estimates
      c<- nrow(df[df$beta > 0 & df$anno == i,])
      d<- nrow(df[df$beta > 0 & df$anno != i,])
      
      #Generate contingency table - cols are class, rows are direction
      c_table<- data.frame(negY = c(a, b),
                           negN = c(c, d),
                           row.names = c(paste(i, "Y"), paste(i,"N")))
      
      df_list[[length(df_list)+1]] = c_table
    } 
  } else if (var_opt=="age") {
    
    #filter model for negative sex estimates
    df<- df %>%
      filter(fdr_age < 0.05)
    
    #create list of contingency tables for sex
    for(i in unique(df$anno)){
      
      #Counts for negative estimates
      a<- nrow(df[df$beta_age < 0 & df$anno == i,])
      b<- nrow(df[df$beta_age < 0 & df$anno != i,])
      
      #Counts positive estimates
      c<- nrow(df[df$beta_age > 0 & df$anno == i,])
      d<- nrow(df[df$beta_age > 0 & df$anno != i,])
      
      #Generate contingency table - cols are class, rows are direction
      c_table<- data.frame(negY = c(a, b),
                           negN = c(c, d),
                           row.names = c(paste(i, "Y"), paste(i,"N")))
      print(c_table)
      
      df_list[[length(df_list)+1]] = c_table
     
    }
  } 
  
  #name table list
  names(df_list)<- chmm
  
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
}

#Age
age_chmm<- directional_enrichment(pqlseq_chmm, "age")
age_chmm$annotation<- factor(age_chmm$annotation, levels = rev(annotations_ordered[1:15]))
age_chmm<- age_chmm %>%
  filter(!estimate == Inf)

age_chmm %>%
  ggplot(aes(x=annotation, y=log_or, fill=padj<0.05)) +
  geom_col() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(ymin = age_chmm$log_ci.lo, ymax = age_chmm$log_ci.hi, width = 0.3) +
  theme_classic(base_size = 32) +
  theme(legend.position = "none") +
  ylim(c(-1, 6)) +
  ylab("Log Odds") +
  xlab("Annotation") +
  coord_flip()

#Sex
sex_chmm<- directional_enrichment(sex_pqlseq_chmm, "sex")
sex_chmm$annotation<- factor(sex_chmm$annotation, levels = rev(annotations_ordered[1:15]))

sex_chmm %>%
  ggplot(aes(x=annotation, y=log_or, fill=padj<0.05)) +
  geom_col() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(ymin = sex_chmm$log_ci.lo, ymax = sex_chmm$log_ci.hi, width = 0.3) +
  theme_classic(base_size = 32) +
  theme(legend.position = "none") +
  ylim(c(-4, 4)) +
  ylab("Log Odds") +
  xlab("Annotation") +
  coord_flip()

#TE enrichment------------------------------------------------------------------
te_enrichment<- function(model_df, var_opt){
  
  df_list<- list()
  
  df<- model_df
  
  te<- unique(df$repClass)
  
  if (var_opt=="sex") {
    
    #filter model for set fdr threshold
    df<- df %>%
      filter(fdr < 0.05)
    
    #generate list of contigency tables for each rep class
    for (i in te) {
      
      #Counts for negative estimates
      a<- nrow(df[df$beta < 0 & df$repClass == i,])
      b<- nrow(df[df$beta < 0 & df$repClass != i,])
      
      #Counts positive estimates
      c<- nrow(df[df$beta > 0 & df$repClass == i,])
      d<- nrow(df[df$beta > 0 & df$repClass != i,])
      
      #Generate contingency table - cols are class, rows are direction
      c_table<- data.frame(negY = c(a, b),
                           negN = c(c, d),
                           row.names = c(paste(i, "Y"), paste(i,"N")))
      print(c_table)
      
      df_list[[length(df_list)+1]] = c_table
    } 
  } else if (var_opt=="age") {
    
    #filter model for set fdr threshold
    df<- df %>%
      filter(fdr_age < 0.05)
    
    #generate list of contigency tables for each rep class
    for(i in te){
      
      #Counts for negative estimates
      a<- nrow(df[df$beta_age < 0 & df$repClass == i,])
      b<- nrow(df[df$beta_age < 0 & df$repClass != i,])
      
      #Counts positive estimates
      c<- nrow(df[df$beta_age > 0 & df$repClass == i,])
      d<- nrow(df[df$beta_age > 0 & df$repClass != i,])
      
      #Generate contingency table - cols are class, rows are direction
      c_table<- data.frame(negY = c(a, b),
                           negN = c(c, d),
                           row.names = c(paste(i, "Y"), paste(i,"N")))
      print(c_table)
      
      df_list[[length(df_list)+1]] = c_table
      
    }
  } 
  
  #name table list
  names(df_list)<- te
  
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
}

#Age
age_te<- te_enrichment(pqlseq_re, "age")

age_te %>%
  filter(!estimate == Inf) %>%
  ggplot(aes(x=reorder(annotation, estimate), y=log_or, fill=padj<0.05)) +
  geom_col(colour = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(ymin = age_te[!age_te$estimate == Inf,]$log_ci.lo, ymax = age_te[!age_te$estimate == Inf,]$log_ci.hi, width = 0.3) +
  scale_fill_brewer(palette = "Set3") +
  theme_classic(base_size = 32) +
  theme(legend.position = "none") +
  ylim(c(-4, 4)) +
  ylab("Log Odds") +
  xlab("Annotation") +
  coord_flip() +
  ggtitle("Age")

#Sex
sex_te<- te_enrichment(sex_pqlseq_re, "sex")
sex_te<- sex_te %>%
  filter(!conf.high == Inf)

sex_te %>%
  ggplot(aes(x=reorder(annotation, estimate), y=log_or, fill=padj<0.05)) +
  geom_col(colour = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(ymin = sex_te[!sex_te$estimate == Inf,]$log_ci.lo, ymax = sex_te[!sex_te$estimate == Inf,]$log_ci.hi, width = 0.3) +
  #scale_fill_brewer(palette = "Set3") +
  theme_classic(base_size = 32) +
  theme(legend.position = "none") +
  ylim(c(-3, 6)) +
  ylab("Log Odds") +
  xlab("Annotation") +
  coord_flip()

#Sex no X
no_x<- sex_pqlseq_re %>%
  filter(!chr == "X")
no_x_te<- te_enrichment(no_x, "sex")
no_x_te<- no_x_te %>%
  filter(!conf.high == Inf)

no_x_te %>%
  ggplot(aes(x=reorder(annotation, estimate), y=log_or, fill=padj<0.05)) +
  geom_col(colour = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(ymin = no_x_te[!no_x_te$estimate == Inf,]$log_ci.lo, ymax = no_x_te[!no_x_te$estimate == Inf,]$log_ci.hi, width = 0.3) +
  #scale_fill_brewer(palette = "Set3") +
  theme_classic(base_size = 32) +
  theme(legend.position = "none") +
  ylim(c(-3, 6)) +
  ylab("Log Odds") +
  xlab("Annotation") +
  coord_flip()


######################################
###        ADDITIVE EFFECTS        ###   
######################################
#Additive chmm------------------------------------------------------------------
add_chmm<- pqlseq_chmm %>%
  filter(!chr == "Y")
df<- sex_pqlseq_chmm[, 11:16]
colnames(df)<- paste(colnames(df), "sex" , sep = "_")
add_chmm<- cbind(add_chmm, df)
rm(df)

add_chmm<- add_chmm %>%
  filter(fdr_age < 0.05 & fdr_sex < 0.05 & !chr == "X") 

add_chmm<- add_chmm %>%
  dplyr::select(c(anno, cpg_loc, beta_age, beta_sex))

add_chmm<- add_chmm %>%
  pivot_longer(cols = c("beta_age", "beta_sex"), names_to = "type", values_to = "beta")

add_chmm %>%
  #filter(anno == "7_Enh") %>%
  ggplot(aes(type, beta, colour = beta < 0, group = cpg_loc)) +
  geom_point(alpha = 0.1) +
  geom_path(alpha = 0.1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic(base_size = 24) +
  scale_x_discrete(labels = c("Age", "Sex"))
  #facet_wrap(vars(anno), nrow=3)

add_chmm_no_x %>%
  ggplot(aes(beta, fill = type)) +
  geom_histogram(colour = "black", bins = 200, alpha = 0.7, position = "identity") +
  theme_classic(base_size = 24) +
  ylab("Count") +
  xlab("Estimate")

#Additive chmm------------------------------------------------------------------
add_re<- pqlseq_re %>%
  filter(!chr == "Y")
df<- sex_pqlseq_re[, 14:19]
colnames(df)<- paste(colnames(df), "sex" , sep = "_")
add_re<- cbind(add_re, df)
rm(df)

add_re<- add_re %>%
  filter(fdr_age < 0.05 & fdr_sex < 0.05)
add_re<- add_re %>%
  mutate(additive = beta_age + beta_sex)

add_re_no_x<- add_re %>%
  filter(!chr == "X") %>%
  dplyr::select(c(repClass, beta_age, beta_sex, additive))

add_re_no_x<- add_re_no_x %>%
  pivot_longer(cols = c("beta_age", "beta_sex", "additive"), names_to = "type", values_to = "beta")

add_re_no_x %>%
  ggplot(aes(beta, fill = type)) +
  geom_density(alpha = 0.7) +
  #geom_histogram(colour = "black", bins = 200, alpha = 0.7, position = "identity") +
  theme_classic(base_size = 24) +
  ylab("Count") +
  xlab("Estimate")

#Save workspace image-----------------------------------------------------------
save.image("wb_age_analysis.RData")


