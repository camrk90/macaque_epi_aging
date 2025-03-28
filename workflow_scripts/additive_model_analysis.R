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
long_data %>%
  mutate(pid = as.factor(pid))
long_data$pid<- reorder(long_data$pid, long_data$pid, FUN = length)

length(unique(long_data$monkey_id))

long_data %>%
  distinct(monkey_id, .keep_all = T) %>%
  ggplot(aes(individual_sex, fill = individual_sex)) +
  geom_bar(colour="black", alpha = 0.8) + 
  scale_fill_manual(values = c("royalblue1", "orangered1"), name = "Sex") +
  geom_text(stat = "count", aes(label = after_stat(count)),
            position = position_dodge(),
            color = "black",
            size = 8,
            vjust = -0.2) +
  ylim(0, 100) +
  theme_classic(base_size = 24)

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

#Generate %methylation matrix (m/cov)
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

av_meth_sum<- av_meth %>%
  group_by(monkey_id) %>%
  summarise(meth = mean(perc_meth),
            age = mean(age))
df<- av_meth %>%
  distinct(monkey_id, .keep_all = T) %>%
  dplyr::select(sex)

av_meth_sum<- cbind(av_meth_sum, df$sex)

av_meth %>%
  ggplot(aes(age, perc_meth)) +
  geom_point(aes(colour=sex)) +
  geom_smooth(aes(colour=sex),method = "loess") +
  geom_smooth(method = "loess", colour='black', linetype = "dashed", se = F) +
  scale_colour_manual(values = c("royalblue2", "orangered1"), name = "Sex") +
  scale_x_continuous(breaks = seq(0, 30, by=5)) +
  ylab("Average % Methylation") +
  xlab("Age") +
  theme_classic(base_size = 24)

av_meth_sum %>%
  ggplot(aes(age, meth)) +
  geom_point(aes(colour=df$sex)) +
  geom_smooth(aes(colour=df$sex),method = "glm") +
  geom_smooth(method = "glm", colour='black', linetype = "dashed", se = F) +
  scale_colour_manual(values = c("royalblue2", "orangered1"), name = "Sex") +
  scale_x_continuous(breaks = seq(0, 30, by=5)) +
  ylab("Average % Methylation") +
  xlab("Mean Age") +
  theme_classic(base_size = 24)

av_meth_diff<- av_meth %>% 
  group_by(monkey_id) %>%
  mutate(diff = max(perc_meth) - min(perc_meth)) %>%
  dplyr::select(c(monkey_id, sex, diff)) %>%
  distinct()

av_meth_diff %>%
  ggplot(aes(sex, diff)) +
  geom_dotplot(aes(fill = sex), binaxis = "y", stackdir = "center", dotsize = 0.7) +
  geom_boxplot(width=0.05) +
  scale_fill_manual(values = c("royalblue2", "orangered1"), name = "Sex") +
  ylab("%Methylation Change") +
  xlab("Sex") +
  theme_classic(base_size = 24)

t.test(diff ~ sex, data = av_meth_diff)

#Import regions info
regions<- readRDS("regions_full")
rownames(regions)<- str_split_i(rownames(regions), "\\.", 2)
regions<- regions_df[rownames(regions) %in% rownames(ratio),]

#Import genes-------------------------------------------------------------------
mm_genes<- rtracklayer::import('/scratch/ckelsey4/Cayo_meth/Macaca_mulatta.Mmul_10.110.chr.gtf')
mm_genes=as.data.frame(mm_genes)

######################################
###  Plot metadata distributions   ###
######################################
#Mean age by sex----------------------------------------------------------------
long_data %>% ggplot(aes(x=mean.age, fill = sex)) +
  geom_histogram(alpha=0.7, position = position_dodge(width = 1.5), bins=10, colour="black") +
  scale_fill_manual(values = c("royalblue2", "orangered1"), name = "Sex") +
  scale_x_continuous(breaks = seq(0, 30, by=5)) +
  theme_classic(base_size = 32) +
  labs(y = "Count", x = "Mean Age")

#N samples by sex---------------------------------------------------------------
long_data %>% ggplot(aes(x=n, fill = sex)) +
  geom_bar(alpha = 0.7, position = position_dodge(width = 0.5), colour = "black") +
  scale_fill_manual(values = c("royalblue2", "orangered1"), name = "Sex") +
  theme_classic(base_size = 32) +
  labs(y = "Count", x = "N")

#N samples per individual-------------------------------------------------------
long_data %>%
  ggplot(aes(x=age_at_sampling, y=reorder(monkey_id, mean.age), colour=sex)) +
  geom_path(linewidth = 1.5, alpha = 0.8) +
  geom_point(colour="black") +
  scale_x_continuous(breaks = seq(0, 30, by=2)) +
  scale_colour_manual(values = c("royalblue1", "orangered1"), name = "Sex") +
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
### Plot X-chromosome methylation  ###
######################################
x_meth<- ratio %>%
  filter(grepl("X", rownames(ratio))) %>%
  t() %>%
  as.data.frame()
colnames(x_meth)<- str_split_i(colnames(x_meth), "\\.", 4)
males<- long_data$lid_pid[long_data$sex == "M"]
females<- long_data$lid_pid[long_data$sex == "F"]
x_meth$sex<- "M"
x_meth$sex[rownames(x_meth) %in% females]<- "F"
x_meth<- x_meth %>%
  relocate(sex, .before = X_199_677)
x_meth_sum<- x_meth %>%
  group_by(sex) %>%
  summarize(across(where(is.numeric), median))
x_meth_sum<- as.data.frame(t(x_meth_sum))
colnames(x_meth_sum)<- c("female", "male")
x_meth_sum<- x_meth_sum[-1,]
x_meth_sum$region<- rownames(x_meth_sum)
x_meth_sum$region<- str_split_i(x_meth_sum$region, "\\_", 2)
x_meth_sum<- x_meth_sum %>%
  mutate_if(is.character, as.numeric)
test<- x_meth_sum %>%
  pivot_longer(!region, names_to = "sex", values_to = "median_methylation")

test %>%
  ggplot(aes(median_methylation, fill=sex)) +
  geom_histogram(bins=100, position = "identity", alpha=0.5, colour="black") +
  theme_classic(base_size = 30)

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

pcs<- cbind(pcs[1:10], dplyr::select(long_data, c("within.age", "mean.age", "sex", "pid")))

pcs %>% 
  ggplot(aes(PC1, PC2, colour=pid)) + 
  geom_point() +
  theme_classic(base_size=24)

pc_lm<- lm(PC1 ~ within.age + mean.age + sex + pid, data = pcs)
predicted_pcs<- predict.lm(pc_lm)
summary(pc_lm)

pcs<- cbind(pcs, predicted_pcs)

pcs %>%
  ggplot(aes(within.age, predicted_pcs)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_classic(base_size = 24)

pc.matrix<- model.matrix(~ PC1 + PC2 + within.age + sex + mean.age + pid, data = pcs)
pc.matrix[, 2:14] %>% 
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
import_pqlseq<- function(x){
  
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
    separate_wider_delim(outcome2, names=c("chrom", "chromStart", "chromEnd"), delim = "_") %>%
    relocate(c(chrom, chromStart, chromEnd), .after = outcome)
  
  #Add length col and filter by length
  model<- model %>%
    mutate(length = 1+(as.numeric(chromEnd) - as.numeric(chromStart))) %>%
    relocate(length, .after=outcome)
  model<- model %>%
    filter(!length > 5000)
  
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

#Age Within
age_w_files<- 'wb_pqlseq2_within_age'
age_w_pqlseq<- import_pqlseq(age_w_files)

#Mean Age
age_m_files<- 'wb_pqlseq2_mean_age'
age_m_pqlseq<- import_pqlseq(age_m_files)

#Sex
sex_files<- 'wb_pqlseq2_sex'
sex_pqlseq<- import_pqlseq(sex_files)

#Rename cols for each df to indicate variable
colnames(age_w_pqlseq)<- c("outcome", "length", "chrom", "chromStart", "chromEnd", "n", paste(names(age_w_pqlseq[,7:14]), "age", sep = "_"))
colnames(age_m_pqlseq)<- c(paste(names(age_m_pqlseq), "mean_age", sep = "_"))
colnames(sex_pqlseq)<- c(paste(names(sex_pqlseq), "sex", sep = "_"))

#Cbind cols for age and sex dfs
pqlseq_model<- cbind(age_w_pqlseq, age_m_pqlseq[,7:14], sex_pqlseq[,7:14])

#Sort chromosome factors
sorted_labels<- str_sort(unique(pqlseq_model$chrom), numeric=T)
pqlseq_model<- pqlseq_model %>%
  mutate(chrom = factor(chrom, levels = sorted_labels))

rm(age_w_pqlseq);rm(age_m_pqlseq);rm(sex_pqlseq)

pqlseq_model %>% 
  distinct(outcome, .keep_all = T) %>% 
  ggplot(aes(length)) + 
  geom_histogram(bins=100) +
  theme_classic() +
  scale_x_continuous(breaks = seq(0, max(pqlseq_full$length), by=1000))

#Plot pqlseq features------------------------------------------------------------
age.w.count<- nrow(pqlseq_model[pqlseq_model$fdr_age < 0.05,])
age.m.count<- nrow(pqlseq_model[pqlseq_model$fdr_mean_age < 0.05,])
sex.count<- nrow(pqlseq_model[pqlseq_model$fdr_sex < 0.05,])
counts<- data.frame(count = c(age.w.count, age.m.count, sex.count),
                    predictor = as.factor(c('Age Within', 'Age Between', 'Sex')))
counts<- counts %>%
  mutate(predictor = as.factor(predictor)) %>%
  mutate(predictor=fct_reorder(predictor, count, .desc=T))

rm(age.w.count);rm(age.m.count);rm(sex.count)

counts %>%
  ggplot(aes(predictor, count, fill = predictor)) +
  geom_bar(stat = 'identity', colour="black") +
  geom_text(label=counts$count, vjust=-1, size =5) +
  theme_classic(base_size = 32) +
  xlab("Predictor") +
  ylab("Count") +
  scale_fill_brewer(palette = "Set2")

## Significant regions Venn diagram
age.w<- pqlseq_model$outcome[pqlseq_model$fdr_age < 0.05]
age.m<- pqlseq_model$outcome[pqlseq_model$fdr_mean_age < 0.05]
sex<- pqlseq_model$outcome[pqlseq_model$fdr_sex < 0.05]

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

#Distribution of effect sizes
#Age
pqlseq_model %>%
  filter(fdr_age < 0.05) %>%
  ggplot(aes(beta_age, fill = beta_age>0)) +
  geom_histogram(bins=100, colour="black") +
  geom_vline(xintercept=0, linetype="dashed") +
  scale_fill_manual(values = c("goldenrod2", "orangered2")) +
  xlab("Age Estimate") +
  ylab("Count") +
  theme_classic(base_size=32)

#Age
test %>%
  filter(fdr_age < 0.05) %>%
  ggplot(aes(add, fill = add>0)) +
  geom_histogram(bins=100, colour="black") +
  geom_vline(xintercept=0, linetype="dashed") +
  scale_fill_manual(values = c("goldenrod2", "orangered2")) +
  xlab("Age Estimate") +
  ylab("Count") +
  theme_classic(base_size=32)
  
#Sex
pqlseq_model %>%
  filter(fdr_sex < 0.05) %>%
  ggplot(aes(beta_sex, fill=beta_sex>0)) +
  geom_histogram(bins=100, colour="black") +
  geom_vline(xintercept=0, linetype="dashed") +
  scale_fill_manual(values = c("darkolivegreen", "darkmagenta")) +
  xlab("Sex Estimate") +
  ylab("Count") +
  theme_classic(base_size=32)
  
#Distribution of coefficients for age and sex combined
age_sex<- pqlseq_model %>% 
  filter(fdr_age < 0.05 & fdr_sex < 0.05)

age_sex$direction[age_sex$beta_age > 0 & age_sex$beta_sex > 0]<- "Both Positive"
age_sex$direction[age_sex$beta_age > 0 & age_sex$beta_sex < 0]<- "Age Positive, Sex Negative"
age_sex$direction[age_sex$beta_age < 0 & age_sex$beta_sex > 0]<- "Age Negative, Sex Positive"
age_sex$direction[age_sex$beta_age < 0 & age_sex$beta_sex < 0]<- "Both Negative"

age_sex<- age_sex %>%
  mutate(direction = as.factor(direction)) %>%
  mutate(direction=reorder(direction, direction, FUN=length))

age_sex %>%
  ggplot(aes(beta_age, beta_sex)) +
  geom_point(aes(colour=direction)) +
  geom_abline() +
  geom_smooth(method = "lm") +
  theme_classic(base_size=32) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_colour_brewer(palette = "Set2")

age_sex %>%
  ggplot(aes(direction, fill=direction)) +
  geom_bar(colour="black") +
  geom_text(aes(label=after_stat(count)), stat="count", vjust=-1) +
  theme_classic(base_size = 32) +
  theme(axis.text.x=element_blank()) +
  scale_fill_brewer(palette = "Set2")

#Distribution of sex estimates by chromosome
#Age
pqlseq_model %>%
  filter(fdr_age < 0.05) %>%
  ggplot(aes(chrom, beta_age, colour=beta_age>0)) +
  geom_jitter(alpha=0.5) +
  geom_hline(yintercept = 0, linetype="dashed") +
  scale_colour_manual(values = c("goldenrod2", "orangered2")) +
  theme_classic(base_size = 32) +
  theme(legend.position = "none")

#Sex
pqlseq_model %>%
  filter(fdr_sex < 0.05) %>%
  ggplot(aes(chrom, beta_sex, colour=beta_sex>0)) +
  geom_jitter(alpha=0.5) +
  geom_hline(yintercept = 0, linetype="dashed") +
  scale_colour_manual(values = c("darkolivegreen", "darkmagenta")) +
  theme_classic(base_size = 32) +
  theme(legend.position = "none")

#Import glm random effects summaries--------------------------------------------
#Generate list of file names
glm_rfx_list<- list.files('/scratch/ckelsey4/Cayo_meth/glmer_model_compare', 
                          pattern = 'glmer_ranefx', 
                          full.names = T)

#Import glm models as list
glm_rfx<- lapply(glm_rfx_list, readRDS)
rfx_order<- gsub("/scratch/ckelsey4/Cayo_meth/glmer_model_compare/glmer_ranefx_", "", glm_rfx_list)
names(glm_rfx)<- rfx_order

#Generate list of file names
glm_sum_list<- list.files('/scratch/ckelsey4/Cayo_meth/glmer_model_compare', 
                          pattern = 'glmer_summaries', 
                          full.names = T)

#Import glm models as list
glm_sum<- lapply(glm_sum_list, readRDS)
sum_order<- gsub("/scratch/ckelsey4/Cayo_meth/glmer_model_compare/glmer_summaries_", "", glm_sum_list)
names(glm_sum)<- sum_order

sum_test<-glm_sum[[1]]

#TESTING WAYS TO COMBINE SLOPES INTO ONE DF
df1<- glm_rfx[[1]][[6]]
df1<- df1 %>% pivot_wider(names_from = term, values_from = c(condval, condsd))

#Remove intercept random effects
df1<- glm_rfx[[1]]
df1<- lapply(df1, function(x) filter(x, term == "age.w"))


################################################################################
##### JOIN INTERSECT FILES WITH REGIONS FROM MODELS
################################################################################
setwd("/scratch/ckelsey4/Cayo_meth")
##Import annotation files-------------------------------------------------------
re_anno<- read_csv("re_annotations.csv")
chmm_intersect<- read_csv("chmm_annotations.csv")

#CHMM---------------------------------------------------------------------------
#Join pqlseq model and chmm 
pqlseq_chmm<- left_join(chmm_intersect, pqlseq_model, by = "chromStart")
#Filter out regions with models that didn't converge resulting in NAs in the annotation join
pqlseq_chmm<- pqlseq_chmm %>%
  drop_na()

#Select out unnecessary cols and rename
pqlseq_chmm<- pqlseq_chmm %>%
  dplyr::select(-c(region_chr, chrom, chromEnd)) %>%
  dplyr::rename(region_start = chromStart)

#REPEAT ELEMENTS----------------------------------------------------------------
#Join repeats annotations and glm_models df
pqlseq_re<- left_join(re_anno, pqlseq_model, by = "chromStart")
pqlseq_re<- pqlseq_re %>%
  drop_na()
#Select out unnecessary cols and rename
pqlseq_re<- pqlseq_re %>%
  dplyr::select(-c(region_chr, range, chrom, chromEnd, repName)) %>%
  relocate(repClass, .after = anno_end) %>%
  dplyr::rename(region_start = chromStart,
                anno = repClass)

## Bind annotation dfs together-------------------------------------------------
#glm_full<- rbind(glm_chmm, glm_re)
pqlseq_full<- rbind(pqlseq_chmm, pqlseq_re)

clean_pqlseq<- function(x){
  
  #Remove non-sensical annotations (NA, Unknown etc)
  x<- x[!x$anno == "Unknown",]
  x<- x[!x$anno == "DNA?",]
  x<- x[!x$anno == "LTR?",]
  x<- x[!x$anno == "RC?",]
  x<- x[!x$anno == "Unspecified",]
  x<- x[!is.na(x$anno),]
  
  #Set annotations as factor and reorder
  annotations_ordered<- str_sort(unique(x$anno), numeric = TRUE)
  x$anno<- factor(x$anno, levels = annotations_ordered)
  
  #Create column of broad categories for annotations
  x$class<- "A"
  x$class[x$anno %in% annotations_ordered[16:length(annotations_ordered)]]<- "Repeat Elements"
  x$class[x$anno %in% annotations_ordered[1:2]]<- "Transcription Start Sites"
  x$class[x$anno %in% annotations_ordered[3:5]]<- "Active Transcription"
  x$class[x$anno %in% annotations_ordered[6:8]]<- "Enhancer Regions"
  x$class[x$anno %in% annotations_ordered[9:15]]<- "Quiescent States"
  
  #Set classes as factors
  x$class<- as.factor(x$class)
  
  #Add and plot region length col
  x<- x %>%
    mutate(length = 1+(region_end-region_start)) %>%
    relocate(length, .before=outcome)
}

pqlseq_full<- clean_pqlseq(pqlseq_full)
pqlseq_full<- pqlseq_full %>%
  group_by(outcome) %>%
  distinct(anno, .keep_all=T)

pqlseq_full %>% 
  #distinct(outcome, .keep_all = T) %>% 
  ggplot(aes(length)) + 
  geom_histogram(bins=100) +
  theme_classic() +
  scale_x_continuous(breaks = seq(0, max(pqlseq_full$length), by=1000))

## Plot class basics------------------------------------------------------------
#Age
pqlseq_full %>%
  filter(fdr_age < .05) %>%
  ggplot(aes(abs(beta_age), fill=class)) +
  geom_density(alpha = 0.5) +
  #geom_boxplot(width = 0.5) +
  #geom_violin(alpha = 0.5) +
  geom_hline(yintercept =0, linetype = 'dashed') +
  theme_classic(base_size=32) +
  scale_fill_brewer(palette = "Set3")

#Sex
pqlseq_full %>%
  filter(fdr_sex < .05) %>%
  ggplot(aes(beta_sex, fill=class)) +
  geom_density(alpha = 0.5) +
  #geom_boxplot(width = 0.5) +
  #geom_violin(alpha = 0.5) +
  geom_hline(yintercept =0, linetype = 'dashed') +
  theme_classic(base_size=32) +
  scale_fill_brewer(palette = "Set3")

### ENRICHMENT #################################################################
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
chmm_enrichment<- function(model_df, var_opt){
  
  df_list<- list()
  
  df<- model_df[,1:20]
  
  if (var_opt=="sex") {
    
    #filter model for negative sex estimates
    df<- df %>%
      filter(fdr_sex < 0.05)
    
    #create list of contingency tables for sex
    for(i in chmm){
      
      #Counts for negative estimates
      a<- nrow(df[df$beta_sex < 0 & df[,i] == 1,])
      b<- nrow(df[df$beta_sex < 0 & df[,i] == 0,])
      
      #Counts positive estimates
      c<- nrow(df[df$beta_sex > 0 & df[,i] == 1,])
      d<- nrow(df[df$beta_sex > 0 & df[,i] == 0,])
      
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
    for(i in chmm){
      
      #Counts for negative estimates
      a<- nrow(df[df$beta_age < 0 & df[,i] == 1,])
      b<- nrow(df[df$beta_age < 0 & df[,i] == 0,])
      
      #Counts positive estimates
      c<- nrow(df[df$beta_age > 0 & df[,i] == 1,])
      d<- nrow(df[df$beta_age > 0 & df[,i] == 0,])
      
      #Generate contingency table - cols are class, rows are direction
      c_table<- data.frame(negY = c(a, b),
                           negN = c(c, d),
                           row.names = c(paste(i, "Y"), paste(i,"N")))
      
      df_list[[length(df_list)+1]] = c_table
    }
  } else if (var_opt=="both") {
    
    #filter model for negative sex estimates
    df<- df %>%
      filter(fdr_sex < 0.05 & fdr_sex < 0.05)
    
    #create list of contingency tables for sex
    for(i in chmm){
      
      #Counts for negative estimates
      a<- nrow(df[df$beta_sex < 0 & df$beta_age < 0 & df[,i] == 1,])
      b<- nrow(df[df$beta_sex < 0 & df$beta_age < 0 & df[,i] == 0,])
      
      #Counts positive estimates
      c<- nrow(df[!df$beta_sex < 0 & !df$beta_age < 0 & df[,i] == 1,])
      d<- nrow(df[!df$beta_sex < 0 & !df$beta_age < 0 & df[,i] == 0,])
      
      #Generate contingency table - cols are class, rows are direction
      c_table<- data.frame(negY = c(a, b),
                           negN = c(c, d),
                           row.names = c("annoY", "annoN"))
      
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
age_chmm<- chmm_enrichment(enrich_df, "age")
age_chmm$annotation<- factor(age_chmm$annotation, levels = rev(annotations_ordered[1:15]))

age_chmm %>%
  ggplot(aes(x=annotation, y=log_or, colour=log_or<0, shape=padj<.05)) +
  geom_point(size=3) +
  #geom_line(aes(group=log_or<0)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(ymin = age_chmm$log_ci.lo, ymax = age_chmm$log_ci.hi, width = 0.3) +
  #scale_shape_manual(values = c("Signif" = 1, "Non-signif"=13)) +
  theme_classic(base_size = 32) +
  theme(legend.position = "none") +
  ylim(c(-2, 8)) +
  ylab("Log Odds") +
  xlab("Annotation") +
  coord_flip()

#Sex
sex_chmm<- chmm_enrichment(enrich_df, "sex")
sex_chmm$annotation<- factor(sex_chmm$annotation, levels = rev(annotations_ordered[1:15]))

sex_chmm %>%
  ggplot(aes(x=annotation, y=log_or, colour=log_or<0, shape=padj<.05, alpha = padj)) +
  geom_point(size=3) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(ymin = sex_chmm$log_ci.lo, ymax = sex_chmm$log_ci.hi, width = 0.3) +
  scale_fill_brewer(palette = "Set3") +
  theme_classic(base_size = 32) +
  theme(legend.position = "none") +
  ylim(c(-5, 5)) +
  ylab("Log Odds") +
  xlab("Annotation") +
  coord_flip()

#Age and Sex
age_sex_chmm<- chmm_enrichment(enrich_df, "both")
age_sex_chmm$annotation<- factor(age_sex_chmm$annotation, levels = rev(annotations_ordered[1:15]))

age_sex_chmm %>%
  filter(!estimate == Inf) %>%
  ggplot(aes(x=annotation, y=log_or, colour=log_or<0, shape=padj<.05)) +
  geom_point(size=3) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(ymin = age_sex_chmm[!age_sex_chmm$estimate == Inf,]$log_ci.lo, ymax = age_sex_chmm[!age_sex_chmm$estimate == Inf,]$log_ci.hi, width = 0.3) +
  scale_fill_brewer(palette = "Set3") +
  theme_classic(base_size = 32) +
  theme(legend.position = "none") +
  ylim(c(-6, 6)) +
  ylab("Log Odds") +
  xlab("Annotation") +
  coord_flip()

#TE enrichment------------------------------------------------------------------
te<- as.factor(annotations_ordered[16:29])
te_enrichment<- function(model_df, var_opt){
  
  df_list<- list()
  
  df<- model_df[,c(1:5, 21:34)]
  
  if (var_opt=="sex") {
    
    #filter model for negative sex estimates
    df<- df %>%
      filter(fdr_sex < 0.05)
    
    #create list of contingency tables for sex
    for(i in te){
      
      #Counts for negative estimates
      a<- nrow(df[df$beta_sex < 0 & df[,i] == 1,])
      b<- nrow(df[df$beta_sex < 0 & df[,i] == 0,])
      
      #Counts positive estimates
      c<- nrow(df[df$beta_sex > 0 & df[,i] == 1,])
      d<- nrow(df[df$beta_sex > 0 & df[,i] == 0,])
      
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
    for(i in te){
      
      #Counts for negative estimates
      a<- nrow(df[df$beta_age < 0 & df[,i] == 1,])
      b<- nrow(df[df$beta_age < 0 & df[,i] == 0,])
      
      #Counts positive estimates
      c<- nrow(df[df$beta_age > 0 & df[,i] == 1,])
      d<- nrow(df[df$beta_age > 0 & df[,i] == 0,])
      
      #Generate contingency table - cols are class, rows are direction
      c_table<- data.frame(negY = c(a, b),
                           negN = c(c, d),
                           row.names = c(paste(i, "Y"), paste(i,"N")))
      
      df_list[[length(df_list)+1]] = c_table
    }
  } else if (var_opt=="both") {
    
    #filter model for negative sex estimates
    df<- df %>%
      filter(fdr_sex < 0.05 & fdr_sex < 0.05)
    
    #create list of contingency tables for sex
    for(i in te){
      
      #Counts for negative estimates
      a<- nrow(df[df$beta_sex < 0 & df$beta_age < 0 & df[,i] == 1,])
      b<- nrow(df[df$beta_sex < 0 & df$beta_age < 0 & df[,i] == 0,])
      
      #Counts positive estimates
      c<- nrow(df[!df$beta_sex < 0 & !df$beta_age < 0 & df[,i] == 1,])
      d<- nrow(df[!df$beta_sex < 0 & !df$beta_age < 0 & df[,i] == 0,])
      
      #Generate contingency table - cols are class, rows are direction
      c_table<- data.frame(negY = c(a, b),
                           negN = c(c, d),
                           row.names = c("annoY", "annoN"))
      
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
age_te<- te_enrichment(enrich_df, "age")

age_te %>%
  filter(!estimate == Inf) %>%
  ggplot(aes(x=reorder(annotation, estimate), y=log_or, fill=annotation)) +
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
sex_te<- te_enrichment(enrich_df, "sex")

sex_te %>%
  filter(!estimate == Inf) %>%
  ggplot(aes(x=reorder(annotation, estimate), y=log_or, fill=annotation)) +
  geom_col(colour = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(ymin = sex_te[!sex_te$estimate == Inf,]$log_ci.lo, ymax = sex_te[!sex_te$estimate == Inf,]$log_ci.hi, width = 0.3) +
  scale_fill_brewer(palette = "Set3") +
  theme_classic(base_size = 32) +
  theme(legend.position = "none") +
  ylim(c(-4, 4)) +
  ylab("Log Odds") +
  xlab("Annotation") +
  coord_flip()

#Age and Sex
age_sex_te<- te_enrichment(enrich_df, "both")

age_sex_te %>%
  filter(!estimate == Inf) %>%
  ggplot(aes(x=annotation, y=log_or, colour=log_or<0, shape=padj<.05)) +
  geom_point(size=3) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(ymin = age_sex_te[!age_sex_te$estimate == Inf,]$log_ci.lo, ymax = age_sex_te[!age_sex_te$estimate == Inf,]$log_ci.hi, width = 0.3) +
  scale_fill_brewer(palette = "Set3") +
  theme_classic(base_size = 32) +
  theme(legend.position = "none") +
  ylim(c(-7, 4)) +
  ylab("Log Odds") +
  xlab("Annotation") +
  coord_flip()


#Save workspace image-----------------------------------------------------------
save.image("wb_age_analysis.RData")


