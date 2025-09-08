library(tidyverse)
library(ggplot2)
library(ggcorrplot)
library(variancePartition)
library(ggvenn)

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

long_data %>%
  ggplot(aes(within.age)) +
  geom_density() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_classic(base_size=24)

long_data %>%
  ggplot(aes(mean.age)) +
  geom_density() +
  geom_vline(xintercept = median(long_data$mean.age), linetype = "dashed") +
  theme_classic(base_size=24)

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
  arrange(chr) %>%
  mutate(region_range = paste(chromStart, "-", chromEnd, sep = " ")) %>%
  relocate(region_range, .after = chromEnd)

#rm(age_w_pqlseq);rm(age_m_pqlseq);rm(chron_age_pqlseq)

######################################
###       Descriptive Stats        ###
######################################

age_full %>%
  select(c(beta_within_age, beta_chron_age, beta_long_cross, beta_short_cross)) %>%
  cor(use="pairwise.complete.obs") %>%
  ggcorrplot(show.diag=FALSE, type="lower", lab=TRUE, lab_size=5, sig.level = 0.05, insig = "blank")

#Import m/cov rds
regions_cov<- readRDS("/scratch/ckelsey4/Cayo_meth/regions_cov_filtered.rds")
regions_m<- readRDS("/scratch/ckelsey4/Cayo_meth/regions_m_filtered.rds")

regions_m<- do.call(rbind, regions_m)
regions_cov<- do.call(rbind, regions_cov)

p_meth<- as.data.frame(t(regions_m/regions_cov))

######################################
###       Plot Distributions       ###
######################################

compare_plot<- function(df, fdr1, fdr2, var1, var2, c1, c2, lab1, lab2) {
  
  df<- df %>%
    #filter({{fdr1}} < .05 & {{fdr2}} < .05) %>%
    mutate(diff = abs({{var1}}) - abs({{var2}}))
  
  eval(substitute(df_lm<- lm(var1 ~ var2, data=df)))
  
  #print(summary(df_lm))
    
  df %>%
    ggplot(aes({{var1}}, {{var2}}, colour = diff)) +
    geom_point() +
    geom_abline() +
    geom_abline(slope = df_lm[["coefficients"]][[2]], 
                intercept = df_lm[["coefficients"]][[1]],
                colour = "red") +
    geom_vline(xintercept=0, linetype="dashed") +
    geom_hline(yintercept=0, linetype="dashed") +
    scale_color_gradient2(low = c1, mid = "white", high = c2, midpoint = 0, name = "") +
    theme_classic(base_size=32) +
    theme(legend.key.height= unit(2, 'cm')) +
    xlab(lab1) +
    ylab(lab2) 

}

#Full CS set plots
compare_plot(age_full, fdr_long_cross, fdr_chron_age,
             beta_long_cross, beta_chron_age,
             "darkgoldenrod2", "steelblue2",
             "Full Cross Age", "Chronological Age")

compare_plot(age_full, fdr_long_cross, fdr_mean_age,
             beta_long_cross, beta_mean_age,
             "chartreuse3", "steelblue2",
             "Full Cross Age", "Mean Age")

compare_plot(age_full, fdr_long_cross, fdr_within_age,
             beta_long_cross, beta_within_age,
             "purple", "steelblue2",
             "Full Cross Age", "Within Age")

compare_plot(age_full, fdr_chron_age, fdr_within_age,
             beta_chron_age, beta_within_age,
             "purple", "darkgoldenrod2",
             "Chronological Age", "Within Age")

#CS Subset Plots
compare_plot(age_full, fdr_short_cross, fdr_within_age,
             beta_short_cross, beta_within_age,
             "purple", "hotpink3",
             "Short Cross Age", "Within Age")

compare_plot(age_full, fdr_short_cross, fdr_chron_age,
             beta_short_cross, beta_chron_age,
             "darkgoldenrod2", "hotpink3",
             "Short Cross Age", "Chronological Age")

compare_plot(age_full, fdr_short_cross, fdr_mean_age,
             beta_short_cross, beta_mean_age,
             "chartreuse3", "hotpink3",
             "Short Cross Age", "Mean Age")

age.w.count<- nrow(age_full[age_full$fdr_within_age < 0.05,])
age.chron.count<- nrow(age_full[age_full$fdr_chron_age < 0.05,])
age.short.count<- nrow(age_full[age_full$fdr_short_cross < 0.05,])
age.long.count<- nrow(age_full[age_full$fdr_long_cross < 0.05,])
counts<- data.frame(count = c(age.w.count, age.chron.count, age.short.count, age.long.count),
                    predictor = as.factor(c('Age Within', 'Age Chron', 'Age Short', 'Age Full')))

counts<- counts %>%
  mutate(predictor = as.factor(predictor)) %>%
  mutate(predictor=fct_reorder(predictor, count, .desc=T))
counts$predictor

rm(age.w.count);rm(age.m.count);rm(sex.count)

counts %>%
  ggplot(aes(predictor, count, fill = predictor)) +
  geom_bar(stat = 'identity', colour="black") +
  geom_text(label=counts$count, vjust=-1, size =5) +
  theme_classic(base_size = 32) +
  xlab("Predictor") +
  ylab("Count") +
  scale_fill_manual(values = c("steelblue2", "purple", "darkgoldenrod2", "hotpink3"))

## Significant regions Venn diagram
age.w<- age_full$outcome[age_full$fdr_within_age < 0.05]
age.chron<- age_full$outcome[age_full$fdr_chron_age < 0.05]
age.cross<- age_full$outcome[age_full$fdr_short_cross < 0.05]

venn_list<- list(age.w, age.chron, age.cross)

rm(age.w);rm(age.m);rm(sex)

ggvenn(venn_list,
       fill_color = c("purple", "hotpink3", "darkgoldenrod2"),
       text_size = 8,
       show_percentage = F)

#Proportions
agew<- pqlseq_anno %>% 
  group_by(agew_signif) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count))

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
pqlseq_prom<- left_join(promoters, age_full, by = c("region_range", "chr"))
pqlseq_prom<- pqlseq_prom %>%
  drop_na() %>%
  distinct(anno, .keep_all = T)
pqlseq_prom$anno_class<- "Promoter"
pqlseq_prom<- pqlseq_prom %>%
  dplyr::select(-c(chromStart, chromEnd, outcome))

#CHMM---------------------------------------------------------------------------
#Join pqlseq model and chmm 
pqlseq_chmm<- left_join(chmm_intersect, age_full, by = c("region_range", "chr"))

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
pqlseq_chmm$anno_class<- "TSSs"
pqlseq_chmm$anno_class[pqlseq_chmm$anno %in% annotations_ordered[3:5]]<- "Active Tr."
pqlseq_chmm$anno_class[pqlseq_chmm$anno %in% annotations_ordered[6:8]]<- "Enhancers"
pqlseq_chmm$anno_class[pqlseq_chmm$anno %in% annotations_ordered[9:15]]<- "Quiescent"

#Set classes as factors
class_factors<- c("TSSs", "Active Tr.", "Enhancers", "Quiescent")
pqlseq_chmm$anno_class<- factor(pqlseq_chmm$anno_class, levels = rev(class_factors))

#REPEAT ELEMENTS----------------------------------------------------------------
#Join repeats annotations and glm_models df
pqlseq_re<- left_join(re_anno, age_full, by = c("region_range", "chr"))

pqlseq_re<- pqlseq_re %>%
  drop_na()

#Remove non-sensical annotations (NA, Unknown etc)
pqlseq_re<- pqlseq_re[!pqlseq_re$repClass == "Unknown",]
pqlseq_re<- pqlseq_re[!pqlseq_re$repClass == "DNA?",]
pqlseq_re<- pqlseq_re[!pqlseq_re$repClass == "LTR?",]
pqlseq_re<- pqlseq_re[!pqlseq_re$repClass == "RC?",]
pqlseq_re<- pqlseq_re[!pqlseq_re$repClass == "Unspecified",]

pqlseq_re$anno_class<- "Simple Repeats"
pqlseq_re$anno_class[pqlseq_re$repClass %in% c("SINE", "LINE", "LTR", "Retroposon")]<- "TE Class I"
pqlseq_re$anno_class[pqlseq_re$repClass %in% "DNA"]<- "TE Class II"
pqlseq_re$anno_class[pqlseq_re$repClass %in% c("rRNA", "snRNA", "tRNA", "srpRNA", "scRNA")]<- "Str. RNAs"

pqlseq_re<- pqlseq_re %>%
  dplyr::rename(anno = repClass) %>%
  dplyr::select(-c(repName, chromStart, chromEnd)) %>%
  dplyr::relocate(anno, .after = anno_end)

#Select out unnecessary cols and rename
pqlseq_re<- pqlseq_re %>%
  dplyr::select(-c(range, outcome))

#Bind annotation dfs together---------------------------------------------------
pqlseq_anno<- rbind(pqlseq_chmm, pqlseq_re, pqlseq_prom)

annos<- c("Promoter", "TSSs", "Active Tr.", "Enhancers", "Quiescent", "Simple Repeats",
          "TE Class I", "TE Class II", "Str. RNAs")

pqlseq_anno$anno_source<- "Repeat Elements"
pqlseq_anno$anno_source[pqlseq_anno$anno_class == "Transcription Start Sites" | pqlseq_anno$anno_class == "Active Transcription" |
                          pqlseq_anno$anno_class == "Enhancer Regions" | pqlseq_anno$anno_class == "Quiescent States" | 
                          pqlseq_anno$anno_class == "Promoter"]<- "Transcription"

pqlseq_anno<- pqlseq_anno %>%
  arrange(anno_source, anno_class)

#Rearrange factors to sort by type then log_or
pqlseq_anno$anno_class<- factor(pqlseq_anno$anno_class, levels = rev(annos))

pqlseq_anno$unique_cpg<- paste(pqlseq_anno$chr, pqlseq_anno$cpg_loc, sep="_")

pqlseq_anno<- pqlseq_anno %>%
  dplyr::relocate(c(anno_class, anno_source, unique_cpg), .after=anno)


test<- pqlseq_anno %>%
  #filter(fdr_within_age < 0.05 & fdr_short_cross < 0.05) %>%
  dplyr::select(c(anno_class, beta_within_age, beta_short_cross)) %>%
  pivot_longer(cols = starts_with("beta"),
               names_to = "age_var",
               values_to = "beta")

test %>%
  ggplot(aes(beta, fill = age_var)) +
  geom_density(alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("hotpink2", "purple")) +
  theme_classic(base_size=24)

test %>%
  ggplot(aes(anno_class, beta, fill = age_var)) +
  geom_violin(position = position_dodge(width = 0.5), width = 2) +
  #geom_boxplot(, fill="white") + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("hotpink2", "purple")) +
  theme_classic(base_size=24) +
  theme(axis.text.x = element_text(angle = 45, hjust=1))


######################################
###    Cross-sectional Regions     ###
######################################

cs_anno<- pqlseq_anno[, c(1:13, 32:37)]

cs_anno$signif<- "Non-Significant"
cs_anno$signif[cs_anno$fdr_long_cross < 0.05 & cs_anno$beta_long_cross < 0]<- "Age-Hypomethylated"
cs_anno$signif[cs_anno$fdr_long_cross < 0.05 & cs_anno$beta_long_cross > 0]<- "Age-Hypermethylated"

cs_anno$signif<- factor(cs_anno$signif, 
                                 levels = c("Age-Hypermethylated", "Non-Significant", "Age-Hypomethylated"))

d1<- cs_anno %>% 
  distinct(unique_cpg, .keep_all = T) %>%
  group_by(anno_class, signif) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count))

d2<- cs_anno %>%
  mutate(unique_cpg = paste(chr, cpg_loc, sep="_")) %>%
  distinct(unique_cpg, .keep_all = T) %>%
  group_by(signif) %>%
  summarise(count = n()) %>%
  mutate(perc = count/sum(count))

d2$anno_class<- "All"

d3<- rbind(d1, d2)
annos2<- unique(d3$anno_class)
d3$anno_class<- factor(d3$anno_class, levels = annos2)

d3 %>%
  ggplot(aes(x = perc*100, y=anno_class, fill = factor(signif))) +
  geom_bar(stat="identity", width = 0.7, colour="black") +
  #geom_text(label=df$count, hjust=-5) +
  geom_vline(xintercept = 50, linetype = 'dashed') +
  geom_vline(xintercept = 41.4, linetype = 'dashed') +
  theme_classic(base_size=32) +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("hotpink", "gray90", "hotpink3")) +
  ylab("Annotation") +
  xlab("Percentage")






