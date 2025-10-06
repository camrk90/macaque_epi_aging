library(tidyverse)
library(fgsea)
library(broom)
library(ggpubr)
library(ggcorrplot)
library(ggborderline)
library(ggrepel)
library(ggVennDiagram)
library(effectsize)
library(fgsea)
library(msigdbr)
options(scipen = 999)
setwd('/scratch/ckelsey4/Cayo_meth')
load("nested_analysis.RData")

######################################
###      Import PQLseq Model       ###
######################################
#Import nested pqlseq model-----------------------------------------------------
#Define import function
import_nested_models<- function(model_path, file_string){
  
  #model_path="/scratch/ckelsey4/Cayo_meth/glmer_model_compare"
  #file_string = "nested_model_f"
  print(model_path)
  print(file_string)
  
  file_list<- list.files(path = model_path, pattern = file_string)
  print(file_list)
  file_order<- str_split_i(file_list, "_", 4)
  
  #Import glm models as list
  model_list<- lapply(paste(model_path, file_list, sep = "/"), readRDS)
  
  #Rename list elements
  names(model_list)<- file_list
  
  #Bind model list of models to dataframe
  nested_pqlseq<- do.call(rbind, model_list)
  
  #Split outcome col to remove repeated chr #'s
  nested_pqlseq$outcome<- str_split_i(nested_pqlseq$outcome, "\\.", 3)
  
  #Generate sex col from rownames
  nested_pqlseq$sex<- rownames(nested_pqlseq)
  nested_pqlseq$sex<- str_split_i(nested_pqlseq$sex, "\\_", 3)
  
  #Separate region coordinates into start and end, delete the chr col, and move region col to front
  nested_pqlseq$outcome2<- nested_pqlseq$outcome
  nested_pqlseq<- nested_pqlseq %>%
    separate_wider_delim(outcome2, "_", names = c("chr", "chromStart", "chromEnd")) %>%
    relocate(c(sex, chr, chromStart, chromEnd), .before = n) %>%
    mutate(region_range = paste(as.character(chromStart), "-", as.character(chromEnd)))
  
  #Filter for true convergences
  nested_pqlseq<- nested_pqlseq %>%
    filter(converged == "TRUE")
  
  #Generate df of adjusted pvalues
  nested_pqlseq_fdr<- p.adjust(nested_pqlseq$pvalue, method = "fdr")
  
  #Bind padj cols to model df and relocate
  nested_pqlseq<- cbind(nested_pqlseq, nested_pqlseq_fdr)
  nested_pqlseq<- nested_pqlseq %>%
    dplyr::rename(fdr = nested_pqlseq_fdr) %>%
    relocate(fdr, .after = pvalue) %>%
    dplyr::select(-c(elapsed_time, converged))
  rm(nested_pqlseq_fdr)
  
  return(nested_pqlseq)
  
}

#Import female and male pqlseq output
f_nested<- import_nested_models(model_path = "/scratch/ckelsey4/Cayo_meth/glmer_model_compare",
                                file_string = "nested_model_f")

m_nested<- import_nested_models(model_path = "/scratch/ckelsey4/Cayo_meth/glmer_model_compare",
                                file_string = "nested_model_m_")

nested_pqlseq<- rbind(f_nested, m_nested)
nested_pqlseq<- nested_pqlseq %>%
  mutate(std_beta = beta/se_beta) %>%
  relocate(std_beta, .before = beta)

nested_pqlseq$type<- "autosomes"
nested_pqlseq$type[nested_pqlseq$chr == "X"]<- "X"

#Import genes-------------------------------------------------------------------
mm_genes<- rtracklayer::import('/scratch/ckelsey4/Cayo_meth/Macaca_mulatta.Mmul_10.110.chr.gtf')
mm_genes=as.data.frame(mm_genes)
mm_genes<- mm_genes %>%
  filter(type == "gene")

#Plot nested pqlseq model-------------------------------------------------------
nested_pqlseq$direction<- "Hyper"
nested_pqlseq$direction[nested_pqlseq$beta < 0]<- "Hypo"

f.count_hyper<- nrow(f_nested[f_nested$beta > 0 & f_nested$fdr < 0.05,])
f.count_hypo<- nrow(f_nested[f_nested$beta < 0 & f_nested$fdr < 0.05,])
m.count_hyper<- nrow(m_nested[m_nested$beta > 0 & m_nested$fdr < 0.05,])
m.count_hypo<- nrow(m_nested[m_nested$beta < 0 & m_nested$fdr < 0.05,])
counts<- data.frame(count = c(f.count_hyper, m.count_hyper, f.count_hypo, m.count_hypo),
                    predictor = as.factor(c('Female Hyper', 'Male Hyper', 'Female Hypo', 'Male Hypo')),
                    sex = as.factor(c("Female", "Male", "Female", "Male")))

rm(f.count_hyper);rm(f.count_hypo);rm(m.count_hyper);rm(m.count_hypo)

nested_pqlseq %>%
  filter(fdr < 0.05) %>%
  ggplot(aes(x=direction, fill=sex)) +
  geom_bar(alpha = 0.7, position = position_dodge(width = 0.7),  colour = "black") +
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -1, position = position_dodge(width = 0.7)) +
  scale_fill_manual(values = c("royalblue1", "orangered1"), name = "Sex") +
  theme_classic(base_size = 32) +
  xlab("Predictor") +
  ylab("Count")

#Plot significant estimates
nested_pqlseq %>%
  filter(fdr < 0.05) %>%
  filter(chr == "X") %>%
  ggplot(aes(beta, fill=sex)) +
  geom_histogram(alpha = 0.7, colour="black", bins = 100, position = "identity") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "red") +
  scale_fill_manual(values = c("darkolivegreen", "darkmagenta"), name = "Sex") +
  ylab("Count") +
  xlab("Estimate") +
  theme_classic(base_size=24)
  #facet_wrap(vars(type), nrow = 2)

#Plot significant hypomethylated estimates
nested_pqlseq %>%
  filter(beta < 0 & fdr < 0.05) %>%
  ggplot(aes(beta, fill=sex)) +
  geom_histogram(alpha = 0.5, colour="black", bins = 100, position = "identity") + 
  scale_fill_manual(values = c("royalblue2", "orangered1"), name = "Sex") +
  ylab("Count") +
  xlab("Estimate") +
  theme_classic(base_size = 30)

#Plot significant hypermethylated estimates
nested_pqlseq %>%
  filter(beta > 0 & fdr < 0.05) %>%
  ggplot(aes(beta, fill=sex)) +
  geom_histogram(alpha = 0.5, colour="black", bins = 100, position = "identity") + 
  scale_fill_manual(values = c("darkolivegreen", "darkmagenta"), name = "Sex") +
  ylab("Count") +
  xlab("Estimate") +
  theme_classic(base_size = 30)

#Boxplot of hypermethylated significant estimates
nested_pqlseq %>%
  filter(beta > 0 & fdr < 0.05) %>%
  ggplot(aes(sex, beta)) +
  #geom_jitter(aes(colour = sex), alpha = 0.7)+
  geom_boxplot(aes(fill = sex), width=0.3) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("darkolivegreen", "darkmagenta"), name = "Sex") +
  ylab("Beta (FDR < 0.05)") +
  xlab("Sex") +
  theme_classic(base_size = 30) +
  theme(legend.position = "none")

#Boxplot of hyp0methylated significant estimates
nested_pqlseq %>%
  filter(beta < 0 & fdr < 0.05) %>%
  ggplot(aes(sex, beta)) +
  #geom_jitter(aes(colour = sex), alpha = 0.7)+
  geom_boxplot(aes(fill = sex), width=0.3) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("darkolivegreen", "darkmagenta"), name = "Sex") +
  ylab("Beta (FDR < 0.05)") +
  xlab("Sex") +
  theme_classic(base_size = 30) +
  theme(legend.position = "none")

#Plot estimates by chromosome by sex
df<- nested_pqlseq %>%
  dplyr::select(c(chrom, beta, fdr, sex))
df$type<- "autosomes"
df$type[df$chrom == "X"]<- "X"

df %>%
  filter(fdr < 0.05) %>%
  ggplot(aes(type, beta, colour=sex)) +
  geom_jitter(alpha=0.5) +
  geom_hline(yintercept = 0, linetype="dashed") +
  scale_colour_manual(values = c("royalblue2", "orangered1"), name = "Sex") +
  theme_classic(base_size = 32) +
  ylab("Estimate (FDR < .05)") +
  xlab("")

#Generate and plot paired nested estimates--------------------------------------
sex<- nested_pqlseq %>%
  dplyr::select(c(outcome, region_range, chr, std_beta, beta, se_beta, fdr, sex, type)) %>%
  pivot_wider(names_from = sex, values_from = c(std_beta, beta, fdr, se_beta), names_sep = "_")

sex<- sex %>%
  mutate(abs_diff = abs(beta_f) - abs(beta_m),
         std_diff = abs(std_beta_f) - abs(std_beta_m))

sex_hypo<- sex %>%
  filter(beta_m < 0 & beta_f < 0)
sex_hyper<- sex %>%
  filter(beta_m > 0 & beta_f > 0)

#Effect size distribution t-test
hypo_t<- sex_hypo %>%
  filter(fdr_f < 0.05 | fdr_m < 0.05)
t.test(hypo_t$beta_m, hypo_t$beta_f, paired = T)
cohens_d(hypo_t$beta_f, hypo_t$beta_m)
mean(hypo_t$beta_m)
mean(hypo_t$beta_f)

hyper_t<- sex_hyper %>%
  filter(fdr_f < 0.05 | fdr_m < 0.05)
t.test(hyper_t$beta_m, hyper_t$beta_f, paired = T)
cohens_d(hyper_t$beta_f, hyper_t$beta_m)
mean(hyper_t$beta_m)
mean(hyper_t$beta_f)

#Distribution of difference in abs estimates for hypomethylated regions with x-chrom
sex_hypo %>%
  filter(fdr_f < 0.05 | fdr_m < 0.05) %>%
  ggplot(aes(abs_diff, fill = after_stat(x))) +
  geom_histogram(bins = 50, colour = "black", position = 'identity') +
  geom_vline(xintercept = 0, linetype = "dashed", colour = 'red') +
  scale_fill_gradient2(low = "darkmagenta", mid = "white", high = "darkolivegreen", midpoint = 0, name = "Beta Difference") +
  xlab("|Estimate (F)| - |Estimate (M)|") +
  ylab("Hypomethylated Regions") +
  theme_classic(base_size = 36) +
  theme(legend.position = "none")

#Distribution of difference in abs estimates for hypomethylated regions without x-chrom
sex_hypo %>%
  filter(fdr_f < 0.05 | fdr_m < 0.05) %>%
  filter(chr != "X") %>%
  ggplot(aes(abs_diff, fill = after_stat(x))) +
  geom_histogram(bins = 100, colour = "black", position = 'identity') +
  geom_vline(xintercept = 0, linetype = "dashed", colour = 'red') +
  scale_fill_gradient2(low = "darkmagenta", mid = "white", high = "darkolivegreen", midpoint = 0, name = "Beta Difference") +
  xlab("|Estimate (F)| - |Estimate (M)|") +
  ylab("Hypomethylated Regions No X") +
  theme_classic(base_size = 30) +
  theme(legend.position = "none")

#Distribution of difference in abs estimates with X-chrom
sex_hyper %>%
  filter(fdr_f < 0.05 | fdr_m < 0.05) %>%
  #filter(abs_diff > -0.2 & abs_diff < 0.2) %>%
  ggplot(aes(abs_diff, fill = after_stat(x))) +
  geom_histogram(bins = 50, colour = "black", position = 'identity') +
  geom_vline(xintercept = 0, linetype = "dashed", colour = 'red') +
  scale_fill_gradient2(low = "darkmagenta", mid = "white", high = "darkolivegreen", midpoint = 0, name = "Beta Difference") +
  xlab("|Estimate (F)| - |Estimate (M)|") +
  ylab("Hypermethylated Regions") +
  theme_classic(base_size = 36) +
  theme(legend.position = "none")

#Distribution of difference in abs estimates without X-chrom
sex_hyper %>%
  filter(fdr_f < 0.05 | fdr_m < 0.05) %>%
  filter(chr != "X") %>%
  ggplot(aes(abs_diff, fill = after_stat(x))) +
  geom_histogram(bins = 100, colour = "black", position = 'identity') +
  geom_vline(xintercept = 0, linetype = "dashed", colour = 'red') +
  scale_fill_gradient2(low = "darkmagenta", mid = "white", high = "darkolivegreen", midpoint = 0, name = "Beta Difference") +
  xlab("|Estimate (F)| - |Estimate (M)|") +
  ylab("Hypermethylated Regions No X") +
  theme_classic(base_size = 30) +
  theme(legend.position = "none")

#Scatterplot of male vs female estimates
sex %>%
  filter(fdr_f < 0.05 | fdr_m < 0.05) %>%
  ggplot(aes(beta_m, beta_f, colour = abs_diff)) +
  geom_point() +
  geom_abline() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_gradient2(low = "darkmagenta", mid = "white", high = "darkolivegreen", midpoint = 0, name = "") +
  geom_smooth(method = "lm") +
  #ylim(-0.5, 0.5) +
  #xlim(-0.5, 0.5) +
  xlab("Estimate (Male)") +
  ylab("Estimate (Female)") +
  theme_classic(base_size=24) +
  theme(legend.key.height= unit(2, 'cm')) +
  facet_wrap(vars(type))

sex %>%
  filter(fdr_f < 0.05 | fdr_m < 0.05) %>%
  ggplot(aes(beta_m, beta_f, colour = abs_diff)) +
  geom_point() +
  geom_abline() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_gradient2(low = "darkmagenta", mid = "white", high = "darkolivegreen", midpoint = 0, name = "") +
  geom_smooth(method = "lm") +
  #ylim(-0.5, 0.5) +
  #xlim(-0.5, 0.5) +
  xlab("Estimate (Male)") +
  ylab("Estimate (Female)") +
  theme_classic(base_size=24) +
  theme(legend.key.height= unit(2, 'cm')) +
  facet_wrap(vars(type))

df3<- sex %>%
  filter(fdr_f < 0.05 | fdr_m < 0.05)

cor.test(df3$beta_f, df3$beta_m)

######################################
###        Add Annotations         ###
######################################
#Import annotation files--------------------------------------------------------
re_anno<- read_csv("re_annotations.csv")
chmm_intersect<- read_csv("chmm_annotations.csv")
promoters<- read_csv("promoters.csv")

annotation_labels<- c("Active TSS", "Flanking Active TSS", "Transcr. at 5' & 3'", "Strong Transcription", "Weak Transcription", "Enhancers", "Genic Enhancers",
                      "ZNF Genes & Repeats", "Heterochromatin", "Bivalent/Poised TSS", "Flanking Bivalent TSS/Enh", "Bivalent Enhancer", "Repressed Polycomb",
                      "Weak Repressed Polycomb", "Quiescent/Low")
#Promoters----------------------------------------------------------------------
m_prom<- left_join(promoters, m_nested, by = c("region_range", "chr"))
m_prom<- m_prom %>%
  drop_na() %>%
  distinct(anno, .keep_all = T)
m_prom$anno_class<- "Promoter"

f_prom<- left_join(promoters, f_nested, by = c("region_range", "chr"))
f_prom<- f_prom %>%
  drop_na() %>%
  distinct(anno, .keep_all = T)
f_prom$anno_class<- "Promoter"

#Annotations--------------------------------------------------------------------
#Define join function
anno_join<- function(output_df, annotation_df, annotation_type){
  
  if (annotation_type == "chmm"){
    
    #Join pqlseq model and chmm 
    model_df<- left_join(annotation_df, output_df, by = c("region_range", "chr"))
    
    #Filter out regions with models that didn't converge resulting in NAs in the annotation join
    model_df<- model_df %>%
      drop_na()
    
    #Set annotations as factor and reorder
    annotations_ordered<- str_sort(unique(model_df$anno), numeric = TRUE)
    annotation_labels<- c("Active TSS", "Flanking Active TSS", "Transcr. at 5' & 3'", "Strong Transcription", "Weak Transcription", "Enhancers", "Genic Enhancers",
                            "ZNF Genes & Repeats", "Heterochromatin", "Bivalent/Poised TSS", "Flanking Bivalent TSS/Enh", "Bivalent Enhancer", "Repressed Polycomb",
                            "Weak Repressed Polycomb", "Quiescent/Low")
    model_df$anno<- factor(model_df$anno, levels = annotations_ordered, labels = annotation_labels)
    
    #Create column of broad categories for annotations
    model_df$anno_class<- "A"
    model_df$anno_class[model_df$anno %in% annotation_labels[1:2]]<- "Transcription Start Sites"
    model_df$anno_class[model_df$anno %in% annotation_labels[3:5]]<- "Active Transcription"
    model_df$anno_class[model_df$anno %in% annotation_labels[6:8]]<- "Enhancer Regions"
    model_df$anno_class[model_df$anno %in% annotation_labels[9:15]]<- "Quiescent States"
    
    #Set classes as factors
    class_factors<- c("Transcription Start Sites", "Active Transcription", "Enhancer Regions", "Quiescent States")
    model_df$anno_class<- factor(model_df$anno_class, levels = class_factors)
    
    return(model_df)
    
  } else if (annotation_type == "re") {
    
    #Join repeats annotations and glm_models df
   model_df<- left_join(annotation_df, output_df, by = c("region_range", "chr"))
    
   model_df<-model_df %>%
      drop_na()
    
    #Select out unnecessary cols and rename
   model_df<-model_df %>%
      dplyr::select(-c(range))
    
    #Remove non-sensical annotations (NA, Unknown etc)
   model_df<-model_df[!model_df$repClass == "Unknown",]
   model_df<-model_df[!model_df$repClass == "DNA?",]
   model_df<-model_df[!model_df$repClass == "LTR?",]
   model_df<-model_df[!model_df$repClass == "RC?",]
   model_df<-model_df[!model_df$repClass == "Unspecified",]
   
   #Set annotations as factor and reorder
   model_df$repClass<- as.factor(model_df$repClass)
   
   #Create column of broad categories for annotations
   model_df$anno_class<- "Simple Repeats"
   model_df$anno_class[model_df$repClass %in% c("SINE", "LINE", "LTR", "Retroposon")]<- "Transposable Elements Class I"
   model_df$anno_class[model_df$repClass %in% "DNA"]<- "Transposable Elements Class II"
   model_df$anno_class[model_df$repClass %in% c("rRNA", "snRNA", "tRNA", "srpRNA", "scRNA")]<- "Structural RNAs"
   
   #Set classes as factors
   class_factors<- c("Simple Repeats", "Transposable Elements Class I", "Transposable Elements Class II", "Structural RNAs")
   model_df$anno_class<- factor(model_df$anno_class, levels = class_factors)
   
   #Arrange cols to match chmm df
   model_df<- model_df %>%
     dplyr::select(-repName) %>%
     dplyr::rename(anno = repClass) %>%
     relocate(c(anno_start, anno_end), .before = anno)
   
   return(model_df)
  }
}

m_chmm<- anno_join(m_nested, chmm_intersect, "chmm")
m_re<- anno_join(m_nested, re_anno, "re")
m_full<- rbind(m_chmm, m_re, m_prom)

annos<- c("Promoter", "Transcription Start Sites", "Active Transcription", "Enhancer Regions", "Quiescent States", "Simple Repeats",
          "Transposable Elements Class I", "Transposable Elements Class II", "Structural RNAs")

m_full$anno_source<- "Repeat Elements"
m_full$anno_source[m_full$anno_class == "Transcription Start Sites" | m_full$anno_class == "Active Transcription" |
                          m_full$anno_class == "Enhancer Regions" | m_full$anno_class == "Quiescent States" | 
                          m_full$anno_class == "Promoter"]<- "Transcription"

m_full<- m_full %>%
  arrange(anno_source, anno_class)

#Rearrange factors to sort by type then log_or
m_full$anno_class<- factor(m_full$anno_class, levels = rev(annos))

m_full$signif<- "Non-Significant"
m_full$signif[m_full$fdr < 0.05 & m_full$beta < 0]<- "Age-Hypomethylated"
m_full$signif[m_full$fdr < 0.05 & m_full$beta > 0]<- "Age-Hypermethylated"

m_full$signif<- factor(m_full$signif, levels=c("Age-Hypermethylated", "Non-Significant", "Age-Hypomethylated"))

m_full<- m_full %>%
  mutate(unique_cpg = paste(chr, cpg_loc, sep="_"))

f_chmm<- anno_join(f_nested, chmm_intersect, "chmm")
f_re<- anno_join(f_nested, re_anno, "re")
f_full<- rbind(f_chmm, f_re, f_prom)

f_full$anno_source<- "Repeat Elements"
f_full$anno_source[f_full$anno_class == "Transcription Start Sites" | f_full$anno_class == "Active Transcription" |
                     f_full$anno_class == "Enhancer Regions" | f_full$anno_class == "Quiescent States" | 
                     f_full$anno_class == "Promoter"]<- "Transcription"

f_full<- f_full %>%
  arrange(anno_source, anno_class)

#Rearrange factors to sort by type then log_or
f_full$anno_class<- factor(f_full$anno_class, levels = rev(annos))

f_full$signif<- "Non-Significant"
f_full$signif[f_full$fdr < 0.05 & f_full$beta < 0]<- "Age-Hypomethylated"
f_full$signif[f_full$fdr < 0.05 & f_full$beta > 0]<- "Age-Hypermethylated"

f_full$signif<- factor(f_full$signif, levels=c("Age-Hypermethylated", "Non-Significant", "Age-Hypomethylated"))

f_full<- f_full %>%
  mutate(unique_cpg = paste(chr, cpg_loc, sep="_"))

#Plot basics--------------------------------------------------------------------
m_proportions<- m_full %>% 
  group_by(anno_class, signif) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count))

d3<- m_full %>%
  mutate(unique_cpg = paste(chr, cpg_loc, sep="_")) %>%
  distinct(unique_cpg, .keep_all = T) %>%
  group_by(signif) %>%
  summarise(count = n()) %>%
  mutate(perc = count/sum(count))

d3$anno_class<- "All"

m_proportions<- rbind(m_proportions, d3)
annos2<- unique(m_proportions$anno_class)
m_proportions$anno_class<- factor(m_proportions$anno_class, levels = annos2)

m_proportions %>%
  ggplot(aes(x = perc*100, y=anno_class, fill = factor(signif))) +
  geom_bar(stat="identity", width = 0.7, colour="black") +
  theme_classic(base_size=32) +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("magenta1", "gray90", "darkmagenta")) +
  ylab("Annotation") +
  xlab("Percentage")

f_proportions<- f_full %>% 
  group_by(anno_class, signif) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count))

d3<- f_full %>%
  mutate(unique_cpg = paste(chr, cpg_loc, sep="_")) %>%
  distinct(unique_cpg, .keep_all = T) %>%
  group_by(signif) %>%
  summarise(count = n()) %>%
  mutate(perc = count/sum(count))

d3$anno_class<- "All"

f_proportions<- rbind(f_proportions, d3)
annos2<- unique(f_proportions$anno_class)
f_proportions$anno_class<- factor(f_proportions$anno_class, levels = annos2)
f_proportions$sex<- "F"

f_proportions %>%
  ggplot(aes(x = perc*100, y=anno_class, fill = factor(signif))) +
  geom_bar(stat="identity", width = 0.7, colour="black") +
  #geom_text(label=d2$count, hjust=-2) +
  theme_classic(base_size=32) +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("darkolivegreen1", "gray90", "darkolivegreen")) +
  ylab("Annotation") +
  xlab("Percentage")

total_proportions<- 

######################################
###           ENRICHMENT           ###   
######################################
#Hypo Enrichment---------------------------------------------------------------------
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
      a<- nrow(df2[df2$fdr < 0.05 & df2$beta < 0,])
      b<- nrow(df3[df3$fdr < 0.05 & df3$beta < 0,])
      
      #Counts for NOT fdr < 0.05 & beta < 0
      c<- nrow(df2[!(df2$fdr < 0.05 & df2$beta < 0),])
      d<- nrow(df3[!(df3$fdr < 0.05 & df3$beta < 0),])
      
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
      a<- nrow(df2[df2$fdr < 0.05 & df2$beta > 0,])
      b<- nrow(df3[df3$fdr < 0.05 & df3$beta > 0,])
      
      #Counts for NOT fdr < 0.05 & beta < 0
      c<- nrow(df2[!(df2$fdr < 0.05 & df2$beta > 0),])
      d<- nrow(df3[!(df3$fdr < 0.05 & df3$beta > 0),])
      
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
  
  if (unique(model_df$sex) == "m") {
    ft$type<- "M"
  } else if (unique(model_df$sex) == "f"){
    ft$type<- "F"
  }
  
  return(ft)
}

#HYPOMETHYLATION
m_hypo<- enrichment(m_full, "hypo")
f_hypo<- enrichment(f_full, "hypo")

full_hypo<- rbind(f_hypo, m_hypo)

full_hypo %>%
  ggplot(aes(x=annotation, y=estimate, fill = type, alpha = padj < 0.05)) +
  geom_col(position = position_dodge(0.5)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_errorbar(ymin = full_hypo$conf.low, ymax = full_hypo$conf.high, width = 0.3, position = position_dodge(0.7)) +
  scale_fill_manual(values = c("darkolivegreen", "darkmagenta"), name = "Sex") +
  theme_classic(base_size = 32) +
  theme(legend.position = "none") +
  #ylim(c(-1, 7)) +
  ylab("Odds Ratio") +
  xlab("Annotation") +
  coord_flip()


#FEMALES
m_hyper<- enrichment(m_full,"hyper")
f_hyper<- enrichment(f_full,"hyper")

full_hyper<- rbind(f_hyper, m_hyper)

full_hyper %>%
  ggplot(aes(x=annotation, y=estimate, fill = type, alpha = padj < 0.05)) +
  geom_col(position = position_dodge(0.5), colour="black") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_errorbar(ymin = full_hyper$conf.low, ymax = full_hyper$conf.high, width = 0.3, position = position_dodge(0.7)) +
  scale_fill_manual(values = c("darkolivegreen1", "magenta1"), name = "Sex") +
  theme_classic(base_size = 32) +
  theme(legend.position = "none") +
  #ylim(c(-1, 7)) +
  ylab("Odds Ratio") +
  xlab("Annotation") +
  coord_flip()

#GSEA---------------------------------------------------------------------------
#Generate hallmark gene set
hallmark.msigdb = msigdbr(species = "Macaca mulatta", category = "H")
hallmark_list = split(x = hallmark.msigdb$ensembl_gene, f = hallmark.msigdb$gs_name)

#Generate gene ontology set
go_set = msigdbr(species = "Macaca mulatta", category = "C5", subcategory = "GO:BP")
go_set = split(x = go_set$ensembl_gene, f = go_set$gs_name)

#Generate rank-ordered vector by pqlseq coefficient
proms<- hypo_prom %>% 
  dplyr::select(anno, std_beta_f) %>% 
  arrange(desc(std_beta_f))

proms2<- proms$std_beta_f
names(proms2) = proms$anno

#Enrichment for Hallmark set
prom_gsea2<- fgsea(pathways = hallmark_list, 
                   stats = proms2,
                   minSize = 15,
                   maxSize = 500,
                   eps = 0.0,
                   scoreType = "neg")


#This version has the abs_diff code so keeping here for now 
hyper_enrichment<- function(model_df, var_opt){
  
  df_list<- list()
  
  df<- model_df
  
  if (var_opt=="class") {
    
    annotation<- unique(df$anno_class)
    
    #create list of contingency tables for sex
    for (i in annotation) {
      
      #Counts for negative estimates
      a<- nrow(df[df$abs_diff < 0 & df$anno_class == i,])
      b<- nrow(df[df$abs_diff < 0 & df$anno_class != i,])
      
      #Counts positive estimates
      c<- nrow(df[df$abs_diff > 0 & df$anno_class == i,])
      d<- nrow(df[df$abs_diff > 0 & df$anno_class != i,])
      
      #Generate contingency table - cols are class, rows are direction
      c_table<- data.frame("abs_diff < 0" = c(a, b),
                           "abs_diff > 0" = c(c, d),
                           row.names = c(paste(i, "Y"), paste(i,"N")))
      
      print(c_table)
      
      df_list[[length(df_list)+1]] = c_table
    } 
  } else if (var_opt=="annotation") {
    
    #filter model for negative sex estimates
    df<- df %>%
      filter(fdr_f < 0.05 | fdr_m < 0.05)
    
    annotation<- unique(df$anno)
    
    #create list of contingency tables for sex
    for(i in annotation){
      
      #Counts for negative estimates
      a<- nrow(df[df$abs_diff < 0 & df$anno == i,])
      b<- nrow(df[df$abs_diff < 0 & df$anno != i,])
      
      #Counts positive estimates
      c<- nrow(df[df$abs_diff > 0 & df$anno == i,])
      d<- nrow(df[df$abs_diff > 0 & df$anno != i,])
      
      #Generate contingency table - cols are class, rows are direction
      c_table<- data.frame("abs_diff < 0" = c(a, b),
                           "abs_diff > 0" = c(c, d),
                           row.names = c(paste(i, "Y"), paste(i,"N")))
      print(c_table)
      
      df_list[[length(df_list)+1]] = c_table
      
    }
  } 
  
  #name table list
  names(df_list)<- annotation
  
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
#Save workspace image-----------------------------------------------------------
save.image("nested_analysis.RData")
