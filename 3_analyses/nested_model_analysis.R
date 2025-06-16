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
  ggplot(aes(beta, fill=sex)) +
  geom_histogram(alpha = 0.7, colour="black", bins = 100, position = "identity") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "red") +
  scale_fill_manual(values = c("darkolivegreen", "darkmagenta"), name = "Sex") +
  ylab("Count") +
  xlab("Estimate") +
  theme_classic(base_size = 30)

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

#Plot significant estimates
nested_pqlseq %>%
  filter(beta > 0 & fdr < 0.05) %>%
  ggplot(aes(sex, beta)) +
  #geom_jitter(aes(colour = sex), alpha = 0.7)+
  geom_boxplot(aes(fill = sex), width=0.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("royalblue2", "orangered1"), name = "Sex") +
  ylab("Estimate (FDR < 0.05)") +
  xlab("Sex") +
  theme_classic(base_size = 30)

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
  dplyr::select(c(outcome, region_range, chr, std_beta, beta, se_beta, fdr, sex)) %>%
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

nrow(hyper_t[hyper_t$fdr_m < .05,])
nrow(hyper_t[hyper_t$fdr_f < .05,])

#Plot t-test distributions
df<- hyper_t %>%
  pivot_longer(cols = c(beta_f, beta_m), names_to = "sex", values_to = "beta")

df %>%
  ggplot(aes(beta, fill=sex)) +
  geom_density(alpha = 0.5, colour="black", position = "identity") + 
  scale_fill_manual(values = c("royalblue2", "orangered1"), name = "Sex") +
  ylab("Count") +
  xlab("Estimate (fdr_f < 0.05 | fdr_m < 0.05)") +
  theme_classic(base_size = 30)

df %>%
  ggplot(aes(sex, beta, fill=sex)) +
  geom_boxplot() +
  scale_fill_manual(values = c("royalblue2", "orangered1"), name = "Sex") +
  ylab("Count") +
  xlab("Estimate (fdr_f < 0.05 | fdr_m < 0.05)") +
  theme_classic(base_size = 30)

df2<- hypo_t %>%
  pivot_longer(cols = c(std_beta_f, std_beta_m), names_to = "sex", values_to = "beta")

df2 %>%
  ggplot(aes(beta, fill=sex)) +
  geom_histogram(alpha = 0.5, colour="black", bins = 100, position = "identity") + 
  scale_fill_manual(values = c("royalblue2", "orangered1"), name = "Sex") +
  ylab("Count") +
  xlab("Estimate (fdr_f < 0.05 | fdr_m < 0.05)") +
  theme_classic(base_size = 30)

sex_hyper %>%
  filter(fdr_f < 0.05 | fdr_m < 0.05) %>%
  filter(!chr == "X") %>%
  ggplot(aes(abs_diff, fill = abs_diff > 0)) +
  geom_histogram(bins = 100, colour = "black", position = 'identity', alpha = 0.8) +
  scale_fill_manual(values = c("royalblue2", "orangered1"), name = "Sex") +
  xlab("|Estimate (F)| - |Estimate (M)|") +
  ylab("Count") +
  theme_classic(base_size = 30) +
  theme(legend.position = "none")

#Distribution of difference in abs estimates for hypomethylated regions with x-chrom
sex_hypo %>%
  filter(fdr_f < 0.05 | fdr_m < 0.05) %>%
  ggplot(aes(abs_diff, fill = after_stat(x))) +
  geom_histogram(bins = 100, colour = "black", position = 'identity') +
  geom_vline(xintercept = 0, linetype = "dashed", colour = 'red') +
  scale_fill_gradient2(low = "darkmagenta", mid = "white", high = "darkolivegreen", midpoint = 0, name = "Beta Difference") +
  xlab("|Estimate (F)| - |Estimate (M)|") +
  ylab("Hypomethylated Regions") +
  theme_classic(base_size = 30) +
  theme(legend.position = "none")

#Distribution of difference in abs estimates for hypomethylated regions without x-chrom
sex_hypo %>%
  filter(fdr_f < 0.05 | fdr_m < 0.05) %>%
  filter(chr != "X") %>%
  ggplot(aes(abs_diff, fill = after_stat(x))) +
  geom_histogram(bins = 100, colour = "black", position = 'identity') +
  geom_vline(xintercept = 0, linetype = "dashed", colour = 'red') +
  scale_fill_gradient2(low = "orangered2", mid = "white", high = "royalblue2", midpoint = 0, name = "Beta Difference") +
  xlab("|Estimate (F)| - |Estimate (M)|") +
  ylab("Hypomethylated Regions No X") +
  theme_classic(base_size = 30) +
  theme(legend.position = "none")

#Distribution of difference in abs estimates with X-chrom
sex_hyper %>%
  filter(fdr_f < 0.05 | fdr_m < 0.05) %>%
  #filter(abs_diff > -0.2 & abs_diff < 0.2) %>%
  ggplot(aes(abs_diff, fill = after_stat(x))) +
  geom_histogram(bins = 100, colour = "black", position = 'identity') +
  geom_vline(xintercept = 0, linetype = "dashed", colour = 'red') +
  scale_fill_gradient2(low = "darkmagenta", mid = "white", high = "darkolivegreen", midpoint = 0, name = "Beta Difference") +
  xlab("|Estimate (F)| - |Estimate (M)|") +
  ylab("Hypermethylated Regions") +
  theme_classic(base_size = 30) +
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
  #filter(fdr_f < 0.05 | fdr_m < 0.05) %>%
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
  theme_classic(base_size = 30) +
  theme(legend.key.height= unit(2, 'cm'))

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
hypo_prom<- left_join(promoters, sex_hypo, by = c("region_range", "chr"))
hypo_prom<- hypo_prom %>%
  drop_na() %>%
  distinct(anno, .keep_all = T)
hypo_prom$anno_class<- "Promoter"
hyper_prom<- left_join(promoters, sex_hyper, by = c("region_range", "chr"))
hyper_prom<- hyper_prom %>%
  drop_na() %>%
  distinct(anno, .keep_all = T)
hyper_prom$anno_class<- "Promoter"

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

hyper_chmm<- anno_join(sex_hyper, chmm_intersect, "chmm")
hyper_re<- anno_join(sex_hyper, re_anno, "re")
hyper_full<- rbind(hyper_chmm, hyper_re, hyper_prom)

hypo_chmm<- anno_join(sex_hypo, chmm_intersect, "chmm")
hypo_re<- anno_join(sex_hypo, re_anno, "re")
hypo_full<- rbind(hypo_chmm, hypo_re, hypo_prom)

#Plot basics--------------------------------------------------------------------
hyper_chmm %>%
  filter(fdr_f < 0.05 | fdr_m < 0.05) %>%
  ggplot(aes(beta_m, beta_f, colour = abs_diff)) +
  geom_point() +
  #geom_abline() +
  geom_smooth(method = "glm") +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  scale_color_gradient2(low = "orangered2", high = "royalblue2", midpoint = 0, name = "Beta Difference") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  xlab("Estimate M") +
  ylab("Estimate F") +
  facet_wrap(vars(anno), nrow=3)

hypo_full %>%
  filter(fdr_f < 0.05 | fdr_m < 0.05) %>%
  ggplot(aes(beta_m, beta_f, colour = abs_diff)) +
  geom_point() +
  geom_abline() +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  scale_color_gradient2(low = "orangered2", high = "royalblue2", midpoint = 0, name = "Beta Difference") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  xlab("Estimate M") +
  ylab("Estimate F")

######################################
###           ENRICHMENT           ###   
######################################
#CHMM Enrichment----------------------------------------------------------------
enrichment<- function(model_df, var_opt){
  
  df_list<- list()
  
  df<- model_df
  
  if (var_opt=="class") {
    
    #filter model for negative sex estimates
    df<- df %>%
      filter(fdr_f < 0.05 | fdr_m < 0.05)
    
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

hyper_class<- enrichment(hyper_full, "class")

annos<- c("Promoter", "Transcription Start Sites", "Active Transcription", "Enhancer Regions", "Quiescent States", "Simple Repeats",
          "Transposable Elements Class I", "Transposable Elements Class II", "Structural RNAs")

hyper_class$anno_source<- "Repeat Elements"
hyper_class$anno_source[hyper_class$annotation == "Transcription Start Sites" | hyper_class$annotation == "Active Transcription" |
                         hyper_class$annotation == "Enhancer Regions" | hyper_class$annotation == "Quiescent States" | 
                          hyper_class$annotation == "Promoter"]<- "Transcription"
hyper_class<- hyper_class %>%
  arrange(anno_source, log_or)
hyper_levels<- as.character(hyper_class$annotation)

#Rearrange factors to sort by type then log_or
hyper_class$annotation<- factor(hyper_class$annotation, levels = rev(annos))

hyper_class %>%
  ggplot(aes(x=annotation, y=log_or, fill = log_or > 0, alpha = padj < 0.05)) +
  geom_col(position = position_dodge(0.7)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(ymin = hyper_class$log_ci.lo, ymax = hyper_class$log_ci.hi, width = 0.3, position = position_dodge(0.7)) +
  scale_fill_manual(values = c("darkolivegreen", "darkmagenta"), name = "Sex") +
  theme_classic(base_size = 32) +
  theme(legend.position = "none") +
  ylim(c(-1, 7)) +
  ylab("Log Odds") +
  xlab("Annotation") +
  coord_flip()

#Hypomethylated region enrichment
hypo_class<- enrichment(hypo_full, "class")

hypo_class$anno_source<- "Repeat Elements"
hypo_class$anno_source[hypo_class$annotation == "Transcription Start Sites" | hypo_class$annotation == "Active Transcription" |
                       hypo_class$annotation == "Enhancer Regions" | hypo_class$annotation == "Quiescent States" | 
                         hypo_class$annotation == "Promoter"]<- "Transcription"
hypo_class<- hypo_class %>%
  arrange(anno_source)
hypo_levels<- as.character(hypo_class$annotation)

#Rearrange factors to sort by type then log_or
hypo_class$annotation<- factor(hypo_class$annotation, levels = rev(annos))

hypo_class %>%
  ggplot(aes(x=annotation, y=log_or, fill = log_or > 0)) +
  geom_col(position = position_dodge(0.7)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(ymin = hypo_class$log_ci.lo, ymax = hypo_class$log_ci.hi, width = 0.3, position = position_dodge(0.7)) +
  scale_fill_manual(values = c("darkolivegreen", "darkmagenta"), name = "Sex") +
  theme_classic(base_size = 32) +
  theme(legend.position = "none") +
  ylim(c(-2, 2)) +
  ylab("Log Odds") +
  xlab("Annotation") +
  coord_flip() +
  ggtitle("Hypomethylated CpGs")

#Hypo Annotation
test<- hypo_full %>%
  filter(!anno_class == "Promoter")

hypo_enrich<- enrichment(test, "annotation")
hypo_re_enrich<- hypo_re_enrich[16:nrow(hypo_re_enrich),]

hypo_re_enrich$anno_class<- "Simple Repeats"
hypo_re_enrich$anno_class[hypo_re_enrich$annotation %in% c("SINE", "LINE", "LTR", "Retroposon")]<- "Transposable Elements Class I"
hypo_re_enrich$anno_class[hypo_re_enrich$annotation %in% "DNA"]<- "Transposable Elements Class II"
hypo_re_enrich$anno_class[hypo_re_enrich$annotation %in% c("rRNA", "snRNA", "tRNA", "srpRNA", "scRNA")]<- "Structural RNAs"

hypo_re_enrich<- hypo_re_enrich %>%
  mutate(anno_class = as.factor(anno_class))
hypo_re_enrich$anno_class<- factor(hypo_re_enrich$anno_class, levels = c("Simple Repeats", "Transposable Elements Class I", 
                                                                         "Transposable Elements Class II", "Structural RNAs"))

hypo_re_enrich<- hypo_re_enrich %>%
  arrange(anno_class, log_or)
hypo_re_enrich$annotation<- factor(hypo_re_enrich$annotation, levels = rev(hypo_re_enrich$annotation))

hypo_re_enrich %>%
  ggplot(aes(x=annotation, y=log_or, fill = log_or > 0, alpha = padj < 0.05)) +
  geom_col(position = position_dodge(0.7)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(ymin = hypo_re_enrich$log_ci.lo, ymax = hypo_re_enrich$log_ci.hi, width = 0.3, position = position_dodge(0.7)) +
  scale_fill_manual(values = c("darkolivegreen", "darkmagenta"), name = "Sex") +
  theme_classic(base_size = 32) +
  theme(legend.position = "none") +
  ylim(c(-3, 3)) +
  ylab("Log Odds") +
  xlab("Annotation") +
  coord_flip() +
  ggtitle("Hypomethylated CpGs")

#Hyper re enrich
test2<- hyper_full %>%
  filter(!anno_class == "Promoter")

hyper_enrich<- enrichment(test2, "annotation")
hyper_chmm_enrich<- hyper_enrich[1:14,]
hyper_chmm_enrich$annotation<- factor(hyper_chmm_enrich$annotation, levels = rev(annotation_labels))
hyper_chmm_enrich<- hyper_chmm_enrich %>%
  filter(!estimate == Inf) %>%
  arrange(annotation)

hyper_chmm_enrich %>%
  ggplot(aes(x=annotation, y=log_or, fill = log_or > 0, alpha = padj < 0.05)) +
  geom_col(position = position_dodge(0.7)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(ymin = hyper_chmm_enrich$log_ci.lo, ymax = hyper_chmm_enrich$log_ci.hi, width = 0.3, position = position_dodge(0.7)) +
  scale_fill_manual(values = c("darkolivegreen", "darkmagenta"), name = "Sex") +
  theme_classic(base_size = 32) +
  theme(legend.position = "none") +
  ylim(c(-3, 3)) +
  ylab("Log Odds") +
  xlab("Annotation") +
  coord_flip()

hyper_re_enrich<- hyper_enrich[15:nrow(hyper_enrich),]

hyper_re_enrich$anno_class<- "Simple Repeats"
hyper_re_enrich$anno_class[hyper_re_enrich$annotation %in% c("SINE", "LINE", "LTR", "Retroposon")]<- "Transposable Elements Class I"
hyper_re_enrich$anno_class[hyper_re_enrich$annotation %in% "DNA"]<- "Transposable Elements Class II"
hyper_re_enrich$anno_class[hyper_re_enrich$annotation %in% c("rRNA", "snRNA", "tRNA", "srpRNA", "scRNA")]<- "Structural RNAs"

hyper_re_enrich<- hyper_re_enrich %>%
  mutate(anno_class = as.factor(anno_class))
hyper_re_enrich$anno_class<- factor(hyper_re_enrich$anno_class, levels = c("Simple Repeats", "Transposable Elements Class I", 
                                                                         "Transposable Elements Class II", "Structural RNAs"))

hyper_re_enrich<- hyper_re_enrich %>%
  arrange(anno_class, log_or)
hyper_re_enrich$annotation<- factor(hyper_re_enrich$annotation, levels = rev(hyper_re_enrich$annotation))

hyper_re_enrich<- hyper_re_enrich %>%
  filter(!estimate == Inf)

hyper_re_enrich %>%
  ggplot(aes(x=annotation, y=log_or, fill = log_or > 0, alpha = padj < 0.05)) +
  geom_col(position = position_dodge(0.7)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(ymin = hyper_re_enrich$log_ci.lo, ymax = hyper_re_enrich$log_ci.hi, width = 0.3, position = position_dodge(0.7)) +
  scale_fill_manual(values = c("darkolivegreen", "darkmagenta"), name = "Sex") +
  theme_classic(base_size = 32) +
  theme(legend.position = "none") +
  ylim(c(-3, 3)) +
  ylab("Log Odds") +
  xlab("Annotation") +
  coord_flip()

#Save workspace image-----------------------------------------------------------
save.image("nested_analysis.RData")
