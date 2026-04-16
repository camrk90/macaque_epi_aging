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
library(ggridges)
options(scipen = 999)
#setwd('/scratch/ckelsey4/Cayo_meth')
load("nested_analysis.RData")

######################################
###      Import PQLseq Model       ###
######################################
#Import nested pqlseq model-----------------------------------------------------
#Define import function
import_models<- function(model_path, file_string, mod_type){
  
  file_list<- list.files(path = model_path, pattern = file_string)
  file_order<- str_split_i(file_list, "_", 4)
  
  #Import glm models as list
  model_list<- lapply(paste(model_path, file_list, sep = "/"), readRDS)
  
  #Rename list elements
  names(model_list)<- file_order
  
  if (mod_type == "interaction") {
    
    df_nms<- c("age", "sexM", "age*sex")
    
  } else if (mod_type == "nested") {
    
    df_nms<- c("ageF", "ageM")
    
  }
  
  mod2<- lapply(model_list, function(mod){
    
    names(mod)<- df_nms
    
    mod<- lapply(names(mod), function(y){
      
      df<- mod[[y]]
      
      df<- df %>%
        dplyr::select(-c(h2, sigma2))
      
      colnames(df)<- c("outcome", "n", paste(colnames(df[3:length(df)]), y, sep = "_"))
      
      df
    })
    
    names(mod)<- df_nms
    
    mod
    
  })
  
  mod2 <- lapply(df_nms, function(nm) {
    do.call(rbind, lapply(mod2, function(x) x[[nm]]))
  })
  
  names(mod2)<- df_nms
  
  mod2<- lapply(mod2, function(df){
    
    df$outcome<- str_split_i(df$outcome, "\\.", 3)
    
    df$outcome2<- df$outcome
    
    df<- df %>%
      separate_wider_delim(outcome2, "_", names = c("chr", "chromStart", "chromEnd")) %>%
      mutate(region_range = paste(as.character(chromStart), "-", as.character(chromEnd))) %>%
      relocate(c(chr, chromStart, chromEnd, region_range), .before = n)
      
  })
  
  return(mod2)
  
}

#Import female and male pqlseq output
nested_list<- import_models(model_path = "/scratch/ckelsey4/Cayo_meth/glmer_model_compare",
                                file_string = "dnam_nested_model", mod_type = "nested")

nested_df<- left_join(nested_list[[1]], nested_list[[2]][,c(1,7:12)], by = "outcome")

nested_df <- nested_df %>%
  mutate(across(6:18, as.numeric))

nested_df<- nested_df %>% 
  mutate(diff = abs(beta_ageF) - abs(beta_ageM))

nested_df$type<- "autosomes"
nested_df$type[nested_df$chr == "X"]<- "X"

#Import additive pqlseq model---------------------------------------------------
additive_list<- import_models(model_path = "/scratch/ckelsey4/Cayo_meth/glmer_model_compare",
                                   file_string = "dnam_interaction_model", mod_type = "interaction")

additive_df<- left_join(additive_list[[1]], additive_list[[2]][,c(1,7:12)], by = "outcome")
additive_df<- left_join(additive_df, additive_list[[3]][,c(1,7:12)], by = "outcome")

additive_df <- additive_df %>%
  mutate(across(6:24, as.numeric))

additive_df$type<- "autosomes"
additive_df$type[additive_df$chr == "X"]<- "X"

#Import genes-------------------------------------------------------------------
mm_genes<- rtracklayer::import('/scratch/ckelsey4/Cayo_meth/Macaca_mulatta.Mmul_10.110.chr.gtf')
mm_genes=as.data.frame(mm_genes)
mm_genes<- mm_genes %>%
  filter(type == "gene")

#Plot nested pqlseq model-------------------------------------------------------
#Plot significant estimates
nested_df %>%
  filter(chr != "X") %>%
  dplyr::select(c(beta_ageF, beta_ageM, fdr_ageF, fdr_ageM)) %>%
  pivot_longer(cols = c(beta_ageF, beta_ageM, fdr_ageF, fdr_ageM),
               names_to = c(".value", "sex"),
               names_sep = "_") %>%
  filter(fdr < .05) %>%
  ggplot(aes(beta, fill=sex)) +
  geom_histogram(alpha = 0.7, colour="black", bins = 100, position = "identity") +
  #geom_density(alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "red") +
  scale_fill_manual(values = c("darkolivegreen", "darkmagenta"), name = "Sex") +
  ylab("Count") +
  xlab("Estimate") +
  theme_classic(base_size=24)

nested_df %>% 
  filter(fdr_ageF < .05 | fdr_ageM < .05) %>%
  filter(chr != "X") %>%
  ggplot(aes(beta_ageF, beta_ageM, colour = diff)) +
  geom_point(alpha = 0.5) +
  geom_vline(xintercept=0, linetype="dashed") +
  geom_hline(yintercept=0, linetype="dashed") +
  geom_abline() +
  scale_color_gradient2(low = "darkolivegreen", mid = "grey70", high = "darkmagenta", midpoint = 0, name = "") +
  theme_classic(base_size = 6) +
  theme(panel.background = element_rect(colour = "black", linewidth=1),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        plot.margin = margin(1, 1, 1, 1, "pt"),
        aspect.ratio = 1,
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey98", linewidth = 0.5)) +
  stat_cor()

nested_df %>% 
  filter(chr != "X") %>%
  ggplot(aes(diff, fill = after_stat(x))) +
  geom_histogram(bins = 50) +
  geom_vline(xintercept=0, linetype="dashed") +
  geom_vline(xintercept=median(nested_df$diff), linetype="dashed", colour = 'red') +
  theme_classic(base_size = 6) +
  scale_fill_gradient2(low = "darkolivegreen", mid = "grey70", high = "darkmagenta", midpoint = 0, name = "") +
  theme(legend.key.width = unit(5, 'mm'), 
        legend.key.height = unit(1, 'mm'),
        legend.position = "top") +
  theme(panel.background = element_rect(colour = "black", linewidth=1),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        plot.margin = margin(1, 1, 1, 1, "pt"),
        aspect.ratio = 1,
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey98", linewidth = 0.5)) +
  scale_x_continuous(breaks = seq(-.10, 0.10, 0.05), limits = c(-0.10, 0.10))

#Counts and Beta Distributions
count_signif_regions<- function(x) {
  
  df<- x %>%
    select(starts_with("fdr_"))
  vars<- gsub("fdr_", "", colnames(df))
    
  counts <- colSums(df < 0.05, na.rm = TRUE)
  
  counts<- data.frame(predictor = vars,
                      count = counts)
  counts$predictor<- factor(counts$predictor, levels = vars)
  counts<- counts %>%
    mutate(perc_signif = count/nrow(df))
  
  counts_plot<- counts %>%
    ggplot(aes(predictor, count, fill = predictor)) +
    geom_bar(stat = 'identity') +
    geom_text(label=counts$count, vjust=-0.25, size = 1.5) +
    theme_classic(base_size = 12) +
    theme(legend.position = "none",
          panel.background = element_rect(colour = "black", linewidth=1),
          axis.line = element_line(colour = "black", linewidth = 0.5),
          axis.title.x = element_blank(),
          plot.margin = margin(1, 1, 1, 1, "pt")) +
    xlab("Predictor") +
    ylab("Count")
  
  return(list(plot = counts_plot, df = counts))
  
}

nested_counts<- count_signif_regions(nested_df)

#Plot additive pqlseq model-----------------------------------------------------
#Plot significant estimates
additive_df %>%
  filter(chr != "X") %>%
  dplyr::select(c(beta_age, beta_sexM, `beta_age*sex`, fdr_age, fdr_sexM, `fdr_age*sex`)) %>%
  pivot_longer(cols = c(beta_age, beta_sexM, `beta_age*sex`,  fdr_age, fdr_sexM, `fdr_age*sex`),
               names_to = c(".value", "var"),
               names_sep = "_") %>%
  filter(fdr <.05) %>%
  ggplot(aes(beta, fill=var)) +
  geom_histogram(alpha = 0.7, bins = 100, position = "identity") +
  #geom_density(alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "red") +
  ylab("Count") +
  xlab("Estimate") +
  theme_classic(base_size=24)

additive_df %>%
  ggplot(aes(`beta_age*sex`)) +
  geom_histogram(bins=50)

additive_df %>%
  ggplot(aes(beta_age)) +
  geom_histogram(bins=50)

additive_df %>%
  ggplot(aes(beta_sexM)) +
  geom_histogram(bins=50)

additive_df %>% 
  #filter(chr != "X") %>%
  filter(fdr_age < .05 & fdr_sexM < .05) %>%
  ggplot(aes(beta_age, beta_sexM)) +
  geom_point(alpha = 0.5) +
  geom_vline(xintercept=0, linetype="dashed") +
  geom_hline(yintercept=0, linetype="dashed") +
  geom_abline() +
  #scale_color_gradient2(low = "darkolivegreen", mid = "grey70", high = "darkmagenta", midpoint = 0, name = "") +
  theme_classic(base_size = 6) +
  theme(panel.background = element_rect(colour = "black", linewidth=1),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        plot.margin = margin(1, 1, 1, 1, "pt"),
        aspect.ratio = 1,
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey98", linewidth = 0.5)) +
  stat_cor()

#Counts and Beta Distributions
additive_counts<- count_signif_regions(additive_df)

######################################
###        Add Annotations         ###
######################################
#Import annotation files--------------------------------------------------------
re_anno<- read_csv("/scratch/ckelsey4/Cayo_meth/re_annotations.csv")
re_anno<- re_anno %>%
  filter(chr != "Y")

re_labels<- c("Simple_repeat", "Satellite", "SINE", "LINE", "LTR", "Retroposon", "DNA")

chmm_intersect<- read_csv("/scratch/ckelsey4/Cayo_meth/chmm_annotations.csv")
chmm_intersect<- chmm_intersect %>%
  filter(chr != "Y")

chmm_labels<- c("Active TSS", "Flanking Active TSS", "Transcr. at 5' & 3'", "Strong Transcription", "Weak Transcription", "Enhancers", "Genic Enhancers",
                "ZNF Genes & Repeats", "Heterochromatin", "Bivalent/Poised TSS", "Flanking Bivalent TSS/Enh", "Bivalent Enhancer", "Repressed Polycomb",
                "Weak Repressed Polycomb", "Quiescent/Low")

promoters<- read_csv("/scratch/ckelsey4/Cayo_meth/promoters.csv")
promoters<- promoters %>%
  filter(chr != "Y")


#Promoters----------------------------------------------------------------------
nested_proms<- inner_join(nested_df, promoters, by = c("region_range", "chr"))
nested_proms<- nested_proms %>%
  drop_na() %>%
  distinct(anno, .keep_all = T)
nested_df$anno_class<- "Promoter"

#Annotations--------------------------------------------------------------------
#Define join function
anno_join<- function(output_df, annotation_df, annotation_type){
  
  if (annotation_type == "chmm"){
    
    #Join pqlseq model and chmm 
    model_df<- left_join(annotation_df, output_df, by = c("region_range", "chr"))
    
    #Filter out regions with models that didn't converge resulting in NAs in the annotation join
    model_df<- model_df %>%
      drop_na() %>%
      dplyr::select(-c(chromStart, chromEnd))
    
    #Set annotations as factor and reorder
    annotations_ordered<- str_sort(unique(model_df$anno), numeric = TRUE)
    model_df$anno<- factor(model_df$anno, levels = annotations_ordered, labels = chmm_labels)
    
    #Create column of broad categories for annotations
    model_df$anno_class<- "Transcription Start Sites"
    model_df$anno_class[model_df$anno %in% chmm_labels[3:5]]<- "Active Transcription"
    model_df$anno_class[model_df$anno %in% chmm_labels[6:8]]<- "Enhancer Regions"
    model_df$anno_class[model_df$anno %in% chmm_labels[9:15]]<- "Quiescent States"
    
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
      dplyr::select(-c(range, chromStart, chromEnd))
    
    #Remove non-sensical annotations (NA, Unknown etc)
   model_df<-model_df[!model_df$repClass == "Unknown",]
   model_df<-model_df[!model_df$repClass == "DNA?",]
   model_df<-model_df[!model_df$repClass == "LTR?",]
   model_df<-model_df[!model_df$repClass == "RC?",]
   model_df<-model_df[!model_df$repClass == "Unspecified",]
   
   #Set annotations as factor and reorder
   model_df$repClass<- as.factor(model_df$repClass)
   
   #Create column of broad categories for annotations
   model_df<- model_df[!model_df$repClass %in% c("RC", "rRNA", "snRNA", "tRNA", "srpRNA", "scRNA", "Low_complexity"),]
   model_df$anno_class<- "Simple Repeats"
   model_df$anno_class[model_df$repClass %in% c("SINE", "LINE", "LTR", "Retroposon")]<- "Transposable Elements Class I"
   model_df$anno_class[model_df$repClass %in% "DNA"]<- "Transposable Elements Class II"
   
   #Set classes as factors
   class_factors<- c("Simple Repeats", "Transposable Elements Class I", "Transposable Elements Class II")
   model_df$anno_class<- factor(model_df$anno_class, levels = class_factors)
   
   #Arrange cols to match chmm df
   model_df<- model_df %>%
     dplyr::select(-repName) %>%
     dplyr::rename(anno = repClass) %>%
     relocate(c(anno_start, anno_end), .before = anno)
   model_df$anno<- factor(model_df$anno, levels = re_labels)
   
   return(model_df)
  }
}

nested_chmm<- anno_join(nested_df, chmm_intersect, "chmm")
nested_re<- anno_join(nested_df, re_anno, "re")
nested_full<- rbind(nested_chmm, nested_re)

annos<- c("Promoter", "Transcription Start Sites", "Active Transcription", "Enhancer Regions", "Quiescent States", "Simple Repeats",
          "Transposable Elements Class I", "Transposable Elements Class II", "Structural RNAs")

nested_full$anno_source<- "Repeat Elements"
nested_full$anno_source[nested_full$anno_class == "Transcription Start Sites" | nested_full$anno_class == "Active Transcription" |
                          nested_full$anno_class == "Enhancer Regions" | nested_full$anno_class == "Quiescent States" | 
                          nested_full$anno_class == "Promoter"]<- "Transcription"

nested_full<- nested_full %>%
  arrange(anno_source, anno_class)

#Rearrange factors to sort by type then log_or
nested_full$anno_class<- factor(nested_full$anno_class, levels = rev(annos))

nested_full$f_signif<- "Non-Significant"
nested_full$f_signif[nested_full$fdr_ageF < 0.05 & nested_full$beta_ageF < 0]<- "Age-Hypomethylated"
nested_full$f_signif[nested_full$fdr_ageF < 0.05 & nested_full$beta_ageF > 0]<- "Age-Hypermethylated"

nested_full$f_signif<- factor(nested_full$f_signif, levels=c("Age-Hypermethylated", "Non-Significant", "Age-Hypomethylated"))

nested_full$m_signif<- "Non-Significant"
nested_full$m_signif[nested_full$fdr_ageM < 0.05 & nested_full$beta_ageM < 0]<- "Age-Hypomethylated"
nested_full$m_signif[nested_full$fdr_ageM < 0.05 & nested_full$beta_ageM > 0]<- "Age-Hypermethylated"

nested_full$m_signif<- factor(nested_full$m_signif, levels=c("Age-Hypermethylated", "Non-Significant", "Age-Hypomethylated"))

nested_full<- nested_full %>%
  mutate(unique_cpg = paste(chr, cpg_loc, sep="_")) %>%
  relocate(unique_cpg, .after = outcome)

#Plot basics--------------------------------------------------------------------
nested_full %>%
  dplyr::select(c(beta_ageF, beta_ageM, fdr_ageF, fdr_ageM, anno, anno_source)) %>%
  pivot_longer(cols = c(beta_ageF, beta_ageM, fdr_ageF, fdr_ageM),
               names_to = c(".value", "sex"),
               names_sep = "_") %>%
  filter(anno_source == "Repeat Elements") %>%
  filter(fdr < .05) %>%
  ggplot(aes(beta, anno, fill = sex)) +
  geom_density_ridges(alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = 'red')

nested_full %>%
  dplyr::select(c(beta_ageF, beta_ageM, fdr_ageF, fdr_ageM, anno, anno_source)) %>%
  pivot_longer(cols = c(beta_ageF, beta_ageM, fdr_ageF, fdr_ageM),
               names_to = c(".value", "sex"),
               names_sep = "_") %>%
  filter(anno_source == "Transcription") %>%
  filter(fdr < .05) %>%
  ggplot(aes(beta, anno, fill = sex)) +
  geom_density_ridges(alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = 'red') +
  theme_classic(base_size = 12) +
  theme(panel.background = element_rect(colour = "black", linewidth=1),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        plot.margin = margin(1, 1, 1, 1, "pt"),
        aspect.ratio = 1,
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey98", linewidth = 0.5))

generate_proportion<- function(df, x){
  
  d1<- df %>% 
    distinct(unique_cpg, .keep_all = T) %>%
    group_by(anno, {{x}}) %>% 
    summarise(count = n()) %>% 
    mutate(perc = count/sum(count))
  
  d2<- df %>%
    distinct(unique_cpg, .keep_all = T) %>%
    group_by({{x}}) %>%
    summarise(count = n()) %>%
    mutate(perc = count/sum(count))
  
  d2$anno<- "All"
  
  d3<- rbind(d1, d2)
  annos2<- unique(d3$anno)
  d3$anno<- factor(d3$anno, levels = annos2)
  
  d3 %>%
    filter(is.na(anno)) %>%
    dplyr::select(count)
  
  col1 <- eval(substitute(x), d2)
  
  d3$percent<- d3$perc*100
  
  pro_plot<- d3 %>%
    arrange(anno) %>%
    ggplot(aes(x = percent, y=anno, fill = factor({{x}}))) +
    geom_bar(stat="identity", width = 0.7, colour="black") +
    theme_classic(base_size=24) +
    geom_vline(xintercept = (1-d2$perc[col1 == "Age-Hypermethylated"])*100, linetype = 'dashed') +
    geom_vline(xintercept = d2$perc[col1 == "Age-Hypomethylated"]*100, linetype = 'dashed') +
    theme(legend.position = "top") +
    #scale_fill_manual(values = c(c1, c2, c3), name = "") +
    ylab("Annotation") +
    xlab("Percentage") 
  
  return(list(data = d3, plot = pro_plot))
  
}

m_proportions<- generate_proportion(nested_full, m_signif)

f_proportions<- generate_proportion(nested_full, f_signif)

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
