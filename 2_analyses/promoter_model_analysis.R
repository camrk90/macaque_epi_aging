library(tidyverse)
library(GenomicFeatures)
library(GenomicRanges)
library(biomaRt)
library(bsseq)
library(fgsea)
library(broom)
library(ggpubr)
library(ggcorrplot)
library(ggborderline)
library(ggrepel)
library(ggVennDiagram)
library(msigdbr)
library(effectsize)
options(scipen = 999)
setwd('/scratch/ckelsey4/Cayo_meth')
load("promoter_analysis.RData")

######################################
###      Import PQLseq Models      ###
######################################
#Import metadata----------------------------------------------------------------
long_data<- read.table("/scratch/ckelsey4/Cayo_meth/long_data_adjusted.txt")

######################################
###      Import PQLseq Model       ###
######################################
#Import nested pqlseq model-----------------------------------------------------
#Define import function
import_models<- function(model_path, file_string, mod_type){
  
  file_list<- list.files(path = model_path, pattern = file_string)
  file_order<- str_split_i(file_list, "_", 5)
  
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
    
    df<- df %>%
      separate_wider_delim(outcome, ".", names = c("chr", "ensembl_gene_name"))
    

  })
  
  return(mod2)
  
}

#Import female and male pqlseq output
nested_list<- import_models(model_path = "/scratch/ckelsey4/Cayo_meth/glmer_model_compare",
                            file_string = "dnam_prom_nested_model", mod_type = "nested")

nested_df<- left_join(nested_list[[1]], nested_list[[2]][,c(2,4:9)], by = "ensembl_gene_name")

nested_df <- nested_df %>%
  mutate(across(3:ncol(nested_df), as.numeric))

nested_df<- nested_df %>% 
  mutate(diff = abs(beta_ageF) - abs(beta_ageM))

nested_df$type<- "autosomes"
nested_df$type[nested_df$chr == "X"]<- "X"

#Load promoters-----------------------------------------------------------------
prom_cov<- readRDS("/scratch/ckelsey4/Cayo_meth/prom_cov_filtered")
prom_cov<- prom_cov[c(1:21)]
prom_m<- readRDS("/scratch/ckelsey4/Cayo_meth/prom_m_filtered")
prom_m<- prom_m[c(1:21)]

#Filter metadata to lids in regions list
long_data<- long_data[long_data$lid_pid %in% colnames(prom_cov[[1]]),]

prom_cov<- lapply(names(prom_cov), function(x){
  prom_cov<- subset(prom_cov[[x]], select=long_data$lid_pid)
  return(prom_cov)
})

prom_m<- lapply(names(prom_m), function(x){
  prom_m<- subset(prom_m[[x]], select=long_data$lid_pid)
  return(prom_m)
})

prom_cov<- do.call(rbind, prom_cov)
prom_m<- do.call(rbind, prom_m)

perc_meth<- prom_m/prom_cov

#Import genes-------------------------------------------------------------------
mm_genes<- rtracklayer::import('/scratch/ckelsey4/Cayo_meth/Macaca_mulatta.Mmul_10.110.chr.gtf')
mm_genes=as.data.frame(mm_genes)
mm_genes<- mm_genes %>%
  filter(type == "gene") %>%
  dplyr::select(c(gene_id, gene_name)) %>%
  dplyr::rename(ensembl_gene_name = gene_id) %>%
  arrange(ensembl_gene_name)

nested_df<- left_join(nested_df, mm_genes, by = "ensembl_gene_name")
nested_df<- nested_df %>%
  relocate(gene_name, .after = ensembl_gene_name)
no_x<- nested_df %>%
  filter(chr != "X")

#Plot top 10 and bottom 10 different hypo/hypermethylated proms
prom_paired %>%
  filter(fdr_f < 0.05 | fdr_m < 0.05) %>%
  filter(gene_name != "Metazoa_SRP") %>%
  filter(std_beta_f < 0 & std_beta_m < 0) %>%
  arrange(desc(abs_diff)) %>%
  dplyr::slice(1:10, (n()-10):n()) %>%
  ggplot(aes(std_diff, reorder(gene_name, std_diff), fill = std_diff)) +
  #geom_point(aes(size=3)) +
  geom_col(colour='black')+
  geom_vline(xintercept = 0, linetype = "dashed", colour = "red") +
  #scale_fill_discrete(values=c("darkmagenta", "darkolivegreen")) +
  scale_fill_gradient2(low = "darkmagenta", mid = "white", high = "darkolivegreen", midpoint = 0, name = "") +
  theme_classic(base_size = 24)

prom_paired %>%
  filter(fdr_f < 0.05 | fdr_m < 0.05) %>%
  filter(gene_name != "Metazoa_SRP") %>%
  filter(beta_f > 0 & beta_m > 0) %>%
  arrange(desc(abs_diff)) %>%
  dplyr::slice(1:10, (n()-10):n()) %>%
  ggplot(aes(abs_diff, reorder(gene_name, abs_diff))) +
  geom_point() +
  theme_classic(base_size = 24)

sig_proms<- prom_paired %>%
  filter(fdr_f < 0.05 | fdr_m < 0.05) %>%
  filter(gene_name != "Metazoa_SRP")

sig_proms %>%
  ggplot(aes(beta_m, beta_f, colour = abs_diff, 
             label=ifelse(abs_diff > 0.1 | abs_diff< -0.1, as.character(gene_name),''))) +
  geom_point() +
  geom_abline() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_gradient2(low = "darkmagenta", mid = "white", high = "darkolivegreen", midpoint = 0, name = "") +
  geom_smooth(method = "lm") +
  xlab("Beta (Male)") +
  ylab("Beta (Female)") +
  geom_text_repel(colour="black", max.overlaps = Inf, min.segment.length = 0, size = 8, force = 2) +
  theme_classic(base_size = 16) +
  theme(legend.key.height= unit(2, 'cm'))

cor.test(sig_proms$beta_f, sig_proms$beta_m)

#Plot promoter estimates--------------------------------------------------------
nested_df %>%
  filter(chr != "X") %>%
  dplyr::select(c(beta_ageF, beta_ageM, fdr_ageF, fdr_ageM)) %>%
  pivot_longer(cols = c(beta_ageF, beta_ageM, fdr_ageF, fdr_ageM),
               names_to = c(".value", "sex"),
               names_sep = "_") %>%
  filter(fdr < .20) %>%
  ggplot(aes(beta, fill=sex)) +
  geom_histogram(alpha = 0.7, colour="black", bins = 100, position = "identity") +
  #geom_density(alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "red") +
  scale_fill_manual(values = c("darkolivegreen", "darkmagenta"), name = "Sex") +
  ylab("Count") +
  xlab("Estimate") +
  theme_classic(base_size=24)

nested_df %>% 
  filter(fdr_ageF < .20 | fdr_ageM < .20) %>%
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

#Counts and Beta Distributions
f.count<- nrow(nested_df[nested_df$fdr_ageF < 0.05,])
m.count<- nrow(nested_df[nested_df$fdr_ageM < 0.05,])
counts<- data.frame(predictor = c('F', 'M'),
                    count = c(f.count, m.count))
counts$predictor<- factor(counts$predictor, levels = c('F', 'M'))
counts<- counts %>%
  mutate(perc_signif = count/nrow(nested_df))

counts %>%
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

################################################################################
#Categorical enrichment
################################################################################
#Generate hallmark gene set
hallmark.msigdb = msigdbr(species = "Macaca mulatta", category = "H")
hallmark.msigdb = split(x = hallmark.msigdb$ensembl_gene, f = hallmark.msigdb$gs_name)

paths<- hallmark.msigdb %>%
  dplyr::select(c(gs_name, gene_symbol))
paths<- paths %>%
  dplyr::rename(gene_name = gene_symbol)
test<- left_join(prom_pqlseq2, paths, by = "gene_name")
test<- test %>%
  drop_na()

#Generate gene ontology set
go_set = msigdbr(species = "Macaca mulatta", category = "C5", subcategory = "GO:MF")
go_set = split(x = go_set$ensembl_gene, f = go_set$gs_name)

#Generate rank-ordered vector by pqlseq coefficient
paired_genes<- prom_paired %>% 
  dplyr::select(outcome, std_diff) %>% 
  arrange(desc(std_diff))

paired_genes2<- paired_genes$std_diff
names(paired_genes2) = paired_genes$outcome

f_genes<- prom_paired %>% 
  dplyr::select(outcome, std_beta_f) %>% 
  arrange(desc(std_beta_f))

f_genes2<- f_genes$std_beta_f
names(f_genes2) = f_genes$outcome

m_genes<- prom_paired %>% 
  dplyr::select(outcome, std_beta_m) %>% 
  arrange(desc(std_beta_m))

m_genes2<- m_genes$std_beta_m
names(m_genes2) = m_genes$outcome

#Enrichment for Hallmark set
m_gsea<- fgsea(pathways = hallmark.msigdb, 
                    stats    = m_genes2,
                    minSize  = 15,
                    maxSize  = 500,
                    eps = 0.0)

test <- gseGO(geneList     = paired_genes2,
              OrgDb        = org.Mmu.eg.db,
              ont          = "MF",
              minGSSize    = 5,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE)

m_gsea$sex<- "M"
f_gsea$sex<- "F"
sex_gsea<- rbind(m_gsea, f_gsea)

#Plot hallmark enrichment
ggplot(sex_gsea, aes(x = reorder(pathway, NES), y = NES, fill = sex)) +
  geom_col(colour="black", position = "dodge2", alpha = sex_gsea$padj > 0.10) +
  scale_fill_manual(values=c("darkolivegreen", "darkmagenta")) +
  theme_classic(base_size=16) +
  coord_flip()

####By hallmark pathway
hallmark_inflam<- hallmark.msigdb %>%
  filter(gs_name == "HALLMARK_OXIDATIVE_PHOSPHORYLATION")
hallmark_inflam<- hallmark_inflam %>%
  dplyr::select(c(gene_symbol, ensembl_gene))
prom_paired<- prom_paired %>%
  dplyr::rename(ensembl_gene = outcome)

df<- left_join(hallmark_inflam, prom_paired, by = "ensembl_gene")
df<- df %>% drop_na()

df %>%
  ggplot(aes(beta_m, beta_f, colour = abs_diff, label=gene_name)) +
  geom_point() +
  geom_abline() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  #scale_color_gradient2(low = "darkmagenta", mid = "white", high = "darkolivegreen", midpoint = 0, name = "") +
  geom_smooth(method = "lm") +
  geom_text_repel(alpha=df$fdr_m<0.05, box.padding = 0.5) +
  xlab("Estimate (Male)") +
  ylab("Estimate (Female)") +
  theme_classic(base_size = 30) +
  theme(legend.key.height= unit(2, 'cm'))

#Save workspace image-----------------------------------------------------------
save.image("promoter_analysis.RData")

