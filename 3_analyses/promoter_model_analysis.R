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
library(clusterProfiler)
library(org.Mmu.eg.db)
options(scipen = 999)
setwd('/scratch/ckelsey4/Cayo_meth')
load("promoter_analysis.RData")

######################################
###      Import PQLseq Models      ###
######################################
#Import metadata----------------------------------------------------------------
long_data<- readRDS("/scratch/ckelsey4/Cayo_meth/long_data_adjusted")
long_data<- long_data %>%
  arrange(lid_pid) %>%
  filter(n > 1)

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
  dplyr::rename(outcome = gene_id) %>%
  arrange(outcome)

#Import promoter pqlseq files---------------------------------------------------
setwd('/scratch/ckelsey4/Cayo_meth/glmer_model_compare')
import_prom_pqlseq<- function(x){
  
  #Generate list of file names
  file_list<- list.files(pattern = x)
  file_order<- str_split_i(file_list, "_", 4)
  
  #Import glm models as list
  model_list<- lapply(file_list, readRDS)
  
  #Rename list elements
  names(model_list)<- file_order
  
  #Bind model list to df and add rownames
  model<- do.call(rbind, model_list)
  model$outcome2<- model$outcome
  model$outcome<- str_split_i(model$outcome, "\\.", 2)
  rownames(model)<- model$region
  
  #Separate region coordinates into start and end, delete the chr col, and move region col to front
  #model$outcome2<- model$outcome
  model<- model %>%
    separate_wider_delim(outcome2, ".", names = c("chrom", "gene")) %>%
    dplyr::select(-c(gene))
  model<- model %>%
    relocate(c(chrom), .before = n)
  #model<- model %>%
  #mutate(chromStart = as.numeric(chromStart))
  
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
prom_agew_files<- 'prom_pqlseq2_agew_'
prom_agew_pqlseq<- import_prom_pqlseq(prom_agew_files)

#Mean Age
prom_agem_files<- 'prom_pqlseq2_agem_'
prom_agem_pqlseq<- import_prom_pqlseq(prom_agem_files)

#Sex
prom_sex_files<- 'prom_pqlseq2_sex_'
prom_sex_pqlseq<- import_prom_pqlseq(prom_sex_files)

#Rename cols for each df to indicate age/sex
colnames(prom_agew_pqlseq)<- c(names(prom_agew_pqlseq[,1:3]), paste(names(prom_agew_pqlseq[,4:11]), "agew", sep = "_"))
colnames(prom_agem_pqlseq)<- c(names(prom_agem_pqlseq[,1:3]), paste(names(prom_agem_pqlseq[,4:11]), "agem", sep = "_"))
colnames(prom_sex_pqlseq)<- c(names(prom_sex_pqlseq[,1:3]), paste(names(prom_sex_pqlseq[,4:11]), "sex", sep = "_"))

#Cbind cols for age and sex dfs
prom_pqlseq<- cbind(prom_agew_pqlseq, prom_sex_pqlseq[,4:11], prom_agem_pqlseq[,4:11])

#Sort chromosome factors
sorted_labels<- str_sort(unique(prom_pqlseq$chrom), numeric=T)
prom_pqlseq<- prom_pqlseq %>%
  mutate(chrom = factor(chrom, levels = sorted_labels))

rm(prom_agew_pqlseq);rm(prom_agem_pqlseq);rm(prom_sex_pqlseq)

#Add gene names
prom_pqlseq<- inner_join(mm_genes, prom_pqlseq)

#Import nested promoter model output--------------------------------------------
f_files<- 'prom_nested_f_'
f_nested<- import_prom_pqlseq(f_files)
f_nested$sex<- "f"

m_files<- 'prom_nested_m_'
m_nested<- import_prom_pqlseq(m_files)
m_nested$sex<- "m"

prom_nested<- rbind(f_nested, m_nested)
prom_nested<- inner_join(mm_genes, prom_nested)

f_nested<- f_nested %>%
  mutate(std_beta = beta/se_beta) %>%
  dplyr::select(c(outcome, chrom, intercept, se_intercept, beta, se_beta, std_beta, pvalue, fdr))

colnames(f_nested)<- c(colnames(f_nested[, 1:2]), paste(colnames(f_nested[, 3:9]), "f", sep = "_"))

m_nested<- m_nested %>%
  mutate(std_beta = beta/se_beta) %>%
  dplyr::select(c(intercept, se_intercept, beta, se_beta, std_beta, pvalue, fdr))
  
colnames(m_nested)<- paste(colnames(m_nested), "m", sep = "_")


prom_paired<- cbind(f_nested, m_nested)

prom_paired<- prom_paired %>%
  mutate(abs_diff = abs(beta_f) - abs(beta_m),
         std_diff = abs(std_beta_f) - abs(std_beta_m))

prom_paired<- inner_join(mm_genes, prom_paired)
prom_paired<- prom_paired %>%
  arrange(chrom)

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
prom_pqlseq %>%
  filter(! chrom == "X") %>%
  #filter(gene_name %in% hallmark.msigdb$gene_symbol) %>%
  filter(fdr_f < 0.05 | fdr_m < 0.05) %>%
  filter(!outcome == "ENSMMUG00000006211") %>%
  filter(!is.na(gene_name) & !gene_name == "Metazoa_SRP") %>%
  ggplot(aes(beta_m, beta_f)) +
  geom_point(aes(colour=fdr_f<0.05),alpha=0.5, size=2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_smooth(method = "lm") +
  geom_abline() +
  #geom_text(aes(label=ifelse(beta_f>0.15 | beta_f< -0.15 | beta_m>0.15 | beta_m< -0.15, as.character(gene_name),'')), hjust=-0.1) +
  #geom_label_repel() +
  theme_classic()

x_prom %>%
  filter(fdr_f < 0.05 | fdr_m < 0.05) %>%
  filter(!gene_name == "Metazoa_SRP") %>%
  ggplot(aes(beta_m, beta_f, colour=fdr_f<0.05, shape = fdr_f<0.05 & fdr_m<0.05, label=gene_name)) +
  geom_point(size=2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  #geom_smooth(method = "lm") +
  geom_text_repel() +
  theme_classic()

prom_pqlseq %>%
  filter(fdr_f < 0.05 & fdr_m < 0.05) %>%
  pivot_longer(cols=c(beta_f, beta_m), names_to = "sex", values_to = "beta") %>%
  ggplot(aes(beta, fill=sex)) +
  geom_histogram(position = "identity", bins = 200, colour="black", alpha = 0.5) +
  theme_classic(base_size = 24)

prom_pqlseq %>%
  filter(fdr_f < 0.05 & fdr_m < 0.05) %>%
  pivot_longer(cols=c(beta_f, beta_m), names_to = "sex", values_to = "beta") %>%
  ggplot(aes(sex, beta, group = outcome)) +
  geom_point(alpha = 0.7) +
  geom_path(alpha = 0.3) +
  theme_classic(base_size = 24) +
  ylab("Beta") + xlab("") +
  scale_x_discrete(labels = c("Female", "Male"))

prom_pqlseq %>%
  filter(!gene_name == "Metazoa_SRP" | gene_name == "NA") %>%
  arrange(beta_m) %>%
  filter(fdr_m < 0.05) %>%
  dplyr::slice_head(n=10) %>%
  ggplot(aes(beta_m, gene_name)) +
  geom_col()

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
test<- left_join(sig_proms, paths, by = "gene_name")
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

