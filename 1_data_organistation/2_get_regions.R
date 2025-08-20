#!/usr/bin/env /packages/apps/spack/18/opt/spack/gcc-11.2.0/r-4.2.2-kpl/bin/Rscript

#load packages
library(bsseq)
library(comethyl)
library(tidyverse)
library(Biostrings)
library(GenomicFeatures)
library(GenomicRanges)

#set working directory
setwd("/scratch/ckelsey4/Cayo_meth")

#### IMPORT DATA ###############################################################
#Import metadata
blood_metadata<- read.table("blood_metadata_full.txt")

#Import cayo bsseq object and edit lid_pid colnames
#Autosomes/ChrX
cayo<- readRDS("/scratch/nsnyderm/cayo_rrbs/bismarkBSseq.rds")
colnames(cayo)=gsub(".CpG_report.merged_CpG_evidence.cov.gz","",str_split_i(colnames(cayo),"\\/",6))

#ChrY
cayo_Y<- readRDS("bismarkBSseq_Y.rds")
colnames(cayo_Y)=gsub(".CpG_report.merged_CpG_evidence.cov.gz","",colnames(cayo_Y))

#### GENERATE PROMOTER & GENE REGIONS ##########################################
#make txdb from macaque gff
macaque_txdb<- makeTxDbFromGFF("Macaca_mulatta.Mmul_10.110.chr.gtf")

#extract unique chromosomes
unique_chrs <- unique(seqlevels(macaque_txdb))

# Select the first 22(1-20 + X/Y) unique chromosome names as a redundancy against missing chromosome levels in the txdb object
# This will not match to the bsseq object if the chromosomes levels in txdb don't match for whatever reason
selected_chrs <- unique_chrs[1:21]
macaque_txdb <- keepSeqlevels(macaque_txdb, selected_chrs)

#Generate genes and promoters
macaque_genes = genes(macaque_txdb)
macaque_promoters=promoters(macaque_genes,upstream=2000,downstream=200)

#### GENERATE CPG REGIONS ######################################################
#Arrange metadata by lid
blood_metadata<- blood_metadata %>%
  arrange(lid_pid)

#Subset cayo bsseq with lids from longitudinal metadata
cayo_blood<- cayo[, cayo@colData@rownames %in% blood_metadata$lid_pid]
cayo_Y<- cayo_Y[, cayo_Y@colData@rownames %in% blood_metadata$lid_pid]
all.equal(cayo_blood@colData@rownames, blood_metadata$lid_pid)

#Split bsseq by chromosome
cayo_blood_list<- parallel::mclapply(selected_chrs,function(x){
  chrSelectBSseq(cayo_blood, seqnames = x, order = TRUE)
}, mc.cores=21)

names(cayo_blood_list)<- 1:21

#### PRE-FILTER CpG SITES ######################################################
#Set coverage and sample minimums
mincov<-5
per_sample<-0.25

#Filter cayo_list(chrs) object
cayo_filtered_list<- parallel::mclapply(names(cayo_blood_list), function(x){
  
  cayo_chr<- comethyl::filterCpGs(cayo_blood_list[[x]], cov = mincov, perSample = per_sample, verbose = FALSE,
                                  save = FALSE, file = NULL)
  return(cayo_chr)
  
}, mc.cores = 20)

Y_filtered<- comethyl::filterCpGs(cayo_Y, cov = mincov, perSample = per_sample, verbose = FALSE,
                                  save = FALSE, file = NULL)

cayo_filtered_list<- append(cayo_filtered_list, Y_filtered)
names(cayo_filtered_list)<- c(1:20, "X", "Y")

#### GENERATE REGIONS ##########################################################
#Generate regions 
regions<- parallel::mclapply(cayo_filtered_list,function(x){
  bsseq:::regionFinder3(x = as.integer(rep(1,length(x))), 
                        chr = as.character(GenomeInfoDb::seqnames(x)),
                        positions = BiocGenerics::start(x), maxGap = 1000, verbose = FALSE)[["up"]]
},mc.cores=20)

#Paste region coordinates in front of variables
regions<- lapply(regions, function(x){
  rownames(x)<- paste(x$chr, x$start, x$end, sep = "_");
  x
})

#Convert list to df
do.call(rbind, regions)->regions_df

regions_df<- regions_df %>%
  mutate(length = 1+(end - start)) %>%
  relocate(length, .after = end)

#Removes pseudo-duplicate regions (consecutive large regions that start 1bp apart)
regions_df<- regions_df %>% distinct(., cluster, length, .keep_all=T)

#Select region start, end and chrom coordinates
regions_cov<- regions_df %>%
  dplyr::select(chr, start, end)

#Make chrs numeric (X = chr 21)
#regions_cov$chr[regions_cov$chr == "X"]<- 21
#regions_cov$chr[regions_cov$chr == "Y"]<- 22
#regions_cov<- regions_cov %>%
  #mutate_at(vars(chr), as.integer)

#### WRITE RDS FILES ###########################################################
#write cayo_filtered_list to rds for get_coverage.R script
write_rds(cayo_filtered_list, "cayo_filtered_list")
write_rds(cayo_blood_list, "cayo_blood_list")

#Write gene/promoter GRanges files for get_coverage.R script
write_rds(macaque_genes, "macaque_genes")
write_rds(macaque_promoters, "macaque_promoters")

#write regions df to .rds file for get_coverage.R script
write_rds(regions_cov, "regions_filtered")
write_rds(regions_df, "regions_full")



