#!/usr/bin/env /packages/apps/spack/18/opt/spack/gcc-11.2.0/r-4.2.2-kpl/bin/Rscript

#SBATCH -c 21
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ckelsey4@asu.edu

#load packages
library(bsseq)
library(comethyl)
library(tidyverse)
library(Biostrings)
library(GenomicFeatures)
library(GenomicRanges)

#set working directory
setwd("/scratch/ckelsey4/Cayo_meth/")

#### IMPORT DATA ###############################################################
#Import metadata
blood_metadata<- read.table("blood_metadata_full.txt")
cayo_blood_list<- readRDS("cayo_blood_full.rds")

#### GENERATE PROMOTER & GENE REGIONS ##########################################
#make txdb from macaque gff
macaque_txdb<- makeTxDbFromGFF("Macaca_mulatta.Mmul_10.110.chr.gtf")

#extract unique chromosomes
unique_chrs <- unique(seqlevels(macaque_txdb))

# Select the first 22(1-20 + X/Y) unique chromosome names as a redundancy against missing chromosome levels in the txdb object
# This will not match to the bsseq object if the chromosomes levels in txdb don't match for whatever reason
selected_chrs <- unique_chrs[1:22]
macaque_txdb <- keepSeqlevels(macaque_txdb, selected_chrs)

#Generate genes and promoters
macaque_genes = genes(macaque_txdb)
macaque_promoters=promoters(macaque_genes,upstream=2000,downstream=200)

#### PRE-FILTER CpG SITES ######################################################
#Set coverage and sample minimums
mincov<-5
per_sample<-0.25

#Filter cayo_list(chrs) object
cayo_filtered_list<- parallel::mclapply(names(cayo_blood_list), function(x){
  
  cayo_chr<- comethyl::filterCpGs(cayo_blood_list[[x]], cov = mincov, perSample = per_sample, verbose = FALSE,
                                  save = FALSE, file = NULL)
  return(cayo_chr)
  
}, mc.cores = 21)

names(cayo_filtered_list)<- c(1:20, "X")

#### GENERATE REGIONS ##########################################################
#Generate regions 
regions<- parallel::mclapply(cayo_filtered_list,function(x){
  bsseq:::regionFinder3(x = as.integer(rep(1,length(x))), 
                        chr = as.character(GenomeInfoDb::seqnames(x)),
                        positions = BiocGenerics::start(x), maxGap = 1000, verbose = FALSE)[["up"]]
},mc.cores=21)

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

#### WRITE RDS FILES ###########################################################
#write cayo_filtered_list to rds for get_coverage.R script
saveRDS(cayo_filtered_list, "cayo_filtered_list.rds")

#Write gene/promoter GRanges files for get_coverage.R script
saveRDS(macaque_genes, "macaque_genes.rds")
saveRDS(macaque_promoters, "macaque_promoters.rds")

#write regions df to .rds file for get_coverage.R script
saveRDS(regions_cov, "regions_filtered.rds")
saveRDS(regions_df, "regions_full.rds")



