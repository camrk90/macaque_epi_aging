library(bsseq)
library(tidyverse)
library(comethyl)
setwd("/scratch/ckelsey4/Cayo_meth")

#Load in data-------------------------------------------------------------------
long_data<- readRDS("long_data_adjusted")
long_data<- long_data %>%
  arrange(lid_pid)

#Generate vector of chromosomes
chrs<- c(1:20, "X", "Y")

#Regions
regions<- readRDS("regions_filtered")

#Genes/Promoters
genes<- readRDS("macaque_genes")
promoters<- readRDS("macaque_promoters")

#Filtered Cayo bsseq for regions
cayo_filtered_list<- readRDS("cayo_filtered_list")
names(cayo_filtered_list)<- chrs

###################################
##### Generate M/Cov Matrices #####
###################################
#Promoters----------------------------------------------------------------------
prom_m_list<- parallel::mclapply(names(cayo_filtered_list),function(x){
  dd=getCoverage(cayo_filtered_list[[x]], regions = promoters[promoters@seqnames==x,], type = "M", 
                 what = "perRegionTotal")
  rownames(dd)=promoters[promoters@seqnames==x,]$gene_id
  return(as.data.frame(dd))
},mc.cores =20)

names(prom_m_list)<- chrs

prom_cov_list<- parallel::mclapply(names(cayo_filtered_list),function(x){
  dd=getCoverage(cayo_filtered_list[[x]], regions = promoters[promoters@seqnames==x,], type = "Cov", 
                 what = "perRegionTotal")
  rownames(dd)=promoters[promoters@seqnames==x,]$gene_id
  return(as.data.frame(dd))
},mc.cores =20)

names(prom_cov_list)<- chrs

#Genes--------------------------------------------------------------------------
genes_m_list<- parallel::mclapply(names(cayo_filtered_list),function(x){
  dd=getCoverage(cayo_filtered_list[[x]], regions = genes[genes@seqnames==x,], type = "M", 
                 what = "perRegionTotal")
  rownames(dd)=genes[genes@seqnames==x,]$gene_id
  return(as.data.frame(dd))
},mc.cores=20)

names(genes_m_list)<- chrs

genes_cov_list<- parallel::mclapply(names(cayo_filtered_list),function(x){
  dd=getCoverage(cayo_filtered_list[[x]], regions = genes[genes@seqnames == x,], type = "Cov", 
                 what = "perRegionTotal")
  rownames(dd) = genes[genes@seqnames == x,]$gene_id
  return(as.data.frame(dd))
},mc.cores = 20)

names(genes_cov_list)<- chrs

#Regions------------------------------------------------------------------------
regions_m_list<- parallel::mclapply(names(cayo_filtered_list),function(x){
  dd=getCoverage(cayo_filtered_list[[x]], regions = regions[regions$chr == x,], type = "M", 
                 what = "perRegionTotal")
  rownames(dd)<- rownames(regions[regions$chr == x,])
  return(as.data.frame(dd))
},mc.cores=20)

names(regions_m_list)<- chrs

regions_cov_list<- parallel::mclapply(names(cayo_filtered_list),function(x){
  dd=getCoverage(cayo_filtered_list[[x]], regions = regions[regions$chr == x,], type = "Cov", 
                 what = "perRegionTotal")
  rownames(dd)<- rownames(regions[regions$chr == x,])
  return(as.data.frame(dd))
},mc.cores=20)

names(regions_cov_list)<- chrs

###################################
#####  Filter M/Cov Matrices  #####
###################################
#Script filtering functions-----------------------------------------------------
#Filters for constitutively/non-constitutively methylated regions/proms/genes
filter_tails<- function(x, y){
  
  #Generate percent methylation matrix
  ratio<- mapply('/', x, y, SIMPLIFY = F)
  
  #Filter out methylated genes that are constituvively hyper/hypomethylated
  ratio<- sapply(names(ratio), function(i){
    
    df<- ratio[[i]]
    #Keep rows where the number of na's does not equal the number of columns
    df<- df[rowSums(is.na(df)) != ncol(df),]
    #Filter for %meth below 0.10 and above 0.90
    df<- df[!rowMeans(df, na.rm=T) > 0.90,] 
    df<- df[!rowMeans(df, na.rm=T) < 0.10,]
    
  },simplify = F, USE.NAMES = T)
  
  #Use filtered ratio df to filter coverage matrix
  cov<- sapply(names(y), function(i){
    
  cov_chr<- y[[i]]
  sub_chr<- ratio[[i]]
  cov_chr<- cov_chr[rownames(cov_chr) %in% rownames(sub_chr),]
  
  },simplify = F, USE.NAMES = T)
  
  return(cov)
}

#Filters samples with less than 70% of regions covered and regions with less than 25% coverageX
filter_regions<- function(region_list){
  
  #Unlist chrs to full matrix of regions across the genome
  df<- do.call(rbind, region_list[c(1:20, "X")])
  df_y<- region_list[["Y"]]
  
  #Remove NAs
  df<- df %>% drop_na()
  df_y<- df_y %>% drop_na()
  
  #Filter samples with less than 70% of regions covered and regions with less than 25% coverage
  df<- df[,colSums(df >= 1) >= 0.70*nrow(df)]
  df<- df[rowSums(df >= 1) >= 0.25*ncol(df),]
  
  df_y<- df_y[,colSums(df_y >= 1) >= 0.70*nrow(df_y)]
  df_y<- df_y[rowSums(df_y >= 1) >= 0.25*ncol(df_y),]
  
  #Generate col with chromosome number
  df$chr<- rownames(df)
  df$chr<- str_split_i(df$chr,"\\.",1)
  
  df_y$chr<- rownames(df_y)
  
  #Re-list regions
  r_list<- split(df[, -ncol(df)], df$chr)
  y_list<- list(df_y)
  names(y_list)<- "Y"
  r_list<- append(r_list, y_list)
  return(r_list)
}

#Filter promoter matrices-------------------------------------------------------
#Filter constitutively methylated regions
#prom_cov_filtered<- filter_tails(x=prom_m_list, y=prom_cov_list) #NOT FILTERING CONSTITUTIVELY METHYLATED PROMS/GENES AT THIS STAGE

#Filter samples with less than 70% of regions covered and regions with less than 25% coverage
prom_cov_filtered<- filter_regions(region_list = prom_cov_list)
prom_cov_filtered<- prom_cov_filtered[c(1:20, "X", "Y")]

#Subset M list with filtered cov list
prom_m_filtered<- sapply(names(prom_m_list), function(x){
  
  if (x == "Y") {
    df<- prom_m_list[[x]]
    sub<- prom_cov_filtered[[x]]
    df<- df[rownames(df) %in% rownames(sub),]
    df<- df[,colnames(df) %in% colnames(sub),]
  } else {
    df<- prom_m_list[[x]]
    rownames(df)<- paste(x, rownames(df), sep=".")
    sub<- prom_cov_filtered[[x]]
    df<- df[rownames(df) %in% rownames(sub),]
    df<- df[,colnames(df) %in% colnames(sub),]
  }
}, simplify = F, USE.NAMES = T)

#Filter gene matrices-----------------------------------------------------------
#Filter constitutively methylated regions
#genes_cov_filtered<- filter_tails(x=prom_m_list, y=genes_cov_list) #NOT FILTERING CONSTITUTIVELY METHYLATED PROMS/GENES AT THIS STAGE

#Filter samples with less than 70% of regions covered and regions with less than 25% coverage
genes_cov_filtered<- filter_regions(region_list = genes_cov_list)
genes_cov_filtered<- genes_cov_filtered[c(1:20, "X", "Y")]

#Subset M list with filtered cov list
genes_m_filtered<- sapply(names(genes_m_list), function(x){
  
  if (x == "Y") {
    df<- genes_m_list[[x]]
    sub<- genes_cov_filtered[[x]]
    df<- df[rownames(df) %in% rownames(sub),]
    df<- df[,colnames(df) %in% colnames(sub),]
  } else {
    df<- genes_m_list[[x]]
    rownames(df)<- paste(x, rownames(df), sep=".")
    sub<- genes_cov_filtered[[x]]
    df<- df[rownames(df) %in% rownames(sub),]
    df<- df[,colnames(df) %in% colnames(sub),]
  }
}, simplify = F, USE.NAMES = T)

#Filter region matrices---------------------------------------------------------
#Filter constitutively methylated regions
regions_cov_filtered<- filter_tails(x=regions_m_list, y=regions_cov_list)

#Filter samples with less than 70% of regions covered and regions with less than 25% coverage
regions_cov_filtered<- filter_regions(region_list = regions_cov_filtered)
regions_cov_filtered<- regions_cov_filtered[c(1:20, "X", "Y")]

#Subset M list with filtered cov list
regions_m_filtered<- sapply(names(regions_m_list), function(x){
  
  if (x == "Y") {
    df<- regions_m_list[[x]]
    sub<- regions_cov_filtered[[x]]
    df<- df[rownames(df) %in% rownames(sub),]
    df<- df[,colnames(df) %in% colnames(sub),]
  } else {
    df<- regions_m_list[[x]]
    rownames(df)<- paste(x, rownames(df), sep=".")
    sub<- regions_cov_filtered[[x]]
    df<- df[rownames(df) %in% rownames(sub),]
    df<- df[,colnames(df) %in% colnames(sub),]
  }
}, simplify = F, USE.NAMES = T)

#Save rds files-----------------------------------------------------------------
#Promoters
saveRDS(prom_m_filtered, "prom_m_filtered")
saveRDS(prom_cov_filtered, "prom_cov_filtered")

#Genes
saveRDS(genes_m_filtered, "genes_m_filtered")
saveRDS(genes_cov_filtered, "genes_cov_filtered")

#Regions
saveRDS(regions_m_filtered, "regions_m_filtered")
saveRDS(regions_cov_filtered, "regions_cov_filtered")


