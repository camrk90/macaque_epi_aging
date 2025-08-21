library(bsseq)
library(GenomicRanges)
library(GenomicFeatures)

cayo<- readRDS("/scratch/nsnyderm/cayo_rrbs/bismarkBSseq.rds")
cayo_Y<- readRDS("/scratch/ckelsey4/Cayo_meth/bismarkBSseq_Y.rds")
cayo_full<- combine(cayo, cayo_Y)

#Import genome annotation-------------------------------------------------------
#make txdb from macaque gff
macaque_txdb<- makeTxDbFromGFF("/scratch/ckelsey4/Cayo_meth/Macaca_mulatta.Mmul_10.110.chr.gtf")

gtf_mm<- rtracklayer::import('/scratch/ckelsey4/Cayo_meth/Macaca_mulatta.Mmul_10.110.chr.gtf')
gtf_mm_df=as.data.frame(gtf_mm)

#extract unique chromosomes
unique_chrs <- unique(seqlevels(macaque_txdb))

# Select the first 22(1-20 + X/Y) unique chromosome names as a redundancy against missing chromosome levels in the txdb object
# This will not match to the bsseq object if the chromosomes levels in txdb don't match for whatever reason
selected_chrs <- unique_chrs[1:22]
macaque_txdb <- keepSeqlevels(macaque_txdb, selected_chrs)

mincov<-5
per_sample<-0.25

Y_filtered<- comethyl::filterCpGs(cayo_Y, cov = mincov, perSample = per_sample, verbose = FALSE,
                                save = FALSE, file = NULL)

cayo_filtered_list<- readRDS("cayo_filtered_list")

cayo_filtered_list<- append(cayo_filtered_list, Y_filtered)
names(cayo_filtered_list)<- 1:22

## Generate regions for Y ######################################################
regions<- parallel::mclapply(cayo_filtered_list,function(x){
  bsseq:::regionFinder3(x = as.integer(rep(1,length(x))), 
                        chr = as.character(GenomeInfoDb::seqnames(x)),
                        positions = BiocGenerics::start(x), maxGap = 1000, verbose = FALSE)[["up"]]
},mc.cores =20)

#Paste region coordinates in front of variables
regions<- lapply(regions, function(x){
  rownames(x)<- paste(x$chr, x$start, x$end, sep = "_");
  x
})

#Convert list to df and filter for regions with minimum 5 cpg sites
do.call(rbind, regions)->regions_df
regions_df<- regions_df[regions_df$n >= 5, ]

regions_df<- regions_df %>%
  mutate(length = end - start) %>%
  relocate(length, .after = end)

regions_df %>% 
  ggplot(aes(reorder(as.factor(chr), length), length)) +
  geom_bar(stat = "summary", fun = "mean")

regions_cov<- regions_df %>%
  dplyr::select(chr, start, end)

#Make chrs numeric (X = chr 21)
regions_cov$chr[regions_cov$chr == "X"]<- 21
regions_cov<- regions_cov %>%
  mutate_at(vars(chr), as.integer)

save.image("y_chrom.RData")



#function to filter your chrom bsseq list
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

