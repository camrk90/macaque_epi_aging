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

