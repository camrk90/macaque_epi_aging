#THIS SCRIPT IS FOR GENERATING .BED FILES FOR GLM MODEL REGIONS AND REPEAT MASK 
#COORDINATES FROM THE UCSC BROWSER
#IT ALSO LIFTS OVER THE CHROMATIN STATE MAP COORDINATES FROM THE HUMAN GENOME TO
#THE MACAQUE GENOME SO THAT ALL THESE FILES CAN BE INTERSECTED WITH BEDTOOLS

library(tidyverse)
library(GenomicFeatures)
library(GenomicRanges)

#Import genome annotation-------------------------------------------------------
#make txdb from macaque gff
macaque_txdb<- makeTxDbFromGFF("/scratch/ckelsey4/Cayo_meth/Macaca_mulatta.Mmul_10.110.chr.gtf")

#extract unique chromosomes
unique_chrs <- unique(seqlevels(macaque_txdb))

# Select the first 22(1-20 + X/Y) unique chromosome names as a redundancy against missing chromosome levels in the txdb object
# This will not match to the bsseq object if the chromosomes levels in txdb don't match for whatever reason
selected_chrs <- unique_chrs[1:22]
macaque_txdb <- keepSeqlevels(macaque_txdb, selected_chrs)

#Generate genes and promoters
macaque_genes = genes(macaque_txdb)
macaque_promoters=promoters(macaque_genes,upstream=2000,downstream=200)

macaque_promoters<- as.data.frame(macaque_promoters)

macaque_promoters<- macaque_promoters %>%
  dplyr::select(c(seqnames, start, end))

macaque_promoters<- macaque_promoters %>%
  arrange(seqnames, start)

macaque_promoters$seqnames<- paste("chr", macaque_promoters$seqnames, sep = "")

#removes one entry with a negative start coordinate
macaque_promoters<- macaque_promoters %>%
  filter(!start < 0)

######################################
### Generate .bed file for regions ###
######################################
#Import metadata and regions rds------------------------------------------------
long_data<- readRDS("/scratch/ckelsey4/Cayo_meth/long_data_adjusted")
long_data<- long_data %>%
  arrange(lid_pid) %>%
  filter(n > 1)

regions_cov<- readRDS("/scratch/ckelsey4/Cayo_meth/regions_cov_filtered")

#Filter metadata to lids in regions list
long_data<- long_data[long_data$lid_pid %in% colnames(regions_cov[[1]]),]

#Subset regions list for lids n > 1
regions_cov<- lapply(names(regions_cov), function(x){
  regions_cov<- subset(regions_cov[[x]], select=long_data$lid_pid)
  return(regions_cov)
})

#Check metadata lids match the lids (cols) of a random chromosome
all.equal(long_data$lid_pid, colnames(regions_cov[[runif(1, 1, 21)]]))

cov<- do.call(rbind, regions_cov)
regions<- rownames(cov)
regions<- str_split_i(regions, "\\.", 3)
regions<- as.data.frame(regions)
regions<- as.data.frame(regions) %>%
  separate_wider_delim(regions, "_", names = c("chrom", "chromStart", "chromEnd"))

#Arrange df by chromosome
regions$chrom[regions$chrom == "X"]<- 21
regions<- regions %>%
  arrange(as.integer(chrom))
regions$chrom[regions$chrom == 21]<- "X"

#Add 'chr' string to chr col
regions$chrom<- gsub(" ", "", paste("chr", regions$chrom, ""))

#Generate .bed file
regions_bed<- GenomicRanges::makeGRangesFromDataFrame(regions)


######################################
###   Generate .bed file for TEs   ###
######################################
#Import repeatmasker repeat coordinates-----------------------------------------
repeats<- read.table("rmsk.txt", header = F)
repeat_header<- cbind("bin", "swScore",	"milliDiv",	"milliDel",	"milliIns",	"genoName",	"genoStart", "genoEnd", "genoLeft",
                      "strand",	"repName", "repClass", "repFamily", "repStart",	"repEnd",	"repLeft", "id")
colnames(repeats)<- repeat_header

repeats$genoName<- gsub("_.*", "", repeats$genoName)

repeats<- repeats %>%
  filter(!genoName == "chrUn") %>%
  filter(!genoName == "chrY") %>%
  filter(!genoName == 'chrM')

#Export repeats as bed file for intersect
repeats_bed<- repeats %>%
  dplyr::select(c(genoName, genoStart, genoEnd, repName, repClass)) 
colnames(repeats_bed)<- c("chrom", "chromStart", "chromEnd", "repName", "repClass")

#Add 1 to coordinates to account for the fact that exporting to bed will subtract 1 to make it 0-based
#the repeats file is already 0-based
repeats_bed<- repeats_bed %>%
  mutate(chromStart = chromStart + 1, chromEnd = chromEnd + 1)

#Generate range col to more easily match joins
repeats_bed<- repeats_bed %>%
  mutate(range = paste(as.character(chromStart), "-", as.character(chromEnd)))
repeats_bed$chrom<- gsub("chr", "", repeats_bed$chrom)
repeats_bed$chrom[repeats_bed$chrom == "X"]<- 21
repeats_bed<- repeats_bed %>%
  arrange(as.integer(chrom))
repeats_bed$chrom[repeats_bed$chrom == 21]<- "X"
repeats_bed$chrom<- gsub(" ", "", paste("chr", repeats_bed$chrom, ""))


######################################
###  Liftover Macaque CHMM states  ###
######################################
#Import bed file----------------------------------------------------------------
bed_file<- ("E062_15_coreMarks_hg38lift_dense.bed")
chmm_bed<- rtracklayer::import(bed_file, format = "bed")

#Import ucsc chain file
hg38_mmul10<- rtracklayer::import.chain("hg38ToRheMac10.over.chain")

#Liftover hg38 coords to mmul10 assembly and sort
seqlevelsStyle(chmm_bed) = "UCSC"
chmm_mmul<- rtracklayer::liftOver(chmm_bed, hg38_mmul10)
chmm_mmul<- unlist(chmm_mmul)
genome(chmm_mmul) = "mmul_10"
chmm_mmul<- sortSeqlevels(chmm_mmul)
selected_chrs<- gsub(" ", "", paste("chr", selected_chrs))
chmm_mmul2<- keepSeqlevels(chmm_mmul, selected_chrs, pruning.mode = "coarse")
chmm_mmul2<- sort(chmm_mmul2)

#Export files as .bed for intersect---------------------------------------------
setwd('/scratch/ckelsey4/Cayo_meth/intersect_files')
rtracklayer::export.bed(macaque_promoters, con = "promoters.bed")
rtracklayer::export.bed(regions_bed, con = "glm_regions.bed")
rtracklayer::export.bed(repeats_bed, con = "repeats.bed")
rtracklayer::export.bed(chmm_mmul2, con = "chmm_mmul.bed")

