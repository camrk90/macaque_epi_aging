library(tidyverse)

setwd("/scratch/ckelsey4/Cayo_meth")
######################################
###        Chromatin States        ###
######################################
#Import chmm state intersect files
chmm_intersect<- read.table("./intersect_files/region_chmm_intersect.txt", header = F)
chmm_intersect<- chmm_intersect %>%
  dplyr::select(c(V1, V2, V3, V4, V11, V17, V18))
colnames(chmm_intersect)<- c("chr", "anno_start", "anno_end", "anno", "cpg_loc", "region_start", "region_end")

#Factor chrom levels
chmm_intersect$chr<- gsub("chr", "", chmm_intersect$chr)
chmm_intersect$chr<- factor(chmm_intersect$chr, levels = c(1:22, "X", "Y"))

#Add range column
chmm_intersect<- chmm_intersect %>%
  mutate(region_range = paste(as.character(region_start), "-", as.character(region_end)))

######################################
###        Repeat Elements         ###
######################################
#Import and edit RE track intersect 
repeat_intersect<- read.table("./intersect_files/region_repeat_intersect.txt", header = F)
repeat_intersect<- repeat_intersect %>%
  dplyr::select(c(V1, V2, V3, V8, V14, V15))
colnames(repeat_intersect)<- c("chr", "anno_start", "anno_end", "cpg_loc", "region_start", "region_end")

#Factor chrom levels and remove "chr" string
repeat_intersect$chr<- gsub("chr", "", repeat_intersect$chr)
repeat_intersect$chr<- factor(repeat_intersect$chr, levels = c(1:22, "X", "Y"))

#Generate range col to facilitate join to repeat annotation df 
repeat_intersect<- repeat_intersect %>%
  mutate(range = paste(as.character(anno_start), "-", as.character(anno_end))) %>%
  mutate(region_range = paste(as.character(region_start), "-", as.character(region_end)))
  relocate(range, .before = chr)

#Generate range col to more easily match joins
#Repeat bed file
repeats_bed<- read.table("rmsk.txt")

#Select coordinate and TE name/class cols
repeats_bed<- repeats_bed[, c(6:8, 11:12)]
colnames(repeats_bed)<- c("chr", "chromStart", "chromEnd", "repName", "repClass")

#Remove nonsensical chromosome names eg "chr10_NW_021160243v1_random", mitchondrial and unknown chrs
repeats_bed$chr<- gsub("chr", "", repeats_bed$chr)
repeats_bed<- repeats_bed %>%
  filter(chr %in% repeat_intersect$chr)

#Refactor chr levels
repeats_bed$chr<- factor(repeats_bed$chr, levels = c(1:21, 'X', 'Y'))
repeats_bed<- repeats_bed %>%
  arrange(chr)

#Generate range col and subset to range and annotation cols
repeats_bed<- repeats_bed %>%
  mutate(range = paste(as.character(chromStart), "-", as.character(chromEnd)))
repeats_bed<- repeats_bed %>%
  dplyr::select(c(chr, repName, repClass, range))

#Join intersect file and original bed file so intersect coordinates have annotation names
re_anno<- right_join(repeats_bed, repeat_intersect, by = c("range", "chr"), keep = F)

write_csv(re_anno, file = "re_annotations.csv")
write_csv(chmm_intersect, file = "chmm_annotations.csv")



