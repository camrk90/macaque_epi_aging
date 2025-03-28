setwd("/scratch/ckelsey4/Cayo_meth")
######################################
###        Chromatin States        ###
######################################
#Import chmm state intersect files
chmm_intersect<- read.table("./intersect_files/region_chmm_intersect.txt", header = F)
chmm_intersect<- chmm_intersect %>%
  dplyr::select(-c(V5, V6, V7, V8, V9, V13, V14, V15))
colnames(chmm_intersect)<- c("chr", "anno_start", "anno_end", "anno", "region_chr", "chromStart", "region_end", "overlap")
chmm_intersect<- chmm_intersect %>%
  mutate(anno_start = anno_start+1, chromStart = chromStart+1)

#Factor chrom levels
chmm_intersect$chr<- gsub("chr", "", chmm_intersect$chr)
chmm_intersect$chr<- factor(chmm_intersect$chr, levels = c(1:21))

######################################
###        Repeat Elements         ###
######################################
#Import and edit RE track intersect 
repeat_intersect<- read.table("./intersect_files/region_repeat_intersect.txt", header = F)
repeat_intersect<- repeat_intersect %>%
  dplyr::select(-c(V4, V5, V6, V10, V11, V12))
colnames(repeat_intersect)<- c("chr", "anno_start", "anno_end", "region_chr", "region_start", "region_end", "overlap")
repeat_intersect<- repeat_intersect %>%
  mutate(anno_end = anno_end-1)

#Generate range col to more easily match joins
#Repeat bed file
repeats_bed<- read.table("rmsk.txt")

#Select coordinate and TE name/class cols
repeats_bed<- repeats_bed[, c(6:8, 11:12)]
colnames(repeats_bed)<- c("chr", "chromStart", "chromEnd", "repName", "repClass")

#Remove nonsensical chromosome names eg "chr10_NW_021160243v1_random"
repeats_bed$chr<- gsub("_.*", "", repeats_bed$genoName)

#Remove mitochondria, unknown chromosome entries and "chr" string
repeats_bed<- repeats_bed %>%
  filter(!chr == "chrUn") %>%
  filter(!chr == 'chrM')

#Refactor chr levels
repeats_bed$chr<- gsub("chr", "", repeats_bed$chr)
repeats_bed$chr<- factor(repeats_bed$chr, levels = c(1:21, 'X'))

#Generate range col
repeats_bed<- repeats_bed %>%
  mutate(range = paste(as.character(chromStart), "-", as.character(chromEnd)))
repeats_bed<- repeats_bed %>%
  dplyr::select(c(repName, repClass, range))

#Generate range col, remove "chr" string, and factor chr col for intersect file 
repeat_intersect<- repeat_intersect %>%
  mutate(range = paste(as.character(anno_start), "-", as.character(anno_end)))
repeat_intersect$chr<- gsub("chr", "", repeat_intersect$chr)
repeat_intersect$chr<- factor(repeat_intersect$chr, levels = c(1:21, 'X'))

#Join intersect file and original bed file so intersect coordinates have annotation names
re_anno<- inner_join(repeat_intersect, repeats_bed, by = "range")

#Rename region start col to easily join with model dfs
re_anno<- re_anno %>%
  dplyr::rename(chromStart = region_start)
re_anno<- re_anno %>%
  mutate(chromStart = chromStart+1)

write_csv(re_anno, file = "re_annotations.csv")
write_csv(chmm_intersect, file = "chmm_annotations.csv")



