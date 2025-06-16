#!/usr/bin/env /packages/apps/spack/18/opt/spack/gcc-11.2.0/r-4.2.2-kpl/bin/Rscript

library(Rsamtools)
library(stringr)
options(scipen=9999999999999)

lid_pid <- read.table("/scratch/nsnyderm/cayo_rrbs/sort_bam/lid_pid")[,1]
bam_file <- paste0("/scratch/nsnyderm/cayo_rrbs/sort_bam/",lid_pid,".sort.bam")  # assuming filePaths is a vector of BAM file paths

# Extract all required fields and the XM tag
param <- ScanBamParam(
  what = c("qname","rname", "pos", "cigar", "qwidth"),  # V3, V4, V6, actual.read.length
  tag = "XM"
)

i=slurm_task_id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
bam_data <- scanBam(bam_file[i], param = param)[[1]]

SAM <- data.frame(
  CHROM = bam_data$rname,
  POS = bam_data$pos,
  CIGAR = bam_data$cigar,
  actual.read.length = bam_data$qwidth,
  plain.string = bam_data$tag$XM,
  query_name=bam_data$qname
)

SAM$zero.based.start <- SAM$POS - 1
SAM$zero.based.end <- SAM$zero.based.start + SAM$actual.read.length

## count number of CpGs in string
SAM$number.of.unmeth.cpg <- str_count(SAM$plain.string, "z")
SAM$number.of.meth.cpg <- str_count(SAM$plain.string, "Z")
SAM$number.of.cpg <- (SAM$number.of.meth.cpg + SAM$number.of.unmeth.cpg)
SAM$actual.read.length <- str_count(SAM$plain.string, "") ##should be same as CIGAR specification if no mismatches
## extract just the CpGs in the methylation string
SAM$cpg.meth.string <- gsub("[^zZ]", "", SAM$plain.string)
## calculate number of neighbors (i.e. observations of pairs that could be concordant/discordant)
SAM$num.neighbors <- (((SAM$number.of.cpg - 2)*2)+(2))
## set negative values (to cytosines with no cytosines) to zero
SAM$num.neighbors[SAM$num.neighbors < 0] <- 0
## count discordant occcurances of a string
SAM$zZ.occur <- str_count(SAM$cpg.meth.string,"zZ")
SAM$Zz.occur <- str_count(SAM$cpg.meth.string,"Zz")
## for each occurance of zZ or Zz, read lose 2 concordance points from max concordance
SAM$discordant.occur <- (SAM$zZ.occur + SAM$Zz.occur)
SAM$discordance.deduction <- (SAM$discordant.occur*2)
SAM$concordance.score <- (SAM$num.neighbors - SAM$discordance.deduction)
SAM$concordance.perc <- (SAM$concordance.score/SAM$num.neighbors)
SAM$read.discordance.perc <- (1 - SAM$concordance.perc) ##perc of discordance per read 

cleaned.SAM <- na.omit(SAM)

cleaned.SAM <- cleaned.SAM[,c("CHROM","zero.based.start","zero.based.end","query_name","read.discordance.perc")]

#Remove scientific notation from coordinates
#cleaned.SAM<- format(cleaned.SAM, scientific=F)

write.table(cleaned.SAM,row.names=F,quote=F, 
            file=paste0("/scratch/ckelsey4/Cayo_meth/epigenetic_drift/pdn_out/", lid_pid[i],".PDN_per_read.bed"),
            sep ="\t", col.names=F)

print(paste0("Completed PDN calc for number ",i," sample: ", lid_pid[i]))
