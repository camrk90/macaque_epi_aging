#!/usr/bin/env /packages/apps/spack/18/opt/spack/gcc-11.2.0/r-4.2.2-kpl/bin/Rscript

SAMP <- Sys.getenv("SLURM_ARRAY_TASK_ID")
SAMP <- as.integer(SAMP)

##just need to change working directory, file names, and file paths and run in directory with files in it. 
##input is SAM format resulting from Bismark alignment (with no header)

library(dplyr)
library(stringr)
library(readr)

##set unique file ID's
unique_fileID<- list.files(path = "/scratch/ckelsey4/Cayo_meth/epigenetic_drift/sam_files")[SAMP]
unique_fileID<- str_split_i(unique_fileID, "\\.", 1)

##set file paths
filePath<- paste0("/scratch/ckelsey4/Cayo_meth/epigenetic_drift/sam_files/", unique_fileID, ".sort.sam")

##make forloop to calculate PDN
##read in array position for job
##read in SAM file written as a plain text file with headers removed
SAM <- read.table(filePath, header = FALSE, fill = TRUE)
## remove XM:Z: string tag from methylation string (column V14)
SAM$plain.string <- str_remove_all(SAM$V14, "XM:Z:")
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
## print only read inforation as a text file for each sample
SAM$zero.based.start <- (SAM$V4 - 1)
SAM$zero.based.end <- (SAM$zero.based.start + SAM$actual.read.length)
SAM$POS <- SAM$V4
SAM$CIGAR <- SAM$V6
SAM$CHROM <- SAM$V3
## remove NAs so you can calculate average and remove any reads with 1 or fewer CpGs
cleaned.SAM <- na.omit(SAM)
## write.table(cleaned.SAM, file=paste0(opt$outputDir,"/", opt$fileID, "_cleaned_SAM.txt"))
## now make txt file in bed format for each with the 0 based start and end positions after removing cytosines without a neighbor on the nearest read (0 or 1 CpG)
perRead_bed <- cleaned.SAM[,c("CHROM","zero.based.start","zero.based.end","V1","read.discordance.perc")]
## write the result table in an individual text file
write.table(perRead_bed, file=paste0("/scratch/ckelsey4/Cayo_meth/epigenetic_drift/pdn_out/", unique_fileID, "_PDN_per_read.txt"))



