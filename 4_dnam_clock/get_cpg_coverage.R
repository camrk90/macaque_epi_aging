#!/usr/bin/env /packages/apps/spack/18/opt/spack/gcc-11.2.0/r-4.2.2-kpl/bin/Rscript

SAMP <- Sys.getenv("SLURM_ARRAY_TASK_ID")
SAMP <- as.integer(SAMP)

library(bsseq)

#Load in data-------------------------------------------------------------------
cayo_filtered_list<- readRDS("/scratch/ckelsey4/Cayo_meth/cayo_filtered_list")
cayo_filtered<- cayo_filtered_list[[SAMP]]

#Generate coverage for individual CpGs
#M count
m=getCoverage(cayo_filtered, regions = NULL, type = "M", what = "perBase")
rownames(m)<- cayo_filtered@rowRanges@ranges@start
m<- as.data.frame(m)

#Coverage
cov=getCoverage(cayo_filtered, regions = NULL, type = "Cov", what = "perBase")
rownames(cov)<- cayo_filtered@rowRanges@ranges@start
cov<- as.data.frame(cov)

file_path<- paste("/scratch/ckelsey4/Cayo_meth/dnam_clock/", SAMP, sep="")

#Write files as txt
write.table(m, paste(file_path, "cpg_m", "txt", sep="."), quote = F)
write.table(cov, paste(file_path, "cpg_cov", "txt", sep="."), quote = F)

