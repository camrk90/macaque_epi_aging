#!/usr/bin/env /packages/apps/spack/18/opt/spack/gcc-11.2.0/r-4.2.2-kpl/bin/Rscript

library(bsseq)

#Load in data-------------------------------------------------------------------
cayo_filtered_list<- readRDS("/scratch/ckelsey4/Cayo_meth/cayo_filtered_list")
cayo_filtered_list<- cayo_filtered_list[1:21]

#Generate coverage for individual CpGs
cpg_m_list<- parallel::mclapply(names(cayo_filtered_list),function(x){
  dd=getCoverage(cayo_filtered_list[[x]], regions = NULL, type = "M", 
                 what = "perBase")
  rownames(dd)<- cayo_filtered_list[[x]]@rowRanges@ranges@start
  return(as.data.frame(dd))
}, mc.cores=21)

names(cpg_m_list)<- 1:21
cpg_m<- do.call(rbind, cpg_m_list)
rm(cpg_m_list)

cpg_cov_list<- parallel::mclapply(names(cayo_filtered_list),function(x){
  dd=getCoverage(cayo_filtered_list[[x]], regions = NULL, type = "Cov", 
                 what = "perBase")
  rownames(dd)<- cayo_filtered_list[[x]]@rowRanges@ranges@start
  return(as.data.frame(dd))
}, mc.cores=21)

names(cpg_cov_list)<- 1:21
cpg_cov<- do.call(rbind, cpg_cov_list)
rm(cpg_cov_list)

write.table(cpg_m, "/scratch/ckelsey4/Cayo_meth/dnam_clock/cpg_m.txt", quote = F)

write.table(cpg_cov, "/scratch/ckelsey4/Cayo_meth/dnam_clock/cpg_cov.txt", quote = F)

