library(tidyverse)

regions_cov<- readRDS("/Users/cameronkelsey/Documents/smack_lab/cayo_data/regions_cov_filtered.rds")
regions_cov<- do.call(rbind, regions_cov)
write.table(regions_cov, "/Users/cameronkelsey/Documents/smack_lab/cayo_data/regions_cov_filtered.txt", 
            quote=F)

regions_m<- readRDS("/Users/cameronkelsey/Documents/smack_lab/cayo_data/regions_m_filtered.rds")
regions_m<- do.call(rbind, regions_m)
write.table(regions_m, "/Users/cameronkelsey/Documents/smack_lab/cayo_data/regions_m_filtered.txt", 
            quote=F)

head(regions_cov)
