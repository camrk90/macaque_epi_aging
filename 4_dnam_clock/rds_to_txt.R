
#### Cov Counts ####
regions_cov<- readRDS("/Users/cameronkelsey/Documents/smack_lab/cayo_data/dnam_clock/regions_cov_filtered")
regions_cov<- regions_cov[1:21]

#Bind list elements to df
regions_cov<- do.call(rbind, regions_cov)

#Reduce rownames and add chr column
rownames(regions_cov)<- str_split_i(rownames(regions_cov), '\\.', 4)

write.table(regions_cov, "/Users/cameronkelsey/Documents/smack_lab/cayo_data/dnam_clock/regions_cov.txt", quote=F)

#### M counts ####
regions_m<- readRDS("/Users/cameronkelsey/Documents/smack_lab/cayo_data/dnam_clock/regions_m_filtered")
regions_m<- regions_m[1:21]

#Bind list elements to df
regions_m<- do.call(rbind, regions_m)

#Reduce rownames and add chr column
rownames(regions_m)<- rownames(regions_cov)

write.table(regions_m, "/Users/cameronkelsey/Documents/smack_lab/cayo_data/dnam_clock/regions_m.txt", quote=F)
