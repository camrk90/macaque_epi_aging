library(tidyverse)
setwd("/scratch/ckelsey4/Cayo_meth")

#THIS SCRIPT IS FOR SELECTING LIDS WITH THE HIGHEST COV/M FOR THOSE THAT ARE DUPLICATED
#AT THE SAME AGE

long_data<- readRDS("wb_long_metadata")
meta_short<- long_data %>%
  group_by(monkey_id) %>%
  distinct(age_at_sampling, .keep_all = T) %>%
  mutate(n = n())
meta_short<- meta_short[!meta_short$n == 1,]
lids<- long_data[! long_data$lid_pid %in% meta_short$lid_pid,]

#Import m/cov rds------------------------------------------------------------
# load region lists that have been filtered for 5x coverage in 90% of samples
regions_cov_list<- readRDS("regions_cov_list")
regions_m_list<- readRDS("regions_m_list")

#Subset regions list for lids with repeated measures
regions_cov<- lapply(names(regions_cov_list), function(x){
  regions_cov<- subset(regions_cov_list[[x]], select=lids$lid_pid)
  return(regions_cov)
})

regions_m<- lapply(names(regions_m_list), function(x){
  regions_m<- subset(regions_m_list[[x]], select=lids$lid_pid)
  return(regions_m)
})

regions_cov<- do.call(rbind, regions_cov)
regions_m<- do.call(rbind, regions_m)
ratio<- as.matrix(regions_m/regions_cov)
ratio[is.nan(ratio)]<- 0

m_sum<- as.data.frame(colMeans(regions_m))
all.equal(rownames(m_sum), lids$lid_pid)
lids<- cbind(lids, m_sum)

cov_sum<- as.data.frame(colMeans(regions_cov))
all.equal(rownames(cov_sum), lids$lid_pid)
lids<- cbind(lids, cov_sum)

perc_meth<- as.data.frame(colMeans(ratio))
all.equal(rownames(perc_meth), lids$lid_pid)
lids<- cbind(lids, perc_meth)

lids_min<- lids %>%
  group_by(monkey_id) %>%
  slice_min(`colMeans(regions_cov)`, n = 1)

#Subset out lids with multiple entries at the same age
long_data<- long_data[!long_data$lid_pid %in% lids_min$lid_pid,]

saveRDS(long_data, 'long_data_adjusted')





