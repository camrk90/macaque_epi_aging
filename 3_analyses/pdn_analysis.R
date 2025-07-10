library(tidyverse)
library(ggplot2)
library(lme4)

#Import metadata----------------------------------------------------------------
long_data<- readRDS("/scratch/ckelsey4/Cayo_meth/long_data_adjusted")
long_data<- long_data %>%
  dplyr::rename(sex = individual_sex) %>%
  filter(n > 1) %>%
  arrange(lid_pid)

#Import PDN and COV lists
avg_pdn<- readRDS("/scratch/ckelsey4/Cayo_meth/epigenetic_drift/avg_pdn/avg_pdn_list")

clean_list<- function(x){
  x<- x %>%
    mutate_if(is.character, as.numeric) %>%
    select(-chrom)
  x<- x[!rowSums(is.na(x)) >= 0.25*ncol(x),]
}

avg_pdn<- lapply(avg_pdn, clean_list)

read_cov<- readRDS("/scratch/ckelsey4/Cayo_meth/epigenetic_drift/avg_pdn/read_cov_list")
subset_cov<- function(x){
  df<- read_cov[[x]]
  df<- df %>%
    mutate_if(is.character, as.numeric) %>%
    select(-chrom)
  df<- df[rownames(df) %in% rownames(avg_pdn[[x]]),]
}

read_cov<- lapply(names(read_cov), subset_cov)

avg_pdn<- do.call(rbind, avg_pdn)
read_cov<- do.call(rbind, read_cov)

mean_pdn<- colMeans(avg_pdn, na.rm=T)
mean_cov<- colMeans(read_cov)

long_data$mean_pdn<- mean_pdn
long_data$mean_cov<- mean_cov

global_pdn=lmer(mean_pdn ~ within.age + sex + (1 + within.age|monkey_id),
          data=long_data)
?lmer
summary(global_pdn)

#Import avg_pdn per region RDS--------------------------------------------------
avg_pdn_f<- readRDS("/scratch/ckelsey4/Cayo_meth/epigenetic_drift/avg_pdn/_pdn_f")
avg_pdn_f<- avg_pdn_f %>%
  filter(converged == TRUE)
padj<- p.adjust(avg_pdn_f$pvalue, method = "fdr")
avg_pdn_f$padj<- padj
avg_pdn_f$sex<- "f"

avg_pdn_m<- readRDS("/scratch/ckelsey4/Cayo_meth/epigenetic_drift/avg_pdn/_pdn_m")
padj<- p.adjust(avg_pdn_m$pvalue, method = "fdr")
avg_pdn_m$padj<- padj
avg_pdn_m<- avg_pdn_m %>%
  filter(converged == TRUE)
avg_pdn_m$sex<- "m"

pdn_full<- rbind(avg_pdn_f, avg_pdn_m)
pdn_full<- pdn_full %>%
  group_by(outcome) %>%
  arrange(max(beta))

pdn_full %>%
  filter(padj < 0.05) %>%
  ggplot(aes(beta, fill=sex)) +
  geom_histogram(alpha = 0.5, colour="black", bins = 100, position = "identity") + 
  scale_fill_manual(values = c("darkolivegreen", "darkmagenta"), name = "Sex")

pdn_full %>%
  ggplot(aes(padj, fill=sex)) +
  geom_histogram(alpha = 0.5, colour="black", bins = 100, position = "identity") + 
  scale_fill_manual(values = c("darkolivegreen", "darkmagenta"), name = "Sex")
  

m<- avg_pdn_m[, c("beta", "padj")]
colnames(m)<- paste(colnames(m), "m", sep="_")
f<- avg_pdn_f[, c("outcome", "beta", "padj")]
colnames(f)<- paste(colnames(f), "f", sep="_")

m_f<- cbind(f, m)
m_f<- m_f %>%
  mutate(diff = abs(beta_f) - abs(beta_m))

m_f %>%
  filter(padj_f < 0.05 | padj_m < 0.05) %>%
  ggplot(aes(beta_m, beta_f, colour = diff)) +
  geom_point() +
  geom_abline() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_gradient2(low = "darkmagenta", mid = "white", high = "darkolivegreen", midpoint = 0, name = "") +
  geom_smooth(method = "lm") +
  #ylim(-0.3, 0.3) +
  #xlim(-0.3, 0.3) +
  xlab("Estimate (Male)") +
  ylab("Estimate (Female)") +
  theme_classic(base_size = 36) +
  theme(legend.key.height= unit(2, 'cm')) +
  ggtitle("padj_f < 0.05 | padj_m < 0.05")

m_f %>%
  filter(padj_f < 0.05 | padj_m < 0.05) %>%
  ggplot(aes(diff, fill = after_stat(x))) +
  geom_histogram(bins = 100, colour = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = 'red') +
  scale_fill_gradient2(low = "darkmagenta", mid = "white", high = "darkolivegreen", midpoint = 0, name = "Beta Difference") +
  xlab("|Estimate (F)| - |Estimate (M)|") +
  ggtitle("padj_f < 0.05 | padj_m < 0.05")

long_data %>%
  ggplot(aes(within.age, mean_pdn, colour=sex)) +
  geom_point() +
  geom_smooth(method = "lm") +
  geom_vline(xintercept = 0, linetype = "dashed")
  
  


