library(tidyverse)
library(purrr)
library(ggplot2)
library(PQLseq2)

#Import metadata----------------------------------------------------------------
long_data<- readRDS("/scratch/ckelsey4/Cayo_meth/long_data_adjusted")
long_data<- long_data %>%
  dplyr::rename(sex = individual_sex) %>%
  filter(n > 1) %>%
  arrange(lid_pid)

#Import avg_pdn per region RDS--------------------------------------------------
avg_pdn_f<- readRDS("/scratch/ckelsey4/Cayo_meth/epigenetic_drift/avg_pdn/avg_pdn_f")
padj<- p.adjust(avg_pdn_f$pvalue, method = "fdr")
avg_pdn_f$padj<- padj
avg_pdn_f<- avg_pdn_f %>%
  filter(converged == TRUE)

avg_pdn_m<- readRDS("/scratch/ckelsey4/Cayo_meth/epigenetic_drift/avg_pdn/avg_pdn_m")
padj<- p.adjust(avg_pdn_m$pvalue, method = "fdr")
avg_pdn_m$padj<- padj
avg_pdn_m<- avg_pdn_m %>%
  filter(converged == TRUE)

m<- avg_pdn_m[, c("beta", "padj")]
colnames(m)<- paste(colnames(m), "m", sep="_")
f<- avg_pdn_f[, c("outcome", "beta", "padj")]
colnames(f)<- paste(colnames(f), "f", sep="_")

m_f<- cbind(f, m)
m_f<- m_f %>%
  mutate(diff = abs(beta_f) - abs(beta_m))

m_f %>%
  #filter(fdr_f < 0.05 | fdr_m < 0.05) %>%
  ggplot(aes(beta_m, beta_f, colour = diff)) +
  geom_point() +
  geom_abline() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_gradient2(low = "darkmagenta", mid = "white", high = "darkolivegreen", midpoint = 0, name = "") +
  geom_smooth(method = "lm") +
  #ylim(-0.5, 0.5) +
  #xlim(-0.5, 0.5) +
  xlab("Estimate (Male)") +
  ylab("Estimate (Female)") +
  theme_classic(base_size = 36) +
  theme(legend.key.height= unit(2, 'cm'))

m_f %>%
  ggplot(aes(diff, fill = after_stat(x))) +
  geom_histogram(bins = 100, colour = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = 'red') +
  scale_fill_gradient2(low = "darkmagenta", mid = "white", high = "darkolivegreen", midpoint = 0, name = "Beta Difference") +
  xlab("|Estimate (F)| - |Estimate (M)|")


