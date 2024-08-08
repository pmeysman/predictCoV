library(ggplot2)
library(dplyr)
library(tidyr)  
library(ggsignif)
library(patchwork)

#Set working directory to source
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

dir= '../results/'

cov <- read.csv(paste0(dir,'allfirst_cov.csv'),header=T)

cov$tcrid <- paste0(cov$TCR,cov$sample)

cov_nodupl <- cov %>% group_by(tcrid,sample,Epitope_df2) %>% summarise(cloneFraction_df1 = sum(cloneFraction_df1),
                cloneFraction_df2=sum(cloneFraction_df2),
                corrected_p_value=mean(corrected_p_value))

print(cov_nodupl)

# Reshape the data frame to long format
cov_long <- cov_nodupl %>%
  pivot_longer(cols = c(cloneFraction_df1, cloneFraction_df2), names_to = "Time", values_to = "Value")

cov_long$cd8 <- lapply(cov_long$Epitope_df2,nchar) < 11
cov_long$exp <- cov_long$corrected_p_value < 0.05

# Create the line plot
p_cd8 <- ggplot(cov_long[cov_long$cd8 == TRUE,], aes(x = Time, y = Value, group = tcrid, color = sample, alpha = exp)) +
  scale_alpha_discrete(range = c(0.3, 1),name = "Expanded") + 
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(title = "MHC class I SARS-CoV-2 Spike-specific T-cells", x = "Time Point", y = "T-cell clonotype frequency",color='Individual') +
  scale_x_discrete(labels=c("cloneFraction_df1" = "D0 (pre-vaccination)","cloneFraction_df2" = "D2_28 (post-vaccination)"))

p_cd4 <- ggplot(cov_long[cov_long$cd8 == FALSE,], aes(x = Time, y = Value, group = tcrid, color = sample, alpha = exp)) +
  scale_alpha_discrete(range = c(0.3, 1),name = "Expanded") + 
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(title = "MHC class II  SARS-CoV-2 Spike-specific T-cells", x = "Time Point", y = "T-cell clonotype frequency",color='Individual') +
  scale_x_discrete(labels=c("cloneFraction_df1" = "D0 (pre-vaccination)","cloneFraction_df2" = "D2_28 (post-vaccination)"))

cov_byepitope <- cov_nodupl %>% group_by(Epitope_df2) %>% summarise(cloneFraction_df1= median(cloneFraction_df1),cloneFraction_df2= median(cloneFraction_df2))
cov_byepitope$logfold <- log((cov_byepitope$cloneFraction_df2 + 1e-8) / (cov_byepitope$cloneFraction_df1 + 1e-8))

p_epi<-ggplot(cov_byepitope[cov_byepitope$logfold > 0,], aes(x=Epitope_df2, y=logfold)) +
  geom_point() + 
  geom_segment( aes(x=Epitope_df2, xend=Epitope_df2, y=0, yend=logfold)) + theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  labs(title = "SARS-CoV-2 epitopes with expanded T-cells after vaccination",x="",y="Log-fold median change")

p_fig1 <- p_epi / (p_cd8 + p_cd4) + plot_annotation(tag_levels = "A") 
ggsave('../results/figure1.pdf',p_fig1, width = 30, height = 30,units = 'cm')

       