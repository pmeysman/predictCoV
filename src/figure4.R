library(ggplot2)
library(dplyr)
library(tidyr)  
library(ggsignif)
library(patchwork)

#Set working directory to source
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

dir= '../results/'

cd8 <- read.table(paste0(dir,'cov_enrich_CD8.txt'),header=T)

print(cd8)

# Reshape the data frame to long format
cd8_long <- cd8 %>%
  pivot_longer(cols = c(D1D2, D3, D4), names_to = "Time", values_to = "Value")

# Convert Time to factor with specific order
cd8_long$Time <- factor(cd8_long$Time, levels = c("D1D2", "D3", "D4"))

# Create the line plot
p_cd8 <- ggplot(cd8_long, aes(x = Time, y = Value, group = Sample)) +
  geom_line(aes(color = Sample)) +
  geom_point(aes(color = Sample)) +
  theme_minimal() +
  labs(title = "Expanding spike-specific CD8+ T-cells", x = "Time Point", y = "Unique CD8+ T Cells")


# Perform the paired non-parametric statistical test (Wilcoxon signed-rank test)
# Compare D1D2 vs D3, D1D2 vs D4, and D3 vs D4
wilcox_test_D1D2_D3 <- wilcox.test(cd8$D1D2, cd8$D3, paired = TRUE,alternative = c("less"))
wilcox_test_D1D2_D4 <- wilcox.test(cd8$D1D2, cd8$D4, paired = TRUE,alternative = c("less"))
wilcox_test_D3_D4 <- wilcox.test(cd8$D3, cd8$D4, paired = TRUE,alternative = c("less"))

# Print test results
print(wilcox_test_D1D2_D3)
print(wilcox_test_D1D2_D4)
print(wilcox_test_D3_D4)

# Add significance bars
p_cd8  <- p_cd8 + 
  geom_signif(comparisons = list(c("D1D2", "D3")), 
              annotations = paste0("p = ", signif(wilcox_test_D1D2_D3$p.value, digits = 3)), 
              y_position = max(cd8_long$Value) * 1.05) +
  geom_signif(comparisons = list(c("D1D2", "D4")), 
              annotations = paste0("p = ", signif(wilcox_test_D1D2_D4$p.value, digits = 3)), 
              y_position = max(cd8_long$Value) * 1.10) +
  geom_signif(comparisons = list(c("D3", "D4")), 
              annotations = paste0("p = ", signif(wilcox_test_D3_D4$p.value, digits = 3)), 
              y_position = max(cd8_long$Value) * 1.18)

ggsave(paste0(dir,'cov_enrich_CD8.pdf'))


cd4 <- read.table(paste0(dir,'cov_enrich_CD4.txt'),header=T)

print(cd4)

# Reshape the data frame to long format
cd4_long <- cd4 %>%
  pivot_longer(cols = c(D1D2, D3, D4), names_to = "Time", values_to = "Value")

# Convert Time to factor with specific order
cd4_long$Time <- factor(cd4_long$Time, levels = c("D1D2", "D3", "D4"))

# Create the line plot
p_cd4 <- ggplot(cd4_long, aes(x = Time, y = Value, group = Sample)) +
  geom_line(aes(color = Sample)) +
  geom_point(aes(color = Sample)) +
  theme_minimal() +
  labs(title = "Expanding spike-specific CD4+ T-cells", x = "Time Point", y = "Unique CD4+ T Cells")


# Perform the paired non-parametric statistical test (Wilcoxon signed-rank test)
# Compare D1D2 vs D3, D1D2 vs D4, and D3 vs D4
wilcox_test_D1D2_D3_cd4 <- wilcox.test(cd4$D1D2, cd4$D3, paired = TRUE,alternative = c("less"))
wilcox_test_D1D2_D4_cd4 <- wilcox.test(cd4$D1D2, cd4$D4, paired = TRUE,alternative = c("less"))
wilcox_test_D3_D4_cd4 <- wilcox.test(cd4$D3, cd4$D4, paired = TRUE,alternative = c("less"))

# Add significance bars
p_cd4 <- p_cd4 + 
  geom_signif(comparisons = list(c("D1D2", "D3")),
              annotations = paste0("p = ", signif(wilcox_test_D1D2_D3_cd4$p.value, digits = 3)), 
              y_position = max(cd4_long$Value) * 1.05) +
  geom_signif(comparisons = list(c("D1D2", "D4")), 
              annotations = paste0("p = ", signif(wilcox_test_D1D2_D4_cd4$p.value, digits = 3)), 
              y_position = max(cd4_long$Value) * 1.10) +
  geom_signif(comparisons = list(c("D3", "D4")), 
              annotations = paste0("p = ", signif(wilcox_test_D3_D4_cd4$p.value, digits = 3)), 
              y_position = max(cd4_long$Value) * 1.18)

ggsave(paste0(dir,'cov_enrich_CD4.pdf'))

p_fig4 <- p_cd8 + p_cd4 + plot_annotation(tag_levels = "A") 
ggsave('../results/figure4.pdf',p_fig4, width = 28, height = 13,units = 'cm')
