
library(ggplot2)
library(cowplot)

x <- read.csv("results/CD8TEM_IFNgamma.csv")

x$PFS_6M[x$PFS_6M == 0] <- "NR"
x$PFS_6M[x$PFS_6M == 1] <- "R"

x$PFS_6M <- as.factor(x$PFS_6M)
x$Timepoint <- paste0("Timepoint ", x$Timepoint)

ggplot(x, aes(PFS_6M, HALLMARK_INTERFERON_GAMMA_RESPONSE)) + ylab("IFN-\u03B3 response scores") + xlab("Phenotype") + ylim(0, 17.5) + 
  geom_boxplot(notch = T, varwidth = F, outlier.shape = NA) + theme_cowplot() + 
  facet_grid(~Timepoint)
ggsave("results/ifn_gamma_response.png", width = 3.92, height = 3.92, dpi = 600)
