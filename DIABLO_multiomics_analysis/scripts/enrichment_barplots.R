
library(clusterProfiler)
library(dplyr)
library(enrichplot)
library(DOSE)
library(ggplot2)
library(patchwork)
library(cowplot)

x <- readRDS("results/tumour.rds")


top5 <- head(x, 5)
top5$score <- -log10(top5$pvalue) 
top5 <- arrange(top5, desc(score))

tumour <- ggplot(top5, aes(score, y = reorder(Description, score))) + 
  geom_col(fill = "#FFA500", color = "black") + labs(x = "-log10pval", y = "go terms") + ggtitle("Enrichment: tumour") + theme_cowplot() + ylab("")


x <- readRDS("results/blood.rds")


top5 <- head(x, 5)
top5$score <- -log10(top5$pvalue) 
top5 <- arrange(top5, desc(score))

blood <- ggplot(top5, aes(score, y = reorder(Description, score))) + 
  geom_col(fill = "#a020f0", color = "black") + labs(x = "-log10pval", y = "go terms") + ggtitle("Enrichment: blood") + theme_cowplot() + ylab("")

blood
ggsave("results/barplots_blood.png", width = 7.10, height = 2.75, dpi = 600)

tumour
ggsave("results/barplots_tumour.png", width = 7.10, height = 2.75, dpi = 600)
