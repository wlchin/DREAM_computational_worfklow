
library(clusterProfiler)
library(dplyr)
library(enrichplot)
library(DOSE)
library(ggplot2)
library(patchwork)
library(cowplot)

#x <- readRDS("results/tumour.rds")


#top5 <- head(x, 5)
#top5$score <- -log10(top5$pvalue) 
#top5 <- arrange(top5, desc(score))

#tumour <- ggplot(top5, aes(score, y = reorder(Description, score))) + 
#  geom_col(fill = "#FFA500", color = "black") + labs(x = "-log10pval", y = "go terms") + ggtitle("Enrichment: tumour") + theme_cowplot() + ylab("")


x25 <- readRDS("results/blood_25_genes.rds")
x50 <- readRDS("results/blood_50_genes.rds")
x100 <- readRDS("results/blood_100_genes.rds")

plot_diagrams <- function(x, fill, savefile){
  top5 <- head(x, 5)
  top5$score <- -log10(top5$pvalue) 
  top5 <- arrange(top5, desc(score))
  blood <- ggplot(top5, aes(score, y = reorder(Description, score))) + 
    geom_col(fill = fill, color = "black") + labs(x = "-log10pval", y = "go terms") + ggtitle("Enrichment: blood") + theme_cowplot() + ylab("")
  ggsave(savefile, blood, width = 7.10, height = 2.75, dpi = 600)
}

plot_diagrams(x25, "#a020f0", "results/barplots_blood_samp25.png")
plot_diagrams(x50, "#a020f0", "results/barplots_blood_samp50.png")
plot_diagrams(x100, "#a020f0", "results/barplots_blood_samp100.png")

x25 <- readRDS("results/tumour_25_genes.rds")
x50 <- readRDS("results/tumour_50_genes.rds")
x100 <- readRDS("results/tumour_100_genes.rds")

plot_diagrams(x25, "#FFA500", "results/barplots_tumour_samp25.png")
plot_diagrams(x50, "#FFA500", "results/barplots_tumour_samp50.png")
plot_diagrams(x100, "#FFA500", "results/barplots_tumour_samp100.png")

#tumour
#ggsave("results/barplots_tumour.png", width = 7.10, height = 2.75, dpi = 600)
