
library(clusterProfiler)
library(enrichplot)
library(pheatmap)
library(org.Hs.eg.db)
library(ggplot2)
library(stringr)
library(cowplot)

top <- read.csv("DE_results/Base_C2_interaction.csv")
top$FeatureID <- rownames(top)
top  <- top %>% dplyr::arrange(desc(logFC))

query <- top$logFC
names(query) <- top$X
geneList2 <- sort(query, decreasing = T)

ego <- gseGO(geneList     = geneList2,
             OrgDb        = org.Hs.eg.db,
             ont          = "BP",
             minGSSize    = 100,
             maxGSSize    = 500,
             pvalueCutoff = 0.05,
             nPermSimple = 100000,
             seed = 42,
             eps = 0,
             verbose      = TRUE)


gseaplot2(ego, geneSetID = c("GO:0050852"), pvalue_table = T)
ggsave("results/GSEAplot_antigen.png", width = 7.96)

gseaplot2(ego, geneSetID = c("GO:0050851"), pvalue_table = T)
ggsave("results/GSEAplot_TCR.png", width = 7.96)

