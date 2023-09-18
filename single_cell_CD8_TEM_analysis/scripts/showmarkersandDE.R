# from pbmc

library(Seurat)
library(dplyr)
library(ggplot2)

x <- readRDS(snakemake@input[[1]])

#recently_activated <- c("KLRG1", "EOMES", "TBX21", "CX3CR1", "FGFBP2", "TCF7", "CXCR5")
#terminally_exhausted <- c("GZMB", "NKG7", "PRF1", "S1PR5")

#memory <- c("TCF7", "IL7R", "CXCR3")
#effector <- c("KLRG1", "GZMA", "CX3CR1")

seurat_obj_smol <- subset(x, idents = c("3","6"), invert = T)
#FeaturePlot(pbmc, recently_activated)
#FeaturePlot(pbmc, c(recently_activated, terminally_exhausted))
#FeaturePlot(pbmc, others)
#ggsave(snakemake@output[[1]], width = 5.67, height = 6.76)


seurat_obj_smol <- RenameIdents(object = seurat_obj_smol, 
                                `0` = "CD8_TEM_1",
                                `1` = "CD8_TEM_2",
                                `2` = "CD8_TEM_3",
                                `4` = "CD8_TEM_4",
                                `5` = "CD8_TEM_5",
                                `7` = "CD8_TEM_6",
                                `8` = "CD8_TEM_7"
)

pbmc.markers <- FindAllMarkers(seurat_obj_smol, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

#pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

#pbmc.markers %>%
#  group_by(cluster) %>%
#  slice_max(n = 2, order_by = avg_log2FC)

pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10


DoHeatmap(seurat_obj_smol, features = top10$gene)
ggsave(snakemake@output[[1]], width = 8.75, height = 10)


