# from pbmc

library(Seurat)
library(dplyr)
library(ggplot2)

x <- readRDS(snakemake@input[[1]])

seurat_obj_smol <- subset(x, idents = c("3","6"), invert = T)

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

pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10


DoHeatmap(seurat_obj_smol, features = top10$gene)
ggsave(snakemake@output[[1]], width = 8.75, height = 10)


