

library("Seurat")
library("ggsci")
library("ggplot2")
library("gridExtra")
library('cowplot')
library("dplyr")

seurat_obj <- readRDS(snakemake@input[[1]])

seurat_obj_smol <- subset(seurat_obj, idents = c("3","6"), invert = T)

seurat_obj_smol <- RenameIdents(object = seurat_obj_smol, 
                     `0` = "CD8_TEM_1",
                     `1` = "CD8_TEM_2",
                     `2` = "CD8_TEM_3",
                     `4` = "CD8_TEM_4",
                     `5` = "CD8_TEM_5",
                     `7` = "CD8_TEM_6",
                     `8` = "CD8_TEM_7"
                     )

DimPlot(seurat_obj_smol, label = T) + NoLegend()

seurat_obj_smol$new_labels <- Idents(seurat_obj_smol)

DimPlot(seurat_obj_smol)
ggsave("results/UMAP.png")

df <- seurat_obj_smol@meta.data

df$PFS_6M[df$PFS_6M == 0] <- "NR"
df$PFS_6M[df$PFS_6M == 1] <- "R"

df$Timepoint[df$Timepoint == 0] <- "Timepoint0"
df$Timepoint[df$Timepoint == 1] <- "Timepoint1"
df$Timepoint[df$Timepoint == 2] <- "Timepoint2"

df$phenotype <- paste0(df$Timepoint, "_", df$PFS_6M)

tally_df <- group_by(df, phenotype, new_labels) %>% count()

ggplot(tally_df) + geom_col(aes(new_labels, n, fill = phenotype)) + theme_cowplot() + ylab("number of cells") + xlab("clusters") +
  theme(axis.text.x = element_text(size = 12, angle = 90), axis.text.y = element_text(size = 12), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_fill_discrete(limits = c("Timepoint0_NR", "Timepoint1_NR", "Timepoint2_NR", "Timepoint0_R", "Timepoint1_R", "Timepoint2_R"), name = "Condition") +
  scale_fill_npg()
ggsave("results/barplot_props.png")
