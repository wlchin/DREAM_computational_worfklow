
library("Seurat")
library("flextable")
library("ComplexHeatmap")
library("tidyr")
library("dplyr")


large_seurat <- readRDS("data/CD8Tcells_seurat.rds")

seurat_obj_filt <- subset(large_seurat, idents = c("3","6"), invert = T)

seurat_obj_filt <- RenameIdents(object = seurat_obj_filt, 
                                `0` = "CD8_TEM_1",
                                `1` = "CD8_TEM_2",
                                `2` = "CD8_TEM_3",
                                `4` = "CD8_TEM_4",
                                `5` = "CD8_TEM_5",
                                `7` = "CD8_TEM_6",
                                `8` = "CD8_TEM_7"
)

seurat_obj_filt[['clusters_renamed']] <- Idents(seurat_obj_filt)

df <- seurat_obj_filt@meta.data
tot <- df %>% group_by(clusters_renamed, Patient) %>% count()
totmat <- spread(tot, Patient, n)
totmat[is.na(totmat)] <- 0
mat <- round(totmat[,-1])
rownames(mat) <- totmat$clusters_renamed



pheno <- df[,c("Patient", "PFS_6M")]

phenotype <- pheno[!duplicated(pheno),] %>% arrange(Patient)
phenotype$PFS_6M[phenotype$PFS_6M == 0] <- "Non-responder" 
phenotype$PFS_6M[phenotype$PFS_6M == 1] <- "Responder"

phenotype$PFS_6M <- as.factor(phenotype$PFS_6M)

ha <- HeatmapAnnotation(Phenotype = phenotype$PFS_6M, col = list(Phenotype = c("Responder" = "red", "Non-responder" = "blue")))

#, col = list(Response = c("Responder" = "Red", "Non-Responder" = "Blue")

ht <- Heatmap(mat, top_annotation = ha, cluster_rows = F, cluster_columns = F, name = "Cell \ncounts", cell_fun = function(j, i, x, y, width, height, fill){
  grid.text(mat[i, j], x, y, gp = gpar(fontsize = 5))}, 
  show_column_names = F, column_title = "Patients", 
  row_title = "Clusters")

png("results/composition_heatmap_cd8_TEM.png", width = 10.71, height = 2.03, units = "in", res = 600)
draw(ht)
dev.off()
