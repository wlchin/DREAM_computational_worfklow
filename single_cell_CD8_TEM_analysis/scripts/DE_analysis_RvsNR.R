
library(Seurat)
library(cowplot)

# plan is do DE on each cluster, then do clustrprofiler

CD8_TEMs <- readRDS("data/CD8Tcells_seurat.rds")

#table(CD8_TEMs@meta.data$PFS_6M, CD8_TEMs@meta.data$seurat_clusters) # there are enough except in cluster 3

CD8_TEMs$cluster_response <- paste0(CD8_TEMs$seurat_clusters, "_", CD8_TEMs$PFS_6M, "_", CD8_TEMs$Timepoint)

#unique(CD8_TEMs$cluster_response)

Idents(CD8_TEMs) <- "cluster_response"

#for(i in unique(CD8_TEMs$seurat_clusters)){
#    NR = paste0(i, "_", 0)
#    R = paste0(i, "_", 1)
#    cluster.markers <- FindMarkers(CD8_TEMs, ident.1 = NR , ident.2 = R, min.pct = 0.10, logfc.threshold = 0.1)
#    filename_csv <- paste0("cluster_", i, ".csv")
#    write.csv(cluster.markers, filename_csv)
#}


threshold = 0.001


cluster.markers0 <- FindMarkers(CD8_TEMs, ident.1 = "2_1_0" , ident.2 = "2_0_0", min.pct = 0.10, logfc.threshold = threshold)
saveRDS(cluster.markers0, "results/cluster2_timepoint0.rds")
cluster.markers1 <- FindMarkers(CD8_TEMs, ident.1 = "2_1_1" , ident.2 = "2_0_1", min.pct = 0.10, logfc.threshold = threshold)
saveRDS(cluster.markers1, "results/cluster2_timepoint1.rds")
cluster.markers2 <- FindMarkers(CD8_TEMs, ident.1 = "2_1_2" , ident.2 = "2_0_2", min.pct = 0.10, logfc.threshold = threshold)
saveRDS(cluster.markers2, "results/cluster2_timepoint2.rds")


cluster.markers0 <- FindMarkers(CD8_TEMs, ident.1 = "1_1_0" , ident.2 = "1_0_0", min.pct = 0.10, logfc.threshold = threshold)
saveRDS(cluster.markers0, "results/cluster1_timepoint0.rds")
cluster.markers1 <- FindMarkers(CD8_TEMs, ident.1 = "1_1_1" , ident.2 = "1_0_1", min.pct = 0.10, logfc.threshold = threshold)
saveRDS(cluster.markers1, "results/cluster1_timepoint1.rds")
cluster.markers2 <- FindMarkers(CD8_TEMs, ident.1 = "1_1_2" , ident.2 = "1_0_2", min.pct = 0.10, logfc.threshold = threshold)
saveRDS(cluster.markers2, "results/cluster1_timepoint2.rds")


cluster.markers0 <- FindMarkers(CD8_TEMs, ident.1 = "0_1_0" , ident.2 = "0_0_0", min.pct = 0.10, logfc.threshold = threshold)
saveRDS(cluster.markers0, "results/cluster0_timepoint0.rds")
cluster.markers1 <- FindMarkers(CD8_TEMs, ident.1 = "0_1_1" , ident.2 = "0_0_1", min.pct = 0.10, logfc.threshold = threshold)
saveRDS(cluster.markers1, "results/cluster0_timepoint1.rds")
cluster.markers2 <- FindMarkers(CD8_TEMs, ident.1 = "0_1_2" , ident.2 = "0_0_2", min.pct = 0.10, logfc.threshold = threshold)
saveRDS(cluster.markers2, "results/cluster0_timepoint2.rds")
