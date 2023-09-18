# This is to allocate seurat clusters for scanpy

library(Seurat)

CD8 <- readRDS("data/CD8Tcells_seurat.rds")

clusterdf <- data.frame(CD8@meta.data[,c("seurat_clusters")])
rownames(clusterdf) <- rownames(CD8@meta.data)
colnames(clusterdf) <- c("cluster_seurat")

write.csv(clusterdf, "data/CD8_seurat_clusters.csv")
