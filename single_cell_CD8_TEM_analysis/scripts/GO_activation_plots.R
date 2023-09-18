library("ComplexHeatmap")

x <- read.csv("results/res_cluster0.csv")
mat <- x[,-1]
rownames(mat) <- x$X
mat1 <- mat[,c("CD81", "HLA.A", "HLA.DRB1", "CCL5", "CD6", "EOMES", "HLA.E", "AKT1", "HLA.DQA1", "IRF1", "RHOH", "RUNX3", "CD8B", "IFNG", "CD5" )]
rownames(mat1) <- c("Timepoint0_NR", 
                    "Timepoint1_NR", 
                    "Timepoint2_NR", 
                    "Timepoint0_R", 
                    "Timepoint1_R", 
                    "Timepoint2_R")

png("results/markers_acitvation_cluster0.png", width = 4.47, height = 1.95, units = "in", res = 600)
Heatmap(mat1, row_split = c("C", "B", "A", "C", "B", "A"), show_row_dend = F, show_column_dend = F, cluster_rows = F, cluster_columns = F)
dev.off()

x <- read.csv("results/res_cluster1.csv")
mat <- x[,-1]
rownames(mat) <- x$X
mat1 <- mat[,c("CD81", "HLA.A", "HLA.DRB1", "CCL5", "CD6", "EOMES", "HLA.E", "AKT1", "HLA.DQA1", "IRF1", "RHOH", "RUNX3", "CD8B", "IFNG", "CD5" )]
rownames(mat1) <- c("Timepoint0_NR", 
                    "Timepoint1_NR", 
                    "Timepoint2_NR", 
                    "Timepoint0_R", 
                    "Timepoint1_R", 
                    "Timepoint2_R")

png("results/markers_acitvation_cluster1.png", width = 4.47, height = 1.95, units = "in", res = 600)
Heatmap(mat1, row_split = c("C", "B", "A", "C", "B", "A"), show_row_dend = F, show_column_dend = F, cluster_rows = F, cluster_columns = F)
dev.off()

x <- read.csv("results/res_cluster2.csv")
mat <- x[,-1]
rownames(mat) <- x$X
mat1 <- mat[,c("CD81", "HLA.A", "HLA.DRB1", "CCL5", "CD6", "EOMES", "HLA.E", "AKT1", "HLA.DQA1", "IRF1", "RHOH", "RUNX3", "CD8B", "IFNG", "CD5" )]
rownames(mat1) <- c("Timepoint0_NR", 
                    "Timepoint1_NR", 
                    "Timepoint2_NR", 
                    "Timepoint0_R", 
                    "Timepoint1_R", 
                    "Timepoint2_R")

png("results/markers_acitvation_cluster2.png", width = 4.47, height = 1.95, units = "in", res = 600)
Heatmap(mat1, row_split = c("C", "B", "A", "C", "B", "A"), show_row_dend = F, show_column_dend = F, cluster_rows = F, cluster_columns = F)
dev.off()

