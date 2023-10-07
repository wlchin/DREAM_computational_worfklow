
library(ComplexHeatmap)

mat <- read.csv("results/GSEA_DIABLO_CD8.csv", row.names = 1, header = T)
mat <- mat[3:5,] 
colnames(mat) <- c("T0", "T1", "T2")
rownames(mat) <- c("CD8_TEM_0", "CD8_TEM_1", "CD8_TEM_2")

a <- Heatmap(na.omit(mat), 
             row_names_gp = gpar(fontsize = 14), 
             column_names_gp = gpar(fontsize = 14), 
             rect_gp = gpar(col = "black", lwd = 1), cluster_rows = F, cluster_columns = F, row_names_side = "left")

png("results/diablo_sig_CD8_TEM.png", width = 3.41, height = 1.96, res = 600, units = "in")
draw(a)
dev.off()
