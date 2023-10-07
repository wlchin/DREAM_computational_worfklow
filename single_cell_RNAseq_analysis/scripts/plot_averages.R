

### this is the new figure 4

mat <- read.csv("results/average_expression_individual_genes_DIABLO.csv", row.names = 1)

library(ComplexHeatmap)

scaled_mat = t(scale(t(mat)))

#scaled_mat[,is.na(scaled_mat)]
#Heatmap(na.omit(scaled_mat), row_names_gp = gpar(fontsize = 5))

a <- Heatmap(na.omit(mat), 
        row_names_gp = gpar(fontsize = 4), 
        column_names_gp = gpar(fontsize = 4), 
        rect_gp = gpar(col = "black", lwd = 1), show_column_dend = F, show_row_dend = F)
png("results/diablo_sig_individual_genes.png", width = 2.24, height = 7.69, res = 600, units = "in")
draw(a)
dev.off()

mat2 <- read.csv("results/GSEA_DIABLO.csv")
matgsea <- mat2[-c(1,2),-1]
rownames(matgsea) <- mat2$source[-c(1,2)]
colnames(matgsea) <- c("Timepoint0", "T1", "T2")

a <- Heatmap(na.omit(matgsea), 
             row_names_gp = gpar(fontsize = 14), 
             column_names_gp = gpar(fontsize = 14), 
             rect_gp = gpar(col = "black", lwd = 1), show_column_dend = F, show_row_dend = F, cluster_columns = F, row_names_side = "left")
png("results/diablo_GSEA.png", width = 3.10, height = 3.91, res = 600, units = "in")
draw(a)
dev.off()

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
