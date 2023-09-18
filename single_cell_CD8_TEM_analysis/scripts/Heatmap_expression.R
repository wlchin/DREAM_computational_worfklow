library("ComplexHeatmap")

x <- read.csv("results/zscaled_proprotions.csv")

mat <- x[,-1]
rownames(mat) <- x$X

mat1 <- mat[1:3, 1:36]

split = c(rep("Exhaustion", 5), rep("Naive-like", 4), rep("Effector", 6), rep("Naive-like", 5), rep("Term-Ex", 8), rep("KIR+", 8))

#split = c(rep("Exhaustion", 5), rep("Naive-like", 4), rep("Effector", 6), rep("Prog-Ex", 5), rep("Term-Ex", 8), )
#split = c(rep("Exhaustion", 5), rep("Naive-like", 4))

ht <- Heatmap(mat1, heatmap_width = unit(21, "cm"), 
        heatmap_height = unit(4, "cm"), 
        column_names_gp = grid::gpar(fontsize = 8),
        row_names_gp = grid::gpar(fontsize = 8),
        rect_gp = gpar(col = "white", lwd = 2),
        column_title_gp = gpar(fontsize = 12), 
        cluster_rows = F, cluster_columns = F, column_split = split, 
        row_names_side = "left", row_title = "",
        heatmap_legend_param = list(
          title = "",
          legend_height = unit(2.2, "cm"))
        )

png("results/markers_CD8Tcells.png", width = 8.74, height = 3.79, units = "in", res = 600)
draw(ht)
dev.off()

#pdf(snakemake@output[[1]])
#draw(ht)
#dev.off()