
library(ComplexHeatmap)
library(RColorBrewer)
library(viridisLite)


x1 <- read.csv("results/composition_mats/CD8_TEM.csv")
x2 <- read.csv("results/composition_mats/CD4_TEM.csv")
x3 <- read.csv("results/composition_mats/CD14_Mono.csv")
x4 <- read.csv("results/composition_mats/NK.csv")
x5 <- read.csv("results/composition_mats/CD4_Naive.csv")
x6 <- read.csv("results/composition_mats/Treg.csv")
x7 <- read.csv("results/composition_mats/CD16_Mono.csv")
x8 <- read.csv("results/composition_mats/B_naive.csv")
x9 <- read.csv("results/composition_mats/CD4_TCM.csv")
x10 <- read.csv("results/composition_mats/CD8_Naive.csv")

colsplit <- c(rep("CD8 TEM", 2),
  rep("CD4 TEM", 2),
  rep("CD14 Mono", 2),
  rep("NK", 2),
  rep("CD4 Naive", 2),
  rep("Treg", 2),
  rep("CD16 Mono", 2),
  rep("B Naive", 2),
  rep("CD4 TCM", 2),
  rep("CD8 Naive", 2))

filenames_list <- list.files("./results/composition_mats", ".csv", full.names = T)


condense_counts <- function(mat){
  x <- mat
  x$pheno[x$pheno == 0] <- "Patients (NR)"
  x$pheno[x$pheno == 1] <- "Patients (R)"
  counts <- x[,-c(1, 12)]
  rownames(counts) <- x$X
  counts_ <- cbind(rowSums(counts[,1:5]), rowSums(counts[,6:10]))
  counts_
}

count_list <- list(condense_counts(x1),
     condense_counts(x2),
     condense_counts(x3),
     condense_counts(x4),
     condense_counts(x5),
     condense_counts(x6),
     condense_counts(x7),
     condense_counts(x8),
     condense_counts(x9),
     condense_counts(x10))



mat <- as.matrix(do.call(cbind, count_list))

column_ha = HeatmapAnnotation(`DA neighbourhoods` = rep(c("Abundant in NR", "Abundant in R"), 10), height = unit(1, "mm"), col = list(`Abundant in R` = "red", `Abundant in NR` = "blue"),  gp = gpar(col = "black"), annotation_legend_param = list(border = "black"))

ht <- Heatmap(mat, cluster_rows = F, name = "Cells in neighbourhoods",
        top_annotation = column_ha, 
        cluster_columns = F, row_split = x1$pheno, 
        column_split = colsplit, column_title_gp = grid::gpar(fontsize = 10), 
        show_row_names = F, column_title_rot = 90, column_title_side = "bottom", border = TRUE, 
        width = ncol(mat)*unit(2, "mm"), 
        height = nrow(mat)*unit(2, "mm"))


png("results/summarised_composition.png", width = 5.63, height = 5.76, res = 600, units = "in")
draw(ht)
dev.off()
