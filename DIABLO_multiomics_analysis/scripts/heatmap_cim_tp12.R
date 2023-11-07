library("mixOmics")
library("ComplexHeatmap")

fulldf <- read.csv("data/phenodata_all.csv")
rownames(fulldf) <- fulldf$patient
final.diablo.model <- readRDS("results/diagnostics/full_model_all_samples_tp12.rds")

p <- cimDiablo(final.diablo.model, comp = 1, 
               color.Y = c("red", "blue"),
               color.blocks = c("yellow", "orange", "pink"))


mat <- p$mat
saveRDS(mat, "results/cim.rds")

mat <- readRDS("results/cim.rds")

#str(final.diablo.model)

ones <- selectVar(final.diablo.model, block = "blood", comp = 1)
blood <- ones[["blood"]]$name
write(blood, "results/test_blood.txt")

ones <- selectVar(final.diablo.model, block = "blood2", comp = 1)
blood <- ones[["blood2"]]$name
write(blood, "results/test_blood2.txt")

ones <- selectVar(final.diablo.model, block = "tumour")
tumour <- ones[["tumour"]]$name
write(tumour, "results/test_tumour.txt")

#rownames(fulldf) <- fulldf$PT_mod
histo <- fulldf[rownames(p$mat),]$histology
bestresp <- fulldf[rownames(p$mat),'ITMBESTCONRESP']
os <- fulldf[rownames(p$mat),'PFS_6M']

anno_ <- ifelse(p$col.sideColors == "orange", "Tumour", "Blood")
annotation_col = HeatmapAnnotation(assay = anno_,
                                   col = list(assay = c("Tumour" = "Orange", "Blood" = "Purple")),
                                   border = T, col_title = NULL, height = unit(3, "mm"))


anno_ <- ifelse(p$row.sideColors == "red", "Responder", "Non-Responder")
annotation_row = rowAnnotation(Response = anno_, 
                               col = list(Response = c("Responder" = "Red", "Non-Responder" = "Blue")),
                               border = T, row_title = NULL, show_annotation_name = F)




png("results/heatmap.png", width = 9.29, height = 3, res = 600, units = "in")
Heatmap(mat, name = "Correlation", 
        right_annotation = annotation_row, row_title = NULL, column_title = NULL,
        top_annotation = annotation_col, row_km = 2, column_km = 2,
        row_names_gp = grid::gpar(fontsize = 2), show_parent_dend_line = FALSE,
        column_names_gp = grid::gpar(fontsize = 2), border = TRUE, clustering_method_rows = "single", clustering_method_columns = "single",
        show_column_names = F, show_row_names = F, show_column_dend = T, show_row_dend = T)
dev.off()



