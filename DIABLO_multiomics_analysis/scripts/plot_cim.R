library("mixOmics")
library("ComplexHeatmap")

fulldf <- read.csv("data/phenodata_all.csv")
rownames(fulldf) <- fulldf$patient
final.diablo.model <- readRDS("results/diagnostics/full_model_all_samples_tp1.rds")

p <- cimDiablo(final.diablo.model, comp = 1, 
               color.Y = c("red", "blue"),
               color.blocks = c("yellow", "orange"))

mat <- p$mat
saveRDS(mat, "results/cim.rds")

#str(final.diablo.model)

ones <- selectVar(final.diablo.model, block = "blood", comp = 1)
blood <- ones[["blood"]]$name
write(blood, "results/test_blood.txt")

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




png("results/heatmap.png", width = 8.96, height = 3, res = 600, units = "in")
Heatmap(mat, name = "Correlation", 
        right_annotation = annotation_row, row_title = NULL, column_title = NULL,
        top_annotation = annotation_col, row_km = 2, column_km = 2,
        row_names_gp = grid::gpar(fontsize = 2), show_parent_dend_line = FALSE,
        column_names_gp = grid::gpar(fontsize = 2), border = TRUE, 
        show_column_names = F, show_row_names = F, show_column_dend = T, show_row_dend = T)
dev.off()

png("results/blood_tree.png", width = 5, height = 7.5, units = "in", res = 600)
plotLoadings(final.diablo.model, 
             block = "blood", 
             comp = 1, contrib = 'max', 
             method = 'median', ndisplay = 50, 
             title = "", legend.color = c("red", "blue"), legend.title = "Repsonse")
dev.off()

png("results/tumour_tree.png", width = 5, height = 7.5, units = "in", res = 600)
plotLoadings(final.diablo.model, 
             block = "tumour", ndisplay = 50,
             comp = 1, 
             contrib = 'max', 
             method = 'median', 
             title = "", legend.color = c("red", "blue"), legend.title = "Repsonse")
dev.off()

