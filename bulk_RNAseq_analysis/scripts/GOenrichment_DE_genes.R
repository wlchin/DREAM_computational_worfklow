library(clusterProfiler)
library(enrichplot)
library(pheatmap)
library(org.Hs.eg.db)
library(ggplot2)
library(cowplot)


nrB_C2 <- read.csv("DE_results/Base_C2_nr.csv")
rB_C2 <- read.csv("DE_results/Base_C2_r.csv")
rC2_C3 <- read.csv("DE_results/C2_C3_r.csv")

parse_geneset <- function(x, FC, FDR, direction, phenogroup){
  if(direction == "up"){
    data <- x[(x$logFC > FC) & (x$adj.P.Val < FDR),]
  }
  if(direction == "down"){
    data <- x[(x$logFC < -FC) & (x$adj.P.Val < FDR),]
  }
  data$pheno <- phenogroup
  data$direction <- direction
  data
}

NR_BC2_up = parse_geneset(nrB_C2, 0.5, 0.05, "up" , "NR_BC2")
NR_BC2_down = parse_geneset(nrB_C2, 0.5, 0.05, "down" , "NR_BC2")
R_BC2_up = parse_geneset(rB_C2, 0.58, 0.05, "up" , "R_BC2")
R_BC2_down = parse_geneset(rB_C2, 0.58, 0.05, "down" , "R_BC2")
#RC2_C3_up = parse_geneset(rC2_C3, 0.58, 0.05, "up" , "R_C2C3")
#RC2_C3_down = parse_geneset(rC2_C3, 0.5, 0.1, "down" , "R_C2C3") # zero genes

dfs <- list(NR_BC2_up, NR_BC2_down, R_BC2_up, R_BC2_down)#, RC2_C3_up)
whole_df <- do.call(rbind, dfs)

formula_res <- compareCluster(X~pheno+direction, 
                              data=whole_df, 
                              fun="enrichGO", 
                              OrgDb="org.Hs.eg.db", 
                              ont = "BP", universe = rC2_C3$X)

formularessimp <- simplify(formula_res)
formularessimp <- formula_res

## Visualisation with dotplots
dplot <- dotplot(formularessimp, font.size = 10, size = "count", showCategory = 3)
dplot + scale_x_discrete(labels = c("Timepoint 0 NR vs \nTimepoint 1 NR \n(downregulated)", "Timepoint 0 NR vs \nTimepoint 1 NR \n(upregulated)", "Timepoint 0 R vs \nTimepoint 1 R \n(downregulated)", "Timepoint 0 R vs \nTimepoint 1 R \n(upregulated)"))
ggsave("results/dotplot.png", width = 8, height = 5.29, dpi = 600)

## Visualisation with barplots

resdf <- formularessimp@compareClusterResult
top5Rup <- head(resdf[resdf$Cluster == "R_BC2.up",], 5)
top5Rdown <- head(resdf[resdf$Cluster == "R_BC2.down",], 5)

top5Rup$log10p <- -log10(top5Rup$p.adjust)
top5Rdown$log10p <- -log10(top5Rdown$p.adjust)

top5Rup$Description[4] <- "immune response mediated by antimicrobial peptide"

a <- ggplot(top5Rup, aes(x = reorder(Description, log10p), log10p)) + 
  geom_col(fill = "#4dbbd5", color = "black") + 
  xlab("") + 
  ylab("log10p") + ylim(0, 8.0) + 
  coord_flip() + 
  theme_cowplot() + theme(axis.text.y = element_text(size = 20), axis.text.x = element_text(size = 20))

b <- ggplot(top5Rdown, aes(x = reorder(Description, log10p), log10p), ) + 
  geom_col(fill = "#e64b35", color = "black") + 
  xlab("") + 
  ylab("log10p") + 
  coord_flip() + 
  theme_cowplot() + theme(axis.text.y = element_text(size = 20), axis.text.x = element_text(size = 20))

plot_grid(a, b, ncol = 1, align = "v")
ggsave("results/barchart.png", width = 11, height = 6.3, dpi = 600)
