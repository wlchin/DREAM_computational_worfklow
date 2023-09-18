library(nebula)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)

x <- readRDS("data/CD8Tcells_seurat.rds")
x$PFS_6M[x$PFS_6M == 0] <- "NR"
x$PFS_6M[x$PFS_6M == 1] <- "R"

stem_all <- readRDS("data/stem_all_gene_name.rds")
stem_all <- stem_all[,c(2,1)]
stemsmol <- stem_all[stem_all$gene.set.name == "PACE" ,]

do_nebula <- function(h5obj) {
  h5obj@meta.data$cdr <- scale(h5obj@meta.data$nFeature_RNA)
  counts <- h5obj@assays$RNA@counts
  df <- model.matrix(~PFS_6M:Timepoint + cdr + Block, data=h5obj@meta.data)
  sid <- h5obj@meta.data$Patient
  res_re <- group_cell(counts, sid, pred = df, offset = NULL)
  re <- nebula(res_re$count,res_re$id,pred=res_re$pred, cpc = 0.005)
  re
}




cluster1 <- subset(x, idents = 0)
x_smol <- subset(cluster1, Timepoint == 0 | Timepoint == 1 | Timepoint == 2)
res <- do_nebula(x_smol)
df <- res$summary
ranks <- df$`logFC_PFS_6MR:Timepoint` #* -log10(x1$adj.P.Val_1 + 10e-6)
names(ranks) <- df$gene
ranks_sorted <- sort(ranks, decreasing = TRUE)
y <- GSEA(ranks_sorted, TERM2GENE = stemsmol, pvalueCutoff = 1, nPermSimple = 100000)
saveRDS(y@result, "results/cluster0_interaction.rds")
saveRDS(df, "results/cluster0_nebula_DE.rds")

cluster1 <- subset(x, idents = 1)
x_smol <- subset(cluster1, Timepoint == 0 | Timepoint == 1 | Timepoint == 2)
res <- do_nebula(x_smol)
df <- res$summary
ranks <- df$`logFC_PFS_6MR:Timepoint` #* -log10(x1$adj.P.Val_1 + 10e-6)
names(ranks) <- df$gene
ranks_sorted <- sort(ranks, decreasing = TRUE)
y <- GSEA(ranks_sorted, TERM2GENE = stemsmol, pvalueCutoff = 1, nPermSimple = 100000)
saveRDS(y@result, "results/cluster1_interaction.rds")
saveRDS(df, "results/cluster1_nebula_DE.rds")

gseaplot2(y, "PACE", pvalue_table = T)
ggsave("results/cluster1_PACE.png")

cluster1 <- subset(x, idents = 2)
x_smol <- subset(cluster1, Timepoint == 0 | Timepoint == 1 | Timepoint == 2)
res <- do_nebula(x_smol)
df <- res$summary
ranks <- df$`logFC_PFS_6MR:Timepoint` #* -log10(x1$adj.P.Val_1 + 10e-6)
names(ranks) <- df$gene
ranks_sorted <- sort(ranks, decreasing = TRUE)
y <- GSEA(ranks_sorted, TERM2GENE = stemsmol, pvalueCutoff = 1, nPermSimple = 100000)
saveRDS(y@result, "results/cluster2_interaction.rds")
saveRDS(df, "results/cluster2_nebula_DE.rds")
