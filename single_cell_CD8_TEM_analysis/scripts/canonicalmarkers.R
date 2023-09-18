
library(Seurat)
library(ggplot2)
library("clusterProfiler")
library(RColorBrewer)
library(enrichplot)

CD8 <- readRDS("data/CD8Tcells_seurat.rds")
CD8smol <- subset(CD8, idents = c(0,1,2))

res <- FindMarkers(CD8smol, ident.1 = "1", min.pct = 0.1, logfc.threshold = 0.001)
saveRDS(res, "results/canonical_cluster1.rds")
res <- FindMarkers(CD8smol, ident.1 = "2", min.pct = 0.1, logfc.threshold = 0.001)
saveRDS(res, "results/canonical_cluster2.rds")
res <- FindMarkers(CD8smol, ident.1 = "0", min.pct = 0.1, logfc.threshold = 0.001)
saveRDS(res, "results/canonical_cluster0.rds")

stem_all <- readRDS("data/stem_all_gene_name.rds")
stem_all <- stem_all[,c(2,1)]
stem_all <- stem_all[(stem_all$gene.set.name == "PACE") | (stem_all$gene.set.name == "Miller_prog_ex") | (stem_all$gene.set.name == "zhengTCF7PEX") | (stem_all$gene.set.name == "Miller_term_ex"),]

res <- readRDS("results/canonical_cluster1.rds")
ranks <- res$avg_log2FC
names(ranks) <- rownames(res)
ranks_sorted <- sort(ranks, decreasing = TRUE)
y <- GSEA(ranks_sorted, TERM2GENE = stem_all, pvalueCutoff = 1)
gseaplot2(y, "PACE", subplots = 1:2, color = "black", pvalue_table = T) 
ggsave("results/PACE_enrichment_CD8TEM2.png")
cluster1 <- y@result
cluster1$cluster = "CD8_TEM_2"

res <- readRDS("results/canonical_cluster0.rds")
ranks <- res$avg_log2FC
names(ranks) <- rownames(res)
ranks_sorted <- sort(ranks, decreasing = TRUE)
y <- GSEA(ranks_sorted, TERM2GENE = stem_all, pvalueCutoff = 1)
cluster0 <- y@result
cluster0$cluster = "CD8_TEM_1"

res <- readRDS("results/canonical_cluster2.rds")
ranks <- res$avg_log2FC
names(ranks) <- rownames(res)
ranks_sorted <- sort(ranks, decreasing = TRUE)
y <- GSEA(ranks_sorted, TERM2GENE = stem_all, pvalueCutoff = 1)
cluster2 <- y@result
cluster2$cluster = "CD8_TEM_3"

testdf4 <- rbind(cluster0, cluster1, cluster2)
testdf4$log10p <- -log10(testdf4$p.adjust)

smol <- testdf4[testdf4$ID %in% c("PACE", "Miller_prog_ex","zhengTCF7PEX", "Miller_term_ex"),]

hec_sp <- ggplot(smol, aes(x = cluster, y = ID)) + xlab("") + ylab("") + 
  geom_point(aes(size = log10p, fill = NES), shape = 21, colour = "black") + 
  scale_size_binned_area(max_size = 10, guide = "legend") + theme(axis.text.y = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust=1)) +
  theme(strip.background =element_rect(fill="black"), panel.border = element_rect(colour = "black", fill=NA, linewidth =1), panel.background = element_blank()) +
  scale_fill_distiller(palette ="RdBu", direction = -1)

hec_sp
ggsave("results/markers_per_cluster.png", width = 4.05, height = 6.33, dpi = 600)



