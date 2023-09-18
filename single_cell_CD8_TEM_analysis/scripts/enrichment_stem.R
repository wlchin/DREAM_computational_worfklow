library(clusterProfiler)
library(cowplot)
library(ggsci)
library(ggplot2)

stem_all <- readRDS("data/stem_all_gene_name.rds")
stem_all <- stem_all[,c(2,1)]
stem_all <- stem_all[(stem_all$gene.set.name == "PACE") | (stem_all$gene.set.name == "TSTEMhi") | (stem_all$gene.set.name == "TPEXhi") | (stem_all$gene.set.name == "GALLETTI_TSTEM") | (stem_all$gene.set.name == "Miller_prog_ex") | (stem_all$gene.set.name == "GALLETTI_TPEX"),]

make_df <- function(fileload, timepoint_lab, clusterlab){
  x <- readRDS(fileload)
  ranks <- x$avg_log2FC 
  names(ranks) <- rownames(x)
  ranks_sorted <- sort(ranks, decreasing = TRUE)
  y <- GSEA(ranks_sorted, TERM2GENE = stem_all, pvalueCutoff = 1)
  df <- y@result[,c("NES", "pvalue", "p.adjust")]
  df$timepoint <- timepoint_lab
  df$cluster <- clusterlab
  df$signature <- rownames(df)
  df$log10p <- -log10(df$pvalue)
  df
}

tpt2_clust2 <- make_df("results/cluster2_timepoint2.rds", "Timepoint2", "Cluster2")
tpt1_clust2 <- make_df("results/cluster2_timepoint1.rds", "Timepoint1", "Cluster2")
tpt0_clust2 <- make_df("results/cluster2_timepoint0.rds", "Timepoint0", "Cluster2")

tpt2_clust1 <- make_df("results/cluster1_timepoint2.rds", "Timepoint2", "Cluster1")
tpt1_clust1 <- make_df("results/cluster1_timepoint1.rds", "Timepoint1", "Cluster1")
tpt0_clust1 <- make_df("results/cluster1_timepoint0.rds", "Timepoint0", "Cluster1")

tpt2_clust0 <- make_df("results/cluster0_timepoint2.rds", "Timepoint2", "Cluster0")
tpt1_clust0 <- make_df("results/cluster0_timepoint1.rds", "Timepoint1", "Cluster0")
tpt0_clust0 <- make_df("results/cluster0_timepoint0.rds", "Timepoint0", "Cluster0")

testdf <- rbind(tpt2_clust2, tpt1_clust2, tpt0_clust2,
      tpt2_clust1, tpt1_clust1, tpt0_clust1, 
      tpt2_clust0, tpt1_clust0, tpt0_clust0)

testdf4 <- testdf[(testdf$signature == "TPEXhi") | (testdf$signature == "TSTEMhi"),]


testdf4$cluster[testdf4$cluster == "Cluster0"] <- "CD8_TEM_1"
testdf4$cluster[testdf4$cluster == "Cluster1"] <- "CD8_TEM_2"
testdf4$cluster[testdf4$cluster == "Cluster2"] <- "CD8_TEM_3"

testdf4$timepoint <- gsub("Timepoint", "", testdf4$timepoint)
testdf4$NES[testdf4$p.adjust > 0.1] <- 0
testdf4$log10p[testdf4$p.adjust > 0.1] <- 0
#testdf4 <- testdf4[testdf4$cluster == "CD8_TEM_2",]

testdf4$timepoint <- factor(testdf4$timepoint, levels = c("2", "1", "0"))

hec_sp <- ggplot(testdf4, aes(x = cluster, y = timepoint)) + xlab("signature enrichment R vs NR in CD8 TEM clusters") +
  geom_point(aes(size = log10p, fill = NES), shape = 21, colour = "black") + ylab("Timepoints") +
  scale_size_area(max_size = 10, guide = "legend") + theme(axis.text.x = element_text(angle=45, hjust=1, size = 12), axis.text.y = element_text(size = 12)) + 
  #geom_text(aes(label = label), vjust = 0.65, colour = "white", size = 10) +
  facet_grid(~signature) +
  theme(strip.background =element_rect(fill="black"), panel.border = element_rect(colour = "black", fill=NA, linewidth =1)) +
  theme(strip.text = element_text(colour = 'white'), panel.background = element_blank()) + 
  scale_fill_distiller(palette ="RdBu", direction = -1)

hec_sp
ggsave("results/stem_signatures_across_time_per_cluster.png", hec_sp, dpi = 600, width = 6.73, height = 2.69)

