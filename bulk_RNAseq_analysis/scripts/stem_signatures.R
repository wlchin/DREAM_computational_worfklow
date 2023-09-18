
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(ggplot2)
library(cowplot)

blood0 <- read.csv("DE_results/Base_r_vs_nr.csv")
blood1 <- read.csv("DE_results/C2_r_vs_nr.csv")
blood2 <- read.csv("DE_results/C3_r_vs_nr.csv")
blood3 <- read.csv("DE_results/Base_C2_r.csv")
blood4 <- read.csv("DE_results/Base_C2_nr.csv")
blood5 <- read.csv("DE_results/C2_C3_r.csv")
blood6 <- read.csv("DE_results/C3_C2_nr.csv")

stems <- readRDS('data/stem_all_entrez.rds')
stems <- stems[(stems$genes == "TPEXhi") | (stems$genes == "TSTEMhi") | (stems$genes == "PACE"),]
stems <- stems[,c(2,1)]

query <- blood0$logFC 
names(query) <- blood0$X
query_sorted <- sort(query, decreasing = T)
y <- GSEA(query_sorted, TERM2GENE = stems, pvalueCutoff = 1)
tpt0 <- y@result
tpt0$comparison <- "Timepoint0"

query <- blood1$logFC
names(query) <- blood1$X
query_sorted <- sort(query, decreasing = T)
y <- GSEA(query_sorted, TERM2GENE = stems, pvalueCutoff = 1)
tpt1 <- y@result
tpt1$comparison <- "Timepoint1"

query <- blood2$logFC
names(query) <- blood2$X
query_sorted <- sort(query, decreasing = T)
y <- GSEA(query_sorted, TERM2GENE = stems, pvalueCutoff = 1)
tpt2 <- y@result
tpt2$comparison <- "Timepoint2"

tot <- do.call(rbind, list(tpt0, tpt1, tpt2))
tot$log10p <- -log10(tot$pvalue)
tot$label <- "*"
tot$label[tot$p.adjust >= 0.1] <- " "
tot$num <- c(1)

tot$comparison <- gsub("Timepoint", "", tot$comparison)
tot$comparison <- factor(tot$comparison, levels = c("2", "1", "0"))
tot <- tot[!rownames(tot) %in% c("PACE", "PACE1", "PACE2"),]

hec_sp <- ggplot(tot, aes(y = comparison, x = num)) +
  geom_point(aes(size = log10p, fill = NES), shape = 21, colour = "black") + 
  scale_size_area(max_size = 10, guide = "legend") + 
  theme_cowplot() + ylab("") + 
  theme(strip.background =element_rect(fill="black"), axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.border = element_rect(colour = "black", fill=NA, linewidth =1)) +
  theme(strip.text = element_text(colour = 'white')) + xlab("") + 
  scale_fill_distiller(palette ="RdBu", direction = -1) + ylab("Timepoint") + facet_grid(~Description)

hec_sp
ggsave("results/stem_enrichment_per_timepoint.png", hec_sp, dpi = 600, width = 4.66, height = 1.92)


stems <- readRDS('data/stem_all_entrez.rds')
stems <- stems[(stems$genes == "PACE"),]
stems <- stems[,c(2,1)]

bloodx <- read.csv("DE_results/Base_C2_interaction.csv")
query <- bloodx$logFC
names(query) <- bloodx$X
query_sorted <- sort(query, decreasing = T)
y <- GSEA(query_sorted, TERM2GENE = stems, pvalueCutoff = 1, seed = 42, eps = 10000)
interact_y_t0_t1 <- y@result
gseaplot2(y, c("PACE"), pvalue_table = T, subplots = 1:2, ES_geom = "dot", color = "black")
ggsave("results/interaction_tpt0_tpt1_PACE.png")