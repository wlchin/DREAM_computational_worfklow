library('dplyr')
library("patchwork")
library("ggplot2")
library('ggbeeswarm')
library("dplyr")
library("scales")

clusters = read.csv("data/CD8_seurat_clusters.csv")
colnames(clusters) <- c("index_cell", "cluster_seurat")

width = 0.3

alpha = 0.1

da.res <- read.csv("data/contrast_2.csv")

da.res <- inner_join(da.res, clusters)

da.res <- da.res[!((da.res$cluster_seurat == 3) | (da.res$cluster_seurat == 6)),]

da.res$cluster_seurat[da.res$cluster_seurat == 0] <- "CD8_TEM_1"
da.res$cluster_seurat[da.res$cluster_seurat == 1] <- "CD8_TEM_2"
da.res$cluster_seurat[da.res$cluster_seurat == 2] <- "CD8_TEM_3"
da.res$cluster_seurat[da.res$cluster_seurat == 4] <- "CD8_TEM_4"
da.res$cluster_seurat[da.res$cluster_seurat == 5] <- "CD8_TEM_5"
da.res$cluster_seurat[da.res$cluster_seurat == 7] <- "CD8_TEM_6"
da.res$cluster_seurat[da.res$cluster_seurat == 8] <- "CD8_TEM_7"

da.res <- da.res[da.res$SpatialFDR < alpha,]
da.res <- mutate(da.res, group_by = da.res[,"cluster_seurat"])
da.res <- mutate(da.res, factor(group_by, levels=unique(group_by)))

test2 <- da.res %>%
  mutate(is_signif = ifelse(SpatialFDR < alpha, 1, 0)) %>%
  mutate(logFC_color = ifelse(is_signif==1, logFC, NA)) %>%
  arrange(group_by) %>%
  mutate(Nhood=factor(index_cell, levels=unique(index_cell)))

p2 <- ggplot(test2, aes(group_by, logFC, color=logFC_color), size = 0.1) +
  scale_fill_distiller(palette ="RdBu", direction = -1) +
  guides(color="none") +
  xlab("celltypes") + ylab("Log Fold Change") +
  geom_quasirandom(width = width) +
  coord_flip() + 
  theme_bw(base_size=15) +
  theme(strip.text.y =  element_text(angle=0), axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ggtitle("Timepoint 2")

test2$timepoint = "Timepoint 2"

p2

da.res <- read.csv("data/contrast_1.csv")
da.res <- inner_join(da.res, clusters)

da.res <- da.res[!((da.res$cluster_seurat == 3) | (da.res$cluster_seurat == 6)),]

da.res$cluster_seurat[da.res$cluster_seurat == 0] <- "CD8_TEM_1"
da.res$cluster_seurat[da.res$cluster_seurat == 1] <- "CD8_TEM_2"
da.res$cluster_seurat[da.res$cluster_seurat == 2] <- "CD8_TEM_3"
da.res$cluster_seurat[da.res$cluster_seurat == 4] <- "CD8_TEM_4"
da.res$cluster_seurat[da.res$cluster_seurat == 5] <- "CD8_TEM_5"
da.res$cluster_seurat[da.res$cluster_seurat == 7] <- "CD8_TEM_6"
da.res$cluster_seurat[da.res$cluster_seurat == 8] <- "CD8_TEM_7"
#da.res <- da.res[da.res$nhood_annotation %in% consistent,]
da.res <- da.res[da.res$SpatialFDR < alpha,]
da.res <- mutate(da.res, group_by = da.res[,"cluster_seurat"])
da.res <- mutate(da.res, factor(group_by, levels=unique(group_by)))

test1 <- da.res %>%
  mutate(is_signif = ifelse(SpatialFDR < alpha, 1, 0)) %>%
  mutate(logFC_color = ifelse(is_signif==1, logFC, NA)) %>%
  arrange(group_by) %>%
  mutate(Nhood=factor(index_cell, levels=unique(index_cell)))

p1 <- ggplot(test1, aes(group_by, logFC, color=logFC_color), size = 0.1) +
  #scale_fill_distiller(palette ="RdBu", direction = -1, ) +
  scale_color_gradient2(low = muted("blue"), high = muted("red")) +
  guides(color="none") +
  xlab("celltypes") + ylab("Log Fold Change") +
  geom_quasirandom(width = width) +
  coord_flip() + 
  theme_bw(base_size=15) +
  theme(strip.text.y =  element_text(angle=0), axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ggtitle("Timepoint 1")

p1
### timepoint 2

test1$timepoint = "Timepoint 1"

da.res <- read.csv("data/contrast_0.csv")
da.res <- inner_join(da.res, clusters)

da.res <- da.res[!((da.res$cluster_seurat == 3) | (da.res$cluster_seurat == 6)),]

da.res$cluster_seurat[da.res$cluster_seurat == 0] <- "CD8_TEM_1"
da.res$cluster_seurat[da.res$cluster_seurat == 1] <- "CD8_TEM_2"
da.res$cluster_seurat[da.res$cluster_seurat == 2] <- "CD8_TEM_3"
da.res$cluster_seurat[da.res$cluster_seurat == 4] <- "CD8_TEM_4"
da.res$cluster_seurat[da.res$cluster_seurat == 5] <- "CD8_TEM_5"
da.res$cluster_seurat[da.res$cluster_seurat == 7] <- "CD8_TEM_6"
da.res$cluster_seurat[da.res$cluster_seurat == 8] <- "CD8_TEM_7"
#da.res <- da.res[da.res$nhood_annotation %in% consistent,]
da.res <- da.res[da.res$SpatialFDR < alpha,]
da.res <- mutate(da.res, group_by = da.res[,"cluster_seurat"])
da.res <- mutate(da.res, factor(group_by, levels=unique(group_by)))

test0 <- da.res %>%
  mutate(is_signif = ifelse(SpatialFDR < alpha, 1, 0)) %>%
  mutate(logFC_color = ifelse(is_signif==1, logFC, NA)) %>%
  arrange(group_by) %>%
  mutate(Nhood=factor(index_cell, levels=unique(index_cell)))

p0 <- ggplot(test0, aes(group_by, logFC, color=logFC_color), size = 0.1) +
  scale_color_gradient2(low = muted("blue"), high = muted("red")) +
  guides(color="none") +
  xlab("celltypes") + ylab("Log Fold Change") +
  geom_quasirandom(width = width) +
  coord_flip() + 
  theme_bw(base_size=15) +
  theme(strip.text.y =  element_text(angle=0), axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ggtitle("Timepoint 0")

test0$timepoint = "Timepoint 0"

options(repr.plot.width=25, repr.plot.height=10, repr.plot.quality = 200, repr.plot.res = 90)
p0 | p1 | p2

total_df <- do.call("rbind", list(test0, test1, test2))

#total_df$clusters <- total_df

total_df1 <- total_df[total_df$group_by %in% c("CD8_TEM_1", "CD8_TEM_2", "CD8_TEM_3"),]
total_df1$group_by <- factor(total_df1$group_by, levels = c("CD8_TEM_3", "CD8_TEM_2", "CD8_TEM_1"))
total_df1$timepoint <- gsub("Timepoint ", "", total_df1$timepoint)

pt <- ggplot(total_df1, aes(group_by, logFC, color=logFC_color), size = 0.1) +
  scale_color_gradient2(low = muted("blue"), high = muted("red")) +
  guides(color="none") +
  xlab("") + ylab("Log Fold Change") +
  geom_quasirandom(width = 0.3) +
  coord_flip() +
  theme_bw() + 
  theme(strip.text.y =  element_text(angle=0), axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title=element_text(size = 20)) + 
  facet_grid(~timepoint) + theme(strip.background = element_rect(fill = "black"), strip.text = element_text(color = "white", size = 20))

pt 
ggsave("results/beewswarm.png", width = 7.06, height = 4.18, dpi = 600)
