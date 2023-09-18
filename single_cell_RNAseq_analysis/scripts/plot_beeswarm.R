library('dplyr')
library("patchwork")
library("ggplot2")
library('ggbeeswarm')
library("scales")

width = 0.6
alpha = 0.1

da.res <- read.csv(snakemake@input[[3]])
da.res <- da.res[da.res$SpatialFDR < alpha,]
da.res <- mutate(da.res, group_by = da.res[,"nhood_annotation"])
da.res <- mutate(da.res, factor(group_by, levels=unique(group_by)))

test2 <- da.res %>%
  mutate(is_signif = ifelse(SpatialFDR < alpha, 1, 0)) %>%
  mutate(logFC_color = ifelse(is_signif==1, logFC, NA)) %>%
  arrange(group_by) %>%
  mutate(Nhood=factor(index_cell, levels=unique(index_cell)))

test2$Timepoint <- "Timepoint 2"
  
p2 <- ggplot(test2, aes(group_by, logFC, color=logFC_color), size = 0.01) +
  scale_color_gradient2() +
  guides(color="none") +
  xlab("celltypes") + ylab("Log Fold Change") +
  geom_quasirandom(width = width) +
  coord_flip() + 
  theme_bw(base_size=15) +
  theme(strip.text.y =  element_text(angle=0)) + ggtitle("Timepoint 2")

### timepoint 2

da.res <- read.csv(snakemake@input[[2]])
da.res <- da.res[da.res$SpatialFDR < alpha,]
da.res <- mutate(da.res, group_by = da.res[,"nhood_annotation"])
da.res <- mutate(da.res, factor(group_by, levels=unique(group_by)))

test1 <- da.res %>%
  mutate(is_signif = ifelse(SpatialFDR < alpha, 1, 0)) %>%
  mutate(logFC_color = ifelse(is_signif==1, logFC, NA)) %>%
  arrange(group_by) %>%
  mutate(Nhood=factor(index_cell, levels=unique(index_cell)))

test1$Timepoint <- "Timepoint 1"

p1 <- ggplot(test1, aes(group_by, logFC, color=logFC_color), size = 0.01) +
  scale_color_gradient2() +
  guides(color="none") +
  xlab("celltypes") + ylab("Log Fold Change") +
  geom_quasirandom(width = width) +
  coord_flip() + 
  theme_bw(base_size=15) +
  theme(strip.text.y =  element_text(angle=0)) + ggtitle("Timepoint 1")

### timepoint 0

da.res <- read.csv(snakemake@input[[1]])
da.res <- da.res[da.res$SpatialFDR < alpha,]
da.res <- mutate(da.res, group_by = da.res[,"nhood_annotation"])
da.res <- mutate(da.res, factor(group_by, levels=unique(group_by)))

test0 <- da.res %>%
  mutate(is_signif = ifelse(SpatialFDR < alpha, 1, 0)) %>%
  mutate(logFC_color = ifelse(is_signif==1, logFC, NA)) %>%
  arrange(group_by) %>%
  mutate(Nhood=factor(index_cell, levels=unique(index_cell)))

test0$Timepoint <- "Timepoint 0"

p0 <- ggplot(test0, aes(group_by, logFC, color=logFC_color), size = 0.01) +
  scale_color_gradient2() +
  guides(color="none") +
  xlab("celltypes") + ylab("Log Fold Change") +
  geom_quasirandom(width = width) +
  coord_flip() + 
  theme_bw(base_size=15) +
  theme(strip.text.y =  element_text(angle=0)) + ggtitle("Timepoint 0")

total_df <- do.call("rbind", list(test0, test1, test2))

pt <- ggplot(total_df, aes(group_by, logFC, color=logFC_color), size = 0.1) +
  scale_color_gradient2(low = muted("blue"), high = muted("red")) +
  guides(color="none") +
  xlab("") + ylab("Log Fold Change") +
  geom_quasirandom(width = 0.3) + scale_x_discrete(drop=FALSE) + 
  coord_flip() +
  theme_bw(base_size=15) + 
  theme(strip.text.y =  element_text(angle=0), axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title=element_text(size=20)) + 
  facet_grid(~Timepoint)

pt + theme(strip.background = element_rect(fill = "black"), strip.text = element_text(color = "white", size = 20))
ggsave("results/beeswarm.png", width = 9, height = 10, dpi = 600)

