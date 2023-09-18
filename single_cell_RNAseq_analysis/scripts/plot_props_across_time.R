library(ggplot2)
library(ggsci)
library(cowplot)

total_df <- readRDS("results/props.rds")

plot_stuff <- function(total_df, celltype){
  totaldf_subset <- total_df[total_df$predicted.celltype.l2 == celltype,]
  totaldf_subset$timepoint <- c("0", "1", "2", "0", "1", "2")
  totaldf_subset$response <- c("R", "R", "R", "NR", "NR", "NR")
  a <- ggplot(totaldf_subset, aes(timepoint, freq)) + geom_col(aes(fill = response), position = "dodge") + 
    ggtitle(celltype) + theme_cowplot() + theme(plot.title = element_text(hjust = 0.5, face = "plain", size = 12), legend.position ='none') +
    scale_fill_aaas()
  a
}


calculate_stats <- function(total_df, celltype){
  totaldf_subset <- total_df[total_df$predicted.celltype.l2 == celltype,]
  totaldf_subset$timepoint <- c("0", "1", "2", "0", "1", "2")
  totaldf_subset$response <- c("R", "R", "R", "NR", "NR", "NR")
  a <- ggplot(totaldf_subset, aes(timepoint, freq)) + geom_col(aes(fill = response), position = "dodge") + 
    ggtitle(celltype) + theme_cowplot() + theme(plot.title = element_text(hjust = 0.5, face = "plain", size = 12), legend.position ='none') +
    scale_fill_aaas()
  a
}

totaldf_subset <- total_df[total_df$predicted.celltype.l2 == "NK",]

celltypes <- c("CD14 Mono",
               "CD4 TCM",
               "NK",
               "CD8 TEM",
               "CD4 Naive",
               "CD16 Mono",
               "B naive",
               "CD8 Naive",
               "CD4 TEM",
               "Treg")

pval <- list()

for(i in celltypes){
  pval[[i]] <- plot_stuff(total_df, i)
}

plot_grid(plotlist = pval, ncol = 5)
ggsave("results/individual_barplots_per_timepoint.png" , width = 9.30, height = 3.89, dpi = 600)

