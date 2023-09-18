library(ggplot2)
library(dplyr)
library(ggplot2)
library(cowplot)
library(ggsci)

stuff <- function(filename){
  filename_ <- paste0("DE_results/", filename)
  celltype <- read.csv(filename_)
  celltype$cell <- gsub(".csv", "", filename)
  celltype$sigpvalue <- ifelse(celltype$adj.P.Val < 0.05, "p < 0.05","p > 0.05")
  celltype$sig <- ifelse(celltype$adj.P.Val < 0.05 & abs(celltype$logFC) > 0.58, "Significant","Not Significant")
  celltype$avg_log2FC <- ifelse(celltype$logFC > 0, "up" , "down")
  celltype
}

filelist <- list.files("DE_results/", "*r.csv")

startList <- list()

for(i in filelist){
  startList[[i]] <- stuff(i)
}

data <- do.call(rbind, startList)
data$cell <- factor(data$cell)

df_Count <- data %>% group_by(sig, cell, avg_log2FC) %>% count()
df_Count <- data.frame(df_Count)
df_Count  <- filter(df_Count, df_Count$sig == "Significant")

df_Count$cell <- as.character(df_Count$cell)

df <- data.frame(sig  = c("Significant", "Significant", "Significant", "Significant", "Significant", "Significant", "Significant", "Significant", "Significant"),
                cell = c("C2_C3_r", "Timepoint 0 R vs Timepoint 0 NR", "Timepoint 0 R vs Timepoint 0 NR", "Timepoint 1 R vs Timepoint 1 NR", "Timepoint 1 R vs Timepoint 1 NR", "Timepoint 2 R vs Timepoint 2 NR", "Timepoint 2 R vs Timepoint 2 NR", "Timepoint 1 NR vs Timepoint 2 NR", "Timepoint 1 NR vs Timepoint 2 NR"),
                avg_log2FC = c("down", "up", "down", "up", "down", "up", "down", "up", "down"),
                n = as.integer(c(0, 0, 0, 0, 0, 0, 0, 0, 0))
)
ddd <- rbind(df_Count, df)

ddd$cell[ddd$cell == "C2_C3_r"] <- "Timepoint 1 R vs Timepoint 2 R"
ddd$cell[ddd$cell == "Base_C2_r"] <- "Timepoint 0 R vs Timepoint 1 R"
ddd$cell[ddd$cell == "Base_C2_nr"] <- "Timepoint 0 NR vs Timepoint 1 NR"

ddd$cell <- factor(ddd$cell, levels = c("Timepoint 0 R vs Timepoint 0 NR", 
          "Timepoint 1 R vs Timepoint 1 NR", 
          "Timepoint 2 R vs Timepoint 2 NR", 
          "Timepoint 0 R vs Timepoint 1 R",
          "Timepoint 1 R vs Timepoint 2 R",
          "Timepoint 0 NR vs Timepoint 1 NR", 
          "Timepoint 1 NR vs Timepoint 2 NR"))

BC3 <- ggplot(ddd, aes(x = cell, y = n, fill = avg_log2FC)) + geom_col(position = "dodge") + 
  #geom_text(aes(label = n),  vjust=-1) + 
  geom_text(aes(label=n), position=position_dodge(width=0.9), vjust=0.5, hjust = -0.5) + labs(fill="logFC direction") + scale_fill_npg() +
  #scale_fill_manual(values = c("red3", "dodgerblue4")) + 
  scale_x_discrete(drop=FALSE)+ theme_cowplot() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 15), axis.text.y = element_text(size = 15)) + ylab("number of DE genes") + xlab("") + coord_flip()

BC3
ggsave("results/tally_DE_genes.png", width = 7.15, height = 3.82)
