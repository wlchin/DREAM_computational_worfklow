
total_df <- readRDS("results/props.rds")

celltypes <- unique(total_df$predicted.celltype.l2)

calculate_anova <- function(celltype){
  cellspecific <- total_df[total_df$predicted.celltype.l2 == celltype,]
  cellspecific$time <- as.factor(c("0", "1", "2", "0", "1", "2"))
  cellspecific$response <- as.factor(c("R", "R", "R", "NR", "NR", "NR"))
  res.aov2 <- aov(freq ~ time + response, data = cellspecific)
  timepval = summary(res.aov2)[[1]]["Pr(>F)"]["time",]
  responsepval = summary(res.aov2)[[1]]["Pr(>F)"]["response",]
  data.frame(celltype, timepval, responsepval)
}

result_list <- list()

for(i in celltypes){
  x <- calculate_anova(i)
  result_list[[i]] <- x
}

df <- do.call(rbind, result_list)
df$timeFDR <- p.adjust(df$timepval)
df$responseFDR <- p.adjust(df$responsepval)

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

df <- df[df$celltype %in% celltypes,]
write.csv(df, "results/anova_results.csv")
