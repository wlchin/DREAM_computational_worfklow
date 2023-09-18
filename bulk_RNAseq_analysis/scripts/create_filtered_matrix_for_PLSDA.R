library(clusterProfiler)
library(org.Hs.eg.db)

x <- read.csv("results/filtered_count_matrix.csv", row.names = "X")

conv <- bitr(rownames(x), fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Hs.eg.db")
common <- intersect(conv$ENTREZID, rownames(x))

all_genes <- x[rownames(x) %in% common,]
conv_all <- conv[conv$ENTREZID %in% common,]

rownames(conv_all) <- conv_all$ENTREZID
rownames(all_genes) <- conv_all[rownames(all_genes),]$SYMBOL

write.csv(all_genes, "results/count_matrix.csv")

