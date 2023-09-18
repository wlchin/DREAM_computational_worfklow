library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)

blood <- read.csv("results/test_blood.txt", header = F)
tumour <- read.csv("results/test_tumour.txt", header = F)

gene_background_blood <- colnames(read.csv("data/peripheralbloodfiltered_PLS.csv"))[-c(1)]
gene_background_tumour <- colnames(read.csv("data/tumour_PLS.csv"))[-c(1)]
                                   
# conversion to entrezIDs
firstCharacter = substr(blood$V1[1][1],1,1)

if(firstCharacter == "X"){
  blood_query_entrez <- gsub("X", "", blood$V1)
  blood_background_entrez <- gsub("X", "", gene_background_blood)
} else{
  blood_query_entrez = bitr(blood$V1, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")$ENTREZID
  blood_background_entrez = bitr(gene_background_blood, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")$ENTREZID
}

tumour_query_entrez = bitr(tumour$V1, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")$ENTREZID
tumour_background_entrez = bitr(gene_background_tumour, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")$ENTREZID


ego <- enrichGO(gene          = blood_query_entrez,
                universe      = blood_background_entrez,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
blood_df <- ego@result
#head(blood_df, 30)
blood_df[blood_df$p.adjust < 0.1,]
saveRDS(blood_df, "results/blood.rds")

ego <- enrichGO(gene          = tumour_query_entrez,
                universe      = tumour_background_entrez,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.1,
                readable      = TRUE)

tumour_df <- ego@result
head(tumour_df)
saveRDS(tumour_df, "results/tumour.rds")
