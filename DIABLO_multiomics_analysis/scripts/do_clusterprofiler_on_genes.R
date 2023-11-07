library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)

blood1 <- read.csv("results/test_blood.txt", header = F)
blood2 <- read.csv("results/test_blood2.txt", header = F)
tumour <- read.csv("results/test_tumour.txt", header = F)

gene_background_blood <- colnames(read.csv("data/peripheralbloodfiltered_PLS.csv"))[-c(1)]
gene_background_blood2 <- colnames(read.csv("data/peripheralbloodC2D1filtered_PLS.csv"))[-c(1)]
gene_background_tumour <- colnames(read.csv("data/tumour_PLS.csv"))[-c(1)]
                                   
# conversion to entrezIDs

blood_query_entrez = bitr(blood1$V1, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")$ENTREZID
blood_background_entrez = bitr(gene_background_blood, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")$ENTREZID

blood2_query_entrez = bitr(blood2$V1, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")$ENTREZID
blood2_background_entrez = bitr(gene_background_blood2, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")$ENTREZID

tumour_query_entrez = bitr(tumour$V1, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")$ENTREZID
tumour_background_entrez = bitr(gene_background_tumour, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")$ENTREZID

comb_query_blood <- c(blood_query_entrez, blood2_query_entrez)
comb_query_blood_background <- c(blood_background_entrez, blood2_background_entrez)

ego <- enrichGO(gene          = comb_query_blood,
                universe      = comb_query_blood_background,
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
saveRDS(blood_df, "results/blood_timepoint0.rds")

ego <- enrichGO(gene          = blood2_query_entrez,
                universe      = blood2_background_entrez,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
blood_df <- ego@result
#head(blood_df, 30)
blood_df[blood_df$p.adjust < 0.1,]
saveRDS(blood_df, "results/blood_timepoint1.rds")

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
