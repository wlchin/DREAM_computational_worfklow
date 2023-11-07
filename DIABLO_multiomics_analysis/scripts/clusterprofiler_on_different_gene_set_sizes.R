library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)

process_enrichment <- function(query_geneset_file, background_geneset_file, enrichment_result_as_rds){
  blood <- read.csv(query_geneset_file, header = F)
  gene_background_blood <- read.csv(background_geneset_file, header = F)
  blood_query_entrez = bitr(blood$V1, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")$ENTREZID
  blood_background_entrez = bitr(gene_background_blood$V1, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")$ENTREZID
  
  ego <- enrichGO(gene          = blood_query_entrez,
                  universe      = blood_background_entrez,
                  OrgDb         = org.Hs.eg.db,
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
  blood_df <- ego@result
  saveRDS(blood_df, enrichment_result_as_rds)
}

process_enrichment("results/blood_responder_sampled_25genes.txt", "results/gene_blood_background_new_seed.txt", "results/blood_25_genes.rds")
process_enrichment("results/blood_responder_sampled_50genes.txt", "results/gene_blood_background_new_seed.txt", "results/blood_50_genes.rds")
process_enrichment("results/blood_responder_sampled_100genes.txt", "results/gene_blood_background_new_seed.txt", "results/blood_100_genes.rds")


process_enrichment("results/tumour_nonresponder_sampled_25genes.txt", "results/gene_tumour_background_new_seed.txt", "results/tumour_25_genes.rds")
process_enrichment("results/tumour_nonresponder_sampled_50genes.txt", "results/gene_tumour_background_new_seed.txt", "results/tumour_50_genes.rds")
process_enrichment("results/tumour_nonresponder_sampled_100genes.txt", "results/gene_tumour_background_new_seed.txt", "results/tumour_100_genes.rds")

