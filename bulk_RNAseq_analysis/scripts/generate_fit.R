library('variancePartition')
library('edgeR')
library('BiocParallel')
library('ggplot2')

count_mat <- readRDS("data/entrez_counts.rds")
pheno_df <- readRDS("data/phenodata.rds")
pheno_df$PFS_6M_response <- paste0(pheno_df$timepoint, "_", pheno_df$pd_6m)

isexpr = rowSums(cpm(count_mat)>10.0) >= 17 
write.csv(count_mat[isexpr,], "results/filtered_count_matrix.csv")

geneExpr = DGEList(count_mat[isexpr,] )

geneExpr = calcNormFactors(geneExpr)
form <- ~0 + PFS_6M_response + (1|patient) 
vobjDream = voomWithDreamWeights(geneExpr, form, pheno_df)


L2 = makeContrastsDream(form, pheno_df, contrasts = 
                           c(Base_C2_interaction = "(PFS_6M_responseC2D1_no - PFS_6M_responseBaseline_no) - (PFS_6M_responseC2D1_yes - PFS_6M_responseBaseline_yes)",
                             Base_C3_interaction = "(PFS_6M_responseC3D1_no - PFS_6M_responseBaseline_no) - (PFS_6M_responseC3D1_yes - PFS_6M_responseBaseline_yes)",
                             C2_C3_interaction = "(PFS_6M_responseC3D1_no - PFS_6M_responseC2D1_no) - (PFS_6M_responseC3D1_yes - PFS_6M_responseC2D1_yes)",
                             C2_C3_r = "(PFS_6M_responseC3D1_no - PFS_6M_responseC2D1_no)",
                             Base_C2_r = "(PFS_6M_responseC2D1_no - PFS_6M_responseBaseline_no)",
                             C3_C2_nr = "(PFS_6M_responseC3D1_yes - PFS_6M_responseC2D1_yes)",
                             Base_C2_nr = "(PFS_6M_responseC2D1_yes - PFS_6M_responseBaseline_yes)",
                             Base_r_vs_nr = "(PFS_6M_responseBaseline_no - PFS_6M_responseBaseline_yes)",
                             C2_r_vs_nr = "(PFS_6M_responseC2D1_no - PFS_6M_responseC2D1_yes)",
                             C3_r_vs_nr = "(PFS_6M_responseC3D1_no - PFS_6M_responseC3D1_yes)"
                             ))

fit = dream(vobjDream, form, pheno_df, L=L2)
fit = eBayes(fit)

saveRDS(fit, "results/fit_obj.rds")

