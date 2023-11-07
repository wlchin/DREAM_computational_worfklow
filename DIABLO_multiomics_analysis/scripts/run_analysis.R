source("scripts/sample_and_plot_features.R")

set.seed(123)


pheno <- read.csv("data/phenodata_PLS.csv")
nanostring <- read.csv("data/tumour_PLS.csv")
peripheral_blood <- read.csv("data/peripheralbloodfiltered_PLS.csv")
peripheral_blood2 <- read.csv("data/peripheralbloodC2D1filtered_PLS.csv")

X1 <- as.matrix(peripheral_blood[,-1])
rownames(X1) <- peripheral_blood[,1]

X2 <- as.matrix(peripheral_blood2[,-1])
rownames(X2) <- peripheral_blood2[,1]

X4 <- as.matrix(nanostring[,-1])
rownames(X4) <- nanostring[,1]

Y <- pheno$pd_6m

assay_list <- list(X1, X2, X4)
names_of_assays <- c("blood_tpt0", "blood_tpt1", "tumour")
names(assay_list) <- names_of_assays
phenotype <- Y
get_multi_block(100, assay = assay_list, 
                block = "blood_tpt1", 
                name_of_assay = names_of_assays, 
                phenotype = Y, 
                saveBERfile = "results/BER_tp12_res_no_norm_seeds.rds",
                saveMODELfile = "data/tp12_sampled.rds") # data/tp12_no_norm_seeds.rds