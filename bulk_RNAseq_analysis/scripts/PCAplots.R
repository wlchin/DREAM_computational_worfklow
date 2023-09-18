library(ggfortify)
library(patchwork)
library(cowplot)

counts <- readRDS("data/entrez_counts.rds")
pheno <- readRDS("data/phenodata.rds")

pheno$PFS_6M_Timepoint <- paste0(pheno$timepoint, "_", pheno$pd_6m)

baselinesamples = pheno$timepoint == "Baseline"
pca_res <- prcomp(t(counts[,baselinesamples]))
baseline <- autoplot(pca_res, data = pheno[baselinesamples,], colour = "PFS_6M_Timepoint") + theme_bw()

baselinesamples = pheno$timepoint == "C2D1"
pca_res <- prcomp(t(counts[,baselinesamples]))
C2D1 <- autoplot(pca_res, data = pheno[baselinesamples,], colour = "PFS_6M_Timepoint") + theme_bw()

baselinesamples = pheno$timepoint == "C3D1"
pca_res <- prcomp(t(counts[,baselinesamples]))
C3D1 <- autoplot(pca_res, data = pheno[baselinesamples,], colour = "PFS_6M_Timepoint") + theme_bw()

baseline / C2D1 / C3D1
ggsave("results/pca_plot.png", dpi = 600, height = 6, width = 4)
