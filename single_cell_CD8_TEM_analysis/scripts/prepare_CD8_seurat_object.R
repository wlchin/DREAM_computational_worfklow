library(SingleCellExperiment)
library(zellkonverter)
library(readr)
library(monocle3)
library(scran)
library(Seurat)
library(ggplot2)
library(dplyr)
#devtools::install_github('cole-trapnell-lab/monocle3', ref="develop")

#file <- system.file("extdata", "CD8TEM.h5ad", package = "zellkonverter")
sce <- readH5AD(snakemake@input[[1]], reader = "R") # ignore warning message

CD8TEM <- read.csv(snakemake@input[[2]]) # metadata file
rownames(CD8TEM) <- CD8TEM$X

expression_matrix <- assay(sce, "X")
rownames(expression_matrix) <- rowData(sce)$gene
colnames(expression_matrix) <- CD8TEM$X

print("row")
head(rownames(expression_matrix))

print("col")
head(colnames(expression_matrix))

CD8_TEM <- CreateSeuratObject(expression_matrix, project = "SeuratProject", assay = "RNA",
                   min.cells = 0, min.features = 0, meta.data = CD8TEM)

pbmc <- NormalizeData(CD8_TEM, normalization.method = "LogNormalize", scale.factor = 10000)

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

pbmc <- subset(pbmc, subset = percent.mt < 10)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))


pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.3)
pbmc <- RunUMAP(pbmc, dims = 1:10)

DimPlot(pbmc, reduction = "umap", label = T) + NoLegend()
ggsave("results/dimplot.pdf", width = 5.02, height = 3.74)

saveRDS(pbmc, snakemake@output[[1]])

