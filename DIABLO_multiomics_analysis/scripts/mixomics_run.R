library("mixOmics")
library("BiocParallel")

set.seed(123)

pheno <- read.csv("data/phenodata_PLS.csv")
nanostring <- read.csv("data/tumour_PLS.csv")
peripheral_blood <- read.csv("data/peripheralbloodfiltered_PLS.csv")
peripheral_blood2 <- read.csv("data/peripheralbloodC2D1filtered_PLS.csv")

X1 <- as.matrix(peripheral_blood[,-1])
rownames(X1) <- peripheral_blood[,1]

X2 <- as.matrix(nanostring[,-1])
rownames(X2) <- nanostring[,1]

X3 <- as.matrix(peripheral_blood2[,-1])
rownames(X3) <- peripheral_blood2[,1]

Y <- pheno$pd_6m

X <- list(blood = X1, tumour = X2, blood2 = X3)

design = matrix(0.1, ncol = length(X), nrow = length(X), 
                dimnames = list(names(X), names(X)))
diag(design) = 0 # set diagonal to 0s

basic.diablo.model = block.splsda(X = X, Y = Y, ncomp = 5, design = design) 

# run component number tuning with repeated CV
perf.diablo <- perf(basic.diablo.model, validation = 'Mfold', 
                   folds = 5, nrepeat = 10) 

ncomp <- perf.diablo$choice.ncomp$WeightedVote["Overall.BER", "max.dist"] 
# show the optimal choice for ncomp for each dist metric

# set grid of values for each component to test
test.keepX <- list (blood = c(50, 100, 150), 
                   blood2 = c(50, 100, 150),
                   tumour = c(50, 100, 150))

# run the feature selection tuning
tune.TCGA <- tune.block.splsda(X = X, Y = Y, ncomp = ncomp, 
                              test.keepX = test.keepX, design = design,
                              validation = 'Mfold', folds = 5, nrepeat = 10,
                              dist = "max.dist", BPPARAM = MulticoreParam(workers = 2))

list.keepX <- tune.TCGA$choice.keepX # set the optimal values of features to retain

final.diablo.model = block.splsda(X = X, Y = Y, ncomp = ncomp, 
                          keepX = list.keepX, design = design)


saveRDS(final.diablo.model, snakemake@output[[1]])
saveRDS(ncomp, snakemake@output[[2]])
saveRDS(list.keepX, snakemake@output[[3]])




