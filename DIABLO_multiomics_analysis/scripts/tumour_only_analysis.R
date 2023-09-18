library("mixOmics")
pheno <- read.csv("data/phenodata_PLS.csv")
nanostring <- read.csv("data/tumour_PLS.csv")

set.seed(123)

X <- as.matrix(nanostring[,-1])
rownames(X) <- nanostring[,1]

Y <- pheno$pd_6m

srbct.splsda <- splsda(X, Y, ncomp = 10)

perf.splsda.srbct <- perf(srbct.splsda, validation = "Mfold", 
                          folds = 5, nrepeat = 10, # use repeated cross-validation
                          progressBar = TRUE, auc = TRUE) # include AUC values

# grid of possible keepX values that will be tested for each component
list.keepX <- c(50, 100, 150)

# undergo the tuning process to determine the optimal number of variables
tune.splsda.srbct <- tune.splsda(X, Y, ncomp = 5, # calculate for first 4 components
                                 validation = 'Mfold',
                                 folds = 5, nrepeat = 10, # use repeated cross-validation
                                 dist = 'max.dist', # use max.dist measure
                                 measure = "BER", # use balanced error rate of dist measure
                                 test.keepX = list.keepX) # allow for paralleliation to decrease runtime

optimal.ncomp <- tune.splsda.srbct$choice.ncomp$ncomp
optimal.keepX <- tune.splsda.srbct$choice.keepX[1:optimal.ncomp]

final.splsda <- splsda(X, Y, 
                       ncomp = 1, 
                       keepX = optimal.keepX)
saveRDS(final.splsda, snakemake@output[[1]])
saveRDS(optimal.ncomp, snakemake@output[[2]])
saveRDS(optimal.keepX, snakemake@output[[3]])
