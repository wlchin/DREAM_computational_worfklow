set.seed(1234)
#set.seed(2345)

characters <- c(0:9, letters, LETTERS)
string_length <- 10
random_string <- paste(sample(characters, string_length, replace = TRUE), collapse = "")


library("mixOmics")
library("BiocParallel")
library(gtools)

get_test_and_train_sets_ <- function(assay, phenotype){
  vec <- 1:34
  testvec <- sample(vec, size = 10, replace = FALSE)
  trainvec <- vec[!(vec %in% testvec)]
  container <- list()
  container[["train_assay"]] <- assay[trainvec,]
  container[["train_labels"]] <- phenotype[trainvec]
  container[["test_assay"]] <- assay[testvec,]
  container[["test_labels"]] <- phenotype[testvec]
  container[["test_vec"]] <- testvec
  container[["train_vec"]] <- trainvec
  container
}

training_loop_single_ <- function(input_assay, pheno_vec){
  X2 <- input_assay
  Y <- pheno_vec
  srbct.splsda <- splsda(X2, Y, ncomp = 5)
  perf.splsda.srbct <- perf(srbct.splsda, validation = "Mfold", 
                            folds = 5, nrepeat = 10, # use repeated cross-validation
                            progressBar = TRUE, auc = TRUE) # include AUC values
  # grid of possible keepX values that will be tested for each component
  list.keepX <- c(50, 100, 150)
  # undergo the tuning process to determine the optimal number of variables
  tune.splsda.srbct <- tune.splsda(X2, Y, ncomp = 5, # calculate for first 4 components
                                   validation = 'Mfold',
                                   folds = 5, nrepeat = 10, # use repeated cross-validation
                                   dist = 'max.dist', # use max.dist measure
                                   measure = "BER", # use balanced error rate of dist measure
                                   test.keepX = list.keepX) # allow for paralleliation to decrease runtime
  optimal.ncomp <- tune.splsda.srbct$choice.ncomp$ncomp
  optimal.keepX <- tune.splsda.srbct$choice.keepX[1:optimal.ncomp]
  final.splsda <- splsda(X2, Y, 
                         ncomp = optimal.ncomp, 
                         keepX = optimal.keepX)
  final.splsda
}


pheno <- read.csv("data/phenodata_PLS.csv")
nanostring <- read.csv("data/peripheralbloodfiltered_PLS.csv")
X2 <- as.matrix(nanostring[,-1])
rownames(X2) <- nanostring[,1]
result <- training_loop_single_ (X2, pheno$pd_6m)
saveRDS(result, "results/p1.rds")

pheno <- read.csv("data/phenodata_PLS.csv")
nanostring <- read.csv("data/peripheralbloodC2D1filtered_PLS.csv")
X2 <- as.matrix(nanostring[,-1])
rownames(X2) <- nanostring[,1]
result <- training_loop_single_ (X2, pheno$pd_6m)
saveRDS(result, "results/p2.rds")

pheno <- read.csv("data/phenodata_PLS.csv")
nanostring <- read.csv("data/peripheralbloodC3D1filtered_PLS.csv")
X2 <- as.matrix(nanostring[,-1])
rownames(X2) <- nanostring[,1]
result <- training_loop_single_ (X2, pheno$pd_6m)
saveRDS(result, "results/p3.rds")

#res <- calculate_average_BER(5, X2, pheno$pd_6m)
