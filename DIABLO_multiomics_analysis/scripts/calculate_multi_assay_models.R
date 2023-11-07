
#set.seed(1234)
set.seed(2345)

library("mixOmics")
library("BiocParallel")
library(gtools)

create_subset_ <- function(assay_list, prefix, names_of_assays, phenotype, subset_vec){
  container <- list()
  assay_container <- list()
  for(i in 1:length(assay_list)){
    assay.subset <- assay_list[[i]][subset_vec,]
    assay.name = names_of_assays[i]
    assay_container[[assay.name]] <- assay.subset
  }
  assay_call <- paste0("assay", "_", prefix)
  container[[assay_call]] <- assay_container
  phenotype.subset <- phenotype[subset_vec]
  label.name = paste0("label", "_", prefix)
  container[[label.name]] <- phenotype.subset
  container
}

get_test_and_train_sets_ <- function(assay_list, names_of_assays, phenotype){
  vec <- 1:34
  testvec <- sample(vec, size = 10, replace = FALSE)
  trainvec <- vec[!(vec %in% testvec)]
  container <- list()
  train.set <- create_subset_(assay_list, "train", names_of_assays, phenotype, trainvec)
  test.set <- create_subset_(assay_list, "test", names_of_assays, phenotype, testvec)
  container[["train_assay"]] <- train.set[[1]]
  container[["train_labels"]] <- train.set[[2]]
  container[["test_assay"]] <- test.set[[1]]
  container[["test_labels"]] <- test.set[[2]]
  container[["test_vec"]] <- testvec
  container[["train_vec"]] <- trainvec
  container
}

training_loop_single_ <- function(input_assay, pheno_vec){
  X <- input_assay
  Y <- as.factor(pheno_vec)
  design = matrix(0.1, ncol = length(X), nrow = length(X), 
                  dimnames = list(names(X), names(X)))
  diag(design) = 0 # set diagonal to 0s
  print("design done")
  
  basic.diablo.model = block.splsda(X = X, Y = Y, ncomp = 5, design = design) 
  print("completed")
  # run component number tuning with repeated CV
  perf.diablo = perf(basic.diablo.model, validation = 'Mfold', 
                     folds = 5, nrepeat = 10) 
  
  ncomp <- perf.diablo$choice.ncomp$WeightedVote["Overall.BER", "max.dist"] 
  
  test.keepX <- list()
  
  for(i in names(input_assay)){
    test.keepX[[i]] <- c(50, 100, 150)
  }
  
  # run the feature selection tuning
  tune.TCGA = tune.block.splsda(X = X, Y = Y, ncomp = ncomp, 
                                test.keepX = test.keepX, design = design,
                                validation = 'Mfold', folds = 5, nrepeat = 10,
                                dist = "centroids.dist", BPPARAM = MulticoreParam(workers = 12))
  
  list.keepX = tune.TCGA$choice.keepX # set the optimal values of features to retain
  
  final.diablo.model = block.splsda(X = X, Y = Y, ncomp = ncomp, 
                                    keepX = list.keepX, design = design)
  final.diablo.model
}


# load data and preprocess

pheno <- read.csv("data/phenodata_PLS.csv")
nanostring <- read.csv("data/tumour_PLS.csv")
peripheral_blood <- read.csv("data/peripheralbloodfiltered_PLS.csv")
peripheral_blood2 <- read.csv("data/peripheralbloodC2D1filtered_PLS.csv")
peripheral_blood3 <- read.csv("data/peripheralbloodC3D1filtered_PLS.csv")

X1 <- as.matrix(peripheral_blood[,-1])
rownames(X1) <- peripheral_blood[,1]

X2 <- as.matrix(peripheral_blood2[,-1])
rownames(X2) <- peripheral_blood2[,1]

X3 <- as.matrix(peripheral_blood3[,-1])
rownames(X3) <- peripheral_blood3[,1]

X4 <- as.matrix(nanostring[,-1])
rownames(X4) <- nanostring[,1]

Y <- pheno$pd_6m

assay_list <- list(X2, X4)
names(assay_list) <- c("blood_tp2", "tumour")
phenotype <- Y
result <- training_loop_single_(assay_list, phenotype)
saveRDS(result, "results/tp2.rds")

assay_list <- list(X3, X4)
names(assay_list) <- c("blood_tp3", "tumour")
phenotype <- Y
result <- training_loop_single_(assay_list, phenotype)
saveRDS(result, "results/tp3.rds")

assay_list <- list(X1, X3, X4)
names(assay_list) <- c("blood_tpt0", "blood_tpt2", "tumour")
phenotype <- Y
result <- training_loop_single_(assay_list, phenotype)
saveRDS(result, "results/tp13.rds")

assay_list <- list(X2, X3, X4)
names(assay_list) <- c("blood_tpt1", "blood_tpt2", "tumour")
phenotype <- Y
result <- training_loop_single_(assay_list, phenotype)
saveRDS(result, "results/tp23.rds")

assay_list <- list(X2, X3)
names(assay_list) <- c("blood_tpt1", "blood_tpt2")
phenotype <- Y
result <- training_loop_single_(assay_list, phenotype)
saveRDS(result, "results/p23.rds")

assay_list <- list(X1, X2)
names(assay_list) <- c("blood_tpt0", "blood_tpt1")
phenotype <- Y
result <- training_loop_single_(assay_list, phenotype)
saveRDS(result, "results/p12.rds")

assay_list <- list(X1, X3)
names(assay_list) <- c("blood_tpt0", "blood_tpt2")
phenotype <- Y
result <- training_loop_single_(assay_list, phenotype)
saveRDS(result, "results/p13.rds")