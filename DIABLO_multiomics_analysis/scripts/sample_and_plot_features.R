
set.seed(1234)
#set.seed(2345)

characters <- c(0:9, letters, LETTERS)
string_length <- 10
random_string <- paste(sample(characters, string_length, replace = TRUE), collapse = "")

library("mixOmics")
library("BiocParallel")
library("gtools")
library("ComplexHeatmap")

create_subset_m <- function(assay_list, prefix, names_of_assays, phenotype, subset_vec){
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

get_test_and_train_sets_m <- function(assay_list, names_of_assays, phenotype){
  vec <- 1:34
  testvec <- sample(vec, size = 10, replace = FALSE)
  trainvec <- vec[!(vec %in% testvec)]
  container <- list()
  train.set <- create_subset_m(assay_list, "train", names_of_assays, phenotype, trainvec)
  test.set <- create_subset_m(assay_list, "test", names_of_assays, phenotype, testvec)
  container[["train_assay"]] <- train.set[[1]]
  container[["train_labels"]] <- train.set[[2]]
  container[["test_assay"]] <- test.set[[1]]
  container[["test_labels"]] <- test.set[[2]]
  container[["test_vec"]] <- testvec
  container[["train_vec"]] <- trainvec
  container
}

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
                                   test.keepX = list.keepX, cpus = 1, progressBar = F) # allow for paralleliation to decrease runtime
  optimal.ncomp <- tune.splsda.srbct$choice.ncomp$ncomp
  optimal.keepX <- tune.splsda.srbct$choice.keepX[1:optimal.ncomp]
  final.splsda <- splsda(X2, Y, 
                         ncomp = optimal.ncomp, 
                         keepX = optimal.keepX)
  final.splsda
}

testing_loop_single_ <- function(model_res, testassay, testlabels){
  predictions <- predict(model_res, testassay)
  confusion.mat = get.confusion_matrix(truth = testlabels,
                                       predicted = predictions$MajorityVote$centroids.dist[,1])
  get.BER(confusion.mat)
}

testing_loop_single_auc_ <- function(model_res, testassay, testlabels){
  auc <- auroc(model_res, newdata = testassay, outcome.test = testlabels)
  auc[[1]][1] # auc value
}

run_whole_test_train_workflow_singleblock <- function(assay, phenotype, filename_save){
  train_test_split1 <- get_test_and_train_sets_(assay, phenotype)
  train.assay <- train_test_split1[["train_assay"]]
  train.labels <- train_test_split1[["train_labels"]]
  print("training model")
  model_res <- training_loop_single_(train.assay, train.labels)
  saveRDS(model_res, filename_save)
  test.assay <- train_test_split1[["test_assay"]]
  test.labels <- train_test_split1[["test_labels"]]
  print("testing_model")
  BER <- testing_loop_single_(model_res, test.assay, test.labels)
  AUC <- testing_loop_single_auc_(model_res, test.assay, test.labels)
  list("BER" = BER, "AUC" = AUC)
}

run_whole_test_train_workflow_multiblock <- function(assay, names_of_assays, phenotype, filename_save){
  train_test_split1 <- get_test_and_train_sets_m(assay, names_of_assays, phenotype)
  train.assay <- train_test_split1[["train_assay"]]
  train.labels <- train_test_split1[["train_labels"]]
  print("training model")
  model_res <- training_loop_single_multiblock(train.assay, train.labels)
  saveRDS(model_res, filename_save)
  test.assay <- train_test_split1[["test_assay"]]
  test.labels <- train_test_split1[["test_labels"]]
  print("testing_model")
  BER <- testing_loop_single_(model_res, test.assay, test.labels)
  AUC <- testing_loop_single_auc_(model_res, test.assay, test.labels)
  list("BER" = BER, "AUC" = AUC)
}

training_loop_single_multiblock <- function(input_assay, pheno_vec){
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
                                test.keepX = test.keepX, design = design, max.iter = 200,
                                validation = 'Mfold', folds = 5, nrepeat = 10,
                                dist = "max.dist", BPPARAM = MulticoreParam(workers = 1))
  
  list.keepX = tune.TCGA$choice.keepX # set the optimal values of features to retain
  
  final.diablo.model = block.splsda(X = X, Y = Y, ncomp = ncomp, 
                                    keepX = list.keepX, design = design)
  final.diablo.model
}

calculate_average_BER_multiblock <- function(number_reps, assay, phenovec, names_of_assays) {
  start_time <- Sys.time()
  #pb <- txtProgressBar(min = 0, max = number_reps, style = 3)
  BER_vec <- list()
  for (i in 1:number_reps) {
    seed_offset <- 0
    repeat {
      set.seed(i + 1000 + seed_offset)  # Increment the seed with an offset
      filename_save <- paste0("results/model_", i, ".rds")
      
      # Add error handling using tryCatch
      tryCatch({
        res_BER <- run_whole_test_train_workflow_multiblock(assay, names_of_assays, phenovec, filename_save)
        print(paste0("Successful repetition :", i))
        BER_vec[[i]] <- res_BER
        #setTxtProgressBar(pb, i)
        if (i %% 5 == 0) {
          collapsed_list <- paste(names_of_assays, collapse = ".")
          filenamesave <- paste0("results/", collapsed_list, "_", gsub(" |-|:", "_", Sys.time()), "_run_", random_string, ".rds")
          saveRDS(BER_vec, filenamesave)
        }
      }, error = function(e) {
        # Handle the error
        print(paste0("Error occurred on repetition ", i, ": ", conditionMessage(e), "\n"))
        seed_offset <- seed_offset + 1  # Increment the offset for the next seed
      })
      
      if (is.null(BER_vec[[i]])) {
        # The result is still missing; continue with a different seed
        print(paste0("Retrying repetition ", i, " with a different seed.\n"))
      } else {
        # Break out of the repeat loop when a valid result is obtained
        break
      }
    }
  }
  end_time <- Sys.time()
  elapsed_time <- end_time - start_time
  print(paste0("Elapsed time: ", elapsed_time[1], " seconds\n"))
  BER_vec
}



calculate_average_BER_singleblock <- function(number_reps, assay, phenovec, names_of_assays){
  start_time <- Sys.time()
  pb <- txtProgressBar(min = 0, max = number_reps, style = 3)
  BER_vec <- list()
  for(i in 1:number_reps){
    filename_save = paste0("results/model_", i, ".rds")
    res_BER <- run_whole_test_train_workflow_singleblock(assay, phenovec, filename_save)
    BER_vec[[i]] <- res_BER
    setTxtProgressBar(pb, i)
    if(i %% 10 == 0){
      collapsed_list <- paste(names_of_assays, collapse = ".")
      filenamesave <- paste0("results/", collapsed_list, "_", gsub(" |-|:", "_", Sys.time()), "_run_", random_string, ".rds")
      saveRDS(BER_vec, filenamesave)
    }
  }
  close(pb)
  end_time <- Sys.time()
  elapsed_time <- end_time - start_time
  cat("Elapsed time: ", elapsed_time[1], " seconds\n")
  BER_vec
}







get_critical_vars_multiblock <- function(model, block_extract){
  res <- selectVar(model, block = block_extract, comp = 1)
  res[[1]][["name"]]
}

get_critical_vars_multiblock <- function(model, block_extract){
  res <- selectVar(model, block = block_extract, comp = 1)
  res[[1]][["name"]]
}

get_critical_vars_singleblock <- function(model){
  res <- selectVar(model)
  res[["name"]]
}

upset_plot_for_runs_multiblock <- function(list_of_models, block_extract){
  upset_list <- list()
  for(i in 1:length(list_of_models)){
    model <- list_of_models[[i]]
    name_of_index <- as.character(i)
    upset_list[[name_of_index]] <- get_critical_vars_multiblock(model, block_extract)
  }
  #m1 = make_comb_mat(upset_list)
  #UpSet(m1)
  #intersection <- Reduce(intersect, upset_list)
  #intersection
  list(comb_mat = "a", list_of_genesets = upset_list)
}

upset_plot_for_runs_singleblock <- function(list_of_models){
  upset_list <- list()
  for(i in 1:length(list_of_models)){
    model <- list_of_models[[i]]
    name_of_index <- as.character(i)
    upset_list[[name_of_index]] <- get_critical_vars_singleblock(model)
  }
  #m1 = make_comb_mat(upset_list)
  #png("results/plot_upset.png", res = 300, width = 7.18, height = 6.15, units = "in")
  list(comb_mat = "a", list_of_genesets = upset_list)
  #dev.off()
}

read_models_and_plot_upset<- function(reps){
  list_of_models <- list()
  for(i in 1:reps){
    filename <- paste0("results/model_", i, ".rds")
    list_of_models[[i]] <- readRDS(filename)
  }
  list_of_models
}

items_in_at_least_n_sets <- function(list_of_sets, n) {
  item_counts <- table(unlist(list_of_sets))
  result <- names(item_counts[item_counts >= n])
  return(result)
} # this extracts from sets

get_single_block <- function(nrep, assay, phenotype, name_of_assay, saveBERfile, saveMODELfile){
  result <- calculate_average_BER_singleblock(nrep, assay, phenotype, name_of_assay)
  saveRDS(result, saveBERfile)
  res_models <- read_models_and_plot_upset(nrep)
  #fname <- paste0("temp_", gsub(" |-|:", "_", Sys.time()), ".rds")
  saveRDS(res_models, saveMODELfile)
  #res_mat <- upset_plot_for_runs_singleblock(res_models)
  #list_mat <- list()
  #for(i in 1:nrep){
  #  col_ind <- as.numeric(colnames(assay) %in% res_mat[[2]][[i]])
  #  list_mat[[i]] <- col_ind
  #}
  #mat <- do.call(rbind, list_mat)
  #list(indicator_mat = mat, gene_list_per_run = list_mat)
}

get_multi_block <- function(nrep, assay, phenotype, name_of_assay, block, saveBERfile, saveMODELfile){
  result <- calculate_average_BER_multiblock(nrep, assay, phenotype, name_of_assay)
  saveRDS(result, saveBERfile)
  res_models <- read_models_and_plot_upset(nrep)
  #fname <- paste0("temp_", gsub(" |-|:", "_", Sys.time()), ".rds")
  saveRDS(res_models, saveMODELfile)
  #res_mat <- upset_plot_for_runs_multiblock(res_models, block)
  #list_mat <- list()
  #assay_genes <- colnames(assay[[block]])
  #for(i in 1:nrep){
  #  col_ind <- as.numeric(assay_genes %in% res_mat[[2]][[i]])
  #  list_mat[[i]] <- col_ind
  #}
  #mat <- do.call(rbind, list_mat)
  #list(indicator_mat = mat, gene_list_per_run = list_mat)
}


