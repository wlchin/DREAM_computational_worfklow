
library("PCGSE")
library(mixOmics)

stem <- readRDS("data/stem_all_gene_name.rds")
colnames(stem) <- c("genes", "gene.set.name")

generate_gene_sets <- function(data, stem){
  
  firstCharacter = substr(colnames(data)[1],1,1)
  if(firstCharacter == "X"){
    colnames(data) <- gsub("X", "", colnames(data))
  }
  
  PACE <- stem[stem$gene.set.name == "PACE",]$genes
  PACE_bin <- colnames(data) %in% PACE + 0
  print("PACE")
  print(sum(PACE_bin))
  
  TSTEMhi <- stem[stem$gene.set.name == "TSTEMhi",]$genes
  TSTEMhi_bin <- (colnames(data) %in% TSTEMhi) + 0
  print("STEMhi")
  print(sum(TSTEMhi_bin))
  
  TPEXhi <- stem[stem$gene.set.name == "TSTEMhi",]$genes
  TPEXhi_bin <- (colnames(data) %in% TPEXhi) + 0
  print("TPEXhi")
  print(sum(TPEXhi_bin))
  
  TSTEMG <- stem[stem$gene.set.name == "GALLETTI_TSTEM",]$genes
  TSTEMG_bin <- (colnames(data) %in% TSTEMG) + 0
  print("TSTEMg")
  print(sum(TSTEMG_bin))
  
  TPEXG <- stem[stem$gene.set.name == "GALLETTI_TPEX",]$genes
  TPEXG_bin <- (colnames(data) %in% TPEXG) + 0
  print("TPEXg")
  print(sum(TPEXG_bin))
  
  GATTINONI <- stem[stem$gene.set.name == "GATTINONI",]$genes
  GATTINONI_bin <- (colnames(data) %in% GATTINONI) + 0
  print("Gattinoni")
  print(sum(GATTINONI_bin))
  
  MILLER <- stem[stem$gene.set.name == "Miller_prog_ex",]$genes
  MILLER_bin <- (colnames(data) %in% MILLER) + 0
  print("Miller")
  print(sum(MILLER_bin))
  
  zheng <- stem[stem$gene.set.name == "zhengTCF7PEX",]$genes
  zheng_bin <- (colnames(data) %in% zheng) + 0
  print("Zheng")
  print(sum(zheng_bin))
  
  #genesets_mat <- as.matrix(rbind(PACE_bin, TSTEMhi_bin, TPEXhi_bin, TSTEMG_bin, TPEXG_bin, GATTINONI_bin, MILLER_bin, SIDDIQ_bin))
  
  #genesets_mat <- as.matrix(rbind(PACE_bin, TSTEMG_bin))
  #rownames(genesets_mat) <- c("PACE", "TSTEMhi", "TPEXhi", "TSTEMG", "TPEXG", "GATTINONI", "MILLER", "SIDDIQ")
  genesets_mat <- as.matrix(rbind(PACE_bin, TSTEMhi_bin, TPEXhi_bin))
  rownames(genesets_mat) <- c("PACE", "TSTEMhi", "TPEXhi")
  #genesets_mat <- as.matrix(rbind(PACE_bin, TSTEMhi_bin, TPEXhi_bin, MILLER_bin, zheng_bin, GATTINONI_bin, TSTEMG_bin, TPEXG_bin))
  #rownames(genesets_mat) <- c("PACE", "TSTEMhi", "TPEXhi", "Miller", "Zheng", "Gattinoni", "TSTEMc", "TPEXc")
  genesets_mat
}

global.StandAveDiff = function(C.mat,local.stats,...) { 
  gene.set.mat = t(as.matrix(C.mat))
  return (function(gene.statistics=local.stats,gene.sets=gene.set.mat,...) {  
    num.gene.sets = nrow(gene.set.mat)
    gene.set.statistics = rep(0, num.gene.sets)
    for (i in 1:num.gene.sets) {
      indexes.for.gene.set = which(gene.set.mat[i,]==1)
      m1 = length(indexes.for.gene.set)
      not.gene.set.indexes = which(!(1:ncol(gene.set.mat) %in% indexes.for.gene.set))
      m2 = length(not.gene.set.indexes)
      # compute the mean difference of the gene-level statistics
      mean.diff = mean(gene.statistics[indexes.for.gene.set]) - mean(gene.statistics[not.gene.set.indexes])
      # compute the pooled standard deviation
      pooled.sd = sqrt(((m1-1)*var(gene.statistics[indexes.for.gene.set]) + (m2-1)*var(gene.statistics[not.gene.set.indexes]))/(m1+m2-2))      
      # compute the t-statistic
      gene.set.statistics[i] = mean.diff/(pooled.sd*sqrt(1/m1 + 1/m2))
    }
    return (gene.set.statistics)
  })
}

assess_pgse <- function(x = diablo_obj, genesets_df, ind = i, pind = pind, statistic = "mean.diff", trans = "none", perm = 10000, gene.statistic = "z", strategy = "cor.adj.parametric"){
  data = x$X[[ind]]
  genesets_mat <- generate_gene_sets(data, genesets_df)
  prcompres <- list(x = x$variates[[ind]])
  res3 <- pcgse(data = data, prcomp.output = prcompres, pc.indexes = pind,
        gene.sets = genesets_mat, transformation = trans, gene.statistic = gene.statistic,
        gene.set.statistic = statistic, gene.set.test = strategy, nperm = perm)
  res3_df <- data.frame(res = res3$p.values, FDR = p.adjust(res3$p.values))
  res3_df$gene_set_names <- rownames(genesets_mat)
  res3_df$assay <- names(x$X)[ind]
  res3_df
}

# this is required for PGSCE to run

local.GeneStatistics = function(X.mat, y.vec,...){ 
  args.local = list(...)[[1]]
  gene.stat = args.local$gene.statistic  
  trans = args.local$transformation    
  return (function(data=X.mat, vector=y.vec, gene.statistic = gene.stat, transformation = trans,...) {  
    p = length(vector)
    n = ncol(data) 
    gene.statistics = rep(0, p)
    # compute the Pearson correlation between the selected PCs and the data
    gene.statistics = cor(t(data), vector) 
    if (gene.statistic == "z") {
      # use Fisher's Z transformation to convert to Z-statisics
      gene.statistics = sapply(gene.statistics, function(x) {
        return (sqrt(n-3)*atanh(x))})      
    }    
    # Absolute value transformation of the gene-level statistics if requested
    if (transformation == "abs.value") {
      gene.statistics = sapply(gene.statistics, abs)
    }
    return (gene.statistics)
  })
}

analyse_diablo_obj <- function(objrds, gs, nperm = 10000, statistic = "rank.sum", trans = "abs.value", strategy = "permutation"){
  x <- objrds
  df_list <- list()
  for(i in 1:length(x$X)){
    for(j in 1:x$ncomp[1]){
      print("doing assay:")
      print(i)
      print("doing component:")
      print(j)
      x1 <- assess_pgse(x, genesets_df = gs, ind = i, pind = j, statistic = statistic, trans = trans, perm = nperm, strategy = strategy)
      x1$assay_index <- i
      x1$component_index <- j
      inder <- paste0(i, "_", j)
      df_list[[inder]] <- x1
    }
    }
  
  df_all <- do.call(rbind, df_list)
  df_all
}

#genesets_mat <- readRDS("results/genesets_PCGSE.rds")

#n_tot = 1000

tp123 <- readRDS("results/diagnostics/full_model_all_samples_tp123.rds")
a <- analyse_diablo_obj(tp123, stem)
a$combination <- "tp123"
saveRDS(a, "results/tp123_pcgse.rds")
print("tp123 done")

p123 <- readRDS("results/diagnostics/full_model_all_samples_p123.rds")
b <- analyse_diablo_obj(p123, stem)
b$combination <- "p123"
saveRDS(b, "results/p123_pcgse.rds")
print("p123 done")

tp1 <- readRDS("results/diagnostics/full_model_all_samples_tp1.rds")
c <- analyse_diablo_obj(tp1, stem)
c$combination <- "tp1"
saveRDS(c, "results/tp1_pcgse.rds")
print("tp1 done")

tp12 <- readRDS("results/diagnostics/full_model_all_samples_tp12.rds")
d <- analyse_diablo_obj(tp12, stem)
d$combination <- "tp12"
saveRDS(d, "results/tp12_pcgse.rds")
print("tp12 done")

totaldf <- do.call(rbind, list(a,b,c,d))
saveRDS(totaldf, "results/pgcse_result.rds")


