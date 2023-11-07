
library(mixOmics)
library(dplyr)

p123 <- readRDS("results/diagnostics/full_model_all_samples_p123.rds")
tp123 <- readRDS("results/diagnostics/full_model_all_samples_tp123.rds")
tp12 <- readRDS("results/diagnostics/full_model_all_samples_tp12.rds")
tp1 <- readRDS("results/diagnostics/full_model_all_samples_tp1.rds")
tonly <- readRDS("results/diagnostics/full_model_all_samples_t.rds")


tp23 <- readRDS("results/tp23.rds")
tp13 <- readRDS("results/tp13.rds")
tp1 <- readRDS("results/tp1.rds")
tp2 <- readRDS("results/tp2.rds")
tp3 <- readRDS("results/tp3.rds")
p1 <- readRDS("results/p1.rds")
p2 <- readRDS("results/p2.rds")
p3 <- readRDS("results/p3.rds")
p12 <- readRDS("results/p12.rds")
p23 <- readRDS("results/p23.rds")
p13 <- readRDS("results/p13.rds")

extract_metrics <- function(perfobj, desc){
  perfobj = perf(perfobj, validation = 'loo', auc = T) 
  rot = length(perfobj$auc)
  sum_auc = 0
  for(i in 1:rot){
    sum_auc = sum_auc + perfobj$auc[[i]][1]
  }
  auc_av = sum_auc/rot
  x <- data.frame(
    Combination <- desc,
    Overall_BER = min(perfobj$AveragedPredict.error.rate["Overall.BER",]),
    NR_BER = min(perfobj$AveragedPredict.error.rate["yes",]),
    R_BER = min(perfobj$AveragedPredict.error.rate["no",]),

    AUC = auc_av

    )
  x
}

extract_metrics_tumour <- function(perfobj, desc){
  perfobj = perf(perfobj, validation = 'loo', auc = T) 
  x <- data.frame(
    Combination <- desc,
    Overall_BER = min(perfobj$error.rate$BER),
    NR_BER = perfobj$error.rate.class$max.dist[2],
    R_BER = perfobj$error.rate.class$max.dist[1],
    AUC = perfobj$auc[[1]][[1]][1]
    #`p-value` = perfobj$auc.all[[1]][[1]][2]
  )
  x
}

list_res <- list()
list_res[[1]] <- extract_metrics(tp1, "Tumour + 0")
list_res[[2]] <- extract_metrics(tp12, "Tumour + 0 + 1")
list_res[[3]] <- extract_metrics(tp123, "Tumour + 0 + 1 + 2")
list_res[[4]] <- extract_metrics(p123, "0 + 1 + 2")
list_res[[5]] <- extract_metrics_tumour(tonly, "Tumour only")


list_res[[6]] <- extract_metrics_tumour(p1, "0")
list_res[[7]] <- extract_metrics_tumour(p2, "1")
list_res[[8]] <- extract_metrics_tumour(p3, "2")
list_res[[9]] <- extract_metrics(p12, "0+1")
list_res[[10]] <- extract_metrics(p23, "1+2")
list_res[[11]] <- extract_metrics(p13, "0+3")
list_res[[12]] <- extract_metrics(tp2, "Tumour + 1")
list_res[[13]] <- extract_metrics(tp3, "Tumour + 2")
list_res[[14]] <- extract_metrics(tp23, "Tumour + 1 + 2")
list_res[[15]] <- extract_metrics(tp13, "Tumour + 0 + 2")

tabs <- do.call(rbind, list_res)
rownames(tabs) <- tabs$Combination....desc
tabs1 <- round(tabs[,-1], 3)
tabs1$Combinations <- tabs$Combination....desc
tabs1 <- tabs1[,c("Combinations", "R_BER", "NR_BER", "Overall_BER")]

saveRDS(tabs1, "results/allscores.rds")

tabs1 <- readRDS("results/allscores.rds")

#tabs1 <- tabs1[c(5, 4, 1, 2, 3),]

#tabs1

print("making table")

library(gt)

gt_tbl <- gt(data.frame(tabs1), rowname_col = "Combinations")

gt_tbl <- 
  gt_tbl |>
  tab_stubhead(label = "Assay combinations")

#gt_tbl <- 
#  gt_tbl |>
#  tab_spanner(
#    label = "Model AUC",
#    columns = c(AUC, p.value)
#  ) |>
#  tab_spanner(
#    label = "Class-specific BER",
#    columns = c(R_BER, NR_BER)
#  )

gt_tbl <- 
  gt_tbl |>
  tab_spanner(
    label = "Class-specific BER",
    columns = c(R_BER, NR_BER)
  )


gt_tbl

gt_tbl <- 
  gt_tbl |>
  cols_label(
    R_BER = html("Responder<sub>BER</sub>"),
    NR_BER = html("Non-Responder<sub>BER</sub>"),
    Overall_BER = html("Overall<sub>BER</sub>"),
    #AUC = html("AUC<sub>value</sub>"),
    #p.value = html("p-value")
  )



gt_tbl <- 
  gt_tbl |> 
  tab_row_group(
    label = md("**Multi-platform (blood and tumour combined)**"),
    rows = c(1, 2, 3, 12, 13, 14, 15)
  ) |>
  tab_row_group(
    label = md("**Single-platform (blood or tumour only)**"),
    rows = c(5, 6, 7, 8, 9, 10, 11, 4)
  )


gt_tbl

gtsave(gt_tbl, "results/performance_table.html")
gtsave(gt_tbl, "results/performance_table.docx")
