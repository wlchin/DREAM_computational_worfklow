
library(mixOmics)
library(dplyr)

p123 <- readRDS("results/diagnostics/full_model_all_samples_p123.rds")
#auc_p123 <- perf(p123, auc = T, validation = "loo")

tp123 <- readRDS("results/diagnostics/full_model_all_samples_tp123.rds")
#auc_tp123 <- perf(tp123, auc = T, validation = "loo")

tp12 <- readRDS("results/diagnostics/full_model_all_samples_tp12.rds")
#auc_tp12 <- perf(tp12, auc = T, validation = "loo")

tp1 <- readRDS("results/diagnostics/full_model_all_samples_tp1.rds")
#auc_tp1 <- perf(tp1, auc = T)

tonly <- readRDS("results/diagnostics/full_model_all_samples_t.rds")
#auc_t <- perf(tonly, auc = T, validation = "loo")

#tp2 <- readRDS("results/diagnostics/full_model_all_samples_tp2.rds")
#auc_t <- perf(tonly, auc = T, validation = "loo")

#tp3 <- readRDS("results/diagnostics/full_model_all_samples_tp3.rds")
#auc_t <- perf(tonly, auc = T, validation = "loo")

#perf.tp2 = perf(tp2, validation = 'loo', auc = T) 
#perf.tp2$AveragedPredict.error.rate

#perf.tp3 = perf(tp3, validation = 'loo', auc = T) 
#perf.tp3$AveragedPredict.error.rate

perf.t = perf(tonly, validation = 'loo', auc = T) 
perf.t$error.rate$BER

perf.tp1 = perf(tp1, validation = "loo", auc = T) 
perf.tp1$AveragedPredict.error.rate

perf.tp12 = perf(tp12, validation = "loo", auc = T) 
perf.tp12$AveragedPredict.error.rate

perf.tp123 = perf(tp123, validation = "loo", auc = T) 
perf.tp123$AveragedPredict.error.rate

perf.p123 = perf(p123, validation = "loo", auc = T) 
perf.p123$AveragedPredict.error.rate


extract_stuff <- function(perfobj, desc){
  x <- data.frame(
    Combination <- desc,
    Overall_BER = min(perfobj$AveragedPredict.error.rate["Overall.BER",]),
    NR_BER = min(perfobj$AveragedPredict.error.rate["yes",]),
    R_BER = min(perfobj$AveragedPredict.error.rate["no",]),
    AUC = perfobj$auc[length(perfobj$auc)][[1]][1],
    `p-value` = perfobj$auc[length(perfobj$auc)][[1]][2]
    )
  x
}

extract_stuff_tumour <- function(perfobj, desc){
  x <- data.frame(
    Combination <- desc,
    Overall_BER = mean(perfobj$error.rate$BER),
    NR_BER = perfobj$error.rate.class$max.dist[2],
    R_BER = perfobj$error.rate.class$max.dist[1],
    AUC = perfobj$auc.all[[1]][[1]][1],
    `p-value` = perfobj$auc.all[[1]][[1]][2]
  )
  x
}

list_res <- list()
list_res[[1]] <- extract_stuff(perf.tp1, "Tumour + 0")
list_res[[2]] <- extract_stuff(perf.tp12, "Tumour + 0 + 1")
list_res[[3]] <- extract_stuff(perf.tp123, "Tumour + 0 + 1 + 2")
list_res[[4]] <- extract_stuff(perf.p123, "0 + 1 + 2")
list_res[[5]] <- extract_stuff_tumour(perf.t, "Tumour only")
#list_res[[6]] <- extract_stuff(perf.tp2, "Tumour + Tpt2")
#list_res[[7]] <- extract_stuff(perf.tp3, "Tumour + Tpt3")

tabs <- do.call(rbind, list_res)
rownames(tabs) <- tabs$Combination....desc
tabs1 <- round(tabs[,-1], 3)
tabs1$Combinations <- tabs$Combination....desc
tabs1 <- tabs1[,c("Combinations", "AUC", "p.value", "R_BER", "NR_BER", "Overall_BER")]
#tabs1 <- arrange(tabs1, desc(Overall_BER), AUC)

tabs1 <- tabs1[c(5, 4, 1, 2, 3),]

library(gt)

gt_tbl <- gt(data.frame(tabs1), rowname_col = "Combinations")

gt_tbl <- 
  gt_tbl |>
  tab_stubhead(label = "Assay combinations")

gt_tbl <- 
  gt_tbl |>
  tab_spanner(
    label = "Model AUC",
    columns = c(AUC, p.value)
  ) |>
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
    AUC = html("AUC<sub>value</sub>"),
    p.value = html("p-value")
  )

gt_tbl <- 
  gt_tbl |>
    tab_style(
        style = list(
        cell_fill(color = "white"),
        cell_text(weight = "bold")
      ),
      locations = cells_body(
        columns = Overall_BER,
        rows = c(F, F, F, F, F)
      )
    )

gt_tbl <- 
  gt_tbl |> 
  tab_row_group(
    label = md("**Single platform (blood or tumour only)**"),
    rows = 1:2
  ) |>
  tab_row_group(
    label = md("**Multi-platform (blood and tumour combined)**"),
    rows = 3:5
  )


gt_tbl

gtsave(gt_tbl, "results/performance_table.html")
gtsave(gt_tbl, "results/performance_table.docx")


