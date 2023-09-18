

library(gtable)
library(gt)
library("tidyr")
library("dplyr")

tp123_pcgse <- readRDS("results/tp123_pcgse.rds")
tp12_pcgse <- readRDS("results/tp12_pcgse.rds")
tp1_pcgse <- readRDS("results/tp1_pcgse.rds")
p123_pcgse <- readRDS("results/p123_pcgse.rds")

tptest <- tp123_pcgse[,-1]
data_wide <- spread(tptest, gene_set_names, FDR)
data_wide$combination <- "Multi-platform (Tumour + 0 + 1 + 2)"
#data_wide$assay <- c("Timepoint 0", "Timepoint 1", "Timepoint 2")

tptest1 <- tp12_pcgse[,-1]
data_wide1 <- spread(tptest1, gene_set_names, FDR)
data_wide1$combination <- "Multi-platform (Tumour + 0 + 1)"
#data_wide1$assay <- c("Timepoint 0", "Timepoint 1")

tptest2 <- p123_pcgse[,-1]
data_wide2 <- spread(tptest2, gene_set_names, FDR)
data_wide2$combination <- "Single-platform (0 + 1 + 2)"
#data_wide2$assay <- c("Timepoint 0", "Timepoint 1", "Timepoint 2")

tptest3 <- tp1_pcgse[,-1]
data_wide3 <- spread(tptest3, gene_set_names, FDR)
data_wide3$combination <- "Multi-platform (Tumour + 0)"
#data_wide3$assay <- c("Timepoint 0")

df_list <- list(data_wide, data_wide1, data_wide2, data_wide3)

df <- do.call(rbind, df_list)

df <- df[grepl("blood", df$assay, fixed = TRUE) & df$component_index == 1,]

df <- df[,-c(2,3)]

#df <- df[,c(1,2,4,3,5)]

#df <- group_by(df, "combination")

gt_tbl <- df |>
  group_by(combination) |>
  gt(rowname_col = "assay")



gt_tbl


gt_tbl <- 
  gt_tbl |>
  tab_stubhead(label = "Timepoint")

gt_tbl



gt_tbl

gt_tbl <- 
  gt_tbl |>
  tab_spanner(
    label = "Signature enrichment (FDR)",
    columns = c(PACE, TPEXhi, TSTEMhi))


gt_tbl <- gt_tbl |>
  cols_label(
    PACE = md("Pace et al."),
    TPEXhi = md("T<sub>PEX</sub>"),
    TSTEMhi = md("T<sub>STEM</sub>")
  )

gt_tbl <- gt_tbl |>
  tab_source_note(
    source_note = "1. Enrichment columns describe statistical significance (FDR) of association between first PC and gene set"
  ) |>
  tab_source_note(
    source_note = "2. Statistical significance indicates that gene sets have discriminative capacity between phenotypes."
  ) |>
  tab_source_note(
    source_note = "3. Timepoint column describes peripheral blood assay at each timepoint in each assay combination"
  )


gt_tbl

gt_tbl <- 
  gt_tbl |> 
  tab_row_group(
    label = md("**Multi-platform (Tumour + 0 + 1 + 2)**"),
    rows = 1:3
  ) |>
  tab_row_group(
    label = md("**Multi-platform (Tumour + 0 + 1)**"),
    rows = 4:5
  ) |>
  tab_row_group(
    label = md("**Single-platform (0 + 1 + 2)**"),
    rows = 6:8
  ) |>
  tab_row_group(
    label = md("**Multi-platform (Tumour + 0)**"),
    rows = 9
  )

gt_tbl

gtsave(gt_tbl, "results/stem_enrichment_table.html")
gtsave(gt_tbl, "results/stem_enrichment_table.docx")
