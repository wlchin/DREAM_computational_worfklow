library(dplyr)
library(mixOmics)

full_model_all_samples_tp12 <- readRDS("results/diagnostics/full_model_all_samples_tp12.rds")

ones <- selectVar(full_model_all_samples_tp12, block = "blood", comp = 1)
x <- data.frame(ones$blood$value)
x$assay <- "Timepoint0"
write(rownames(x), "results/test_blood_tp12.txt")

ones <- selectVar(full_model_all_samples_tp12, block = "blood2", comp = 1)
x2 <- data.frame(ones$blood$value)
x2$assay <- "Timepoint1"
write(rownames(x2), "results/test_blood2_tp12.txt")

totalx <- rbind(x, x2) # combine two assays

final_loadings_tp12df <- arrange(totalx, value.var)
write.csv(final_loadings_tp12df, "results/top_rank_loadings_df_blood.csv")
#saveRDS(final_loadings_tp12df, "results/top_rank_loadings_df_blood.rds")

write(unique(rownames(totalx[totalx$value.var < 0,])), "results/blood_responder.txt")
write(unique(rownames(totalx[totalx$value.var > 0,])), "results/blood_nonresponder.txt")

ones <- selectVar(full_model_all_samples_tp12, block = "tumour", comp = 1)
x3 <- data.frame(ones$tumour$value)
x3$assay <- "Tumour"
write.csv(x3, "results/top_rank_loadings_df_tumour.csv")
#saveRDS(x3, "results/top_rank_loadings_df_tumour.rds")

write(unique(rownames(x3[x3$value.var < 0,])), "results/tumour_responder.txt")
write(unique(rownames(x3[x3$value.var > 0,])), "results/tumour_nonresponder.txt")

