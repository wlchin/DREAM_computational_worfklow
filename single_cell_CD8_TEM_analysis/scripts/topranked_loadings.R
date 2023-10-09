library(dplyr)

full_model_all_samples_tp12 <- readRDS("results/diagnostics/full_model_all_samples_tp12.rds")

ones <- selectVar(full_model_all_samples_tp12, block = "blood", comp = 1)
x <- data.frame(ones$blood$value)
x$assay <- "Timepoint0"
#write(blood, "results/test_blood_tp12.txt")

ones <- selectVar(full_model_all_samples_tp12, block = "blood2", comp = 1)
x2 <- data.frame(ones$blood$value)
x2$assay <- "Timepoint1"
#write(blood, "results/test_blood2_tp12.txt")

totalx <- rbind(x, x2)

final_loadings_tp12df <- arrange(totalx, value.var)
write.csv(final_loadings_tp12df, "results/top_rank_loadings_df.csv")








tail(sort(full_model_all_samples_tp12$loadings$blood[,1]), 10)
head(sort(full_model_all_samples_tp12$loadings$blood[,1]), 10)
tail(sort(full_model_all_samples_tp12$loadings$blood2[,1]), 10)
head(sort(full_model_all_samples_tp12$loadings$blood2[,1]), 10)

tail(sort(full_model_all_samples_tp12$loadings$blood[,1]), 10)
head(sort(full_model_all_samples_tp12$loadings$blood[,1]), 10)
tail(sort(full_model_all_samples_tp12$loadings$blood2[,1]), 10)

full_model_all_samples_tp1 <- readRDS("results/diagnostics/full_model_all_samples_tp1.rds")

tail(sort(full_model_all_samples_tp1$loadings$blood[,1]), 10)
head(sort(full_model_all_samples_tp1$loadings$blood[,1]), 10)

sum(full_model_all_samples_tp1$loadings$blood[,1] > 0)
sum(full_model_all_samples_tp1$loadings$blood[,1] < 0)


sum(full_model_all_samples_tp12$loadings$blood2[,1] > 0)
sum(full_model_all_samples_tp12$loadings$blood[,1] < 0)
sum(full_model_all_samples_tp12$loadings$blood[,1] > 0)
