
library(dplyr)

read_models_and_plot_upset<- function(reps){
  list_of_models <- list()
  for(i in 1:reps){
    filename <- paste0("results/model_", i, ".rds")
    list_of_models[[i]] <- readRDS(filename)
  }
  list_of_models
}

# I need to plot an average factor loading for each one. 
no_norm <- read_models_and_plot_upset(100)

saveRDS(no_norm, "data/tp12_sampled.rds")

get_multiblock_matrix <- function(nreps, model_list, block){
  loading_mat <- list()
  for(i in 1:nreps){
    col_ind <- model_list[[i]]$loadings[[block]][,1]
    loading_mat[[i]] <- col_ind
  }
  loading_heatmap <- do.call(rbind, loading_mat)
  loading_heatmap
}

mat_tumour <- get_multiblock_matrix(100, no_norm, "tumour")
tumour <- data.frame(colMeans(mat_tumour))
tumour$features <- rownames(tumour)
tumour$timepoint <- "Timepoint0"
colnames(tumour) <- c("weights", "features", "timepoint")
tumour <- tumour %>% arrange(desc(weights))
write.csv(tumour, "results/tumour_weights.csv")
write(unique(tail(tumour, 25)$features), "results/tumour_responder_sampled_25genes.txt")
write(unique(head(tumour, 25)$features), "results/tumour_nonresponder_sampled_25genes.txt")
write(unique(tail(tumour, 50)$features), "results/tumour_responder_sampled_50genes.txt")
write(unique(head(tumour, 50)$features), "results/tumour_nonresponder_sampled_50genes.txt")
write(unique(tail(tumour, 100)$features), "results/tumour_responder_sampled_100genes.txt")
write(unique(head(tumour, 100)$features), "results/tumour_nonresponder_sampled_100genes.txt")
write(unique(tumour$features), "results/gene_tumour_background_new_seed.txt")

mat_blood1 <- get_multiblock_matrix(100, no_norm, "blood_tpt0")
blood1 <- data.frame(colMeans(mat_blood1))
blood1$features <- rownames(blood1)
blood1$timepoint <- "Timepoint0"
colnames(blood1) <- c("weights", "features", "timepoint")
blood1 <- blood1 %>% arrange(desc(weights))

mat_blood2 <- get_multiblock_matrix(100, no_norm, "blood_tpt1")
blood2 <- data.frame(colMeans(mat_blood2))
blood2$features <- rownames(blood2)
blood2$timepoint <- "Timepoint1"
colnames(blood2) <- c("weights", "features", "timepoint")
blood2 <- blood2 %>% arrange(desc(weights))

dftotal <- rbind(blood1, blood2)
dftotal <- dftotal %>% arrange(desc(weights))


write(unique(tail(dftotal, 25)$features), "results/blood_responder_sampled_25genes.txt")
write(unique(head(dftotal, 25)$features), "results/blood_nonresponder_sampled_25genes.txt")

write(unique(tail(dftotal, 50)$features), "results/blood_responder_sampled_50genes.txt")
write(unique(head(dftotal, 50)$features), "results/blood_nonresponder_sampled_50genes.txt")

write(unique(tail(dftotal, 100)$features), "results/blood_responder_sampled_100genes.txt")
write(unique(head(dftotal, 100)$features), "results/blood_nonresponder_sampled_100genes.txt")

write.csv(dftotal, "results/combined_weights.csv")

write(unique(dftotal$features), "results/gene_blood_background_new_seed.txt")

topdf <- rbind(head(dftotal, 100), tail(dftotal, 100))
write.csv(topdf, "results/top_rank_loadings_df_blood.csv")
