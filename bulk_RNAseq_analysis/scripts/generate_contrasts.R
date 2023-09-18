library('variancePartition')
library('edgeR')
library('BiocParallel')
library('ggplot2')

fit <- readRDS('results/fit_obj.rds')


coef_to_test <- colnames(fit)

output_to_csv <- function(fit_obj, test_coeff, fname){
  res <- topTable(fit_obj, coef=test_coeff, number = "inf")
  write.csv(res, fname)
}

for(i in coef_to_test){  
    fname_csv <- paste0("DE_results/", i, ".csv")
    output_to_csv(fit, i, fname_csv)
    
}


