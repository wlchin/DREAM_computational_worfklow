# load data

library("survminer")
require("survival")
library("ggplot2")

filename <- snakemake@input[[1]]

fulldf2 <- read.csv(filename)

#res.cox.pfs <- coxph(Surv(adj_time_ir, 
#                      adjir_status) ~ bin_surv, data = fulldf2)

fit <- survfit(Surv(adj_time_ir, adjir_status) ~ bin_surv, data = fulldf2)
ggsurvplot(fit, data = fulldf2, pval = TRUE) + ggtitle(filename)
ggsave(snakemake@output[[1]], width = 7.11, height = 4.79, dpi = 300)
df <- surv_pvalue(fit)
df$file <- filename
df$type <- "PFS"
write.table(df, "results/pval_summary.csv", append = TRUE, col.names = F, sep = ",")

fit <- survfit(Surv(os_time, os_status) ~ bin_surv, data = fulldf2)
ggsurvplot(fit, data = fulldf2, pval = TRUE) + ggtitle(filename)
ggsave(snakemake@output[[2]], width = 7.11, height = 4.79, dpi = 300)
df <- surv_pvalue(fit)
df$file <- filename
df$type <- "OS"
write.table(df, "results/pval_summary.csv", append = TRUE, col.names = F, sep = ",")
