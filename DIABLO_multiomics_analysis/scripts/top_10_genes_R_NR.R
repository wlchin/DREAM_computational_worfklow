# this is the gene plots

library(dplyr)
library(cowplot)
library(ggplot2)

df <- data.frame(df1)

df_res <- head(arrange(df, comp1), 10)
df_res$genes <- rownames(df_res)
df_res <- df_res %>%
  mutate(genes = reorder(genes, -comp1))

df_nonres <- tail(arrange(df, comp1), 10)
df_nonres$genes <- rownames(df_nonres)
df_nonres <- df_nonres %>%
  mutate(genes = reorder(genes, comp1))

x3 <- tail(arrange(x3, value.var), 10)
df_nonres$genes <- rownames(df_nonres)
df_nonres <- df_nonres %>%
  mutate(genes = reorder(genes, comp1))

x3_res <- head(arrange(x3, value.var), 10)
x3_res$genes <- rownames(x3_res)
x3_res <- x3_res %>%
  mutate(genes = reorder(genes, -value.var))

x3_nonres <- tail(arrange(x3, value.var), 10)
x3_nonres$genes <- rownames(x3_nonres)
x3_nonres <- x3_nonres %>%
  mutate(genes = reorder(genes, value.var))

p<-ggplot(data=df_res, aes(x=-(comp1), y=genes)) +
  geom_bar(stat="identity") + xlim(0, 0.4)
p + theme_cowplot() + xlab("loadings") + ylab("Blood: top 10 in R")

p<-ggplot(data=df_nonres, aes(x=comp1, y=genes)) +
  geom_bar(stat="identity") + xlim(0, 0.4)
p + theme_cowplot() + xlab("loadings") + ylab("Blood: top 10 in NR")

p<-ggplot(data=x3_res, aes(x=-(value.var), y=genes)) +
  geom_bar(stat="identity") + xlim(0, 0.4)
p + theme_cowplot() + xlab("loadings") + ylab("Tumour: top 10 in R")

p<-ggplot(data=x3_nonres, aes(x=value.var, y=genes)) +
  geom_bar(stat="identity") + xlim(0, 0.4)
p + theme_cowplot() + xlab("loadings") + ylab("Tumour: top 10 in NR")

#3.48 x 3.80
