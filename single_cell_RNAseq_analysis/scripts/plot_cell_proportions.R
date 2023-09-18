
library(ggplot2)
library(dplyr)
library(readr)
library(Polychrome)

#mypal <- kelly.colors(26)

set.seed(123)

lastpal <- createPalette(31, c("#010101", "#ff0000"), M=10000)
lastpal

x <- read_csv('results/cell_metadata.csv')

length(unique(x$predicted.celltype.l2))

Tpt0_PFS1 <- x[x$Timepoint == 0 & x$PFS_6M == 1,]
T0P1 <- Tpt0_PFS1 %>% 
  group_by(predicted.celltype.l2) %>%
  summarise(n = n(), .groups = "drop") %>% 
  mutate(freq = n / sum(n))
T0P1$pheno <- "T0P1"

Tpt1_PFS1 <- x[x$Timepoint == 1 & x$PFS_6M == 1,]
T1P1 <- Tpt1_PFS1 %>% 
  group_by(predicted.celltype.l2) %>%
  summarise(n = n(), .groups = "drop") %>% 
  mutate(freq = n / sum(n))
T1P1$pheno <- "T1P1"

Tpt2_PFS1 <- x[x$Timepoint == 2 & x$PFS_6M == 1,]
T2P1 <- Tpt2_PFS1 %>% 
  group_by(predicted.celltype.l2) %>%
  summarise(n = n(), .groups = "drop") %>% 
  mutate(freq = n / sum(n))
T2P1$pheno <- "T2P1"

Tpt0_PFS0 <- x[x$Timepoint == 0 & x$PFS_6M == 0,]
T0P0 <- Tpt0_PFS0 %>% 
  group_by(predicted.celltype.l2) %>%
  summarise(n = n(), .groups = "drop") %>% 
  mutate(freq = n / sum(n))
T0P0$pheno <- "T0P0"

Tpt1_PFS0 <- x[x$Timepoint == 1 & x$PFS_6M == 0,]
T1P0 <- Tpt1_PFS0 %>% 
  group_by(predicted.celltype.l2) %>%
  summarise(n = n(), .groups = "drop") %>% 
  mutate(freq = n / sum(n))
T1P0$pheno <- "T1P0"

Tpt2_PFS0 <- x[x$Timepoint == 2 & x$PFS_6M == 0,]
T2P0 <- Tpt2_PFS0 %>% 
  group_by(predicted.celltype.l2) %>%
  summarise(n = n(), .groups = "drop") %>% 
  mutate(freq = n / sum(n))
T2P0$pheno <- "T2P0"

total_df <- do.call(rbind, list(T0P1, T1P1, T2P1, T0P0, T1P0, T2P0))

total_df$predicted.celltype.l2 <- as.factor(total_df$predicted.celltype.l2)

names(lastpal) <- unique(total_df$predicted.celltype.l2)

ggplot(total_df) + geom_col(aes(pheno, freq, fill = predicted.celltype.l2), color = "black") + 
  scale_fill_manual(values = lastpal, name="Cell types") + theme_classic() + 
  ylab("Cell proportions") + xlab("Phenotype") + theme(axis.title.x = element_text(size = 14), 
                                                       axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
                                                       axis.text.y = element_text(size = 12),
                                                       axis.title.y = element_text(size = 14),
                                                       legend.text=element_text(size=12))
ggsave("results/barplot.png", width = 6, height = 6, dpi = 600)
saveRDS(total_df, "results/props.rds")


colscheme <- data.frame(names(lastpal), lastpal)
colnames(colscheme) <- c("celltype", "colour")
write.csv(colscheme, "results/colscheme.csv", quote = F, row.names = F)


