


library(ggplot2)

data <- read.csv("results/combined_weights.csv")
data <- data[,c(3,2)]
colnames(data) <- c("Category", "Value")

data$Category[duplicated(data$Category)] <- paste0(data$Category[duplicated(data$Category)], "_2")

data <- data[data$Value != 0, ]

data$Category <- factor(data$Category, levels = data$Category[order(-data$Value)])
# Create a bar chart
ggplot(data, aes(x = Category, y = Value, fill = Value > 0)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("positive" = "green", "negative" = "red")) +
  labs(title = "Average feature importance across 100 runs (Timepoint 0 and Timepoint 1)",
       x = "Genes (features)",
       y = "Average loading score across 100 runs") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),  # Remove x-axis ticks
    panel.grid.major.x = element_blank(),  # Remove major x-axis gridlines
    panel.grid.minor.x = element_line(color = "gray", size = 0.2)  # Add minor x-axis gridlines
  ) + 
  scale_x_discrete(breaks = unique(data$Category)[seq(1, length(unique(data$Category)), by = 50)]) + 
  coord_cartesian(ylim = c(-0.2, max(data$Value))) + 
  geom_vline(xintercept = seq(25, length(unique(data$Category)), by = 50),
             linetype = "dashed", color = "gray")  # Add vertical gridlines
ggsave("results/feature_plot_with_dashes_blood.png", dpi = 600, width = 7.84, height = 3.07, units = "in")




data <- read.csv("results/tumour_weights.csv")
data <- data[,c(3,2)]
colnames(data) <- c("Category", "Value")

data$Category[duplicated(data$Category)] <- paste0(data$Category[duplicated(data$Category)], "_2")

data <- data[data$Value != 0, ]

data$Category <- factor(data$Category, levels = data$Category[order(-data$Value)])
# Create a bar chart
ggplot(data, aes(x = Category, y = Value, fill = Value > 0)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("positive" = "green", "negative" = "red")) +
  labs(title = "Average feature importance across 100 runs (Tumour)",
       x = "Genes (features)",
       y = "Average loading score across 100 runs") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),  # Remove x-axis ticks
    panel.grid.major.x = element_blank(),  # Remove major x-axis gridlines
    panel.grid.minor.x = element_line(color = "gray", size = 0.2)  # Add minor x-axis gridlines
  ) + 
  scale_x_discrete(breaks = unique(data$Category)[seq(1, length(unique(data$Category)), by = 50)]) + 
  coord_cartesian(ylim = c(-0.2, max(data$Value))) + 
  geom_vline(xintercept = seq(25, length(unique(data$Category)), by = 50),
             linetype = "dashed", color = "gray")  # Add vertical gridlines
ggsave("results/feature_plot_with_dashes_tumour.png", dpi = 600, width = 7.84, height = 3.07, units = "in")


