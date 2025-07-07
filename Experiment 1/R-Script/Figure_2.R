

library(ggplot2)
library(dplyr)

mean_surprisals <- read_csv("mean_surprisals.csv")

plot_whole <- ggplot(mean_surprisals %>%
                       group_by(Condition, deletion_rate) %>%
                       summarise(Surprisal = mean(Surprisal)) %>%
                       filter(deletion_rate != 0.25),
                     aes(x = 1 - deletion_rate, y = Surprisal, color = as.factor(Condition), group = Condition)) +
  geom_line(size = 1.2) +  
  scale_color_manual(values = c("a" = "#1f77b4", "b" = "#ff7f0e", "c" = "#2ca02c", "d" = "#d62728")) +  # Custom colors
  labs(x = "Context Sizes", y = "Mean Surprisal", color = "Condition") +  # Axis labels and legend title
  scale_x_continuous(
    breaks = c(0.75, 0.5, 0.25),  
    labels = c(15, 10, 5)         
  ) +
  theme_minimal(base_size = 14) +  
  theme(
    legend.position = "top",  
    legend.title = element_text(face = "bold"),  
    legend.text = element_text(size = 12),  
    axis.title = element_text(face = "bold"),  
    axis.text = element_text(size = 12),  
    panel.grid.major = element_line(color = "gray", size = 0.5),  
    panel.grid.minor = element_blank(),  
    strip.background = element_blank(),  
    strip.text = element_text(face = "bold")  
  ) +
  guides(color = guide_legend(title.position = "top", title.hjust = 0.5))  


print(plot_whole)
ggsave("plot_whole.png", plot = plot_whole, width = 8, height = 6, dpi = 300)