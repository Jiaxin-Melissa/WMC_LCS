

# Load necessary libraries
library(dplyr)
library(ggplot2)

# Read the anonymized dataset
human_surprisal_data_anon <- read_csv("human_surprisal_data_anon.csv")


# Compute participant-level summary statistics:
# - mean_accuracy: proportion of trials where the response was "yes"
# - mean_RT: average reading time to the first answer, excluding missing values
summary_data <- human_surprisal_data_anon %>%
    group_by(subject_ID) %>%
    summarise(
        mean_accuracy = mean(Correct == "yes", na.rm = TRUE),  
        mean_RT = mean(Reading.time.to.first.answer, na.rm = TRUE)  
    )

# Create a scatter plot of mean accuracy vs. mean RT for each participant
plot <- ggplot(summary_data, aes(x = mean_accuracy, y = mean_RT)) +
    geom_point(color = "blue", size = 3, alpha = 0.7) +  
    labs(
        x = "Mean Accuracy (Fraction)", 
        y = "Mean RT (ms)", 
        title = "Mean Accuracy vs. Mean RT (Reading.time.to.first.answer) by Participant"
    ) +
    theme_minimal(base_size = 14) +  
    theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),  
        axis.title = element_text(size = 14),  
        axis.text = element_text(size = 12),  
        panel.grid.major = element_line(color = "gray90", size = 0.5),  
        panel.grid.minor = element_blank(),  
        plot.margin = margin(10, 10, 10, 10)  
    ) +
    geom_smooth(method = "lm", color = "red", linetype = "dashed")  

# Display the plot
print(plot)

# save the plot
ggsave("mean_accuracy_vs_mean_RT_plot.png", plot = plot, width = 8, height = 6, dpi = 300)



####desriptive statistics for  coloumn "mean_RT" in dataset "summary_data" that displayed in "Results of A-Maze task"

# Remove NA values from mean_RT column
rt_values <- summary_data$mean_RT[!is.na(summary_data$mean_RT)]

# Calculate summary statistics
mean_rt <- mean(rt_values)
sd_rt <- sd(rt_values)
median_rt <- median(rt_values)

# Calculate 95% confidence interval using normal approximation
n <- length(rt_values)
error_margin <- qnorm(0.975) * sd_rt / sqrt(n)
ci_lower <- mean_rt - error_margin
ci_upper <- mean_rt + error_margin

# Print results
cat("Mean RT:", mean_rt, "\n")
cat("SD RT:", sd_rt, "\n")
cat("Median RT:", median_rt, "\n")
cat("95% CI for RT: [", ci_lower, ",", ci_upper, "]\n")


####Pearson correlation test between "mean_RT" and "mean_accuracy" in dataset "summary_data"

cor.test(summary_data$mean_accuracy, summary_data$mean_RT, method = "pearson")