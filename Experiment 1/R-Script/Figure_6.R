###########combined plot#####

library(dplyr)
library(readr)

# Read the anonymized human surprisal data
human_surprisal_data_anon <- read_csv("human_surprisal_data_anon.csv")

# Extract unique subject_ID and corresponding ConsonantsAccuracy_original
consonants_accuracy_data <- human_surprisal_data_anon %>%
    select(subject_ID, ConsonantsAccuracy_original) %>%
    distinct()

# Extract unique subject_ID and corresponding MeaningfulnessAccuracy_original
Meaningful_accuracy_data <- human_surprisal_data_anon %>%
    select(subject_ID, MeaningfulnessAccuracy_original) %>%
    distinct()




library(gridExtra)

# Create histogram for Consonants Accuracy
plot1 <- ggplot(consonants_accuracy_data, aes(x = ConsonantsAccuracy_original)) +
    geom_histogram(binwidth = 0.05, fill = "skyblue", color = "black", na.rm = TRUE) +
    geom_vline(aes(xintercept = mean(ConsonantsAccuracy_original, na.rm = TRUE)), 
               color = "red", linetype = "dashed", linewidth = 1) +
    geom_vline(aes(xintercept = median(ConsonantsAccuracy_original, na.rm = TRUE)), 
               color = "blue", linetype = "dotted", linewidth = 1) +
    labs(
        title = "Distribution of Consonants Accuracy",
        x = "Consonants Accuracy",
        y = "Frequency"
    ) +
    annotate("text", x = mean(consonants_accuracy_data$ConsonantsAccuracy_original, na.rm = TRUE) + 0.02, 
             y = 40, 
             label = paste("Mean =", round(mean(consonants_accuracy_data$ConsonantsAccuracy_original, na.rm = TRUE), 3)), 
             color = "red", size = 4) +
    annotate("text", x = median(consonants_accuracy_data$ConsonantsAccuracy_original, na.rm = TRUE) - 0.05, 
             y = 35, 
             label = paste("Median =", round(median(consonants_accuracy_data$ConsonantsAccuracy_original, na.rm = TRUE), 3)), 
             color = "blue", size = 4) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
    scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.05))) +
    theme_minimal() +
    theme(
        panel.grid.major = element_line(color = "grey90"),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    )

# Create histogram for Meaningfulness Accuracy
plot2 <- ggplot(Meaningful_accuracy_data, aes(x = MeaningfulnessAccuracy_original)) +
    geom_histogram(binwidth = 0.05, fill = "skyblue", color = "black", na.rm = TRUE) +
    geom_vline(aes(xintercept = mean(MeaningfulnessAccuracy_original, na.rm = TRUE)), 
               color = "red", linetype = "dashed", linewidth = 1) +
    geom_vline(aes(xintercept = median(MeaningfulnessAccuracy_original, na.rm = TRUE)), 
               color = "blue", linetype = "dotted", linewidth = 1) +
    labs(
        title = "Distribution of Meaningfulness Accuracy",
        x = "Meaningfulness Accuracy",
        y = "Frequency"
    ) +
    annotate("text", x = mean(Meaningful_accuracy_data$MeaningfulnessAccuracy_original, na.rm = TRUE) + 0.02, 
             y = 40, 
             label = paste("Mean =", round(mean(Meaningful_accuracy_data$MeaningfulnessAccuracy_original, na.rm = TRUE), 3)), 
             color = "red", size = 4) +
    annotate("text", x = median(Meaningful_accuracy_data$MeaningfulnessAccuracy_original, na.rm = TRUE) - 0.05, 
             y = 35, 
             label = paste("Median =", round(median(Meaningful_accuracy_data$MeaningfulnessAccuracy_original, na.rm = TRUE), 3)), 
             color = "blue", size = 4) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
    scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.05))) +
    theme_minimal() +
    theme(
        panel.grid.major = element_line(color = "grey90"),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    )

# Combine the two plots side by side with proper space for annotations
grid.arrange(plot1, plot2, ncol = 2, widths = c(1, 1))  # Ensure both plots have equal width
# Save the combined plot
ggsave("Combined_Plot.png", plot = grid.arrange(plot1, plot2, ncol = 2, widths = c(1, 1)), 
       width = 14, height = 7, dpi = 300)




# Descriptive statistics for MeaningfulnessAccuracy_original
mean_val   <- mean(Meaningful_accuracy_data$MeaningfulnessAccuracy_original, na.rm = TRUE)
sd_val     <- sd(Meaningful_accuracy_data$MeaningfulnessAccuracy_original, na.rm = TRUE)
n          <- sum(!is.na(Meaningful_accuracy_data$MeaningfulnessAccuracy_original))
se_val     <- sd_val / sqrt(n)
range_vals <- range(Meaningful_accuracy_data$MeaningfulnessAccuracy_original, na.rm = TRUE)
median_val <- median(Meaningful_accuracy_data$MeaningfulnessAccuracy_original, na.rm = TRUE)

# Print results with specific rounding
cat("Mean:", round(mean_val, 3), "\n")
cat("SD:", round(sd_val, 3), "\n")
cat("SE:", round(se_val, 3), "\n")
cat("Range:", paste0("[", format(round(range_vals[1], 3), nsmall = 3), ", ", 
                     format(round(range_vals[2], 3), nsmall = 3), "]"), "\n")
cat("Median:", round(median_val, 3), "\n")


# Descriptive statistics for ConsonantsAccuracy_original
mean_val   <- mean(consonants_accuracy_data$ConsonantsAccuracy_original, na.rm = TRUE)
sd_val     <- sd(consonants_accuracy_data$ConsonantsAccuracy_original, na.rm = TRUE)
n          <- sum(!is.na(consonants_accuracy_data$ConsonantsAccuracy_original))
se_val     <- sd_val / sqrt(n)
range_vals <- range(consonants_accuracy_data$ConsonantsAccuracy_original, na.rm = TRUE)
median_val <- median(consonants_accuracy_data$ConsonantsAccuracy_original, na.rm = TRUE)

# Print results with 3 decimal points, including range
cat("Mean:", round(mean_val, 3), "\n")
cat("SD:", round(sd_val, 3), "\n")
cat("SE:", round(se_val, 3), "\n")
cat("Range:", paste0("[", format(round(range_vals[1], 3), nsmall = 3), ", ", 
                     format(round(range_vals[2], 3), nsmall = 3), "]"), "\n")
cat("Median:", round(median_val, 3), "\n")


# Merge the two datasets by subject_ID
merged_accuracy_data <- merge(
  consonants_accuracy_data,
  Meaningful_accuracy_data,
  by = "subject_ID",
  all = TRUE  
)

#Pearson correlation test between "ConsonantsAccuracy_original" and "MeaningfulnessAccuracy_original" in dataset "merged_accuracy_data"

cor.test(merged_accuracy_data$ConsonantsAccuracy_original, merged_accuracy_data$MeaningfulnessAccuracy_original, method = "pearson")