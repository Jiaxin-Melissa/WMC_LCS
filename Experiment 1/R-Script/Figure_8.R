####Generate the bottom plot of Figure 8####


# Load necessary libraries
library(dplyr)


#read the csv file
human_surprisal_data_anon <- read_csv("human_surprisal_data_anon.csv")


########wmc bin##########
# Extract posterior samples
posterior_samples <- as_draws_df(model1_brm)

# Extract random effects for each participant
random_effects <- ranef(model1_brm)$subject_ID

# Sort WMC values before binning
sorted_wmc <- sort(human_surprisal_data_anon$ConsonantsAccuracy, na.last = TRUE)
wmc_quantiles <- quantile(sorted_wmc, probs = c(0, 0.25, 0.50, 0.75, 1), na.rm = TRUE)

# Function to categorize WMC into bins
bin_wmc <- function(wmc) {
    if (wmc <= wmc_quantiles[2]) {
        return("0-25%")
    } else if (wmc <= wmc_quantiles[3]) {
        return("25-50%")
    } else if (wmc <= wmc_quantiles[4]) {
        return("50-75%")
    } else {
        return("75-100%")
    }
}

# Sort trial_order before selecting the first & last 10%
sorted_trial_order <- sort(human_surprisal_data_anon$trial_order, na.last = TRUE)
trial_order_quantiles <- quantile(sorted_trial_order, probs = c(0.10, 0.90), na.rm = TRUE)
trial_orders <- c(trial_order_quantiles[1], trial_order_quantiles[2])

# Function to categorize trial order
bin_trial_order <- function(trial_order) {
    if (trial_order <= trial_order_quantiles[1]) {
        return("first")
    } else {
        return("last")
    }
}

# Initialize result lists
results_b <- list()
results_c <- list()
results_d <- list()

# Loop through each participant
for (i in 1:nrow(human_surprisal_data_anon)) {
    # Get participant data
    prolific_id <- as.character(human_surprisal_data_anon$subject_ID[i])
    wmc <- human_surprisal_data_anon$ConsonantsAccuracy[i]
    wmc_bin <- bin_wmc(wmc)
    
    # Check if the participant exists in the random effects
    if (prolific_id %in% rownames(random_effects)) {
        # Extract random intercepts and slopes
        random_intercept <- random_effects[prolific_id, "Estimate", "Intercept"]
        random_bVSc <- random_effects[prolific_id, "Estimate", "bVSc"]
        random_cVSd <- random_effects[prolific_id, "Estimate", "cVSd"]
        random_trial_order <- random_effects[prolific_id, "Estimate", "trial_order"]
        random_trial_order_bVSc <- random_effects[prolific_id, "Estimate", "bVSc:trial_order"]
        random_trial_order_cVSd <- random_effects[prolific_id, "Estimate", "cVSd:trial_order"]
    } else {
        stop(paste("subject_ID not found in random_effects:", prolific_id))
    }
    
    # Loop through trial order conditions
    for (trial_order in trial_orders) {
        trial_order_name <- bin_trial_order(trial_order)
        
        # ====================
        # Condition B
        # ====================
        bc <- -1
        result_name_b <- paste0(prolific_id, "_", wmc_bin, "_", trial_order_name, "_b")
        results_b[[result_name_b]] <- exp(
            posterior_samples[["b_Intercept"]] +
                random_intercept +
                bc * (posterior_samples[["b_bVSc"]] + random_bVSc) +
                trial_order * (posterior_samples[["b_trial_order"]] + random_trial_order) +
                trial_order * bc * (posterior_samples[["b_trial_order:bVSc"]] + random_trial_order_bVSc) +
                wmc * posterior_samples[["b_ConsonantsAccuracy"]] +
                wmc * trial_order * (posterior_samples[["b_ConsonantsAccuracy:trial_order"]] + random_trial_order) +
                wmc * bc * trial_order * (posterior_samples[["b_ConsonantsAccuracy:trial_order:bVSc"]] + random_trial_order_bVSc)
        )
        
        # ====================
        # Condition C
        # ====================
        bc <- 1
        cd <- -1
        result_name_c <- paste0(prolific_id, "_", wmc_bin, "_", trial_order_name, "_c")
        results_c[[result_name_c]] <- exp(
            posterior_samples[["b_Intercept"]] +
                random_intercept +
                bc * (posterior_samples[["b_bVSc"]] + random_bVSc) +
                cd * (posterior_samples[["b_cVSd"]] + random_cVSd) +
                trial_order * (posterior_samples[["b_trial_order"]] + random_trial_order) +
                trial_order * bc * (posterior_samples[["b_trial_order:bVSc"]] + random_trial_order_bVSc) +
                trial_order * cd * (posterior_samples[["b_trial_order:cVSd"]] + random_trial_order_cVSd) +
                wmc * posterior_samples[["b_ConsonantsAccuracy"]] +
                wmc * trial_order * (posterior_samples[["b_ConsonantsAccuracy:trial_order"]] + random_trial_order) +
                wmc * bc * trial_order * (posterior_samples[["b_ConsonantsAccuracy:trial_order:bVSc"]] + random_trial_order_bVSc) +
                wmc * cd * trial_order * (posterior_samples[["b_ConsonantsAccuracy:trial_order:cVSd"]] + random_trial_order_cVSd)
        )
        
        # ====================
        # Condition D
        # ====================
        cd <- 1
        result_name_d <- paste0(prolific_id, "_", wmc_bin, "_", trial_order_name, "_d")
        results_d[[result_name_d]] <- exp(
            posterior_samples[["b_Intercept"]] +
                random_intercept +
                cd * (posterior_samples[["b_cVSd"]] + random_cVSd) +
                trial_order * (posterior_samples[["b_trial_order"]] + random_trial_order) +
                trial_order * cd * (posterior_samples[["b_trial_order:cVSd"]] + random_trial_order_cVSd) +
                wmc * posterior_samples[["b_ConsonantsAccuracy"]] +
                wmc * trial_order * (posterior_samples[["b_ConsonantsAccuracy:trial_order"]] + random_trial_order) +
                wmc * cd * trial_order * (posterior_samples[["b_ConsonantsAccuracy:trial_order:cVSd"]] + random_trial_order_cVSd)
        )
    }
}

# Print results
cat("Results for condition b:\n")
print(results_b)
cat("\nResults for condition c:\n")
print(results_c)
cat("\nResults for condition d:\n")
print(results_d)




## Load necessary package
library(tidyr)

# Initialize an empty data frame to store mean RT values
mean_rt_data <- data.frame(
    subject_ID = character(),
    WMC_Bin = character(),  # Use WMC Bin directly from result_name
    trial_order = character(),
    condition = character(),
    Mean_RT = numeric(),
    stringsAsFactors = FALSE
)

# Function to extract details from the result name
extract_details <- function(result_name) {
    parts <- unlist(strsplit(result_name, "_"))
    list(
        subject_ID = parts[1],
        WMC_Bin = parts[2],  # Store bin directly instead of numeric conversion
        trial_order = parts[3],
        condition = parts[4]
    )
}

# Helper function to process results and append them to mean_rt_data
process_results <- function(results_list) {
    for (result_name in names(results_list)) {
        details <- extract_details(result_name)
        
        mean_rt_data <<- rbind(mean_rt_data, data.frame(
            subject_ID = details$subject_ID,
            WMC_Bin = details$WMC_Bin,  # Use extracted bin directly
            trial_order = details$trial_order,
            condition = details$condition,
            Mean_RT = mean(results_list[[result_name]], na.rm = TRUE)
        ))
    }
}

# Process results for each condition
process_results(results_b)
process_results(results_c)
process_results(results_d)

# Print the resulting data frame
print(mean_rt_data)










###########plot with WMC bin  slope ############

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(patchwork)  # For combining plots

# Ensure WMC_Bin is a factor with the correct order
mean_rt_data$WMC_Bin <- factor(mean_rt_data$WMC_Bin, levels = c("0-25%", "25-50%", "50-75%", "75-100%"))

# Function to plot RT data, including individual random slopes
plot_rt <- function(data, trial_order_label) {
    ggplot(data, aes(x = WMC_Bin, y = Mean_RT, group = subject_ID, color = condition)) +
        
        # Add individual random slopes (Thin, semi-transparent lines per participant)
        geom_line(aes(group = subject_ID), color = "black", alpha = 0.2, size = 0.5) +  
        
        # Add individual data points (smaller dots)
        geom_point(size = 1, alpha = 0.6) +
        
        # Add mean trend lines in corresponding condition colors
        geom_line(data = data %>%
                      group_by(WMC_Bin, condition) %>%
                      summarise(Mean_RT = mean(Mean_RT, na.rm = TRUE), .groups = 'drop'),
                  aes(x = WMC_Bin, y = Mean_RT, group = condition, color = condition),
                  size = 1.5) +
        
        scale_color_manual(values = c("b" = "orange", "c" = "green", "d" = "red")) +
        labs(title = paste(trial_order_label, "Trial"),
             x = "WMC Bin",
             y = "Mean RT (ms)") +
        theme_minimal() +
        theme(legend.title = element_blank(),
              axis.text = element_text(size = 12),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 16, face = "bold"),
              legend.position = "top")
}

# Separate data for first and last trials
first_trial_data <- mean_rt_data %>% filter(trial_order == "first")
last_trial_data <- mean_rt_data %>% filter(trial_order == "last")

# Generate plots
first_trial_plot <- plot_rt(first_trial_data, "First 10%")
last_trial_plot <- plot_rt(last_trial_data, "Last 10%")

# Combine the plots side by side
combined_plot <- first_trial_plot + last_trial_plot + plot_layout(ncol = 2)

# Print the combined plot
print(combined_plot)

# Save the combined plot
ggsave("combined_mean_rt_with_random_slopes.png", combined_plot, width = 12, height = 6, dpi = 300)

### combine the top and bottom images to make the complete Figure 8###

# Load required packages
library(magick)

# Read the two images
top_image <- image_read("conceptual figure adaptation.png")
bottom_image <- image_read("combined_mean_rt_with_random_slopes.png")

# Resize the top image (e.g., 1.5x its original dimensions)
top_image_resized <- image_scale(top_image, geometry = "200%")


# Combine vertically
combined_image <- image_append(c(top_image_resized, bottom_image), stack = TRUE)

# Save the combined image
image_write(combined_image, "combined_figure_scaled_top.png")


######statistic results that relevant to the bottom plot of Figure 8 and the results are displayed in "Adaptation effects"


library(ggplot2)

# Function to create a histogram showing the distribution of RT differences
# and marking the mean with a red dashed line
create_difference_histogram <- function(data, title, x_limits) {
    mean_val <- mean(data)  # Calculate mean of RT differences
    max_density <- max(density(data)$y)  # For text positioning

    ggplot(data.frame(value = data), aes(x = value)) +
        geom_histogram(aes(y = after_stat(density)), bins = 30, 
                       fill = "lightblue", color = "black", alpha = 0.7) +
        geom_vline(aes(xintercept = mean_val), color = "red", 
                   linetype = "dashed", size = 1) +
        annotate("text", x = mean_val, y = max_density * 0.85, 
                 label = paste("Mean:", round(mean_val, 2)), 
                 color = "red", size = 3, hjust = -0.1) +
        labs(title = title, x = "RT Difference (ms)", y = "Density") +
        scale_x_continuous(limits = x_limits) +
        theme_minimal() +
        theme(
            plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
            axis.title.x = element_text(size = 10),
            axis.title.y = element_text(size = 10),
            axis.text = element_text(size = 8)
        )
}

# Function to aggregate results into posterior samples across participants
aggregate_results_with_summary <- function(results) {
    aggregated_data <- list()

    for (result_name in names(results)) {
        parts <- unlist(strsplit(result_name, "_"))
        wmc_bin <- parts[2]
        trial_order <- parts[3]
        group_key <- paste(wmc_bin, trial_order, sep = "_")

        if (!group_key %in% names(aggregated_data)) {
            aggregated_data[[group_key]] <- list()
        }

        aggregated_data[[group_key]][[length(aggregated_data[[group_key]]) + 1]] <- results[[result_name]]
    }

    summary_list <- list()

    for (group_key in names(aggregated_data)) {
        mat <- do.call(cbind, aggregated_data[[group_key]])
        combined_draws <- as.vector(mat)
        summary_list[[group_key]] <- combined_draws
    }

    return(summary_list)
}

# Function to compute the difference in posterior samples between two conditions
compute_differences <- function(results_1, results_2) {
    differences <- list()
    for (group_key in names(results_1)) {
        if (group_key %in% names(results_2)) {
            differences[[group_key]] <- results_2[[group_key]] - results_1[[group_key]]
        }
    }
    return(differences)
}

### --- MAIN ANALYSIS WORKFLOW --- ###

# Step 1: Aggregate posterior samples for each condition
aggregated_b_results <- aggregate_results_with_summary(results_b)
aggregated_c_results <- aggregate_results_with_summary(results_c)
aggregated_d_results <- aggregate_results_with_summary(results_d)

# Step 2: Compute differences between conditions
difference_c_b <- compute_differences(aggregated_b_results, aggregated_c_results)  # Condition C - B
difference_d_c <- compute_differences(aggregated_c_results, aggregated_d_results)  # Condition D - C

# Step 3: Determine x-axis limits for consistent histograms
all_diff_data <- unlist(c(difference_c_b, difference_d_c))
x_limits_diff <- range(all_diff_data, na.rm = TRUE)

# Step 4: Create and save histograms for Condition C - B
for (group_key in names(difference_c_b)) {
    plot <- create_difference_histogram(
        data = difference_c_b[[group_key]],
        title = paste("Histogram for", group_key, "- Condition C - B"),
        x_limits = x_limits_diff
    )
    print(plot)
    ggsave(filename = paste0("hist_C_minus_B_", group_key, ".png"), 
           plot = plot, width = 6, height = 4, dpi = 300)
}

# Step 5: Create and save histograms for Condition D - C
for (group_key in names(difference_d_c)) {
    plot <- create_difference_histogram(
        data = difference_d_c[[group_key]],
        title = paste("Histogram for", group_key, "- Condition D - C"),
        x_limits = x_limits_diff
    )
    print(plot)
    ggsave(filename = paste0("hist_D_minus_C_", group_key, ".png"), 
           plot = plot, width = 6, height = 4, dpi = 300)
}

