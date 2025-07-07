# Load necessary libraries
library(dplyr)


#read the csv file
human_surprisal_data_anon <- read_csv("human_surprisal_data_anon.csv")

human_surprisal_data_anon <- na.omit(human_surprisal_data_anon)

bayesian_model <- brm(
    formula = Correct_binary ~ ConsonantsAccuracy * bVSc*trial_order + ConsonantsAccuracy * cVSd*trial_order + 
        (1 | itemID) + (1 + bVSc + cVSd + trial_order + bVSc*trial_order + cVSd*trial_order | subject_ID),
    data = human_surprisal_data_anon,
    family = bernoulli(link = "logit"),  # For binary response
    chains = 4,       # Number of MCMC chains
    cores = 4,        # Number of cores for parallel computing
    iter = 2000,      # Number of iterations per chain
    warmup = 1000,    # Number of warm-up iterations
    control = list(adapt_delta = 0.95)  # Control settings to handle divergent transitions
)

summary(bayesian_model)


# Define sigmoid function
sigmoid <- function(x) {
    return(1 / (1 + exp(-x)))
}

# Extract posterior samples
posterior_samples <- as_draws_df(bayesian_model)

# Extract random effects for each participant
random_effects <- ranef(bayesian_model)$subject_ID

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
        results_b[[result_name_b]] <- sigmoid(
            posterior_samples[["b_Intercept"]] +
                random_intercept +
                bc * (posterior_samples[["b_bVSc"]] + random_bVSc) +
                trial_order * (posterior_samples[["b_trial_order"]] + random_trial_order) +
                trial_order * bc * (posterior_samples[["b_bVSc:trial_order"]] + random_trial_order_bVSc) +
                wmc * posterior_samples[["b_ConsonantsAccuracy"]] +
                wmc * trial_order * (posterior_samples[["b_ConsonantsAccuracy:trial_order"]] + random_trial_order) +
                wmc * bc * trial_order * (posterior_samples[["b_ConsonantsAccuracy:bVSc:trial_order"]] + random_trial_order_bVSc)
        )
        
        # ====================
        # Condition C
        # ====================
        bc <- 1
        cd <- -1
        result_name_c <- paste0(prolific_id, "_", wmc_bin, "_", trial_order_name, "_c")
        results_c[[result_name_c]] <- sigmoid(
            posterior_samples[["b_Intercept"]] +
                random_intercept +
                bc * (posterior_samples[["b_bVSc"]] + random_bVSc) +
                cd * (posterior_samples[["b_cVSd"]] + random_cVSd) +
                trial_order * (posterior_samples[["b_trial_order"]] + random_trial_order) +
                trial_order * bc * (posterior_samples[["b_bVSc:trial_order"]] + random_trial_order_bVSc) +
                trial_order * cd * (posterior_samples[["b_trial_order:cVSd"]] + random_trial_order_cVSd) +
                wmc * posterior_samples[["b_ConsonantsAccuracy"]] +
                wmc * trial_order * (posterior_samples[["b_ConsonantsAccuracy:trial_order"]] + random_trial_order) +
                wmc * bc * trial_order * (posterior_samples[["b_ConsonantsAccuracy:bVSc:trial_order"]] + random_trial_order_bVSc) +
                wmc * cd * trial_order * (posterior_samples[["b_ConsonantsAccuracy:trial_order:cVSd"]] + random_trial_order_cVSd)
        )
        
        # ====================
        # Condition D
        # ====================
        cd <- 1
        result_name_d <- paste0(prolific_id, "_", wmc_bin, "_", trial_order_name, "_d")
        results_d[[result_name_d]] <- sigmoid(
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


# Initialize an empty data frame to store proportion correct values
mean_accuracy_data <- data.frame(
    subject_ID = character(),
    WMC_bin = character(),
    trial_order = character(),
    condition = character(),
    Proportion_Correct = numeric(),
    stringsAsFactors = FALSE
)

# Function to extract details from the result name
extract_details <- function(result_name) {
    parts <- unlist(strsplit(result_name, "_"))
    list(
        subject_ID = parts[1],
        WMC_bin = parts[2],
        trial_order = parts[3],
        condition = parts[4]
    )
}

# Helper function to process results and update data frame
process_results <- function(results, data_frame) {
    for (result_name in names(results)) {
        details <- extract_details(result_name)
        
        data_frame <- rbind(data_frame, data.frame(
            subject_ID = details$subject_ID,
            WMC_bin = details$WMC_bin,
            trial_order = details$trial_order,
            condition = details$condition,
            Proportion_Correct = mean(results[[result_name]], na.rm = TRUE)
        ))
    }
    return(data_frame)
}

# Process results for each condition
mean_accuracy_data <- process_results(results_b, mean_accuracy_data)
mean_accuracy_data <- process_results(results_c, mean_accuracy_data)
mean_accuracy_data <- process_results(results_d, mean_accuracy_data)

# Print the updated data frame
print(mean_accuracy_data)



##########plot#######

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(patchwork)  # For combining plots

# Ensure WMC_Bin is a factor with the correct order
mean_accuracy_data$WMC_bin <- factor(mean_accuracy_data$WMC_bin, levels = c("0-25%", "25-50%", "50-75%", "75-100%"))

# Function to plot RT data, including individual random slopes
plot_rt <- function(data, trial_order_label) {
    ggplot(data, aes(x = WMC_bin, y = Proportion_Correct, group = subject_ID, color = condition)) +
        
        # Add individual random slopes (Thin, semi-transparent lines per participant)
        geom_line(aes(group = subject_ID), color = "black", alpha = 0.2, size = 0.5) +  
        
        # Add individual data points (smaller dots)
        geom_point(size = 1, alpha = 0.6) +
        
        # Add mean trend lines in corresponding condition colors
        geom_line(data = data %>%
                      group_by(WMC_bin, condition) %>%
                      summarise(Proportion_Correct = mean(Proportion_Correct, na.rm = TRUE), .groups = 'drop'),
                  aes(x = WMC_bin, y = Proportion_Correct, group = condition, color = condition),
                  size = 1.5) +
        
        scale_color_manual(values = c("b" = "orange", "c" = "green", "d" = "red")) +
        labs(title = paste(trial_order_label, "Trial"),
             x = "WMC Bin",
             y = "A-Maze Accuracy") +
        theme_minimal() +
        theme(legend.title = element_blank(),
              axis.text = element_text(size = 12),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 16, face = "bold"),
              legend.position = "top")
}

# Separate data for first and last trials
first_trial_data <- mean_accuracy_data %>% filter(trial_order == "first")
last_trial_data <- mean_accuracy_data %>% filter(trial_order == "last")

# Generate plots
first_trial_plot <- plot_rt(first_trial_data, "First 10%")
last_trial_plot <- plot_rt(last_trial_data, "Last 10%")

# Combine the plots side by side
combined_plot <- first_trial_plot + last_trial_plot + plot_layout(ncol = 2)

# Print the combined plot
print(combined_plot)

# Save the combined plot
ggsave("combined_mean_accuracy_with_random_slopes.png", combined_plot, width = 12, height = 6, dpi = 300)

