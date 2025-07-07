



# Read the anonymized dataset
human_surprisal_data_anon <- read_csv("human_surprisal_data_anon.csv")



model_combined <- brm(
    formula = log(Total.time.to.correct.answer) ~ 
        ConsonantsAccuracy * trial_order * Full_Surprisal + 
        ConsonantsAccuracy * trial_order * Surprisal0.8.R + 
        (1 | itemID) + 
        (1 + Full_Surprisal + Surprisal0.8.R || subject_ID),
    data = human_surprisal_data_anon,
    family = gaussian(),
    prior = c(
        set_prior("normal(0, 5)", class = "b"),
        set_prior("normal(0, 5)", class = "Intercept"),
        set_prior("student_t(3, 0, 5)", class = "sd")
    ),
    iter = 4000,
    warmup = 1000,
    cores = 4,
    control = list(adapt_delta = 0.95)
)
summary(model_combined)


#


#Regression Coefficients:
#                                              Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#Intercept                                         6.96      0.02     6.92     7.00 1.00     1463     3098
#ConsonantsAccuracy                                0.11      0.06    -0.01     0.24 1.00     1617     2712
#trial_order                                      -0.00      0.00    -0.00    -0.00 1.00    11877     7299
#Full_Surprisal                                    0.03      0.00     0.03     0.04 1.00     4632     6647
#Surprisal0.8.R                                    0.01      0.01    -0.00     0.02 1.00     3681     5590
#ConsonantsAccuracy:trial_order                   -0.00      0.00    -0.00     0.00 1.00    13007     7672
#ConsonantsAccuracy:Full_Surprisal                -0.00      0.01    -0.02     0.01 1.00    14366    10048
#trial_order:Full_Surprisal                        0.00      0.00    -0.00     0.00 1.00    12241     7520
#ConsonantsAccuracy:Surprisal0.8.R                -0.00      0.01    -0.03     0.02 1.00    14669     9197
#trial_order:Surprisal0.8.R                       -0.00      0.00    -0.00    -0.00 1.00    12410     8329
#ConsonantsAccuracy:trial_order:Full_Surprisal     0.00      0.00    -0.00     0.00 1.00    12343     7408
#ConsonantsAccuracy:trial_order:Surprisal0.8.R     0.00      0.00     0.00     0.00 1.00    11557     6836


library(posterior)
# Extract posterior samples
posterior_samples <- as_draws_df(model_combined)
# List of fixed-effect coefficients (adjust names as necessary based on your model)
fixed_effects <- names(posterior_samples)[grepl("^b_", names(posterior_samples))]

# Calculate proportions of positive and negative samples for each effect
coef_proportions <- sapply(fixed_effects, function(coef) {
    prop_pos <- mean(posterior_samples[[coef]] > 0)
    prop_neg <- mean(posterior_samples[[coef]] < 0)
    c(proportion_positive = prop_pos, proportion_negative = prop_neg)
})

# Convert to data frame for easier interpretation
coef_proportions_df <- as.data.frame(t(coef_proportions))
coef_proportions_df

#coef_proportions_df
#                                                proportion_positive proportion_negative
#b_Intercept                                              1.00000000          0.00000000
#b_ConsonantsAccuracy                                     0.96241667          0.03758333
#b_trial_order                                            0.00000000          1.00000000
#b_Full_Surprisal                                         1.00000000          0.00000000
#b_Surprisal0.8.R                                         0.89816667          0.10183333
#b_ConsonantsAccuracy:trial_order                         0.08208333          0.91791667
#b_ConsonantsAccuracy:Full_Surprisal                      0.31033333          0.68966667
#b_trial_order:Full_Surprisal                             0.74141667          0.25858333
#b_ConsonantsAccuracy:Surprisal0.8.R                      0.40708333          0.59291667
#b_trial_order:Surprisal0.8.R                             0.00000000          1.00000000
#b_ConsonantsAccuracy:trial_order:Full_Surprisal          0.74391667          0.25608333
#b_ConsonantsAccuracy:trial_order:Surprisal0.8.R          0.98475000          0.01525000

##### get the exact number of 95% CI ####

library(posterior)

# Extract posterior samples from the model
posterior_samples <- as_draws_df(model_combined)

# Get the names of the fixed effects (coefficients) that start with "b_"
fixed_effects <- names(posterior_samples)[grepl("^b_", names(posterior_samples))]

# Calculate the posterior means (estimates) and 95% credible intervals (CrI) for each coefficient
coef_summary <- sapply(fixed_effects, function(coef) {
    mean_coef <- mean(posterior_samples[[coef]])  # Posterior mean (Estimate)
    ci_coef <- quantile(posterior_samples[[coef]], probs = c(0.025, 0.975))  # 95% CrI
    
    c(Estimate = mean_coef, `2.5%` = ci_coef[1], `97.5%` = ci_coef[2])  # Combine into one vector
})

# Convert to data frame for easier viewing
coef_summary_df <- as.data.frame(t(coef_summary))


# Rename the columns correctly
colnames(coef_summary_df) <- c("Estimate", "2.5%", "97.5%")

# Print the summary data frame with proper column names
print(coef_summary_df)

#print(coef_summary_df)
#                                                     Estimate          2.5%         97.5%
#b_Intercept                                      6.9636681845  6.924927e+00  7.0018237726
#b_ConsonantsAccuracy                             0.1091012573 -1.051641e-02  0.2371286600
#b_trial_order                                   -0.0031312993 -3.445436e-03 -0.0028184663
#b_Full_Surprisal                                 0.0337315996  2.655007e-02  0.0409161725
#b_Surprisal0.8.R                                 0.0088677990 -4.946368e-03  0.0227156053
#b_ConsonantsAccuracy:trial_order                -0.0008414816 -2.023511e-03  0.0003324715
#b_ConsonantsAccuracy:Full_Surprisal             -0.0032018177 -1.578586e-02  0.0094365221
#b_trial_order:Full_Surprisal                     0.0000336958 -7.228488e-05  0.0001382244
#b_ConsonantsAccuracy:Surprisal0.8.R             -0.0029422785 -2.765704e-02  0.0216870199
#b_trial_order:Surprisal0.8.R                    -0.0003832470 -5.815816e-04 -0.0001834678
#b_ConsonantsAccuracy:trial_order:Full_Surprisal  0.0001293896 -2.584763e-04  0.0005194808
#b_ConsonantsAccuracy:trial_order:Surprisal0.8.R  0.0008117842  7.542282e-05  0.0015380800

library(posterior)
# Extract posterior samples
posterior_samples <- as_draws_df(model_combined)
# List of fixed-effect coefficients (adjust names as necessary based on your model)
fixed_effects <- names(posterior_samples)[grepl("^b_", names(posterior_samples))]

# Calculate proportions of positive and negative samples for each effect
coef_proportions <- sapply(fixed_effects, function(coef) {
    prop_pos <- mean(posterior_samples[[coef]] > 0)
    prop_neg <- mean(posterior_samples[[coef]] < 0)
    c(proportion_positive = prop_pos, proportion_negative = prop_neg)
})

# Convert to data frame for easier interpretation
coef_proportions_df <- as.data.frame(t(coef_proportions))
coef_proportions_df

coefs_model_combined <- fixef(model_combined)
print(coefs_model_combined)


#####

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)  # For combining plots
library(scales)     # For formatting numbers
library(ggdist)     # For stat_halfeye()

# Extract all fixed effects (including interaction terms and b_Intercept)
fixed_effects_all <- colnames(posterior_samples) %>%
    grep(pattern = "^b_", value = TRUE)  # Include all b_ terms

# Cleaned names for plotting (does NOT affect model)
clean_effect_names <- fixed_effects_all %>%
    gsub("^b_", "", .) %>%
    gsub("Full_Surprisal", "Full Surprisal", .) %>%
    gsub("ConsonantsAccuracy", "WMC", .) %>%
    gsub("Surprisal0.8.R", "Forgetful LCS", .) %>%
    gsub("trial_order", "trial order", .)

# Define the desired custom order of effects for plotting
desired_order <- c(
    "Intercept",
    "WMC",
    "Full Surprisal",
    "Forgetful LCS",
    "trial order",
    "WMC:Full Surprisal",
    "WMC:Forgetful LCS",
    "WMC:trial order",
    "trial order:Full Surprisal",
    "trial order:Forgetful LCS",
    "WMC:trial order:Full Surprisal",
    "WMC:trial order:Forgetful LCS"
)

# Match cleaned names to the custom order
ordered_indices <- match(desired_order, clean_effect_names)

# Check for any missing terms
missing_terms <- desired_order[is.na(ordered_indices)]
if(length(missing_terms) > 0) {
    message("Warning: The following terms were not found in posterior samples and will be skipped: ", 
            paste(missing_terms, collapse = ", "))
}

# Remove any NAs in case of mismatches
valid_indices <- which(!is.na(ordered_indices))

# Filter and reorder based on desired order
fixed_effects_all_ordered <- fixed_effects_all[ordered_indices[valid_indices]]
clean_effect_names_ordered <- clean_effect_names[ordered_indices[valid_indices]]

# Function to create individual plots with dynamic x-axis
create_plot_dynamic <- function(effect, clean_name) {
    ggplot(posterior_samples, aes(x = .data[[effect]])) +
        stat_halfeye(
            .width = c(0.95),
            fill = "skyblue",
            color = "black"
        ) +
        geom_vline(
            xintercept = 0, 
            linetype = "dashed", 
            color = "red", 
            linewidth = 1
        ) +
        scale_x_continuous(
            labels = scales::number_format(scale = 1, accuracy = 0.0001),
            breaks = scales::pretty_breaks(n = 5)
        ) +
        labs(
            title = clean_name,
            x = "Parameter Estimate",
            y = "Density"
        ) +
        theme_minimal(base_size = 14) +
        theme(
            plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
            axis.title.x = element_text(size = 14),
            axis.title.y = element_text(size = 14),
            axis.text = element_text(size = 12)
        )
}

# Create plots for each effect in custom order
all_fixed_plots <- mapply(create_plot_dynamic, fixed_effects_all_ordered, clean_effect_names_ordered, SIMPLIFY = FALSE)

# Combine the plots into a grid (adjust columns as needed)
combined_all_fixed_plot <- wrap_plots(all_fixed_plots, ncol = 3) +
    plot_annotation(
        theme = theme(
            plot.title = element_text(size = 18, hjust = 0.5, face = "bold")
        )
    )

# Save the plot
ggsave("posterior_distributions_combined.png", combined_all_fixed_plot, width = 16, height = 16, dpi = 300)

# Display the plot
print(combined_all_fixed_plot)
