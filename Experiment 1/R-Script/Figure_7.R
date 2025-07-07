
human_surprisal_data_anon <- read_csv("human_surprisal_data_anon.csv")

human_surprisal_data_anon <- na.omit(human_surprisal_data_anon)
# Define the model using brms
model1_brm <- brm(
  bf(log(Total.time.to.correct.answer) ~ 
       ConsonantsAccuracy * trial_order * bVSc +
       ConsonantsAccuracy * trial_order * cVSd +
       (1 | itemID) + (1 + bVSc + cVSd + trial_order + trial_order*bVSc + trial_order*cVSd | subject_ID)),
  data = human_surprisal_data_anon,
  family = gaussian(), 
  chains = 4,          
  cores = 4,           
  iter = 4000,         
  warmup = 1000,       
  control = list(adapt_delta = 0.95) 
)

summary(model1_brm)


##### Calculate Proportions of Positive and Negative Posterior Samples for Fixed Effects #####

# Load necessary package
library(posterior)

# Extract posterior samples from the brms model
posterior_samples <- as_draws_df(model1_brm)

# Identify fixed effects (those starting with "b_")
fixed_effects <- grep("^b_", names(posterior_samples), value = TRUE)

# Compute proportion of positive and negative samples for each fixed effect
coef_proportions <- sapply(fixed_effects, function(coef_name) {
  samples <- posterior_samples[[coef_name]]
  c(
    proportion_positive = mean(samples > 0),
    proportion_negative = mean(samples < 0)
  )
})

# Convert to a clean data frame for display
coef_proportions_df <- as.data.frame(t(coef_proportions))

# Display the result
print(coef_proportions_df)


##### Extract 95% Credible Intervals for Fixed Effects #####

# Load necessary package
library(posterior)

# Extract posterior samples from the brms model
posterior_samples <- as_draws_df(model1_brm)

# Identify fixed effects (those starting with "b_")
fixed_effects <- grep("^b_", names(posterior_samples), value = TRUE)

# Compute mean and 95% credible intervals for each fixed effect
coef_summary <- sapply(fixed_effects, function(coef_name) {
  samples <- posterior_samples[[coef_name]]
  c(
    Estimate = mean(samples),
    `2.5%` = quantile(samples, 0.025),
    `97.5%` = quantile(samples, 0.975)
  )
})

# Convert to a clean data frame for display
coef_summary_df <- as.data.frame(t(coef_summary))
colnames(coef_summary_df) <- c("Estimate", "2.5%", "97.5%")

# Display the result
print(coef_summary_df)



#####Figure 7####

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)  
library(scales)     
library(ggdist)     

# ---- Define function to plot individual effects ----
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

#  Desired order of fixed effects (matching exact names in posterior_samples)
desired_order <- c(
    "b_Intercept",
    "b_ConsonantsAccuracy",
    "b_bVSc",
    "b_cVSd",
    "b_trial_order",
    "b_ConsonantsAccuracy:bVSc",  
    "b_ConsonantsAccuracy:cVSd",
    "b_ConsonantsAccuracy:trial_order",  
    "b_trial_order:bVSc", 
    "b_trial_order:cVSd",
    "b_ConsonantsAccuracy:trial_order:bVSc",  
    "b_ConsonantsAccuracy:trial_order:cVSd"  
)
# Clean names to display (matching the corrected desired order)
clean_names_ordered <- c(
    "Intercept",
    "WMC",
    "bVSc",
    "cVSd",
    "trial_order",
    "WMC:bVSc",
    "WMC:cVSd",
    "WMC:trial_order",
    "bVSc:trial_order",
    "cVSd:trial_order",
    "WMC:bVSc:trial_order",
    "WMC:cVSd:trial_order"
)

# Filter to those that exist in the data 
available_effects <- desired_order[desired_order %in% colnames(posterior_samples)]
matched_clean_names <- clean_names_ordered[desired_order %in% colnames(posterior_samples)]

# Create and combine plots for all desired effects 
ordered_fixed_plots <- mapply(create_plot_dynamic, available_effects, matched_clean_names, SIMPLIFY = FALSE)

# Combine the plots in a grid (adjust `ncol` for grid size) 
combined_all_fixed_plot <- wrap_plots(ordered_fixed_plots, ncol = 3) +
    plot_annotation(
        theme = theme(
            plot.title = element_text(size = 18, hjust = 0.5, face = "bold")
        )
    )

#  Display the combined plot 
combined_all_fixed_plot
# Save the combined plot with adjusted dimensions
ggsave("posterior_distributions_all_fixed_effects_dynamic_human_data.png", combined_all_fixed_plot, width = 16, height = 16, dpi = 300)





