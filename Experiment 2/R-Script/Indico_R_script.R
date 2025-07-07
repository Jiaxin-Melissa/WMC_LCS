library(lme4)
library(dplyr)

model1_indico_SS_fprt <- lmer(
    log(fprt) ~ Full_Surprisal + SSmean.C * Surprisal_0.8.R + 
        lex_freq.C + word_length.C + wordId.C +
        (1 | item) + (1 + Surprisal_0.8.R | subjId),
    data = merged_indico_clean %>%
        filter(experiment == "et", lex_freq_unk != "True"),
    REML = FALSE,
    control = lmerControl(optimizer = "bobyqa")
)

summary(model1_indico_SS_fprt)






##############

model1_indico_SS_tft <- lmer(
     log(tft) ~ Full_Surprisal + SSmean.C * Surprisal_0.8.R + 
         lex_freq.C + word_length.C + wordId.C +
         (1 | item) + (1 + Surprisal_0.8.R | subjId),
     data = merged_indico_clean %>% filter(experiment == "et", lex_freq_unk != "True"),
     REML = FALSE,
     control = lmerControl(optimizer = "bobyqa")
 )
 summary(model1_indico_SS_tft)







#######



model1_indico_OS_fprt <- lmer(
     log(fprt) ~ Full_Surprisal + OSmean.C * Surprisal_0.8.R + 
         lex_freq.C + word_length.C + wordId.C +
         (1 | item) + (1 + Surprisal_0.8.R | subjId),
     data = merged_indico_clean %>% filter(experiment == "et",, lex_freq_unk != "True"),
     REML = FALSE,
     control = lmerControl(optimizer = "bobyqa")
 )
 summary(model1_indico_OS_fprt)



######



model1_indico_OS_tft <- lmer(
     log(tft) ~ Full_Surprisal + OSmean.C * Surprisal_0.8.R + 
         lex_freq.C + word_length.C + wordId.C +
         (1 | item) + (1 + Surprisal_0.8.R | subjId),
     data = merged_indico_clean %>% filter(experiment == "et", lex_freq_unk != "True"),
     REML = FALSE,
     control = lmerControl(optimizer = "bobyqa")
 )
 summary(model1_indico_OS_tft)



#####


######





model1_indico_MU_fprt <- lmer(
     log(fprt) ~ Full_Surprisal + MUmean.C * Surprisal_0.8.R + 
         lex_freq.C + word_length.C + wordId.C +
         (1 | item) + (1 + Surprisal_0.8.R | subjId),
     data = merged_indico_clean %>% filter(experiment == "et", lex_freq_unk != "True"),
     REML = FALSE,
     control = lmerControl(optimizer = "bobyqa")
 )
 summary(model1_indico_MU_fprt)




#######
#





##########
model1_indico_MU_tft <- lmer(
     log(tft) ~ Full_Surprisal + MUmean.C * Surprisal_0.8.R + 
         lex_freq.C + word_length.C + wordId.C +
         (1 | item) + (1 + Surprisal_0.8.R | subjId),
     data = merged_indico_clean %>% filter(experiment == "et", lex_freq_unk != "True"),
     REML = FALSE,
     control = lmerControl(optimizer = "bobyqa")
 )
 summary(model1_indico_MU_tft)



#####



####

model1_indico_SSTM_fprt <- lmer(
     log(fprt) ~ Full_Surprisal + SSTMRelScore.C * Surprisal_0.8.R + 
         lex_freq.C + word_length.C + wordId.C +
         (1 | item) + (1 + Surprisal_0.8.R | subjId),
     data = merged_indico_clean %>% filter(experiment == "et", lex_freq_unk != "True"),
     REML = FALSE,
     control = lmerControl(optimizer = "bobyqa")
 )
 summary(model1_indico_SSTM_fprt) 







###########





model1_indico_SSTM_tft <- lmer(
     log(tft) ~ Full_Surprisal + SSTMRelScore.C * Surprisal_0.8.R + 
         lex_freq.C + word_length.C + wordId.C +
         (1 | item) + (1 + Surprisal_0.8.R | subjId),
     data = merged_indico_clean %>% filter(experiment == "et", lex_freq_unk != "True"),
     REML = TRUE,
     control = lmerControl(optimizer = "bobyqa")
 )
 summary(model1_indico_SSTM_tft) 



####


spr_data <- subset(merged_indico_clean, experiment == "spr")


mean_tft <- mean(spr_data$tft, na.rm = TRUE)
sd_tft <- sd(spr_data$tft, na.rm = TRUE)


lower_bound <- mean_tft - 3 * sd_tft
upper_bound <- mean_tft + 3 * sd_tft


spr_data_clean <- subset(spr_data, tft >= lower_bound & tft <= upper_bound)



n_before <- nrow(spr_data)


n_after <- nrow(spr_data_clean)


n_removed <- n_before - n_after


percent_removed <- (n_removed / n_before) * 100


cat(sprintf("Removed %d out of %d trials (%.2f%%) based on ±3 SD cutoff.\n",
            n_removed, n_before, percent_removed))

spr_data_clean["tft"] = spr_data_clean["tft"] * 1000

########

# Subset only SPR data
spr_data <- subset(merged_indico_clean, experiment == "spr")

# Remove rows where lex_freq_unk == "True"
spr_data <- subset(spr_data, lex_freq_unk != "True")

# Compute mean and SD of tft
mean_tft <- mean(spr_data$tft, na.rm = TRUE)
sd_tft <- sd(spr_data$tft, na.rm = TRUE)

# Define cutoff bounds
lower_bound <- mean_tft - 3 * sd_tft
upper_bound <- mean_tft + 3 * sd_tft

# Remove outliers
spr_data_clean <- subset(spr_data, tft >= lower_bound & tft <= upper_bound)

# Count before and after
n_before <- nrow(spr_data)
n_after <- nrow(spr_data_clean)
n_removed <- n_before - n_after
percent_removed <- (n_removed / n_before) * 100

# Print summary
cat(sprintf("Removed %d out of %d trials (%.2f%%) based on ±3 SD cutoff.\n",
            n_removed, n_before, percent_removed))

# Convert tft to milliseconds
spr_data_clean$tft <- spr_data_clean$tft * 1000



#########


# Load libraries
library(lme4)
library(broom.mixed)
library(dplyr)

# Step 1: Scale the necessary predictors manually with suffix ".S"
spr_data_clean <- spr_data_clean %>%
  mutate(
    Full_Surprisal.S = scale(Full_Surprisal)[, 1],
    Surprisal_0.8.R.S = scale(Surprisal_0.8.R)[, 1],
    lex_freq.C.S = scale(lex_freq.C)[, 1],
    word_length.C.S = scale(word_length.C)[, 1],
    wordId.C.S = scale(wordId.C)[, 1],
    MUmean.C.S = scale(MUmean.C)[, 1],
    SSmean.C.S = scale(SSmean.C)[, 1],
    OSmean.C.S = scale(OSmean.C)[, 1],
    SSTMRelScore.C.S = scale(SSTMRelScore.C)[, 1]
  )

# Step 2: Fit scaled models
model1_indico_MU_spr <- lmer(
  log(tft) ~ Full_Surprisal.S + MUmean.C.S * Surprisal_0.8.R.S +
    lex_freq.C.S + word_length.C.S + wordId.C.S +
    (1 | item) + (1 + Surprisal_0.8.R.S | subjId),
  data = spr_data_clean,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

model1_indico_SS_spr <- lmer(
  log(tft) ~ Full_Surprisal.S + SSmean.C.S * Surprisal_0.8.R.S +
    lex_freq.C.S + word_length.C.S + wordId.C.S +
    (1 | item) + (1 + Surprisal_0.8.R.S | subjId),
  data = spr_data_clean,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

model1_indico_OS_spr <- lmer(
  log(tft) ~ Full_Surprisal.S + OSmean.C.S * Surprisal_0.8.R.S +
    lex_freq.C.S + word_length.C.S + wordId.C.S +
    (1 | item) + (1 + Surprisal_0.8.R.S | subjId),
  data = spr_data_clean,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

model1_indico_SSTM_spr <- lmer(
  log(tft) ~ Full_Surprisal.S + SSTMRelScore.C.S * Surprisal_0.8.R.S +
    lex_freq.C.S + word_length.C.S + wordId.C.S +
    (1 | item) + (1 + Surprisal_0.8.R.S | subjId),
  data = spr_data_clean,
  REML = FALSE,
  control = lmerControl(optimizer = "Nelder_Mead")
)

# Step 3: Define scaling information
scaled_vars <- c(
  "Full_Surprisal", "Surprisal_0.8.R", "lex_freq.C",
  "word_length.C", "wordId.C", "MUmean.C", "SSmean.C",
  "OSmean.C", "SSTMRelScore.C"
)

scale_map <- setNames(
  scaled_vars,
  paste0(scaled_vars, ".S")
)

sd_values <- sapply(spr_data_clean[scaled_vars], sd, na.rm = TRUE)

# Step 4: Function to extract back-transformed fixed effects
extract_raw_results <- function(model) {
  broom.mixed::tidy(model, conf.int = TRUE) %>%
    filter(effect == "fixed") %>%
    mutate(
      beta = case_when(
        term == "(Intercept)" ~ estimate,
        term %in% names(scale_map) ~ estimate / sd_values[scale_map[term]],
        TRUE ~ estimate
      ),
      SE = case_when(
        term == "(Intercept)" ~ std.error,
        term %in% names(scale_map) ~ std.error / sd_values[scale_map[term]],
        TRUE ~ std.error
      )
    ) %>%
    select(term, beta, SE, t_value = statistic)
}

# Step 5: Get raw fixed effects for each model
raw_MU <- extract_raw_results(model1_indico_MU_spr)
raw_SS <- extract_raw_results(model1_indico_SS_spr)
raw_OS <- extract_raw_results(model1_indico_OS_spr)
raw_SSTM <- extract_raw_results(model1_indico_SSTM_spr)

# Step 6: Print results
cat("\n=== MU Model ===\n"); print(raw_MU, n = Inf)
cat("\n=== SS Model ===\n"); print(raw_SS, n = Inf)
cat("\n=== OS Model ===\n"); print(raw_OS, n = Inf)
cat("\n=== SSTM Model ===\n"); print(raw_SSTM, n = Inf)
#######

cat("\n=== MU Model ===\n"); print(raw_MU, n = Inf)



library(ggplot2)
library(dplyr)

# Filter to rows with non-missing values
plot_data <- spr_data_clean %>%
    filter(!is.na(LogWordFreq.C), !is.na(lex_freq.C), is.finite(LogWordFreq.C), is.finite(lex_freq.C))

ggplot(plot_data, aes(x = LogWordFreq, y = lex_freq.C)) +
    geom_point(alpha = 0.2, color = "steelblue") +
    geom_smooth(method = "lm", color = "darkred", se = FALSE) +
    labs(
        title = "Correlation between LogWordFreq and lex_freq.C",
        x = "Log Word Frequency (Wikipedia)",
        y = "Lexical Frequency (Lemma-based)"
    ) +
    theme_minimal()


