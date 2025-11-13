################################################################################
# Trigeminal Sensitivity Analysis - Bormann Study
# Author: Jorn Lotsch
# Description: Analysis of trigeminal sensitivity study data, including
#              preprocessing, translation, categorization, descriptive stats,
#              and tabulations per variable category.
################################################################################


# ======================== #
# 1. Load Required Libraries
# ======================== #
library(stringr)
library(dplyr)
library(readr)
library(tidyr)
library(psych)


# ======================== #
# 2. Define some basic functions
# ======================== #

one_hot_encode_var <- function(var, var_name) {
  factor_var <- factor(var, exclude = NULL)
  levels_var <- levels(factor_var)

  one_hot_list <- lapply(levels_var, function(lev) {
    if (is.na(lev)) {
      as.integer(is.na(var))  # 1 if NA in original, else 0
    } else {
      # Use ifelse to convert NA comparisons to 0
      as.integer(ifelse(is.na(var), 0, var == lev))
    }
  })

  names(one_hot_list) <- paste0(var_name, "_", ifelse(is.na(levels_var), "NA", levels_var))

  return(as.data.frame(one_hot_list))
}

convert_jn_to_numeric <- function(x) {
  ifelse(is.na(x), NA, ifelse(x == "j", 1, 0))
}


# =============================== #
# 3. Read all raw data
# =============================== #

trigeminale_daten_corrected_translated <- read.csv("trigeminale_daten_corrected_translated.csv", check.names = FALSE)

character_vars <- names(trigeminale_daten_corrected_translated)[sapply(trigeminale_daten_corrected_translated, is.character)]
print(character_vars)

variable_categories <- read.csv("trigeminale_daten_variable_categories.csv")

# Remove duplicated columns, keeping the first occurrence
trigeminale_daten_corrected_translated <- trigeminale_daten_corrected_translated[ , !duplicated(names(trigeminale_daten_corrected_translated))]

# Load global mappings (e.g., category labels, translation dictionaries).
source("globals.R")

# =============================== #
# 4. Start analysis_dataset dataframe
# =============================== #

analysis_dataset <- cbind.data.frame(ID = trigeminale_daten_corrected_translated$ID)

# =============================== #
# 5. Prepare data category wise
# =============================== #
#
# Category 1: Demographics
#
# =============================== #


print("Demographics")
print(variables_by_categories$Demographics)

analysis_dataset <- cbind.data.frame(
  analysis_dataset,
  trigeminale_daten_corrected_translated[,variables_by_categories$Demographics[variables_by_categories$Demographics %in% c("Age", "Weight", "Height")]]
  )

analysis_dataset <- cbind.data.frame(
  analysis_dataset,
  one_hot_encode_var(
    trigeminale_daten_corrected_translated[,variables_by_categories$Demographics[variables_by_categories$Demographics %in% c("Gender")]],
    var_name = variables_by_categories$Demographics[variables_by_categories$Demographics %in% c("Gender")]
    )
)

head(analysis_dataset)
max(table(analysis_dataset$ID))


# =============================== #
#
# Category 2: Disorders or health complaints
#
# =============================== #

print("Disorders_or_health_complaints")
print(variables_by_categories$Disorders_or_health_complaints)

# Chronic diseases
print("Reading one-hot encoded chronic diseases variables")

onehot_chronic_diseases <- read.csv("onehot_chronic_diseases.csv", check.names = FALSE)
head(onehot_chronic_diseases)
analysis_dataset <- analysis_dataset %>%
  left_join(onehot_chronic_diseases, by = "ID")
head(analysis_dataset)
max(table(analysis_dataset$ID))


# ENT surgery
print("Reading one-hot encoded ENT surgeries variables")

onehot_ENT_surgery <- read.csv("onehot_ENT_surgery.csv", check.names = FALSE)
head(onehot_ENT_surgery)

analysis_dataset <- analysis_dataset %>%
  left_join(onehot_ENT_surgery, by = "ID")
head(analysis_dataset)
max(table(analysis_dataset$ID))

# =============================== #
#
# Category 3: COVID-19 history
#
# =============================== #

print("COVID19_history")
print(variables_by_categories$COVID19_history)

analysis_dataset[variables_by_categories$COVID19_history[1]] <-
  convert_jn_to_numeric(trigeminale_daten_corrected_translated$`Have you had COVID-19?`)

covid_times_in_past <- read.csv("covid_times_in_past.csv", check.names = FALSE)

df_covid_info <- subset(covid_times_in_past,
                        select = c("ID", "How many times have you had COVID-19?", "Is or was there smell reduction after COVID-19?"))
df_covid_info["Longest time since covid"] <-
  apply(covid_times_in_past[,c("months_since_Period_1","months_since_Period_2","months_since_Period_3")],1, function(x) max(x,na.rm = TRUE))
df_covid_info[["Longest time since covid"]][is.infinite(df_covid_info[["Longest time since covid"]]) & df_covid_info[["Longest time since covid"]] < 0] <- NA

df_covid_info["Shortest time since covid"] <-
  apply(covid_times_in_past[,c("months_since_Period_1","months_since_Period_2","months_since_Period_3")],1, function(x) min(x,na.rm = TRUE))
df_covid_info[["Shortest time since covid"]][is.infinite(df_covid_info[["Shortest time since covid"]]) ] <- NA

analysis_dataset <- analysis_dataset %>%
  left_join(df_covid_info, by = "ID")

analysis_dataset["Covid-related change in smell ability"] <- trigeminale_daten_corrected_translated$`Smell ability before COVID-19` -
  trigeminale_daten_corrected_translated$`Smell ability immediately after COVID-19`
head(analysis_dataset)
max(table(analysis_dataset$ID))

# =============================== #
#
# Category 4: Smoking and alcohol use
#
# =============================== #

print("Smoking_and_alcohol_use")
print(variables_by_categories$Smoking_and_alcohol_use)

head(trigeminale_daten_corrected_translated[,variables_by_categories$Smoking_and_alcohol_use])

analysis_dataset[variables_by_categories$Smoking_and_alcohol_use[1]] <-
  convert_jn_to_numeric(trigeminale_daten_corrected_translated$`Do you smoke?`)
analysis_dataset[variables_by_categories$Smoking_and_alcohol_use[4]] <-
  convert_jn_to_numeric(trigeminale_daten_corrected_translated$`Have you ever smoked?`)


print("Reading preprocesses smoking behavior related variables")

smoking_summary <- read.csv("smoking_summary.csv", row.names = 1)
actual_year <- as.numeric(format(d <- as.Date("2023-11-01") + (as.Date("2024-05-01") - as.Date("2023-11-01"))/2, "%Y")) + (as.numeric(format(d, "%j")) - 1) / ifelse(((y <- as.numeric(format(d, "%Y"))) %% 4 == 0 & y %% 100 != 0) | (y %% 400 == 0), 366, 365)

smoking_summary$time_since_quitting_smoking <- round(as.numeric(actual_year) - smoking_summary$bar_end, 1)
table(smoking_summary$smoker_type)
max(table(smoking_summary$ID))

# Convert smoker_type to character for pivoting
smoking_summary <- smoking_summary %>%
  mutate(smoker_type = as.character(smoker_type))

# Get unique ID-smoker_type combinations
smoker_unique <- smoking_summary %>%
  distinct(ID, smoker_type)

# Add presence indicator
smoker_unique <- smoker_unique %>%
  mutate(present = 1)

# Pivot to wide format: one column per smoker type category
smoker_one_hot <- smoker_unique %>%
  pivot_wider(
    id_cols = ID,
    names_from = smoker_type,
    values_from = present,
    values_fill = list(present = 0),
    values_fn = list(present = max)  # Handle duplicates by taking max
  )

# Rename columns with underscores for consistency
colnames(smoker_one_hot) <- gsub("^([a-zA-Z0-9]+)$", "smoker_type_\\1", colnames(smoker_one_hot))
colnames(smoker_one_hot)[colnames(smoker_one_hot) == "smoker_type_ID"] <- "ID"

# View result
print(head(smoker_one_hot))

smoking_summary <- smoking_summary %>%
  left_join(smoker_one_hot, by = "ID")

analysis_dataset <- analysis_dataset %>%
  left_join(smoking_summary[c("ID", "mean_cigs", "pack_years", "time_since_quitting_smoking", "smoker_type_current", "smoker_type_former")], by = "ID")

# Fill missing values for rows without smoking period info:
analysis_dataset <- analysis_dataset %>%
  mutate(
    mean_cigs = if_else(is.na(mean_cigs), 0, mean_cigs),
    pack_years = if_else(is.na(pack_years), 0, pack_years),
    # Use large sentinel value to distinguish never-smokers from recent quitters
    # Value chosen to be beyond plausible maximum (>800 years)
    time_since_quitting_smoking = if_else(is.na(time_since_quitting_smoking), Inf, time_since_quitting_smoking),
    smoker_type_current = if_else(is.na(smoker_type_current), 0, smoker_type_current),
    smoker_type_former = if_else(is.na(smoker_type_former), 0, smoker_type_former),
    # Create never_smoker indicator: 1 if all smoking indicators are zero or missing
    never_smoker = if_else(
      is.na(smoker_type_current) & is.na(smoker_type_former),
      1,
      if_else(smoker_type_current == 0 & smoker_type_former == 0, 1, 0)
    )
  )


analysis_dataset <- analysis_dataset %>%
  left_join(trigeminale_daten_corrected_translated[c("ID", "Do you drink alcohol?")])
head(analysis_dataset)
max(table(analysis_dataset$ID))

# =============================== #
#
# Category 5: Subjective nasal chemosensory perception (TriFunQ and related)
#
# =============================== #

print("Nasal_chemosensory_perception")
print(variables_by_categories$Nasal_chemosensory_perception)

variables_nasal_chemosensory_perception <- cbind.data.frame(ID = trigeminale_daten_corrected_translated$ID,
                                                            trigeminale_daten_corrected_translated[variables_by_categories$Nasal_chemosensory_perception])

head(variables_nasal_chemosensory_perception)

analysis_dataset <- analysis_dataset %>%
  left_join(variables_nasal_chemosensory_perception)
head(analysis_dataset)
max(table(analysis_dataset$ID))

# =============================== #
#
# Category 6: Rated_olfactory_function
#
# =============================== #

print("Rated_olfactory_function")
print(variables_by_categories$Rated_olfactory_function)

variables_rated_olfactory_function <- cbind.data.frame(ID = trigeminale_daten_corrected_translated$ID,
                                                            trigeminale_daten_corrected_translated[variables_by_categories$Rated_olfactory_function])
head(variables_rated_olfactory_function)
apply(variables_rated_olfactory_function,2,table)

one_hot_reduced_smell_and_taste_ability_what <- one_hot_encode_var(variables_rated_olfactory_function$`Reduced smell and taste ability`,
                                                              "reduced_smell_and_taste_ability_what")

one_hot_if_yes_how_did_the_problem_start <- one_hot_encode_var(variables_rated_olfactory_function$`If yes, how did the problem start`,
                                                               "reduced_smell_and_taste_ability_start")

analysis_dataset <- analysis_dataset %>%
  left_join(variables_rated_olfactory_function[,c("ID", "Current smell ability",
                                                  "Reduced smell and taste ability", "If yes, how did the problem start", "How has the problem changed")])
analysis_dataset <- analysis_dataset %>%
  left_join(cbind.data.frame(ID = trigeminale_daten_corrected_translated$ID, one_hot_reduced_smell_and_taste_ability_what))
analysis_dataset <- analysis_dataset %>%
  left_join(cbind.data.frame(ID = trigeminale_daten_corrected_translated$ID, one_hot_if_yes_how_did_the_problem_start))
head(analysis_dataset)
max(table(analysis_dataset$ID))


# =============================== #
#
# Category 7: Ratings of nasal irritation and airflow
#
# =============================== #

print("Nasal_irritation_and_airflow")
print(variables_by_categories$Nasal_irritation_and_airflow)

variables_nasal_irritation_and_airflow <- cbind.data.frame(ID = trigeminale_daten_corrected_translated$ID,
                                                       trigeminale_daten_corrected_translated[variables_by_categories$Nasal_irritation_and_airflow])
head(variables_nasal_irritation_and_airflow)

analysis_dataset <- analysis_dataset %>%
  left_join(variables_nasal_irritation_and_airflow)
head(analysis_dataset)
max(table(analysis_dataset$ID))


# =============================== #
#
# Category 7: Psychophysical measurements: Odor identification
#
# =============================== #

print("Psychophysical_measurements")
print(variables_by_categories$Psychophysical_measurements)

analysis_dataset[variables_by_categories$Psychophysical_measurements[3]] <- trigeminale_daten_corrected_translated[variables_by_categories$Psychophysical_measurements[3]]
head(analysis_dataset)


# =============================== #
# 6. Check and analyze the complete assembled data set
# =============================== #

print(dim(analysis_dataset))
str(analysis_dataset)

# Calculate descriptive stats for numeric variables only
desc_stats <- describe(analysis_dataset[sapply(analysis_dataset, is.numeric)], na.rm = TRUE)

# Calculate IQR strings for numeric variables only
iqr_values <- sapply(analysis_dataset[sapply(analysis_dataset, is.numeric)], function(x) {
  qs <- quantile(x, probs = c(0.25, 0.75), na.rm = TRUE)
  paste0(format(qs[1], digits=3), " - ", format(qs[2], digits=3))
})

# Convert to data frame to combine with desc_stats
iqr_df <- data.frame(IQR = iqr_values, row.names = names(iqr_values))

# Combine describe output and IQR into one data frame by row names (variable names)
desc_table <- cbind(desc_stats, IQR = iqr_df$IQR)

# Calculate column sums for numeric columns in the original data
col_sums <- colSums(analysis_dataset[sapply(analysis_dataset, is.numeric)], na.rm = TRUE)

# Convert to named vector matching desc_table rownames
col_sums_vec <- col_sums[rownames(desc_table)]

# Add sum as a new column to desc_table
desc_table$Sum <- col_sums_vec

# Print updated table
print(desc_table)

# =============================== #
# 7. Save data and statistics
# =============================== #

write.csv(analysis_dataset, "analysis_dataset.csv")
write.csv(desc_table, "descriptive_statistics_analysis_dataset_not_imputed.csv")

# =============================== #
# 8. Handling of missing data
# =============================== #

# Check which variables or cases have missings more than acceptable
print("Checking which variables or cases have missings more than acceptable.")

n_missings_per_variable <- apply(analysis_dataset,2,function(x) sum(is.na(x)))

f_missings <- 0.2
f_missings_per_variable <- n_missings_per_variable / nrow(analysis_dataset)
names(analysis_dataset)[which(f_missings_per_variable > f_missings)]
