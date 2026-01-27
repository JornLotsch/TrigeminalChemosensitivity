################################################################################
# Trigeminal Sensitivity Analysis - Bormann Study
# Author: Jorn Lotsch
# Description: Analysis of trigeminal sensitivity study data, including
#              preprocessing, translation, categorization, descriptive stats,
#              and tabulations per variable category.
################################################################################


# ==============================================================================
# 1. Load Required Libraries
# ==============================================================================

library(stringr)
library(dplyr)
library(readr)
library(tidyr)
library(psych)
library(ggplot2)
library(ggthemes)
library(nortest)
library(missRanger)
library(missForest)
library(opdisDownsampling)
library(twosamples)


source("globals.R")

# ==============================================================================
# 2. Define Helper Functions
# ==============================================================================

#' One-hot encode a nominal vector
#'
#' Converts a factor or character vector to a one-hot encoded data.frame.
#'
#' @param var The vector to encode (factor or character).
#' @param var_name Base name for new columns.
#' @return A data.frame with one column per factor level, values are 0/1.
one_hot_encode_var <- function(var, var_name) {
  factor_var <- factor(var, exclude = NULL)
  levels_var <- levels(factor_var)
  one_hot_list <- lapply(levels_var, function(lev) {
    if (is.na(lev)) {
      as.integer(is.na(var))
    } else {
      as.integer(ifelse(is.na(var), 0, var == lev))
    }
  })
  names(one_hot_list) <- paste0(var_name, "_", ifelse(is.na(levels_var), "NA", levels_var))
  return(as.data.frame(one_hot_list))
}

#' Convert yes/no (j/n) to numeric
#'
#' Recodes "j" to 1, "n" to 0, NA stays NA.
#'
#' @param x Character vector of "j", "n", or NA.
#' @return Numeric vector (1/0/NA).
convert_jn_to_numeric <- function(x) {
  ifelse(is.na(x), NA, ifelse(x == "j", 1, 0))
}

#' Convert TRUE/FALSE logicals to numeric
#'
#' Converts TRUE to 1, FALSE to 0, NA stays NA.
#'
#' @param x Logical vector (TRUE/FALSE/NA).
#' @return Numeric vector (1/0/NA).
convert_truefalse_to_numeric <- function(x) {
  ifelse(is.na(x), NA, ifelse(x == TRUE, 1, 0))
}

#' Get dimension of a data.frame without ID column
#'
#' @param df A data.frame or tibble with optional "ID" column.
#' @return Vector: number of rows and columns (excluding "ID").
dim_wo_ID <- function(df) {
  if ("ID" %in% names(df)) {
    dim(df[, -1])
  } else dim(df)
}

#' Count variables in a data.frame without ID column
#'
#' @param df A data.frame or tibble with optional "ID" column.
#' @return Integer, number of variables not counting "ID".
nVars_wo_ID <- function(df) {
  if ("ID" %in% names(df)) {
    length(names(df)) - 1
  } else length(names(df))
}

#' Get variable names excluding ID column
#'
#' @param df A data.frame or tibble with optional "ID" column.
#' @return Character vector of variable names (excluding "ID").
var_names_wo_ID <- function(df) {
  names(df)[!names(df) %in% "ID"]
}

#' Left join two dataframes and set NA values in new columns to a given value
#'
#' @param df_primary Left (main) dataframe
#' @param df_secondary Right (joined) dataframe
#' @param by Character vector of columns to join by (default: "ID")
#' @param fill_value Value to replace NAs with (default: 0)
#' @return A joined dataframe with NAs in added columns from df_secondary replaced by fill_value
left_join_fill <- function(df_primary, df_secondary, by = "ID", fill_value = 0) {
  # Do join
  df_combined <- left_join(df_primary, df_secondary, by = by)

  # Find columns that were added from secondary dataframe
  added_cols <- setdiff(names(df_combined), names(df_primary))

  # For each added column, replace NA with fill_value
  df_combined <- df_combined %>%
    mutate(across(all_of(added_cols), ~ ifelse(is.na(.), fill_value, .)))

  return(df_combined)
}


# ==============================================================================
# 3. Read Raw Data
# ==============================================================================

trigeminale_daten_corrected_translated <- read.csv("trigeminale_daten_corrected_translated.csv", check.names = FALSE)

character_vars <- names(trigeminale_daten_corrected_translated)[sapply(trigeminale_daten_corrected_translated, is.character)]
print(character_vars)

variable_categories <- read.csv("trigeminale_daten_variable_categories.csv", check.names = FALSE)

# Remove duplicated columns, keeping the first occurrence
trigeminale_daten_corrected_translated <- trigeminale_daten_corrected_translated[, !duplicated(names(trigeminale_daten_corrected_translated))]

# Load global mappings (e.g., category labels, translation dictionaries).
source("globals.R")

# ==============================================================================
# 4. Initialize Analysis Dataset
# ==============================================================================

analysis_dataset_empty <- cbind.data.frame(ID = trigeminale_daten_corrected_translated$ID)
analysis_dataset <- analysis_dataset_empty

analysis_dataset_variables <- data.frame(category = character(0), variable_name = character(0), stringsAsFactors = FALSE)

# ==============================================================================
# 5. Prepare Data by Category
# ==============================================================================
# =============================== #
#
# Category 1: Demographics
#
# =============================== #


print("Demographics")
print(variables_by_categories$Demographics)
cat1_variables_directly_included <- c("Age", "Weight", "Height")

Category_1_data <- analysis_dataset_empty %>%
  left_join(onehot_chronic_diseases, by = "ID")

Category_1_data <- analysis_dataset_empty %>%
  left_join(trigeminale_daten_corrected_translated[, c("ID", cat1_variables_directly_included)], by = "ID")

Category_1_data <- cbind.data.frame(
  Category_1_data,
  one_hot_encode_var(
    trigeminale_daten_corrected_translated[, variables_by_categories$Demographics[variables_by_categories$Demographics %in% c("Gender")]],
    var_name = variables_by_categories$Demographics[variables_by_categories$Demographics %in% c("Gender")]
  )
)

head(Category_1_data, n = 2)
max(table(Category_1_data$ID))

print("Number of variables in category 1: Demographics")
print(nVars_wo_ID(Category_1_data))
print("Variables in category 1: Demographics")
print(var_names_wo_ID(Category_1_data))

# Add variables to the variable_names list
analysis_dataset_variables <- rbind.data.frame(analysis_dataset_variables, cbind.data.frame(category = 1, variable_name = var_names_wo_ID(Category_1_data)))

# =============================== #
#
# Category 2: Disorders or health complaints
#
# =============================== #

print("Disorders_or_health_complaints")
print(variables_by_categories$Disorders_or_health_complaints)

## Chronic diseases
print("Reading one-hot encoded chronic diseases variables")

# One hot encoded
onehot_chronic_diseases <- read.csv("onehot_chronic_diseases.csv", check.names = FALSE)
dim(onehot_chronic_diseases)
head(onehot_chronic_diseases, n = 2)
print(dim_wo_ID(onehot_chronic_diseases))
print(var_names_wo_ID(onehot_chronic_diseases))

Category_2_data <- analysis_dataset_empty %>%
  left_join(onehot_chronic_diseases, by = "ID")
head(Category_2_data, n = 2)
max(table(Category_2_data$ID))

Category_2_data <- Category_2_data %>%
  rename(Unpecified_chronic_disease = Unspecified)

# # Remove migraine due to unreliability
#
# columns_to_keep <- !grepl("igrai", names(Category_2_data))
# Category_2_data <- Category_2_data[, columns_to_keep]

# Add further variables (migraine (2))
Category_2_data["Has migraine changed over the last 10 years"] <- 2 - trigeminale_daten_corrected_translated$`Has migraine changed over the last 10 years`
print(var_names_wo_ID(Category_2_data))
Category_2_data["How often do you have migraine per month"] <- trigeminale_daten_corrected_translated$`How often do you have migraine per month`

## Actual respiratory tract infections
Category_2_data["Upper respiratory infection"] <- convert_jn_to_numeric(trigeminale_daten_corrected_translated$`Upper respiratory infection`)

# Intermediate check for other category 2 variables not yet included
category_2_vars_not_included <- setdiff(variables_by_categories$Disorders_or_health_complaints, names(Category_2_data))
head(trigeminale_daten_corrected_translated[, category_2_vars_not_included], n = 2)

## ENT surgery
print("Reading one-hot encoded ENT surgeries variables")
onehot_ENT_surgery <- read.csv("onehot_ENT_surgery.csv", check.names = FALSE)
dim(onehot_ENT_surgery)
head(onehot_ENT_surgery, n = 2)
print(dim_wo_ID(onehot_ENT_surgery))
print(var_names_wo_ID(onehot_ENT_surgery))

Category_2_data <- Category_2_data %>%
  left_join(onehot_ENT_surgery, by = "ID")
head(Category_2_data, n = 2)
max(table(Category_2_data$ID))

Category_2_data <- Category_2_data %>%
  rename(Unpecified_ENT_surgery = Unspecified)

# Intermediate check for other category 2 variables not yet included
category_2_vars_not_included <- setdiff(variables_by_categories$Disorders_or_health_complaints, names(Category_2_data))
head(trigeminale_daten_corrected_translated[, category_2_vars_not_included], n = 2)

## Nasal breathing variables
nasal_breathing_df_clean <- read.csv("nasal_breating_prepocessed.csv", check.names = FALSE)
dim(nasal_breathing_df_clean)

# Nasal problems summary variables
names(nasal_breathing_df_clean)
Category_2_data["Medical consultation for nasal breathing problems"] <- convert_jn_to_numeric(trigeminale_daten_corrected_translated$`Medical consultation for nasal breathing problems`)
Category_2_data$nasal_breathing_problem <- nasal_breathing_df_clean$nasal_breathing_problem
Category_2_data$nasal_breathing_problem_ongoing <- convert_truefalse_to_numeric(nasal_breathing_df_clean$ongoing)

# Prepare one_hot with
nasal_breating_therpies_one_hot <- read.csv("nasal_breating_therpies_one_hot.csv", check.names = FALSE)
dim(nasal_breating_therpies_one_hot)
names(nasal_breating_therpies_one_hot)

Category_2_data <- Category_2_data %>%
  left_join_fill(nasal_breating_therpies_one_hot, by = "ID", fill_value = 0)
head(Category_2_data, n = 2)
max(table(Category_2_data$ID))

# Intermediate check for other category 2 variables not yet included
category_2_vars_not_included <- setdiff(variables_by_categories$Disorders_or_health_complaints, names(Category_2_data))
head(trigeminale_daten_corrected_translated[, category_2_vars_not_included])

## Facial pain
facial_pain_preprocessed <- read.csv("facial_pain_preprocessed.csv", check.names = FALSE)
head(facial_pain_preprocessed)
dim(facial_pain_preprocessed)
print(dim_wo_ID(facial_pain_preprocessed))
print(var_names_wo_ID(facial_pain_preprocessed))

facial_pain_preprocessed["Has facial pain"] <- ifelse(!is.na(facial_pain_preprocessed$end_year),
                                                       ifelse(facial_pain_preprocessed$actual_facial_pain == 1, 2, 1), 0)
facial_pain_vars <- c("ID", "Has facial pain", "facial_pain_freq_code", "facial_pain_pulling", "facial_pain_stabbing", "facial_pain_pressing", "facial_pain_burning")

Category_2_data <- Category_2_data %>%
  left_join(facial_pain_preprocessed[facial_pain_vars], by = "ID")
head(Category_2_data)
max(table(Category_2_data$ID))

# Intermediate check for other category 2 variables not yet included
category_2_vars_not_included <- setdiff(variables_by_categories$Disorders_or_health_complaints, names(Category_2_data))
head(trigeminale_daten_corrected_translated[, category_2_vars_not_included])

print("Number of variables in category 2: Disorders or health complaints")
print(nVars_wo_ID(Category_2_data))
print("Variables in category 2: Disorders or health complaints")
print(var_names_wo_ID(Category_2_data))

# Add variables to the variable_names list
analysis_dataset_variables <- rbind.data.frame(analysis_dataset_variables, cbind.data.frame(category = 2, variable_name = var_names_wo_ID(Category_2_data)))

# =============================== #
#
# Category 3: COVID-19 history
#
# =============================== #

print("COVID19_history")
print(variables_by_categories$COVID19_history)
covid_times_in_past <- read.csv("covid_times_in_past.csv", check.names = FALSE)
print(names(covid_times_in_past))
dim(covid_times_in_past)

print("Variables in both original data and processed COVID data")
print(union(names(covid_times_in_past), variables_by_categories$COVID19_history))
print("Variables differing between original data and processed COVID data")
print(union(setdiff(names(covid_times_in_past), variables_by_categories$COVID19_history),
              setdiff(variables_by_categories$COVID19_history, names(covid_times_in_past))))

# Convert logical columns to numeric
logical_columns <- sapply(covid_times_in_past, is.logical)
covid_times_in_past[logical_columns] <- lapply(covid_times_in_past[logical_columns], as.numeric)

# =============================== #
# Compute temporal COVID variables
# =============================== #

# Compute row-wise maximum time since any COVID infection, ensuring `Inf` handling
covid_times_in_past[["Longest time since covid"]] <-
  apply(
    covid_times_in_past[, c("months_since_covid_period_1", "months_since_covid_period_2", "months_since_covid_period_3")],
    1,
    function(x) {
      # Remove Inf values from the calculation
      x <- x[!is.infinite(x)]

      # If all remaining values are NA after removing Inf, return NA
      if (all(is.na(x))) {
        return(NA)
      }

      # Otherwise return the max of the remaining values
      return(max(x, na.rm = TRUE))
    }
  )

# Compute row-wise minimum time since most recent COVID infection, ensuring `Inf` handling
covid_times_in_past["Shortest time since covid"] <-
  apply(covid_times_in_past[, c("months_since_covid_period_1", "months_since_covid_period_2", "months_since_covid_period_3")], 1, function(x) min(x, na.rm = TRUE))

# Replace Inf values (from min of all-NA rows) with NA
covid_times_in_past[["Shortest time since covid"]][is.infinite(covid_times_in_past[["Shortest time since covid"]])] <- NA

# =============================== #
# Create derived categorical variable for COVID temporal patterns
# =============================== #
# This variable provides a complete representation for all subjects:
# 0 = Never infected
# 1 = Infected >24 months ago
# 2 = Infected 12-24 months ago
# 3 = Infected 6-12 months ago
# 4 = Infected <6 months ago

covid_times_in_past["covid_temporal_category"] <- apply(
  covid_times_in_past[, c("months_since_covid_period_1", "months_since_covid_period_2", "months_since_covid_period_3")],
  1,
  function(x) {
    # Remove Inf values
    x <- x[!is.infinite(x)]

    # If all are NA, person never had COVID
    if (all(is.na(x))) {
      return(0)
    }

    # Get shortest time since infection (most recent)
    min_time <- min(x, na.rm = TRUE)

    # Categorize based on time since most recent infection
    if (min_time < 6) {
      return(4) # <6 months
    } else if (min_time < 12) {
      return(3) # 6-12 months
    } else if (min_time < 24) {
      return(2) # 12-24 months
    } else {
      return(1) # >24 months
    }
  }
)

dim(covid_times_in_past)

# =============================== #
# Assemble Category_3 data with proper handling for non-COVID patients
# =============================== #

# Variables that should be 0 for non-COVID patients (counts, binary indicators, ratings)
covid_vars_zero_in_no_covid <- c("Have you had COVID-19?", "How many times have you had COVID-19?", "Smell ability immediately after COVID-19",
                                  "Is or was there smell reduction after COVID-19?", "Smell ability before COVID-19",
                                  "Smell reduction 1", "Smell reduction 2", "Smell reduction 3")

Category_3_data <- analysis_dataset_empty %>%
  left_join_fill(covid_times_in_past[, c("ID", covid_vars_zero_in_no_covid)], by = "ID", fill_value = 0)
head(Category_3_data)
dim(Category_3_data)
max(table(Category_3_data$ID))

# Add categorical temporal variable: 0 for never-infected, 1-4 for time categories
# This provides a complete representation that works in all analysis contexts
Category_3_data <- Category_3_data %>%
  left_join_fill(covid_times_in_past[, c("ID", "covid_temporal_category")], by = "ID", fill_value = 0)
head(Category_3_data)
dim(Category_3_data)
max(table(Category_3_data$ID))

# =============================== #
# Data integrity verification
# =============================== #

cat("\n=== COVID-19 Data Summary ===\n")
cat("Total subjects:", nrow(Category_3_data), "\n")
cat("Subjects with COVID history:", sum(Category_3_data$`Have you had COVID-19?` > 0, na.rm = TRUE), "\n")
cat("Subjects without COVID history:", sum(Category_3_data$`Have you had COVID-19?` == 0, na.rm = TRUE), "\n")

cat("\nCOVID temporal category distribution:\n")
print(table(Category_3_data$covid_temporal_category, useNA = "ifany"))

cat("\n=== End of COVID-19 Data Assembly ===\n")

print("Number of variables in category 3: COVID-19 history")
print(nVars_wo_ID(Category_3_data))
print("Variables in category 3: COVID-19 history")
print(var_names_wo_ID(Category_3_data))

# Add variables to the variable_names list
analysis_dataset_variables <- rbind.data.frame(analysis_dataset_variables, cbind.data.frame(category = 3, variable_name = var_names_wo_ID(Category_3_data)))

# =============================== #
#
# Category 4: Smoking and alcohol use
#
# =============================== #

print("Smoking_and_alcohol_use")
print(variables_by_categories$Smoking_and_alcohol_use)

head(trigeminale_daten_corrected_translated[, variables_by_categories$Smoking_and_alcohol_use])

Category_4_data <- analysis_dataset_empty

Category_4_data[variables_by_categories$Smoking_and_alcohol_use[1]] <-
  convert_jn_to_numeric(trigeminale_daten_corrected_translated$`Do you smoke?`)
Category_4_data[variables_by_categories$Smoking_and_alcohol_use[4]] <-
  convert_jn_to_numeric(trigeminale_daten_corrected_translated$`Have you ever smoked?`)


print("Reading preprocessed smoking behavior related variables")

smoking_summary <- read.csv("smoking_summary.csv", check.names = FALSE)
dim(smoking_summary)

actual_year <- as.numeric(format(d <- as.Date("2023-11-01") + (as.Date("2024-05-01") - as.Date("2023-11-01")) / 2, "%Y")) + (as.numeric(format(d, "%j")) - 1) / ifelse(((y <- as.numeric(format(d, "%Y"))) %% 4 == 0 & y %% 100 != 0) | (y %% 400 == 0), 366, 365)

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
    values_fn = list(present = max) # Handle duplicates by taking max
  )

# Rename columns with underscores for consistency
colnames(smoker_one_hot) <- gsub("^([a-zA-Z0-9]+)$", "smoker_type_\\1", colnames(smoker_one_hot))
colnames(smoker_one_hot)[colnames(smoker_one_hot) == "smoker_type_ID"] <- "ID"

# View result
print(head(smoker_one_hot))

smoking_summary <- smoking_summary %>%
  left_join(smoker_one_hot, by = "ID")

# =============================== #
# Handle smoking data similar to COVID pattern
# =============================== #

# Variables that should be 0 for never-smokers (counts, binary indicators)
smoking_vars_zero_for_never_smokers <- c("mean_cigs", "pack_years", "smoker_type_current", "smoker_type_former")

Category_4_data <- Category_4_data %>%
  left_join_fill(
    smoking_summary[c("ID", smoking_vars_zero_for_never_smokers)],
    by = "ID",
    fill_value = 0
  )

# =============================== #
# Create smoking temporal category variable
# =============================== #
# Provides complete representation for all subjects (parallels covid_temporal_category):
# 0 = Never smoked
# 1 = Quit >10 years ago
# 2 = Quit 5-10 years ago
# 3 = Quit 1-5 years ago
# 4 = Quit <1 year ago or current smoker

Category_4_data <- Category_4_data %>%
  mutate(
    smoking_temporal_category = case_when(
# Current smokers
      smoker_type_current == 1 ~ 4L,

# Former smokers - categorize by time since quitting
      smoker_type_former == 1 &
        !is.na(smoking_summary$time_since_quitting_smoking[match(ID, smoking_summary$ID)]) &
        smoking_summary$time_since_quitting_smoking[match(ID, smoking_summary$ID)] < 1 ~ 4L,
      smoker_type_former == 1 &
        !is.na(smoking_summary$time_since_quitting_smoking[match(ID, smoking_summary$ID)]) &
        smoking_summary$time_since_quitting_smoking[match(ID, smoking_summary$ID)] < 5 ~ 3L,
      smoker_type_former == 1 &
        !is.na(smoking_summary$time_since_quitting_smoking[match(ID, smoking_summary$ID)]) &
        smoking_summary$time_since_quitting_smoking[match(ID, smoking_summary$ID)] < 10 ~ 2L,
      smoker_type_former == 1 &
        !is.na(smoking_summary$time_since_quitting_smoking[match(ID, smoking_summary$ID)]) &
        smoking_summary$time_since_quitting_smoking[match(ID, smoking_summary$ID)] >= 10 ~ 1L,

# Never smokers (default)
      TRUE ~ 0L
    )
  )

# =============================== #
# Handle alcohol data
# =============================== #

Category_4_data <- Category_4_data %>%
  left_join_fill(
    trigeminale_daten_corrected_translated[c("ID", "Do you drink alcohol?")],
    by = "ID",
    fill_value = 0
  )

head(Category_4_data)
max(table(Category_4_data$ID))

print("Number of variables in category 4: Smoking and alcohol use")
print(nVars_wo_ID(Category_4_data))
print("Variables in category 4: Smoking and alcohol use")
print(var_names_wo_ID(Category_4_data))

# Add variables to the variable_names list
analysis_dataset_variables <- rbind.data.frame(
  analysis_dataset_variables,
  cbind.data.frame(category = 4, variable_name = var_names_wo_ID(Category_4_data))
)


# =============================== #
#
# Category 5: Subjective nasal chemosensory perception (TriFunQ and related)
#
# =============================== #

print("Nasal_chemosensory_perception")
print(variables_by_categories$Nasal_chemosensory_perception)

Category_5_data <- analysis_dataset_empty %>%
  left_join(trigeminale_daten_corrected_translated[, c("ID", variables_by_categories$Nasal_chemosensory_perception)], by = "ID")
head(Category_5_data)
max(table(Category_5_data$ID))

print("Number of variables in category 5: Subjective nasal chemosensory perception (TriFunQ and related)")
print(nVars_wo_ID(Category_5_data))
print("Variables in category 5: Subjective nasal chemosensory perception (TriFunQ and related)")
print(var_names_wo_ID(Category_5_data))

# Add variables to the variable_names list
analysis_dataset_variables <- rbind.data.frame(analysis_dataset_variables, cbind.data.frame(category = 5, variable_name = var_names_wo_ID(Category_5_data)))

# =============================== #
#
# Category 6: Rated_olfactory_function
#
# =============================== #

print("Rated_olfactory_function")
print(variables_by_categories$Rated_olfactory_function)

variables_rated_olfactory_function <-
  trigeminale_daten_corrected_translated[, c("ID", variables_by_categories$Rated_olfactory_function)]
head(variables_rated_olfactory_function, n = 2)
apply(variables_rated_olfactory_function, 2, table)

# Recoding of how the olfactory problem has started
recode_problem_start <- function(x) {
  # Preserve NA, then recode others as specified
  ifelse(
    is.na(x), 0,
    dplyr::recode(x,
                   `0` = 1L,
                   `1` = 2L,
                   `3` = 3L,
                   `2` = 4L,
                   .default = NA_integer_
    )
  )
}

variables_rated_olfactory_function$`If yes, how did the problem start` <- recode_problem_start(variables_rated_olfactory_function$`If yes, how did the problem start`)

# Recoding of how the olfactory problem has changed
recode_problem_change <- function(x) {
  # Preserve NA, then recode others as specified
  ifelse(
    is.na(x), 0,
    dplyr::recode(x,
                   `0` = 1L,
                   `1` = 2L,
                   `2` = 3L,
                   .default = NA_integer_
    )
  )
}

variables_rated_olfactory_function$`How has the problem changed` <- recode_problem_change(variables_rated_olfactory_function$`How has the problem changed`)

# Correct uninformative "j" to "Unknwon"
variables_rated_olfactory_function$`Reduced smell and taste ability`[variables_rated_olfactory_function$`Reduced smell and taste ability` == "j"] <- "Unknwon"

one_hot_reduced_smell_and_taste_ability_what <- one_hot_encode_var(variables_rated_olfactory_function$`Reduced smell and taste ability`,
                                                                    "reduced_smell_and_taste_ability_what")

Category_6_data <- analysis_dataset_empty %>%
  left_join(variables_rated_olfactory_function[, c("ID", "Current smell ability", "How has the problem changed", "If yes, how did the problem start")], by = "ID")

Category_6_data <- cbind.data.frame(Category_6_data, one_hot_reduced_smell_and_taste_ability_what)


head(Category_6_data)
max(table(Category_6_data$ID))


print("Number of variables in category 6: Rated_olfactory_function")
print(nVars_wo_ID(Category_6_data))
print("Variables in category 6: Rated_olfactory_function")
print(var_names_wo_ID(Category_6_data))

# Add variables to the variable_names list
analysis_dataset_variables <- rbind.data.frame(analysis_dataset_variables, cbind.data.frame(category = 6, variable_name = var_names_wo_ID(Category_6_data)))

# =============================== #
#
# Category 7: Ratings of nasal irritation and airflow
#
# =============================== #

print("Nasal_irritation_and_airflow")
print(variables_by_categories$Nasal_irritation_and_airflow)

Category_7_data <- analysis_dataset_empty %>%
  left_join(trigeminale_daten_corrected_translated[, c("ID", variables_by_categories$Nasal_irritation_and_airflow)], by = "ID")
head(Category_7_data)
max(table(Category_7_data$ID))

print("Number of variables in category 7: Ratings of nasal irritation and airflow")
print(nVars_wo_ID(Category_7_data))
print("Variables in category 7: Ratings of nasal irritation and airflow")
print(var_names_wo_ID(Category_7_data))

# Add variables to the variable_names list
analysis_dataset_variables <- rbind.data.frame(analysis_dataset_variables, cbind.data.frame(category = 7, variable_name = var_names_wo_ID(Category_7_data)))

# =============================== #
#
# Category 8: Psychophysical measurements
#
# =============================== #

print("Psychophysical_measurements")
print(variables_by_categories$Psychophysical_measurements)

Category_8_data <- analysis_dataset_empty %>%
  left_join(trigeminale_daten_corrected_translated[, c("ID", variables_by_categories$Psychophysical_measurements)], by = "ID")
head(Category_8_data)
max(table(Category_8_data$ID))

print("Number of variables in category 8: Psychophysical measurements")
print(nVars_wo_ID(Category_8_data))
print("Variables in category 8: Psychophysical measurements")
print(var_names_wo_ID(Category_8_data))

# Add variables to the variable_names list
analysis_dataset_variables <- rbind.data.frame(analysis_dataset_variables, cbind.data.frame(category = 8, variable_name = var_names_wo_ID(Category_8_data)))

# ==============================================================================
# 6. Summarize Assembled Data per Category
# ==============================================================================

# Print variable names by category
print_names_by_category <- function() {
  for (i in 1:8) {
    df <- get(paste0("Category_", i, "_data"))
    cat("\n", paste0("Number of variables in Category_", i), ": ")
    cat(paste(ncol(df[, !names(df) %in% c("ID")])))
    cat("\n", paste0("Category_", i, "_data"), "\n")
    cat(paste(names(df)[!names(df) %in% c("ID")], collapse = ", "), "\n")
  }
}

print_names_by_category()

print(table(analysis_dataset_variables$category))


# Generate a html table with all variables plus first 2 values
generate_categories_html_table_with_types <- function() {
  html <- '<table border="1" style="border-collapse: collapse; width: 100%;">'
  for (i in 1:8) {
    df <- get(paste0("Category_", i, "_data"))
    # Subheading
    html <- paste0(html,
                    '<tr style="background-color: #cce5ff; font-weight: bold;">',
                    '<td colspan="', ncol(df), '" style="text-align:center;">',
                    'Category_', i, '_data',
                    '</td></tr>')
    # Variable names row
    html <- paste0(html, '<tr>',
                    paste0('<th>', colnames(df), '</th>', collapse = ''),
                    '</tr>')
    # Variable types row, from class of each column
    types <- sapply(df, function(x) paste(class(x), collapse = ","))
    html <- paste0(html, '<tr style="background-color:#e9ecef;">',
                    paste0('<td><i>', types, '</i></td>', collapse = ''),
                    '</tr>')
    # First two data rows
    for (row in 1:min(2, nrow(df))) {
      row_values <- sapply(df[row,], as.character)
      html <- paste0(html, '<tr>',
                      paste0('<td>', row_values, '</td>', collapse = ''),
                      '</tr>')
    }
  }
  html <- paste0(html, '</table>')
  return(html)
}

# Save to an HTML file to open externally
html_content <- generate_categories_html_table_with_types()
writeLines(html_content, "variables_by_categories_summary.html")

# ==============================================================================
# 7. Combine All Category Data for Analysis
# ==============================================================================

combine_all_categories <- function() {
  combined_df <- Category_1_data
  for (i in 2:8) {
    next_df <- get(paste0("Category_", i, "_data"))
    combined_df <- left_join(combined_df, next_df, by = "ID")
  }
  return(combined_df)
}

# Create final data set and write it to csv
analysis_dataset <- combine_all_categories()
write.csv(analysis_dataset, "analysis_dataset.csv", row.names = FALSE)
dim(analysis_dataset[, -1])

# =============================== #
# 8. Remove the CO2 thresholds acquired with breathing
# =============================== #

analysis_dataset_CO2_breath_hold <- analysis_dataset
analysis_dataset_CO2_breath_hold$`CO2 threshold`[1:549] <- NA

# =============================== #
# 9. Handling of missing data
# =============================== #

# Check which variables or cases have missings more than acceptable
print("Checking which variables or cases have missings more than acceptable.")

n_missings_per_variable <- apply(analysis_dataset_CO2_breath_hold, 2, function(x) sum(is.na(x)))
cat("\n", paste0("Missing values per variable (descending order): "))
sort(n_missings_per_variable[which(n_missings_per_variable > 0)], decreasing = TRUE)
sort(n_missings_per_variable[which(n_missings_per_variable > 0)], decreasing = TRUE) / nrow(analysis_dataset_CO2_breath_hold) *100


f_missings <- 0.2
f_missings_per_variable <- n_missings_per_variable / nrow(analysis_dataset_CO2_breath_hold)

vars_to_drop <- names(analysis_dataset_CO2_breath_hold)[which(f_missings_per_variable > f_missings)]
cat("\n", paste0("Variables with more than ", 100 * f_missings, "% missing values: "))
print(vars_to_drop)
cat("\n", paste0("Total variables dropped: ", length(vars_to_drop), "\n"))

f_missings <- 0.8
vars_to_drop <- names(analysis_dataset_CO2_breath_hold)[which(f_missings_per_variable > f_missings)]
cat("\n", paste0("Variables with more than ", 100 * f_missings, "% missing values: "))
print(vars_to_drop)
cat("\n", paste0("Total variables dropped: ", length(vars_to_drop), "\n"))

# Filter dataset to remove high-missing variables
analysis_dataset_remaining_vars <- analysis_dataset_CO2_breath_hold[, !names(analysis_dataset_CO2_breath_hold) %in% vars_to_drop]

n_missings_analysis_dataset_remaining_vars <- sum(is.na(analysis_dataset_remaining_vars))

cat("\n", paste0("Total number of missings in remaining dataset: "))
print(n_missings_analysis_dataset_remaining_vars)
print(n_missings_analysis_dataset_remaining_vars / length(unlist(analysis_dataset_remaining_vars[, -1])) * 100)


# =============================== #
# 10. Train/Validation Split (Before Imputation)
# =============================== #

cat("\n=== Splitting data into training and validation sets ===\n")

# Prepare dataset for splitting (remove ID for split, will rejoin later)
analysis_dataset_to_split <- analysis_dataset_remaining_vars

# Perform stratified downsampling split
set.seed(42)
analysis_dataset_split <- opdisDownsampling::opdisDownsampling(
  analysis_dataset_to_split,
  Size = 0.8,
  Seed = 42,
  nTrials = 100000,
  MaxCores = min(12, parallel::detectCores() - 1)
)

# write split result IDs
write.table(
  as.matrix(analysis_dataset_split$ReducedInstances),
  "analysis_training_instances.csv",
  sep = ",",
  row.names = FALSE,
  col.names = FALSE
)


# Extract training and validation sets
training_data <- analysis_dataset_split$ReducedData
validation_data <- analysis_dataset_split$RemovedData

# Rename the variables to orginal non-R conform versions
names(training_data) <- names(analysis_dataset_to_split)
names(validation_data) <- names(analysis_dataset_to_split)


cat("  Training set size:", nrow(training_data), "observations\n")
cat("  Validation set size:", nrow(validation_data), "observations\n")
cat("  Training set missing values:", sum(is.na(training_data[, -1])), "\n")
cat("  Validation set missing values:", sum(is.na(validation_data[, -1])), "\n")

# =============================== #
# 11. Imputation Functions
# =============================== #

cat("\n=== Defining imputation functions ===\n")

#' Detect variable type for imputation strategy
#'
#' @param x A numeric vector
#' @return Character string: "binary", "ordinal", "continuous"
detect_variable_type <- function(x) {
  x_clean <- na.omit(x)
  unique_vals <- unique(x_clean)
  n_unique <- length(unique_vals)

  # Binary (including one-hot encoded)
  if (n_unique <= 2 && all(unique_vals %in% c(0, 1))) {
    return("binary")
  }

  # Ordinal: small number of integer values
  if (n_unique <= 10 && all(x_clean == floor(x_clean))) {
    return("ordinal")
  }

  # Otherwise continuous
  return("continuous")
}

#' Type-aware fallback imputation
#'
#' @param data Data frame without ID column
#' @return Data frame with imputed values
fallback_impute <- function(data) {
  cat("  Applying type-aware fallback imputation...\n")

  for (col in names(data)) {
    if (any(is.na(data[[col]]))) {
      var_type <- detect_variable_type(data[[col]])

      if (var_type == "binary") {
        # Mode imputation (most common value, typically 0)
        mode_val <- as.numeric(names(sort(table(data[[col]]), decreasing = TRUE)[1]))
        data[[col]][is.na(data[[col]])] <- mode_val
        cat("    ", col, ": binary -> mode (", mode_val, ")\n", sep = "")

      } else if (var_type == "ordinal") {
        # Median imputation + rounding
        median_val <- round(median(data[[col]], na.rm = TRUE))
        data[[col]][is.na(data[[col]])] <- median_val
        data[[col]] <- as.integer(round(data[[col]]))
        cat("    ", col, ": ordinal -> median rounded (", median_val, ")\n", sep = "")

      } else {
        # Continuous: median imputation
        median_val <- median(data[[col]], na.rm = TRUE)
        data[[col]][is.na(data[[col]])] <- median_val
        cat("    ", col, ": continuous -> median (", signif(median_val, 4), ")\n", sep = "")
      }
    }
  }

  return(data)
}

#' #' Impute a dataset using missForest with fallback
#' #'
#' #' @param dataset Data frame with ID column
#' #' @param dataset_name Name for logging
#' #' @param reference_data Optional reference data for type detection (use training data for validation)
#' #' @return Imputed data frame
#' impute_dataset <- function(dataset, dataset_name, reference_data = NULL) {
#'   cat("\n--- Imputing", dataset_name, "---\n")
#'   cat("  Variables:", ncol(dataset) - 1, "\n")
#'   cat("  Observations:", nrow(dataset), "\n")
#'
#'   # Separate ID column
#'   ids <- dataset$ID
#'   data_no_id <- dataset[, -1, drop = FALSE]
#'
#'   # Check if there are any missing values
#'   n_missing <- sum(is.na(data_no_id))
#'   cat("  Missing values:", n_missing, "\n")
#'
#'   if (n_missing == 0) {
#'     cat("  No missing values - skipping imputation\n")
#'     return(dataset)
#'   }
#'
#'   # Use reference data for type detection if provided (for validation set)
#'   if (is.null(reference_data)) {
#'     reference_data <- dataset
#'   }
#'
#'   # Try missForest first
#'   cat("  Attempting missForest imputation...\n")
#'   set.seed(42)
#'
#'   imputation_result <- try(
#'     suppressWarnings(missForest::missForest(training_data, maxiter = 10000, ntree = 1000)),
#'     silent = TRUE
#'   )
#'
#'   if (!inherits(imputation_result, "try-error")) {
#'     cat("  ✓ missForest succeeded\n")
#'     imputed_data <- imputation_result$ximp
#'
#'     # Post-process: ensure proper data types based on reference data
#'     cat("  Post-processing data types...\n")
#'     for (col in names(imputed_data)) {
#'       var_type <- detect_variable_type(reference_data[[col]])
#'
#'       if (var_type == "binary") {
#'         # Round and clip to 0/1
#'         imputed_data[[col]] <- pmin(1, pmax(0, round(imputed_data[[col]])))
#'       } else if (var_type == "ordinal") {
#'         # Round to nearest integer
#'         imputed_data[[col]] <- as.integer(round(imputed_data[[col]]))
#'       }
#'     }
#'
#'   } else {
#'     cat("  ✗ missForest failed - using fallback\n")
#'     cat("  Error:", attr(imputation_result, "condition")$message, "\n")
#'     imputed_data <- fallback_impute(data_no_id)
#'   }
#'
#'   # Rejoin ID column
#'   imputed_dataset <- cbind(ID = ids, imputed_data)
#'
#'   # Verify no missing values remain
#'   remaining_na <- sum(is.na(imputed_dataset[, -1]))
#'   cat("  Remaining missing values:", remaining_na, "\n")
#'
#'   if (remaining_na > 0) {
#'     cat("  ⚠ Warning: Some NAs remain after imputation!\n")
#'   }
#'
#'   return(imputed_dataset)
#' }

#' Impute a dataset using missRanger with fallback
#'
#' @param dataset Data frame with ID column
#' @param dataset_name Name for logging
#' @param reference_data Optional reference data for type detection (use training data for validation)
#' @return Imputed data frame
impute_dataset <- function(dataset, dataset_name, reference_data = NULL) {
  cat("\n--- Imputing", dataset_name, "---\n")
  cat("  Variables:", ncol(dataset) - 1, "\n")
  cat("  Observations:", nrow(dataset), "\n")

  # Separate ID column
  ids <- dataset$ID
  data_no_id <- dataset[, -1, drop = FALSE]

  # Check if there are any missing values
  n_missing <- sum(is.na(data_no_id))
  cat("  Missing values:", n_missing, "\n")

  if (n_missing == 0) {
    cat("  No missing values - skipping imputation\n")
    return(dataset)
  }

  # Use reference data for type detection if provided (for validation set)
  if (is.null(reference_data)) {
    reference_data <- dataset
  }

  # Try missForest first
  cat("  Attempting missForest imputation...\n")
  set.seed(42)

  imputation_result <- try(missRanger(data_no_id,
                                        pmm.k = 5, # Predictive Mean Matching (key fix)
                                        num.trees = 500,
                                        min.node.size = 5),
                            silent = TRUE)

  if (!inherits(imputation_result, "try-error")) {
    cat("  ✓ missForest succeeded\n")
    imputed_data <- imputation_result

    # Post-process: ensure proper data types based on reference data
    cat("  Post-processing data types...\n")
    for (col in names(imputed_data)) {
      var_type <- detect_variable_type(reference_data[[col]])

      if (var_type == "binary") {
        # Round and clip to 0/1
        imputed_data[[col]] <- pmin(1, pmax(0, round(imputed_data[[col]])))
      } else if (var_type == "ordinal") {
        # Round to nearest integer
        imputed_data[[col]] <- as.integer(round(imputed_data[[col]]))
      }
    }

  } else {
    cat("  ✗ missRanger failed - using fallback\n")
    cat("  Error:", attr(imputation_result, "condition")$message, "\n")
    imputed_data <- fallback_impute(data_no_id)
  }

  # Rejoin ID column
  imputed_dataset <- cbind(ID = ids, imputed_data)

  # Verify no missing values remain
  remaining_na <- sum(is.na(imputed_dataset[, -1]))
  cat("  Remaining missing values:", remaining_na, "\n")

  if (remaining_na > 0) {
    cat("  ⚠ Warning: Some NAs remain after imputation!\n")
  }

  return(imputed_dataset)
}

# =============================== #
# 12. Impute Training and Validation Sets
# =============================== #

cat("\n=== Starting Imputation ===\n")

# imputation_result <- missRanger(training_data[,-1],
#                                 pmm.k = 5,     # Predictive Mean Matching (key fix)
#                                 num.trees = 500,
#                                 min.node.size = 5)


# Impute training set
training_data_imputed <- impute_dataset(
  dataset = training_data,
  dataset_name = "Training Set"
)

# Impute validation set (using training data as reference for type detection)
validation_data_imputed <- impute_dataset(
  validation_data,
  "Validation Set",
  reference_data = training_data
)

# =============================== #
# 13. Combine and Save Results
# =============================== #

cat("\n=== Combining Imputed Data ===\n")

# Combine training and validation back together
analysis_dataset_imputed <- rbind(training_data_imputed, validation_data_imputed)

# Sort by ID to maintain original order
analysis_dataset_imputed <- analysis_dataset_imputed[order(analysis_dataset_imputed$ID),]

cat("  Final dimensions (without ID):", paste(dim(analysis_dataset_imputed[, -1]), collapse = " x "), "\n")
cat("  Total missing values:", sum(is.na(analysis_dataset_imputed[, -1])), "\n")

# Save imputed dataset
write.csv(analysis_dataset_imputed, "analysis_dataset_imputed.csv", row.names = FALSE)

# Save training and validation sets separately for modeling
write.csv(training_data_imputed, "analysis_dataset_training_imputed.csv", row.names = FALSE)
write.csv(validation_data_imputed, "analysis_dataset_validation_imputed.csv", row.names = FALSE)

# =============================== #
# 14. Summary Statistics
# =============================== #

cat("\n=== Imputation Summary ===\n")

cat("\nBefore imputation (after dropping high-missing vars):\n")
cat("  Variables:", ncol(analysis_dataset_remaining_vars) - 1, "\n")
cat("  Observations:", nrow(analysis_dataset_remaining_vars), "\n")
cat("  Missing values:", n_missings_analysis_dataset_remaining_vars, "\n")
cat("  Missing %:", signif(n_missings_analysis_dataset_remaining_vars /
                               length(unlist(analysis_dataset_remaining_vars[, -1])) * 100, 3), "%\n")

cat("\nAfter split:\n")
cat("  Training observations:", nrow(training_data), "\n")
cat("  Validation observations:", nrow(validation_data), "\n")
cat("  Training missing values:", sum(is.na(training_data[, -1])), "\n")
cat("  Validation missing values:", sum(is.na(validation_data[, -1])), "\n")

cat("\nAfter imputation:\n")
cat("  Variables:", ncol(analysis_dataset_imputed) - 1, "\n")
cat("  Total observations:", nrow(analysis_dataset_imputed), "\n")
cat("  Training observations:", nrow(training_data_imputed), "\n")
cat("  Validation observations:", nrow(validation_data_imputed), "\n")
cat("  Total missing values:", sum(is.na(analysis_dataset_imputed[, -1])), "\n")
cat("  Training missing values:", sum(is.na(training_data_imputed[, -1])), "\n")
cat("  Validation missing values:", sum(is.na(validation_data_imputed[, -1])), "\n")
cat("  Missing %:", signif(sum(is.na(analysis_dataset_imputed[, -1])) /
                               length(unlist(analysis_dataset_imputed[, -1])) * 100, 3), "%\n")

cat("\n=== Imputation Complete ===\n")
cat("\nOutput files created:\n")
cat("  - analysis_dataset_imputed.csv (combined imputed data)\n")
cat("  - analysis_dataset_training_imputed.csv (training set only)\n")
cat("  - analysis_dataset_validation_imputed.csv (validation set only)\n")


# =============================== #
# 15. Check imputation of trigeminal psychophysics with many missings
# =============================== #

# After imputation
observed_CO2_thresholds <- training_data$`CO2 threshold`[!is.na(training_data$`CO2 threshold`)]
imputed_CO2_thresholds <- training_data_imputed$`CO2 threshold`[is.na(training_data$`CO2 threshold`)]

# 1. Overlay density plot
df_plot_CO2_thresholds <- data.frame(
  value = c(observed_CO2_thresholds, imputed_CO2_thresholds),
  type = rep(c("Observed", "Imputed"), c(length(observed), length(imputed)))
)

p_observed_imputed_CO2_thresholds <- ggplot(df_plot_CO2_thresholds, aes(x = value, fill = type, color = type)) +
  geom_density(alpha = 0.2, size = 1) +
  scale_fill_manual(values = c("cornsilk2", "cornsilk4")) +
  scale_color_manual(values = c("cornsilk2", "cornsilk4")) +
  labs(title = "CO2 Threshold: Observed vs imputed distributions",
        x = "CO2 Threshold (ms)", y = "Density") +
  theme_plot() +
  theme(legend.position.inside = TRUE, legend.position = c(.7, .8), legend.direction = "horizontal")

print(p_observed_imputed_CO2_thresholds)
ggsave("p_observed_imputed_CO2_thresholds.svg", p_observed_imputed_CO2_thresholds,
        width = 8, height = 8)

df_plot_CO2_thresholds_wide <- df_plot_CO2_thresholds %>%
  group_by(type) %>%
  mutate(idx = row_number()) %>%
  ungroup() %>%
  pivot_wider(
    names_from = type,
    values_from = value,
    values_fill = NA_real_
  ) %>%
  select(-idx)


pPDE_observed_imputed_CO2_thresholds <- PDEplotGG(df_plot_CO2_thresholds_wide) +
  theme_plot() +
  labs(title = "CO2 Threshold: Observed vs imputed distributions", x = "CO2 threshold (ms)", y = "PDE", color = "type") +
  scale_color_manual(values = c("grey83", "cornsilk4"), labels = c("Obsvered", "Imputed")) +
  theme(legend.position.inside = TRUE, legend.position = c(.7, .8), legend.direction = "horizontal")

print(pPDE_observed_imputed_CO2_thresholds)


# 2. Anderson-Darling test (H0: same distribution)
ad_result_CO2_thresholds <- twosamples::ad_test(observed_CO2_thresholds, observed_CO2_thresholds)
print(ad_result_CO2_thresholds) # p > 0.05 = distributions similar

observed_Lateralization <- training_data$`Lateralization (x/20)`[!is.na(training_data$`Lateralization (x/20)`)]
imputed_Lateralization <- training_data_imputed$`Lateralization (x/20)`[is.na(training_data$`Lateralization (x/20)`)]

# 1. Overlay density plot
df_plot_Lateralization <- data.frame(
  value = c(observed_Lateralization, imputed_Lateralization),
  type = rep(c("Observed", "Imputed"), c(length(observed), length(imputed)))
)

p_observed_imputed_Lateralization <- ggplot(df_plot_Lateralization, aes(x = value, fill = type, color = type)) +
  geom_density(alpha = 0.2, size = 1) +
  scale_fill_manual(values = c("cornsilk2", "cornsilk4")) +
  scale_color_manual(values = c("cornsilk2", "cornsilk4")) +
  labs(title = "Lateralization: Observed vs imputed distributions",
        x = "Lateralization [n correct]", y = "Density") +
  theme_plot() +
  theme(legend.position.inside = TRUE, legend.position = c(.2, .8), legend.direction = "vertical")

print(p_observed_imputed_Lateralization)
ggsave("p_observed_imputed_Lateralization.svg", p_observed_imputed_Lateralization,
        width = 8, height = 8)

df_plot_Lateralization_wide <- df_plot_Lateralization %>%
  group_by(type) %>%
  mutate(idx = row_number()) %>%
  ungroup() %>%
  pivot_wider(
    names_from = type,
    values_from = value,
    values_fill = NA_real_
  ) %>%
  select(-idx)


pPDE_observed_imputed_Lateralization <- PDEplotGG(df_plot_Lateralization_wide) +
  theme_plot() +
  labs(title = "Lateralization: Observed vs imputed distributions", x = "Lateralization [n correct]", y = "PDE", color = "type") +
  scale_color_manual(values = c("grey83", "cornsilk4"), labels = c("Obsvered", "Imputed")) +
  theme(legend.position.inside = TRUE, legend.position = c(.2, .8), legend.direction = "horizontal")

print(pPDE_observed_imputed_Lateralization)

# 2. Anderson-Darling test (H0: same distribution)
ad_result_Lateralization <- twosamples::ad_test(observed_Lateralization, observed_Lateralization)
print(ad_result_Lateralization) # p > 0.05 = distributions similar


# =============================== #
# 16. Prepare for Python-based check howto transform interval scaled data
# =============================== #

training_data_to_transform <- training_data_imputed[, which(sapply(training_data_imputed, detect_variable_type) == "continuous")][, -1]
names(training_data_to_transform)

# Write data to file for transformation analysis in Python
write.csv(x = training_data_to_transform, file = "training_data_to_transform.csv")


