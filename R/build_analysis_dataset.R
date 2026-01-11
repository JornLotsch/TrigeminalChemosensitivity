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
    mutate(across(all_of(added_cols), ~ifelse(is.na(.), fill_value, .)))

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
trigeminale_daten_corrected_translated <- trigeminale_daten_corrected_translated[ , !duplicated(names(trigeminale_daten_corrected_translated))]

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
  left_join(trigeminale_daten_corrected_translated[,c("ID", cat1_variables_directly_included)], by = "ID")

Category_1_data <- cbind.data.frame(
  Category_1_data,
  one_hot_encode_var(
    trigeminale_daten_corrected_translated[,variables_by_categories$Demographics[variables_by_categories$Demographics %in% c("Gender")]],
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
Category_2_data["Has migraine changed over the last 10 years"] <- 2-trigeminale_daten_corrected_translated$`Has migraine changed over the last 10 years`
print(var_names_wo_ID(Category_2_data))
Category_2_data["How often do you have migraine per month"] <- trigeminale_daten_corrected_translated$`How often do you have migraine per month`

## Actual respiratory tract infections
Category_2_data["Upper respiratory infection"] <- convert_jn_to_numeric(trigeminale_daten_corrected_translated$`Upper respiratory infection`)

# Intermediate check for other category 2 variables not yet included
category_2_vars_not_included <- setdiff(variables_by_categories$Disorders_or_health_complaints, names(Category_2_data))
head(trigeminale_daten_corrected_translated[,category_2_vars_not_included], n = 2)

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
head(trigeminale_daten_corrected_translated[,category_2_vars_not_included], n = 2)

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
head(trigeminale_daten_corrected_translated[,category_2_vars_not_included])

## Facial pain
facial_pain_preprocessed <- read.csv("facial_pain_preprocessed.csv", check.names = FALSE)
head(facial_pain_preprocessed)
dim(facial_pain_preprocessed)
print(dim_wo_ID(facial_pain_preprocessed))
print(var_names_wo_ID(facial_pain_preprocessed))

facial_pain_preprocessed["Has facial pain"] <- ifelse(!is.na(facial_pain_preprocessed$end_year),
                                                      ifelse(facial_pain_preprocessed$actual_facial_pain == 1, 2, 1), 0)
facial_pain_vars <- c("ID", "Has facial pain","facial_pain_freq_code", "facial_pain_pulling", "facial_pain_stabbing", "facial_pain_pressing", "facial_pain_burning")

Category_2_data <- Category_2_data %>%
  left_join(facial_pain_preprocessed[facial_pain_vars], by = "ID")
head(Category_2_data)
max(table(Category_2_data$ID))

# Intermediate check for other category 2 variables not yet included
category_2_vars_not_included <- setdiff(variables_by_categories$Disorders_or_health_complaints, names(Category_2_data))
head(trigeminale_daten_corrected_translated[,category_2_vars_not_included])

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
print(union(names(covid_times_in_past),variables_by_categories$COVID19_history))
print("Variables differing between original data and processed COVID data")
print(union(setdiff(names(covid_times_in_past),variables_by_categories$COVID19_history),
            setdiff(variables_by_categories$COVID19_history, names(covid_times_in_past))))

logical_columns <- sapply(covid_times_in_past, is.logical)
covid_times_in_past[logical_columns] <- lapply(covid_times_in_past[logical_columns], as.numeric)

# Compute row-wise maximum, ensuring `Inf` handling
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

# Compute row-wise minimum, ensuring `Inf` handling

covid_times_in_past["Shortest time since covid"] <-
  apply(covid_times_in_past[,c("months_since_covid_period_1","months_since_covid_period_2","months_since_covid_period_3")],1, function(x) min(x,na.rm = TRUE))

covid_times_in_past[["Shortest time since covid"]][is.infinite(covid_times_in_past[["Shortest time since covid"]]) ] <- NA
dim(covid_times_in_past)

# Assemble category_3 data
covid_vars_zero_in_no_covid <- c("Have you had COVID-19?", "How many times have you had COVID-19?" ,"Smell ability immediately after COVID-19" ,
                                 "Is or was there smell reduction after COVID-19?", "Smell ability before COVID-19",
                                 "Smell reduction 1", "Smell reduction 2", "Smell reduction 3")

Category_3_data <- analysis_dataset_empty %>%
  left_join_fill(covid_times_in_past[,c("ID", covid_vars_zero_in_no_covid)], by = "ID")
head(Category_3_data)
dim(Category_3_data)
max(table(Category_3_data$ID))

covid_vars_Inf_in_no_covid <- c("months_since_covid_period_1", "months_since_covid_period_2", "months_since_covid_period_3", "Longest time since covid", "Shortest time since covid")

Category_3_data <- Category_3_data %>%
  left_join_fill(covid_times_in_past[,c("ID", covid_vars_Inf_in_no_covid)], by = "ID")
head(Category_3_data)
dim(Category_3_data)
max(table(Category_3_data$ID))


print("Number of variables in category 2: COVID-19 history")
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

head(trigeminale_daten_corrected_translated[,variables_by_categories$Smoking_and_alcohol_use])

Category_4_data <- analysis_dataset_empty

Category_4_data[variables_by_categories$Smoking_and_alcohol_use[1]] <-
  convert_jn_to_numeric(trigeminale_daten_corrected_translated$`Do you smoke?`)
Category_4_data[variables_by_categories$Smoking_and_alcohol_use[4]] <-
  convert_jn_to_numeric(trigeminale_daten_corrected_translated$`Have you ever smoked?`)


print("Reading preprocesses smoking behavior related variables")

smoking_summary <- read.csv("smoking_summary.csv", check.names = FALSE)
dim(smoking_summary)

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

Category_4_data <- Category_4_data %>%
  left_join(smoking_summary[c("ID", "mean_cigs", "pack_years", "time_since_quitting_smoking", "smoker_type_current", "smoker_type_former")], by = "ID")

# Fill missing values for rows without smoking period info:
Category_4_data <- Category_4_data %>%
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

Category_4_data <- Category_4_data %>%
  left_join(trigeminale_daten_corrected_translated[c("ID", "Do you drink alcohol?")])
head(Category_4_data)
max(table(Category_4_data$ID))

print("Number of variables in category 4: Smoking and alcohol use")
print(nVars_wo_ID(Category_4_data))
print("Variables in category 4: Smoking and alcohol use")
print(var_names_wo_ID(Category_4_data))

# Add variables to the variable_names list
analysis_dataset_variables <- rbind.data.frame(analysis_dataset_variables, cbind.data.frame(category = 4, variable_name = var_names_wo_ID(Category_4_data)))

# =============================== #
#
# Category 5: Subjective nasal chemosensory perception (TriFunQ and related)
#
# =============================== #

print("Nasal_chemosensory_perception")
print(variables_by_categories$Nasal_chemosensory_perception)

Category_5_data <- analysis_dataset_empty %>%
  left_join(trigeminale_daten_corrected_translated[,c("ID", variables_by_categories$Nasal_chemosensory_perception)], by = "ID")
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
apply(variables_rated_olfactory_function,2,table)

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
  left_join(variables_rated_olfactory_function[,c("ID", "Current smell ability", "How has the problem changed", "If yes, how did the problem start")], by = "ID")

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
  left_join(trigeminale_daten_corrected_translated[,c("ID", variables_by_categories$Nasal_irritation_and_airflow)], by = "ID")
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
  left_join(trigeminale_daten_corrected_translated[,c("ID", variables_by_categories$Psychophysical_measurements)], by = "ID")
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
      row_values <- sapply(df[row, ], as.character)
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
dim(analysis_dataset[,-1])

# =============================== #
# 8. Handling of missing data
# =============================== #

# Check which variables or cases have missings more than acceptable
print("Checking which variables or cases have missings more than acceptable.")

n_missings_per_variable <- apply(analysis_dataset,2,function(x) sum(is.na(x)))

f_missings <- 0.2
f_missings_per_variable <- n_missings_per_variable / nrow(analysis_dataset)
vars_to_drop <- names(analysis_dataset)[which(f_missings_per_variable > f_missings)]
cat("\n", paste0("Variables with more than 20% missing values: "))
print(vars_to_drop)

analysis_dataset_remaining_vars <- analysis_dataset[,!names(analysis_dataset) %in% vars_to_drop]

n_missings_analysis_dataset_remaining_vars <- sum(is.na(analysis_dataset_remaining_vars))

cat("\n", paste0("Total number of missings in remaining dataset: "))
print(n_missings_analysis_dataset_remaining_vars)
print(n_missings_analysis_dataset_remaining_vars / length(unlist(analysis_dataset_remaining_vars)) * 100)

