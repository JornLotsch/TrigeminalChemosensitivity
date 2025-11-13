################################################################################
# Trigeminal Sensitivity Analysis - Bormann Study
# Author: Jorn Lotsch
# Description: Full analysis pipeline for trigeminal sensitivity study data.
#              Covers data import, variable categorization, type inspection,
#              data corrections, translations, and final export.
################################################################################

# ======================== #
# 1. Load Required Libraries
# --------------------------#
# Load external packages needed for reading Excel files, string manipulation,
# data wrangling, and descriptive statistics.
library(readxl)
library(stringr)
library(dplyr)
library(tidyr)
library(psych)
library(purrr)

# =============================== #
# 2. Data Import & Metadata Setup
# -------------------------------#
# Load raw study data from Excel file.
trigeminale_daten_raw_original <- read_excel(
  "/home/joern/Aktuell/TrigeminalSensitivity/09Originale/Bormann Trigeminale Studie Daten.xlsx",
  sheet = "Tabelle1"
)

# Load variable metadata, including category assignments and labels.
variable_categories <- read_excel(
  "/home/joern/Aktuell/TrigeminalSensitivity/09Originale/Bormann Trigeminale Studie Daten_Categories.xlsx",
  sheet = "einteilung_4_R"
)

# Load global mappings (e.g., category labels, translation dictionaries).
source("globals.R")

# Assign human-readable category names using preloaded labels.
variable_categories$Category_name <- category_labels[as.character(variable_categories$Category)]

# Map German variable names to English equivalents.
# If translation is missing, retain original German name.
variable_categories$Variable_english <- translation_map[variable_categories$Variable_german]
variable_categories$Variable_english[is.na(variable_categories$Variable_english)] <-
  variable_categories$Variable_german[is.na(variable_categories$Variable_english)]

# Prepare lists of variables split by category for future processing.
category_names <- make.names(unique(variable_categories$Category_name))
varlist <- split(variable_categories$Variable_english, variable_categories$Category_name)
names(varlist) <- category_names

# =============================== #
# 3. Inspect Variable Types and Content
# -------------------------------#
# Summarize and display current data types of all variables in the dataset.
# Output grouped lists of variables by their class/type to ease data checking.
var_types_1 <- sapply(trigeminale_daten_raw_original, class)
type_list <- split(names(var_types_1), var_types_1)

for (type in names(type_list)) {
  cat(paste0(type, " Variables:\n"))       # Header showing type name
  cat(paste0(type_list[[type]], collapse = "\n"))  # Each variable name on separate line
  cat("\n\n")
}

# =============================== #
# 4. Convert Variables and Data Corrections
# -------------------------------#
# Explicitly convert specific variables originally read as characters to numeric type.
# This ensures subsequent numeric analyses work correctly without coercion issues.
character_to_numeric_vars <- c(
  "Alter",
  "Körpergröße",
  "Riechvermögen unmittelbar nach Covid-19",
  "Wie oft Covid-19?",
  "Trinken Sie Alkohol?",
  "Wenn ich einen frischen Minzkaugummi gekaut habe, habe ich das Gefühl besser Luft durch die Nase zu bekommen",
  "Wenn ja, wie begann das Problem",
  "Wie hat sich das Problem verändert?"
)

trigeminale_daten_raw_original[character_to_numeric_vars] <-
  lapply(trigeminale_daten_raw_original[character_to_numeric_vars], as.numeric)

# Re-inspect data types post-conversion to confirm changes.
var_types_2 <- sapply(trigeminale_daten_raw_original, class)
type_list <- split(names(var_types_2), var_types_2)

for (type in names(type_list)) {
  cat(paste0(type, " Variables:\n"))
  cat(paste0(type_list[[type]], collapse = "\n"))
  cat("\n\n")
}

# Identify and print variables for which data type changed due to conversion.
changed_vars <- names(var_types_1[var_types_1 != var_types_2])

cat("Changed Variable Types:\n")
for (var in changed_vars) {
  cat(paste0(var, ": ", var_types_1[var], " → ", var_types_2[var], "\n"))
}

# =============================== #
# Additional Data Cleaning Steps
# -------------------------------#
# Define a function to set percentage variables to NA wherever original data
# is non-numeric or missing; fixes false zero values imported from Excel.
set_percent_na_if_not_numeric_multi <- function(df, col_pairs) {
  for (pair in col_pairs) {
    original_col <- pair[1]
    percent_col <- pair[2]
    non_numeric_mask <- is.na(as.numeric(df[[original_col]]))
    df[[percent_col]][non_numeric_mask] <- NA
  }
  return(df)
}

pairs_percent_from_other_variable <- list(
  c("Riechvermögen vor Covid",  "R1 in %"),
  c("Riechvermögen unmittelbar nach Covid-19", "R2 in %"),
  c("Derzeitiges Riechvermögen",  "R3 in %"),
  c("Nasenatmung für beide Nasenlöcher", "R23"),
  c("Nasenatmung für das rechte Nasenloch", "R24"),
  c("Nasenatmung für das linke Nasenloch", "R25")
)

trigeminale_daten_raw_original <- set_percent_na_if_not_numeric_multi(trigeminale_daten_raw_original, pairs_percent_from_other_variable)

# Check for implausible or outlier values, correcting or flagging as necessary.
print("Indices with 'Alter' < 18:")
print(which(trigeminale_daten_raw_original$Alter < 18))

print("Indices with 'Gewicht' < 30 and setting to NA:")
low_weight_idx <- which(trigeminale_daten_raw_original$Gewicht < 30)
print(low_weight_idx)
trigeminale_daten_raw_original$Gewicht[low_weight_idx] <- NA

print("Indices with 'Körpergröße' < 90 and correction (+100):")
low_height_idx <- which(trigeminale_daten_raw_original$Körpergröße < 90)
print(low_height_idx)
trigeminale_daten_raw_original$Körpergröße[low_height_idx] <- trigeminale_daten_raw_original$Körpergröße[low_height_idx] + 100

# Address inconsistencies in nasal airflow variables: invalid values replaced with NA.
rows_to_update <- with(trigeminale_daten_raw_original, R23 < 1 & R24 < 1 & R25 < 1 & R28 > 10)
print("Rows where nasal airflow variables will be set to NA due to inconsistency:")
print(which(rows_to_update))
trigeminale_daten_raw_original[rows_to_update, c("R23", "R24", "R25")] <- NA

# Display distribution of lateralization variable for reference.
print("All obsvered lateralization results:")
print(table(trigeminale_daten_raw_original$`Lateralisierung (x/20)`))

# Replace COVID missing with "j" when there is COVID-related information.

# Convert 'Wie oft Covid-19?' to integer and find indices > 0
covid_gt_0 <- which(as.integer(trigeminale_daten_raw_original$`Wie oft Covid-19?`) > 0)

# For those indices, update 'Waren Sie bereits an Covid erkrankt?'
trigeminale_daten_raw_original$`Waren Sie bereits an Covid erkrankt?`[covid_gt_0] <- "j"

trigeminale_daten_raw_original$`Waren Sie bereits an Covid erkrankt?`

# =============================== #
# 5. Rename Variables to English
# -------------------------------#
# Create copy of the data with variable names translated to English,
# replacing German names where translations are available.
trigeminale_daten_raw_original_selected <- trigeminale_daten_raw_original
names(trigeminale_daten_raw_original_selected) <- ifelse(
  names(trigeminale_daten_raw_original_selected) %in% names(translation_map),
  translation_map[names(trigeminale_daten_raw_original_selected)],
  names(trigeminale_daten_raw_original_selected)
)

# Select only variables that have been successfully translated.
translated_vars <- intersect(names(trigeminale_daten_raw_original_selected), translation_map)
selected_vars_df <- trigeminale_daten_raw_original_selected[, translated_vars]

# =============================== #
# 6. List types of exported variables
# -------------------------------#
# Re-re-inspect data types post-conversion to confirm changes.
var_types_3 <- sapply(selected_vars_df, class)
type_list <- split(names(var_types_3), var_types_3)

for (type in names(type_list)) {
  cat(paste0(type, " Variables:\n"))
  cat(paste0(type_list[[type]], collapse = "\n"))
  cat("\n\n")
}

# =============================== #
# 7. Export Cleaned Data and Metadata
# -------------------------------#


# Save the final cleaned, translated dataset and variable metadata files for analysis.
write.csv(cbind.data.frame(ID = trigeminale_daten_raw_original_selected$Probandennummer,
                           selected_vars_df), "trigeminale_daten_corrected_translated.csv", row.names = FALSE)
write.csv(variable_categories, "trigeminale_daten_variable_categories.csv", row.names = FALSE)
