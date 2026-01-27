################################################################################
# COVID-19 Impact Analysis - Trigeminal Sensitivity Study
# Author: Jorn Lotsch
# Description: Comprehensive analysis of COVID-19 effects on olfactory function
#              - Parsing mixed-format date strings (seasons, years, ranges)
#              - Temporal visualization of smell ability changes
#              - Data cleaning and type conversions
#              - Timeline analysis of COVID-19 infections
################################################################################

# ========================================================================== #
# 1. LOAD REQUIRED LIBRARIES
# ========================================================================== #

library(stringr) # String manipulation and regex operations
library(dplyr) # Data manipulation and transformation
library(tidyr) # Data tidying and reshaping (pivot operations)
library(ggplot2) # Data visualization
library(datefixR) # Date fixing utilities
library(psych) # Psychological and statistical methods
library(readxl) # Excel file import

# ========================================================================== #
# 2. DATA IMPORT
# ========================================================================== #

# Read the main corrected and translated dataset
trigeminale_daten_corrected_translated <- read.csv(
  "trigeminale_daten_corrected_translated.csv",
  check.names = FALSE
)

# ========================================================================== #
# 3. DEFINE VARIABLES OF INTEREST
# ========================================================================== #

# Variables related to COVID-19 and olfactory function
# Includes infection status, smell reduction, smell ability measurements,
# and temporal information about infection periods
Covid_vars <- c(
  "Have you had COVID-19?",
  "Is or was there smell reduction after COVID-19?",
  "Smell ability before COVID-19",
  "Smell ability immediately after COVID-19",
  "Current smell ability",
  "How many times have you had COVID-19?",
  "Period 1",
  "Smell reduction 1",
  "Period 2",
  "Smell reduction 2",
  "Period 3",
  "Smell reduction 3"
)

# Subset the dataframe for these variables
Covid_data <- trigeminale_daten_corrected_translated[, Covid_vars]

# ========================================================================== #
# 4. HELPER FUNCTIONS
# ========================================================================== #

#' Convert mixed-format date strings to Date objects
#'
#' @param dates Character vector of date strings in various formats
#' @return Vector of Date objects
#' @details Handles multiple date formats:
#'   - Excel numeric dates
#'   - "Jahr YYYY" (Year YYYY)
#'   - German seasons with year: "Frühling 2020", "Winter 2021"
#'   - Two-digit years: "Winter 21"
#'   - Date ranges: "Zeitraum 2020-21"
#'   - Ambiguous dates: "01.02.2022 o. 02/23"
#'   - Missing/unknown values mapped to NA
convert_mixed_dates <- function(dates) {
  season_month_map <- list(
    "Frühling" = 4, "Fruehling" = 4, "Frühjahr" = 4, "Frühjar" = 4,
    "Frrühjahr" = 4, "Fruhjahr" = 4,
    "Sommer" = 7,
    "Herbst" = 10, "Jerbst" = 10,
    "Winter" = 1
  )

  # Clean and normalize date input strings
  dates_clean <- dates
  dates_clean <- sub("^jahr ?(\\d{4})$", "jahr \\1", dates_clean)
  dates_clean <- gsub("nicht bekannt|n\\.b\\.|unbekannt|na|^$", NA_character_, dates_clean)
  dates_clean <- gsub("zeitraum.*", NA_character_, dates_clean)

  out <- rep(as.Date(NA), length(dates))

  # Handle Excel numeric date values
  is_excel_num <- suppressWarnings(
    !is.na(as.numeric(dates)) &
      !is.na(as.Date(as.numeric(dates), origin = "1899-12-30"))
  )
  excel_indices <- which(is_excel_num)
  if (length(excel_indices) > 0) {
    out[excel_indices] <- as.Date(
      as.numeric(dates[excel_indices]),
      origin = "1899-12-30"
    )
  }

  # Handle "Jahr YYYY" format
  jahr_indices <- grep("^jahr \\d{4}$", dates_clean)
  if (length(jahr_indices) > 0) {
    years <- as.numeric(sub("jahr (\\d{4})", "\\1", dates_clean[jahr_indices]))
    out[jahr_indices] <- as.Date(paste0(years, "-07-01"))
  }

  # Handle German seasons followed by year
  season_indices <- grep("^(Frühling|Fruehling|Frühjahr|Frühjar|Frrühjahr|Fruhjahr|Sommer|Herbst|Jerbst|Winter) \\d{4}$", dates_clean)

  if (length(season_indices) > 0) {
    for (i in season_indices) {
      parts <- unlist(strsplit(dates_clean[i], " "))
      season <- gsub("frrühjahr|fruhjahr|frühjar|fruehling", "Frühjahr", parts[1])
      season <- gsub("jerbst", "Herbst", season)
      year <- as.numeric(parts[2])
      month <- season_month_map[[season]]
      if (!is.null(month)) {
        out[i] <- as.Date(sprintf("%04d-%02d-01", year, month))
      }
    }
  }

  # Handle strings that are only a year (e.g., "2020")
  year_only_indices <- grep("^\\d{4}$", dates_clean)
  if (length(year_only_indices) > 0) {
    years <- as.numeric(dates_clean[year_only_indices])
    out[year_only_indices] <- as.Date(paste0(years, "-07-01"))
  }

  # Handle two-digit years for seasons (e.g., "Winter 21")
  handle_two_digit_year <- function(ind, season_name, month_num) {
    if (length(ind) > 0) {
      for (i in ind) {
        year_suffix <- as.numeric(sub(paste0("^", season_name, " ?"), "", dates_clean[i]))
        year <- ifelse(year_suffix >= 50, 1900 + year_suffix, 2000 + year_suffix)
        out[i] <- as.Date(paste0(year, sprintf("-%02d-01", month_num)))
      }
    }
  }
  handle_two_digit_year(grep("^winter ?\\d{2}$", dates_clean), "winter", 1)
  handle_two_digit_year(grep("^(frühjahr|frühling) ?\\d{2}$", dates_clean), "frühjahr", 4)
  handle_two_digit_year(grep("^sommer ?\\d{2}$", dates_clean), "sommer", 7)

  # Handle date ranges like "Zeitraum YYYY-YY"
  range_indices <- grep("^zeitraum (\\d{4})-(\\d{2,4})$", dates_clean)
  if (length(range_indices) > 0) {
    for (i in range_indices) {
      years_str <- sub("zeitraum (\\d{4})-(\\d{2,4})", "\\1-\\2", dates_clean[i])
      parts <- unlist(strsplit(years_str, "-"))
      start_year <- as.numeric(parts[1])
      end_year_raw <- parts[2]
      if (nchar(end_year_raw) == 2) {
        end_year <- ifelse(
          as.numeric(end_year_raw) > 50,
          1900 + as.numeric(end_year_raw),
          2000 + as.numeric(end_year_raw)
        )
      } else {
        end_year <- as.numeric(end_year_raw)
      }
      start_date <- as.Date(paste0(start_year, "-07-01"))
      end_date <- as.Date(paste0(end_year, "-07-01"))
      mid_date <- start_date + floor((end_date - start_date) / 2)
      out[i] <- mid_date
    }
  }

  # Handle ambiguous date strings like "01.02.2022 o. 02/23"
  ambiguous_indices <- grep("\\d{2}\\.\\d{2}\\.\\d{4} o", dates_clean)
  if (length(ambiguous_indices) > 0) {
    for (i in ambiguous_indices) {
      first_date_str <- sub("^(\\d{2}\\.\\d{2}\\.\\d{4}).*$", "\\1", dates[i])
      first_date <- as.Date(first_date_str, format = "%d.%m.%Y")
      out[i] <- first_date
    }
  }

  return(out)
}

#' Clean COVID-19 data with appropriate type conversions
#'
#' @param df Data frame containing COVID-19 variables
#' @return Cleaned data frame with converted types
#' @details Performs the following conversions:
#'   - Date variables (Period 1/2/3) converted using convert_mixed_dates()
#'   - Yes/no responses converted to logical TRUE/FALSE
#'   - Percentage values converted to numeric
#'   - Count variables converted to numeric
clean_covid_data <- function(df) {
  # Variable groups by expected types
  date_vars <- grep("Period", names(df), value = TRUE)
  yes_no_vars <- c(
    "Have you had COVID-19?",
    "Is or was there smell reduction after COVID-19?",
    "Smell reduction 1",
    "Smell reduction 2",
    "Smell reduction 3"
  )
  percent_vars <- c(
    "Smell ability before COVID-19",
    "Smell ability immediately after COVID-19",
    "Current smell ability"
  )
  count_vars <- c("How many times have you had COVID-19?")

  clean_df <- df

  # Convert Period variables using mixed date converter
  for (v in date_vars) {
    clean_df[[v]] <- convert_mixed_dates(dates = clean_df[[v]])
  }

  # Convert yes/no responses to logical TRUE/FALSE, NA if unknown
  yes_no_map <- function(x) {
    x <- tolower(as.character(x))
    ifelse(x %in% c("j", "ja", "yes"), TRUE,
           ifelse(x %in% c("n", "nein", "no"), FALSE, NA))
  }
  for (v in yes_no_vars) {
    if (v %in% names(clean_df)) {
      clean_df[[v]] <- yes_no_map(clean_df[[v]])
    }
  }

  # Convert percent columns to numeric, suppress warnings
  for (v in percent_vars) {
    if (v %in% names(clean_df)) {
      clean_df[[v]] <- suppressWarnings(as.numeric(clean_df[[v]]))
    }
  }

  # Convert count column to numeric
  for (v in count_vars) {
    if (v %in% names(clean_df)) {
      clean_df[[v]] <- suppressWarnings(as.numeric(clean_df[[v]]))
    }
  }

  return(clean_df)
}

# ========================================================================== #
# 5. DATA CLEANING AND PREPARATION
# ========================================================================== #

# Clean the COVID-19 data
cleaned_covid_data <- clean_covid_data(df = Covid_data)

# Replace zeros with NA in smell ability columns (considered data errors)
cols <- c(
  "Smell ability before COVID-19",
  "Smell ability immediately after COVID-19",
  "Current smell ability"
)
cleaned_covid_data[cols] <- lapply(cleaned_covid_data[cols], function(x) {
  x[x == 0] <- NA
  x
})

# ========================================================================== #
# 6. PREPARE DATA FOR VISUALIZATION
# ========================================================================== #

# Prepare plotting dataframe
df <- cleaned_covid_data

# Add unique ID per case for grouping
df$ID <- seq_len(nrow(df))

# Filter rows: keep cases with confirmed COVID-19, or if unknown status,
# any Period date present (indicating infection despite missing status)
df_filtered <- df %>%
  filter(
    `Have you had COVID-19?` == TRUE |
      (is.na(`Have you had COVID-19?`) &
        (!is.na(`Period 1`) | !is.na(`Period 2`) | !is.na(`Period 3`)))
  )

# ========================================================================== #
# 7. RESHAPE DATA TO LONG FORMAT
# ========================================================================== #

# Manually reshape the data since column structure doesn't have numeric suffixes
# Create long format by binding three separate subsets (one per time period)
plot_data <- bind_rows(
# Time period 1: Before COVID-19
  df_filtered %>%
    transmute(
      ID = ID,
      `Have you had COVID-19?` = `Have you had COVID-19?`,
      `Is or was there smell reduction after COVID-19?` =
        `Is or was there smell reduction after COVID-19?`,
      `How many times have you had COVID-19?` =
        `How many times have you had COVID-19?`,
      Period = "1",
      Date = `Period 1`,
      R_value = `Smell ability before COVID-19`
    ),
# Time period 2: Immediately after COVID-19
  df_filtered %>%
    transmute(
      ID = ID,
      `Have you had COVID-19?` = `Have you had COVID-19?`,
      `Is or was there smell reduction after COVID-19?` =
        `Is or was there smell reduction after COVID-19?`,
      `How many times have you had COVID-19?` =
        `How many times have you had COVID-19?`,
      Period = "2",
      Date = `Period 2`,
      R_value = `Smell ability immediately after COVID-19`
    ),
# Time period 3: Current smell ability
  df_filtered %>%
    transmute(
      ID = ID,
      `Have you had COVID-19?` = `Have you had COVID-19?`,
      `Is or was there smell reduction after COVID-19?` =
        `Is or was there smell reduction after COVID-19?`,
      `How many times have you had COVID-19?` =
        `How many times have you had COVID-19?`,
      Period = "3",
      Date = `Period 3`,
      R_value = `Current smell ability`
    )
) %>%
# Filter out rows with missing dates or smell ability values
filter(!is.na(Date), !is.na(R_value))

# ========================================================================== #
# 8. CALCULATE MINIMUM OLFACTORY FUNCTION PER CASE
# ========================================================================== #

# Calculate minimum olfactory function per case for color mapping
# This identifies subjects with the most severe smell reduction
min_values <- plot_data %>%
  group_by(ID) %>%
  summarise(min_R = min(R_value, na.rm = TRUE))

# Join min_R to plotting data
plot_data <- plot_data %>%
  left_join(min_values, by = "ID")

# ========================================================================== #
# 9. CREATE MAIN VISUALIZATION
# ========================================================================== #

# Plot: lines and points per case, color scaled by minimum olfactory function
# Darker colors indicate more severe smell reduction
p_Covid <- ggplot(
  plot_data,
  aes(x = Date, y = R_value, group = as.factor(ID), color = min_R)
) +
  geom_point() +
  geom_line() +
  scale_color_gradient(
    high = "cornsilk3",
    low = "ivory4",
    na.value = "grey80"
  ) +
  theme_light(base_size = 8) +
  labs(
    x = "Date",
    y = "Olfactory function [%]",
    color = "Olfactory\nFunction [%]",
    title = "Olfactory function in individuals with history of COVID-19"
  ) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1, vjust = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(
      color = "black",
      size = 0.3,
      linetype = "dashed"
    ),
    legend.position.inside = TRUE,
    legend.position = c(0.9, 0.15)
  )

# Display plot
print(p_Covid)

# Save plot as SVG
ggsave(paste0("p_Covid", ".svg"), p_Covid, width = 8, height = 8)

# ========================================================================== #
# 10. TEMPORAL ANALYSIS - TIME SINCE INFECTION
# ========================================================================== #

# Calculate how long ago the COVID-19 infections occurred
# Study period: mid-point between November 2023 and May 2024
actual_year <- as.Date("2023-11-01") +
  (as.Date("2024-05-01") - as.Date("2023-11-01")) / 2

# Calculate months since each infection period
df_times_in_past <- df_filtered %>%
  mutate(
    months_since_covid_period_1 = as.numeric(actual_year - `Period 1`) / 30.44,
    months_since_covid_period_2 = as.numeric(actual_year - `Period 2`) / 30.44,
    months_since_covid_period_3 = as.numeric(actual_year - `Period 3`) / 30.44
  )
names(df_times_in_past)

# Conditional updates
# Correct number of times COVID according to information about olfactory changes
# Load necessary library for date validation
library(lubridate)

# Helper: check if x is a date or not - TRUE if it is a date
is_date <- function(x) {
  inherits(x, "Date") | inherits(x, "POSIXct") | inherits(x, "POSIXlt")
}

df_times_in_past <- df_times_in_past %>%
  rowwise() %>%
  mutate(
    covid_count_calc = sum(
      (!is_date(`Period 1`) | (!is.na(`Smell reduction 1`) & `Smell reduction 1` > 0)),
      (!is_date(`Period 2`) | (!is.na(`Smell reduction 2`) & `Smell reduction 2` > 0)),
      (!is_date(`Period 3`) | (!is.na(`Smell reduction 3`) & `Smell reduction 3` > 0))
    ),
    `How many times have you had COVID-19?` = if_else(
      covid_count_calc > `How many times have you had COVID-19?`,
      covid_count_calc,
      `How many times have you had COVID-19?`
    )
  ) %>%
  ungroup()


# Case 1: If "How many times have you had COVID-19?" == 1
df_times_in_past$months_since_covid_period_2[df_times_in_past$`How many times have you had COVID-19?` == 1] <- Inf
df_times_in_past$months_since_covid_period_3[df_times_in_past$`How many times have you had COVID-19?` == 1] <- Inf
df_times_in_past$`Smell reduction 2`[df_times_in_past$`How many times have you had COVID-19?` == 1] <- 0
df_times_in_past$`Smell reduction 3`[df_times_in_past$`How many times have you had COVID-19?` == 1] <- 0

# Case 2: If "How many times have you had COVID-19?" == 2
df_times_in_past$months_since_covid_period_3[df_times_in_past$`How many times have you had COVID-19?` == 2] <- Inf
df_times_in_past$`Smell reduction 3`[df_times_in_past$`How many times have you had COVID-19?` == 2] <- 0



# Display descriptive statistics
cat("\nDescriptive statistics of time since infection:\n")
psych::describe(df_times_in_past)

# Tabulate COVID-19 status
cat("\nCOVID-19 status distribution:\n")
table(df_times_in_past$`Have you had COVID-19?`)

# Save the temporal analyis
write.csv(df_times_in_past, "covid_times_in_past.csv", row.names = FALSE)

# ========================================================================== #
# END OF SCRIPT
# ========================================================================== #
