# Load necessary libraries
library(readxl)
library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(datefixR)
library(psych)

# Read Excel file with trigeminal sensitivity data
trigeminale_daten_table1 <- read_excel(
  "/home/joern/Aktuell/TrigeminalSensitivity/09Originale/Bormann Trigeminale Studie Daten.xlsx",
  sheet = "Tabelle1"
)

# Correct zeros when percentages were calculated from empty cells
set_percent_na_if_not_numeric_multi <- function(df, col_pairs) {
  # col_pairs: Named list or two-column matrix/data.frame,
  # where each element/pair contains c("original_col", "percent_col")

  for (pair in col_pairs) {
    original_col <- pair[1]
    percent_col <- pair[2]

    # Check non-numeric entries
    non_numeric_mask <- is.na(as.numeric(df[[original_col]]))

    # Set respective percent column entries to NA where original is non-numeric
    df[[percent_col]][non_numeric_mask] <- NA
  }

  return(df)
}

pairs_percent_from_other_variable = list(
  c("Riechvermögen vor Covid",  "R1 in %"),
  c("Riechvermögen unmittelbar nach Covid-19", "R2 in %"),
  c("Derzeitiges Riechvermögen",  "R3 in %")
)

trigeminale_daten_table1 <- set_percent_na_if_not_numeric_multi(trigeminale_daten_table1, pairs_percent_from_other_variable)

# Variables of interest related to Covid and olfactory function
Covid_vars <- c(
  "Waren Sie bereits an Covid erkrankt?",
  "Besteht oder bestand eine Riechminderung nach Covid-19?",
  "R1 in %",
  "R2 in %",
  "R3 in %",
  "Wie oft Covid-19?",
  "Zeitraum 1",
  "Riechminderung?",
  "Zeitraum 2",
  "Riechminderung?2",
  "Zeitraum 3",
  "Riechminderung?3"
)

# Subset the dataframe for these variables
Covid_data <- trigeminale_daten_table1[, Covid_vars]

# Function to convert mixed-format date strings to Date objects
convert_mixed_dates <- function(dates) {
  # Mapping of German seasons to mid-season months
  season_month_map <- list(
    "Frühling" = 4, "Fruehling" = 4, "Frühjahr" = 4, "Frühjar" = 4, "Frrühjahr" = 4, "Fruhjahr" = 4,
    "Sommer" = 7,
    "Herbst" = 10, "Jerbst" = 10,
    "Winter" = 1
  )

  # Clean and normalize date input strings
  dates_clean <- tolower(dates)
  dates_clean <- sub("^jahr ?(\\d{4})$", "jahr \\1", dates_clean)
  dates_clean <- gsub("frrühjahr|fruhjahr|frühjar", "frühjahr", dates_clean)
  dates_clean <- gsub("nicht bekannt|n\\.b\\.|unbekannt|na|^$", NA_character_, dates_clean)
  dates_clean <- gsub("^jahr ?(\\d{3})$", NA_character_, dates_clean)
  dates_clean <- gsub("zeitraum.*", NA_character_, dates_clean)

  out <- rep(as.Date(NA), length(dates))

  # Handle Excel numeric date values
  is_excel_num <- suppressWarnings(!is.na(as.numeric(dates)) & !is.na(as.Date(as.numeric(dates), origin = "1899-12-30")))
  excel_indices <- which(is_excel_num)
  if (length(excel_indices) > 0) {
    out[excel_indices] <- as.Date(as.numeric(dates[excel_indices]), origin = "1899-12-30")
  }

  # Handle "Jahr YYYY"
  jahr_indices <- grep("^jahr \\d{4}$", dates_clean)
  if (length(jahr_indices) > 0) {
    years <- as.numeric(sub("jahr (\\d{4})", "\\1", dates_clean[jahr_indices]))
    out[jahr_indices] <- as.Date(paste0(years, "-07-01"))
  }

  # Handle German seasons followed by year
  season_indices <- grep("^(frühling|frühjahr|sommer|herbst|winter) \\d{4}$", dates_clean)
  if (length(season_indices) > 0) {
    for (i in season_indices) {
      parts <- unlist(strsplit(dates_clean[i], " "))
      season <- gsub("[^a-zA-Z]", "", parts[1])
      year <- as.numeric(parts[2])
      month <- season_month_map[[season]]
      if (!is.null(month)) {
        out[i] <- as.Date(sprintf("%04d-%02d-01", year, month))
      }
    }
  }

  # Handle strings that are only a year
  year_only_indices <- grep("^\\d{4}$", dates_clean)
  if (length(year_only_indices) > 0) {
    years <- as.numeric(dates_clean[year_only_indices])
    out[year_only_indices] <- as.Date(paste0(years, "-07-01"))
  }

  # Handle two-digit years for Winter, Frühjahr, Sommer
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
        end_year <- ifelse(as.numeric(end_year_raw) > 50, 1900 + as.numeric(end_year_raw), 2000 + as.numeric(end_year_raw))
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

# Function to clean Covid data variables with appropriate conversions
clean_covid_data <- function(df) {
  # Variable groups by expected types
  date_vars <- grep("Zeitraum", names(df), value = TRUE)
  yes_no_vars <- c("Waren Sie bereits an Covid erkrankt?", "Besteht oder bestand eine Riechminderung nach Covid-19?",
                   "Riechminderung?", "Riechminderung?2", "Riechminderung?3")
  percent_vars <- c("R1 in %", "R2 in %", "R3 in %")
  count_vars <- c("Wie oft Covid-19?")

  clean_df <- df

  # Convert Zeitraum variables using mixed date converter
  for (v in date_vars) {
    clean_df[[v]] <- convert_mixed_dates(clean_df[[v]])
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

# Clean the data
cleaned_covid_data <- clean_covid_data(Covid_data)

# Replace zeros with NA in R percentage columns (considered errors)
cols <- c("R1 in %", "R2 in %", "R3 in %")
cleaned_covid_data[cols] <- lapply(cleaned_covid_data[cols], function(x) {
  x[x == 0] <- NA
  x
})

# Prepare for plotting
df <- cleaned_covid_data

# Add unique ID per case for grouping
df$id <- seq_len(nrow(df))

# Filter rows: keep cases with Covid TRUE, or if unknown, any Zeitraum date present
df_filtered <- df %>%
  filter(`Waren Sie bereits an Covid erkrankt?` == TRUE |
           (is.na(`Waren Sie bereits an Covid erkrankt?`) &
              (!is.na(`Zeitraum 1`) | !is.na(`Zeitraum 2`) | !is.na(`Zeitraum 3`))))

# Reshape Zeitraum dates and corresponding R values into long format, paired by period
plot_data <- df_filtered %>%
  pivot_longer(
    cols = c(starts_with("Zeitraum"), starts_with("R")),
    names_to = c(".value", "Period"),
    names_pattern = "(Zeitraum|R) ?(\\d)"
  ) %>%
  rename(Date = Zeitraum, R_value = R) %>%
  filter(!is.na(Date), !is.na(R_value))

# Calculate minimum olfactory function per case for color mapping
min_values <- plot_data %>%
  group_by(id) %>%
  summarise(min_R = min(R_value, na.rm = TRUE))

# Join min_R to plotting data
plot_data <- plot_data %>%
  left_join(min_values, by = "id")

# Plot: lines and points per case, color scaled by minimum olfactory function
p_Covid <- ggplot(plot_data, aes(x = Date, y = R_value, group = as.factor(id), color = min_R)) +
  geom_point() +
  geom_line() +
  scale_color_gradient(high = "cornsilk3", low = "ivory4", na.value = "grey80") +
  theme_light(base_size = 8) +
  labs(
    x = "Date",
    y = "Olfactory function [%]",
    color = "Olfactory\nFunction [%]",
    title = "Olfactory function in individuals with history of Covid-19"
  ) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1, vjust = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "black", size = 0.3, linetype = "dashed"),
    legend.position.inside = TRUE,
    legend.position = c(0.9, 0.15)
  )

p_Covid

ggsave(paste0("p_Covid", ".svg"), p_Covid, width = 8, height = 8)


# Calculate how long ago were the CCOVID-19 infections
df_times_in_past <- df_filtered %>%
  mutate(
    months_since_Zeitraum_1 = as.numeric(Sys.Date() - `Zeitraum 1`) / 30.44,
    months_since_Zeitraum_2 = as.numeric(Sys.Date() - `Zeitraum 2`) / 30.44,
    months_since_Zeitraum_3 = as.numeric(Sys.Date() - `Zeitraum 3`) / 30.44
  )

psych::describe(df_times_in_past)
table(df_times_in_past$`Waren Sie bereits an Covid erkrankt?`)
