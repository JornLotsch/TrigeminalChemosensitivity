# Load necessary libraries
library(readxl)
library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(forcats)
library(scales)

# Read Excel file with trigeminal sensitivity data
trigeminale_daten_table1 <- read_excel(
  "/home/joern/Aktuell/TrigeminalSensitivity/09Originale/Bormann Trigeminale Studie Daten.xlsx",
  sheet = "Tabelle1"
)

# Variables of interest related to facial pain
facial_pain_vars <- c("Gesichtsschmerzen", "R26", "R27")

# Subset the dataframe for these variables
facial_pain_data <- trigeminale_daten_table1[, facial_pain_vars]

# Get current year as numeric
current_year <- year(Sys.Date())

# Clean 'Gesichtsschmerzen' by removing quotes and trimming whitespace
clean_pain <- str_trim(str_replace_all(facial_pain_data$Gesichtsschmerzen, '"', ""))

# Function to extract numeric years anywhere in a string
extract_years <- function(x) {
  years <- str_extract_all(x, "\\b\\d{4}\\b")[[1]]
  if (length(years) == 0) {
    return(NA_integer_)
  } else {
    return(as.integer(years))
  }
}

# Find the minimum year across all entries, to be used for "always" cases
all_years <- unlist(lapply(clean_pain, extract_years))
min_year <- min(all_years, na.rm = TRUE)

# Function to parse start and end years from 'Gesichtsschmerzen' string
parse_dates <- function(txt) {
  txt <- tolower(txt)

  # Handle "j" indicating 'yes' with no date -- assign current year as start and end
  if (txt == "j") {
    return(c(start = current_year, end = current_year))
  }

  # Handle "seit vielen jahren" or "schon immer" as very long durations
  if (str_detect(txt, "seit vielen jahren") | str_detect(txt, "schon immer")) {
    return(c(start = min_year - 2, end = current_year))
  }

  # Handle "seit <year>" meaning ongoing pain from year to current
  if (str_detect(txt, "^seit\\s*\\d{4}")) {
    year_extracted <- as.integer(str_extract(txt, "\\d{4}"))
    return(c(start = year_extracted, end = current_year))
  }

  # Handle single year entries
  if (str_detect(txt, "^\\d{4}$")) {
    year_extracted <- as.integer(txt)
    return(c(start = year_extracted, end = year_extracted))
  }

  # Handle entries containing one year somewhere else
  years_found <- extract_years(txt)
  if (length(years_found) == 1) {
    return(c(start = years_found, end = years_found))
  }

  # If no year found, return NA
  return(c(start = NA_integer_, end = NA_integer_))
}


# Apply date parsing function to all rows
date_matrix <- t(sapply(clean_pain, parse_dates))

# Combine parsed years back into the original dataframe along with R26 and R27
result_df <- facial_pain_data %>%
  mutate(
    start_year = date_matrix[, "start"],
    end_year = date_matrix[, "end"]
  ) %>%
  select(Gesichtsschmerzen, start_year, end_year, R26, R27)

result_df$Prb <- trigeminale_daten_table1$Probandennummer

# Split comma-separated values in R27 into multiple columns
max_parts <- max(str_count(result_df$R27, ",") + 1, na.rm = TRUE)
split_pains <- str_split_fixed(result_df$R27 %>% as.character(), ",", max_parts)
colnames(split_pains) <- paste0("painsensation", 1:max_parts)
split_pains <- apply(split_pains, 2, str_trim) %>% as.data.frame(stringsAsFactors = FALSE)

# Bind the split pain sensation columns to result_df
result_df <- bind_cols(result_df, split_pains)

# Translation maps for frequency and pain sensations from German to English
translate_freq <- c(
  "monatlich" = "monthly",
  "wöchentlich" = "weekly",
  "immer" = "always",
  "mehrfach täglich" = "multiple daily",
  "einmal täglich" = "once daily",
  "NA" = "unknown"
)

translate_sensation <- c(
  "ziehend" = "pulling",
  "stechend" = "stabbing",
  "drückend" = "pressing",
  "brennend" = "burning"
)

# Apply translation maps to the relevant columns, handling NAs
result_df <- result_df %>%
  mutate(
    pain_freq_eng = translate_freq[R26] %>% replace_na("unknown"),
    painsensation1_eng = translate_sensation[painsensation1] %>% replace_na(NA_character_),
    painsensation2_eng = translate_sensation[painsensation2] %>% replace_na(NA_character_)
  )

# Combine frequency and sensations into a single descriptive factor variable
result_df <- result_df %>%
  mutate(
    pain_combo = case_when(
# Both painsensation1 and painsensation2 present and non-empty
      !is.na(painsensation1_eng) & painsensation1_eng != "" &
        !is.na(painsensation2_eng) & painsensation2_eng != "" ~
        paste0(pain_freq_eng, ", ", painsensation1_eng, " and ", painsensation2_eng),

# Only painsensation1 present
      !is.na(painsensation1_eng) & painsensation1_eng != "" ~
        paste(pain_freq_eng, painsensation1_eng, sep = ", "),

# Only painsensation2 present
      !is.na(painsensation2_eng) & painsensation2_eng != "" ~
        paste(pain_freq_eng, painsensation2_eng, sep = ", "),

# Neither present, just frequency
      TRUE ~ pain_freq_eng
    ),
    pain_combo = str_squish(pain_combo),
    pain_combo = factor(pain_combo)
  )

# Prepare data for plotting, filtering out missing start years
plot_df <- result_df %>%
  filter(!is.na(start_year)) %>%
  mutate(
    person = factor(Prb),
    bar_start = start_year,
    bar_end = ifelse(is.na(end_year), start_year, end_year),
    pain_combo = factor(pain_combo)
  )

# Plot facial pain duration bars colored by combined pain frequency and sensation
ggplot(plot_df, aes(
  y = fct_rev(person),
  xmin = bar_start,
  xmax = bar_end,
  color = pain_combo
)) +
  geom_errorbarh(height = 0.3, linewidth = 2) +
  scale_color_viridis_d(option = "plasma", na.value = "grey80") +
  labs(
    x = "Year",
    y = "Person",
    color = "Pain frequency and sensation",
    title = "Facial pain duration colored by frequency and sensation"
  ) +
  theme_minimal(base_size = 11) +
  guides(color = guide_legend(ncol = 1))
