# Load necessary libraries
library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(forcats)
library(scales)

# =============================== #
# Read data
# =============================== #

trigeminale_daten_corrected_translated <- read.csv("trigeminale_daten_corrected_translated.csv", check.names = FALSE)

character_vars <- names(trigeminale_daten_corrected_translated)[sapply(trigeminale_daten_corrected_translated, is.character)]
print(character_vars)

# Variables of interest related to facial pain
facial_pain_vars <- c("Facial pain", "How often do you have facial pain", "What is the nature of facial pain")

# Subset the dataframe for these variables
facial_pain_data <- trigeminale_daten_corrected_translated[, facial_pain_vars]

# Get current year as numeric
actual_year <- as.numeric(format(d <- as.Date("2023-11-01") + (as.Date("2024-05-01") - as.Date("2023-11-01")) / 2, "%Y")) + (as.numeric(format(d, "%j")) - 1) / ifelse(((y <- as.numeric(format(d, "%Y"))) %% 4 == 0 & y %% 100 != 0) | (y %% 400 == 0), 366, 365)


# Clean 'Gesichtsschmerzen' by removing quotes and trimming whitespace
clean_pain <- str_trim(str_replace_all(facial_pain_data$`Facial pain`, '"', ""))

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
    return(c(start = actual_year, end = actual_year))
  }

  # Handle "seit vielen jahren" or "schon immer" as very long durations
  if (str_detect(txt, "seit vielen jahren") | str_detect(txt, "schon immer")) {
    return(c(start = min_year - 2, end = actual_year))
  }

  # Handle "seit <year>" meaning ongoing pain from year to current
  if (str_detect(txt, "^seit\\s*\\d{4}")) {
    year_extracted <- as.integer(str_extract(txt, "\\d{4}"))
    return(c(start = year_extracted, end = actual_year))
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

# Combine parsed years back into the original dataframe along with "How often do you have facial pain", "What is the nature of facial pain"
result_df <- facial_pain_data %>%
  mutate(
    start_year = date_matrix[, "start"],
    end_year = date_matrix[, "end"]
  ) %>%
  dplyr::select("Facial pain", start_year, end_year, "How often do you have facial pain", "What is the nature of facial pain")

result_df$ID <- trigeminale_daten_corrected_translated$ID
result_df$actual_facial_pain <- ifelse(result_df$end_year < actual_year, 0, 1)

# Find max length of code strings to determine number of columns needed
max_parts <- max(nchar(as.character(result_df$`What is the nature of facial pain`)), na.rm = TRUE)

# Split each string into single characters, pad with NA for shorter strings
split_pains <- t(sapply(as.character(result_df$`What is the nature of facial pain`), function(x) {
  chars <- str_split(x, "")[[1]]
  length(chars) <- max_parts # pad with NAs if needed
  chars
}))

colnames(split_pains) <- paste0("painsensation", 1:max_parts)

# Convert to data frame with strings as characters (not factors)
split_pains_df <- as.data.frame(split_pains, stringsAsFactors = FALSE)

# Bind the split pain sensation columns to result_df
result_df <- bind_cols(result_df, split_pains)

table(result_df$`How often do you have facial pain`)

# Translation maps for frequency and pain sensations from German to English
translate_freq <- c(
  "beim passivrauchen" = "during passive smoking",
  "et" = "once a day",
  "i" = "always",
  "m" = "monthly",
  "mt" = "several times a day",
  "nur beim apfelessen" = "only when eating apples",
  "w" = "weekly",
  "NA" = "unknown frequency"
)

translate_sensation <- c(
  "z" = "pulling",
  "s" = "stabbing",
  "sz" = "stabbing,pulling",
  "d" = "pressing",
  "sd" = "stabbing,pressing",
  "b" = "burning",
  "db" = "pressing,burning",
  "dz" = "pressing,pulling",
  "zb" = "pulling,burning",
  "sdzb" = "stabbing,pressing,pulling,burning",
  "NA" = "unspecified sensation"
)


# Apply translation maps to the relevant columns, handling NAs
result_df$`How often do you have facial pain` <- tolower(as.character(result_df$`How often do you have facial pain`))
result_df$`What is the nature of facial pain` <- tolower(as.character(result_df$`What is the nature of facial pain`))
result_df <- result_df %>%
  mutate(
    pain_freq_eng = translate_freq[`How often do you have facial pain`] %>% replace_na("unknown frequency"),
    painsensation1_eng = translate_sensation[painsensation1] %>% replace_na(NA_character_),
    painsensation2_eng = translate_sensation[painsensation2] %>% replace_na(NA_character_),
    painsensation3_eng = translate_sensation[painsensation3] %>% replace_na(NA_character_),
    painsensation4_eng = translate_sensation[painsensation4] %>% replace_na(NA_character_)
  )

library(dplyr)
library(stringr)

result_df <- result_df %>%
  rowwise() %>%
  mutate(
# Collect non-empty, non-NA sensations
    pain_sensations_combined = {
  vals <- c(painsensation1_eng, painsensation2_eng, painsensation3_eng, painsensation4_eng)
  vals <- vals[!is.na(vals) & vals != ""]
  if (length(vals) > 0) {
    paste0(pain_freq_eng, ": ", paste(vals, collapse = ", "))
  } else {
    pain_freq_eng
  }
},
    pain_sensations_combined = str_squish(pain_sensations_combined)
  ) %>%
  ungroup()

result_df <- result_df %>%
  rowwise() %>%
  mutate(
    combined_sensations = {
  vals <- c(painsensation1_eng, painsensation2_eng, painsensation3_eng, painsensation4_eng)
  vals <- vals[vals != "" & !is.na(vals)]
  paste(vals, collapse = ", ")
}
  ) %>%
  ungroup()


result_df <- result_df %>%
  mutate(
    facial_pain_freq_code = case_when(
      `Facial pain` == "n" ~ 0L,
      pain_freq_eng == "during passive smoking" ~ 1L,
      pain_freq_eng == "only when eating apples" ~ 1L,
      pain_freq_eng == "monthly" ~ 2L,
      pain_freq_eng == "weekly" ~ 3L,
      pain_freq_eng == "once a day" ~ 4L,
      pain_freq_eng == "several times a day" ~ 5L,
      pain_freq_eng == "always" ~ 6L,
      pain_freq_eng == "unknown frequency" ~ NA_integer_,
      TRUE ~ NA_integer_ # default to NA if no match
    )
  )

# Prepare data for plotting, filtering out missing start years
plot_df <- result_df %>%
  filter(!is.na(start_year)) %>%
  mutate(
    person = factor(ID),
    bar_start = start_year,
    bar_end = ifelse(is.na(end_year), start_year, end_year),
    pain_sensations_combined = factor(pain_sensations_combined),

  )

# Plot facial pain duration bars colored by combined pain frequency and sensation
p_facial_pain <- ggplot(plot_df, aes(
  y = fct_rev(person),
  xmin = bar_start,
  xmax = bar_end,
  color = pain_sensations_combined
)) +
  geom_errorbar(
    aes(xmin = bar_start, xmax = bar_end),
    width = 0.3,
    linewidth = 2,
    orientation = "y"
  ) +
  scale_color_viridis_d(option = "plasma", na.value = "grey80") +
  labs(
    x = "Year",
    y = "Person",
    color = "Pain frequency and sensation",
    title = "Facial pain duration colored by frequency and sensation"
  ) +
  theme_minimal(base_size = 9) +
  theme(
    legend.position.inside = TRUE, legend.position = c(.15, .8),
    legend.title = element_text(size = 5),
    legend.text = element_text(size = 6),
    legend.key.size = unit(0.3, "cm")
  ) +
  guides(color = guide_legend(override.aes = list(size = 1), ncol = 1))

p_facial_pain

ggsave(paste0("p_facial_pain", ".svg"), p_facial_pain, width = 9, height = 9)


# One-hot

one_hot_df <- result_df %>%
  dplyr::select(ID, combined_sensations) %>%
  filter(!is.na(combined_sensations) & combined_sensations != "") %>%
  separate_rows(combined_sensations, sep = ",\\s*") %>%
  mutate(combined_sensations = str_squish(combined_sensations)) %>%
  distinct(ID, combined_sensations) %>%
  mutate(present = 1) %>%
  pivot_wider(
    names_from = combined_sensations,
    values_from = present,
    values_fill = list(present = 0) # fill missing with 0 instead of NA
  ) %>%
  rename_with(~paste0("facial_pain_", .x), - ID) # add prefix except ID

result_df <- result_df %>%
  left_join(one_hot_df, by = "ID") %>%
  mutate(across(starts_with("facial_pain_"), ~ replace_na(., 0)))



write.csv(result_df, "facial_pain_preprocessed.csv", row.names = FALSE)
