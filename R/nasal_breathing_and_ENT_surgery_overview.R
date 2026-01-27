################################################################################
# ENT Surgery History an Nasal Breathing Problems Analysis -
# Trigeminal Sensitivity Study
# Author: Jorn Lotsch
# Description: Comprehensive analysis of ear, nose, and throat surgical
#              procedures over time
#              - Parsing surgical procedure descriptions with years from free text
#              - German to English translation of medical procedures
#              - Temporal visualization of surgical interventions
#              - One-hot encoding for statistical modeling
################################################################################



# ========================================================================== #
# Load required libraries ----------------------------------------------------
# ========================================================================== #

library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(forcats)
library(scales)
library(purrr)

source("globals.R")


# ========================================================================== #
# CONFIGURATION SWITCHES
# ========================================================================== #

# Control whether to order plot by procedure frequency (descending)
# Set to TRUE to sort procedures by count, FALSE for alphabetical
order_plot <- TRUE


# =============================== #
# Read and inspect data
# =============================== #

trigeminale_daten_corrected_translated <- read.csv(
  "trigeminale_daten_corrected_translated.csv",
  check.names = FALSE
)

character_vars <- names(trigeminale_daten_corrected_translated)[
  sapply(trigeminale_daten_corrected_translated, is.character)
]
print(character_vars)

# ========================================================================== #
# A. Nasal Breathing Problems Analysis
# ========================================================================== #


# ========================================================================== #
# 1. Select nasal breathing variables -------------------------------------------
# ========================================================================== #

nasal_breathing_vars <- c(
  "Nasal airflow (both nostrils)",
  "Nasal airflow (right nostril)",
  "Nasal airflow (left nostril)",
  "Medical consultation for nasal breathing problems",
  "If yes, was therapy performed (and which one)"
)

nasal_breathing_df <- trigeminale_daten_corrected_translated[, c("ID", nasal_breathing_vars)] %>%
  rename(
    ID = ID,
    raw_text = `Medical consultation for nasal breathing problems`,
    therapy_text = `If yes, was therapy performed (and which one)`
  ) %>%
  mutate(age = as.integer(trigeminale_daten_corrected_translated$Age))

# Numeric representation of actual year (midpoint between given dates)
actual_year <- as.numeric(format(d <- as.Date("2023-11-01") + (as.Date("2024-05-01") - as.Date("2023-11-01")) / 2, "%Y")) +
  (as.numeric(format(d, "%j")) - 1) / ifelse(((y <- as.numeric(format(d, "%Y"))) %% 4 == 0 & y %% 100 != 0) | (y %% 400 == 0), 366, 365)

# ========================================================================== #
# 2. Robust year extraction function --------------------------------------------
# ========================================================================== #

extract_years_custom <- function(text, age, actual_year) {
  # Return NA for missing or "n" text
  if (is.na(text) || text == "n") return(NA_character_)

  text <- stringr::str_trim(text)
  text_no_j <- stringr::str_replace(text, "^j\\s+", "")

  # Handle "j MM/YY" format e.g. "j 01/23"
  if (stringr::str_detect(text, "^j\\s*\\d{2}/\\d{2}$")) {
    year_suffix <- as.numeric(stringr::str_sub(text_no_j, 4, 5))
    year_full <- ifelse(year_suffix > 30, 1900 + year_suffix, 2000 + year_suffix)
    return(as.character(year_full))
  }

  # Excel serial date (5+ digits)
  if (stringr::str_detect(text_no_j, "^\\d{5,}$")) {
    return(format(as.Date(as.numeric(text_no_j), origin = "1899-12-30"), "%Y"))
  }

  # Year range, e.g. "2014-2016"
  if (stringr::str_detect(text_no_j, "\\d{4}-\\d{4}")) {
    years <- as.numeric(stringr::str_split(text_no_j, "-", simplify = TRUE))
    return(as.character(seq(years[1], years[2])))
  }

  # Slash-separated years: e.g. "2018/19"
  if (stringr::str_detect(text_no_j, "\\d{4}/\\d{2}")) {
    first <- as.numeric(stringr::str_sub(text_no_j, 1, 4))
    second <- as.numeric(paste0(stringr::str_sub(first, 1, 2), stringr::str_sub(text_no_j, 6, 7)))
    return(as.character(c(first, second)))
  }

  # Plus-separated years: e.g. "1997 + 2023"
  if (stringr::str_detect(text_no_j, "\\d{4}\\s*\\+\\s*\\d{4}")) {
    return(stringr::str_extract_all(text_no_j, "\\d{4}")[[1]])
  }

  # Single year or phrases like "ab 2015", "Winter 2023"
  if (stringr::str_detect(text_no_j, "^(ab\\s*)?\\d{4}$|^Winter\\s*\\d{4}$")) {
    year <- stringr::str_extract(text_no_j, "\\d{4}")
    # If "ab YYYY" means ongoing from that year
    if (stringr::str_detect(text_no_j, "^ab\\s*\\d{4}$")) {
      return(c(year, "ongoing"))
    } else {
      return(year)
    }
  }

  # Handle explicit ongoing with year plus keywords such as "regelmäßig"
  if (stringr::str_detect(text, "\\d{4}\\s*regelmäßig")) {
    year <- stringr::str_extract(text, "\\d{4}")
    return(c(year, "ongoing"))
  }

  # "2024 regelmäßig" exact case handled implicitly above (regex matches)

  # Age-based forecasts
  if (text_no_j == "Grundschüler" && !is.na(age)) {
    return(as.character(round(actual_year - age + 8)))
  }
  if (text_no_j == "zur Geburt" && !is.na(age)) {
    return(as.character(round(actual_year - age)))
  }

  # Detect ongoing keywords alone
  if (stringr::str_detect(text_no_j, regex("regelmäßig|halbjährlich|Quartal|Kindesalter", ignore_case = TRUE))) {
    years_found <- stringr::str_extract_all(text_no_j, "\\b\\d{4}\\b")[[1]]
    if (length(years_found) == 0) {
      return("ongoing")
    } else {
      return(c(years_found, "ongoing"))
    }
  }

  # Extract 4-digit years anywhere as fallback
  years_found <- stringr::str_extract_all(text_no_j, "\\b\\d{4}\\b")[[1]]
  if (length(years_found) > 0) return(unique(years_found))

  return(NA_character_)
}


# ========================================================================== #
# 3. Helper function to strip years and flags -----------------------------------
# ========================================================================== #

strip_years_time <- function(text, pattern) {
  if (is.na(text) || is.na(pattern)) return(text)
  cleaned <- str_remove_all(text, pattern)
  cleaned <- str_squish(str_replace_all(cleaned, ",{2,}", ","))
  cleaned <- str_trim(cleaned)
  if (cleaned == "") return(NA_character_) else return(cleaned)
}

# ========================================================================== #
# 4. Apply extraction, ongoing flag and clean therapy ---------------------------
# ========================================================================== #

nasal_breathing_df <- nasal_breathing_df %>%
  rowwise() %>%
  mutate(
    estimated_birth_year = if_else(
      str_detect(raw_text, "Kindesalter"),
      round(actual_year - age),
      NA_real_
    ),
# Remove ongoing_start_year concept, use extracted 'ongoing' flag directly
    years_raw = list(extract_years_custom(raw_text, age, actual_year)),
    years_therapy = list(extract_years_custom(therapy_text, age, actual_year)),

# Combine years from raw and therapy and birth year estimate
    extracted_years = list(unique(na.omit(c(
      unlist(years_raw), unlist(years_therapy),
      if (!is.na(estimated_birth_year)) estimated_birth_year else NULL
    )))),

# Extract ongoing flag explicitly for each row based on both raw and therapy texts
    ongoing = any(c("ongoing") %in% c(years_raw, years_therapy)),

# Prepare years pattern regex for therapy_clean
    years_pattern = if (length(extracted_years) > 0) {
      paste0("\\b(", paste(extracted_years, collapse = "|"), ")\\b")
    } else {
  NA_character_
},

# Clean therapy text by removing years/time info
    therapy_clean = strip_years_time(therapy_text, years_pattern)

  ) %>%
  ungroup()

# Correct therapy yes/no if necessary
nasal_breathing_df$therapy_clean[
  (is.na(nasal_breathing_df$therapy_clean) | nasal_breathing_df$therapy_clean == "") &
    sapply(nasal_breathing_df$extracted_years, function(x) length(x) > 0 && !all(is.na(x)))
] <- "j"


# ========================================================================== #
# 5. Translation dictionary and translation function ----------------------------
# ========================================================================== #

# Get condensed_ENT_dictionary from globals.R

# breathing_and_ENT_surgery_therapy_dict <- c(
#   "n" = "No therapy",
#   "j" = "Unknown",
#   "OP Nasenscheidewand" = "Nasal septum surgery (septoplasty)",
#   "OP NSW" = "Nasal septum surgery (septoplasty)",
#   "OP" = "Operation",
#   "Septumplastik" = "Septoplasty (nasal septum surgery)",
#   "regelmäßige Untersuchungen" = "Regular checkups",
#   "Desensibilisierung" = "Desensitization",
#   "Cortisonspray" = "Cortisone spray",
#   "Verödung" = "Cauterization",
#   "Atemgerät Schlaf" = "Sleep breathing device",
#   "Antibiotika" = "Antibiotics",
#   "Nasenspray" = "Nasal spray",
#   "Bedarfsspray Asthma" = "Asthma rescue spray",
#   "Antiallergika" = "Antiallergics",
#   "Polypen-Op" = "Polyp operation",
#   "Asthmadiagnostik" = "Asthma diagnostics",
#   "Op geplant" = "Operation planned",
#   "Salbe" = "Ointment",
#   "Doprident" = "Doprident",
#   "Meersalzinhalation" = "Sea salt inhalation",
#   "Cortison-NS" = "Cortisone nasal spray",
#   "Creme" = "Cream",
#   "Mometason" = "Mometasone",
#   "AB" = "Antibiotics",
#   "Rotlichthterapie" = "Red light therapy",
#   "Stromwellentherapie" = "Radio wave therapy",
#   "Trommelfellschnitt" = "Eardrum incision",
#   "Akupunktur" = "Acupuncture",
#   "Hyposensibilisierung" = "Hyposensitization",
#   "Lasertherapie" = "Laser therapy",
#   "Begradigung NSW" = "Nasal septum straightening",
#   "Kortison-NS" = "Cortisone nasal spray",
#   "Septumplastik li." = "Septoplasty left side",
#   "NSW OP" = "Nasal septum surgery",
#   "Nasensalbe" = "Nasal ointment",
#   "Cortison-Nasenspray" = "Cortisone nasal spray",
#   "Allergie Nasenspray Mometa" = "Allergy nasal spray mometasone",
#   "Cortison NS" = "Cortisone nasal spray",
#   "Corticoidinhalation" = "Corticoid inhalation",
#   "Laserung Nasenmuscheln" = "Laser therapy nasal concha",
#   "Fluimveil" = "Fluimveil (medication)",
#   "Nasenmuschelverkleinerung" = "Nasal concha reduction",
#   "Stirnhöhlenvereiderung?" = "Unknown",
#   "Cortison NS, Desensibilisierung" = "Cortisone nasal spray and desensitization",
#   "Hyposensibilisierung, NT, Tabletten" = "Hyposensitization, nasal therapy, tablets",
#   "Cetirizin" = "Cetirizine"
# )

translate_therapy <- function(therapy_text) {
  if (is.na(therapy_text) || therapy_text == "") return(NA_character_)

  # Split therapy_text by commas into individual therapies
  therapies <- str_split(therapy_text, ",")[[1]]
  therapies <- str_trim(therapies)

  translated <- map_chr(therapies, function(single_therapy) {
    therapy_text_lower <- tolower(single_therapy)

    match_key <- purrr::map_lgl(names(condensed_ENT_dictionary), ~ grepl(tolower(.x), therapy_text_lower, fixed = TRUE))
    matched_keys <- names(condensed_ENT_dictionary)[match_key]

    if (length(matched_keys) == 0) {
      return(stringr::str_to_sentence(single_therapy))
    }

    # Choose the longest match for best specificity
    best_key <- matched_keys[which.max(nchar(matched_keys))]
    return(condensed_ENT_dictionary[best_key])
  })

  # Concatenate translated therapies back with commas
  return(paste(translated, collapse = ", "))
}

# Apply this function on therapy_clean column
nasal_breathing_df <- nasal_breathing_df %>%
  mutate(
    therapy_english = map_chr(therapy_clean, translate_therapy)
  )

# Specify generic operation
nasal_breathing_df$therapy_english[nasal_breathing_df$therapy_english == "Operation"] <- "Unspecified surgery for nasal breathing problems"


# ========================================================================== #
# 6. Flag nasal breathing problem ------------------------------------------------
# ========================================================================== #

nasal_breathing_df <- nasal_breathing_df %>%
  mutate(
    nasal_breathing_problem = if_else(
      raw_text != "n" | !is.na(therapy_text),
      1L,
      0L
    )
  )

nasal_breathing_df <- nasal_breathing_df %>%
  rowwise() %>%
  mutate(
# Filter years to numeric only
    numeric_years = list(as.character(Filter(function(x) grepl("^\\d{4}$", x), extracted_years))),
# Create years_clean as comma-separated years or "No year" if only ongoing present
    years_clean = if (length(numeric_years) == 0 && any("ongoing" %in% extracted_years)) {
      "No year"
    } else {
  paste(numeric_years, collapse = ",")
}
  ) %>%
  ungroup() %>%
  dplyr::select(-numeric_years)


# ========================================================================== #
# 7. Plot therapy occurrences ---------------------------------------------------
# ========================================================================== #

# Safely extract numeric years from extracted_years list-column
numeric_years <- nasal_breathing_df$extracted_years %>%
  unlist() %>%
  as.character() %>%
  keep(~grepl("^\\d{4}$", .)) %>% # keep only 4-digit year strings
as.numeric()

# Remove NA and provide fallback if empty
if (length(numeric_years) == 0 || all(is.na(numeric_years))) {
  min_year_value <- 2000 # fallback year if no valid years
} else {
  min_year_value <- min(numeric_years, na.rm = TRUE)
}

df_nasal_breathing_therapy_long <- nasal_breathing_df %>%
  filter(!is.na(therapy_english)) %>%
  dplyr::select(ID, therapy_english, years_clean, ongoing) %>%
  separate_rows(therapy_english, sep = ",\\s*") %>%
  separate_rows(years_clean, sep = ",\\s*") %>%
  mutate(
    therapy_english = str_trim(therapy_english),
    years_clean = str_trim(years_clean),
    year_num = suppressWarnings(as.numeric(years_clean)),
    pseudoyear_pos = min_year_value - 2,
    year_plot = if_else(is.na(year_num), pseudoyear_pos, year_num),
    therapy_english = factor(therapy_english, levels = sort(unique(therapy_english))),
    case_id = ID
  )

count_df <- df_nasal_breathing_therapy_long %>%
  group_by(therapy_english) %>%
  summarise(count = n()) %>%
  mutate(
    xpos = as.numeric(therapy_english),
    ypos = max(df_nasal_breathing_therapy_long$year_plot, na.rm = TRUE) + 3
  )

# Create alternating background stripe data for therapies
therapy_levels <- levels(unique(df_nasal_breathing_therapy_long$therapy_english))


# Reorder therapy_english factor by descending counts for x-axis
if (order_plot) {
  # Order by descending count
  therapy_levels_ordered <- count_df %>%
    arrange(desc(count)) %>%
    pull(therapy_english)
} else {
  # Keep original order in df_nasal_breathing_therapy_long
  therapy_levels_ordered <- levels(df_nasal_breathing_therapy_long$therapy_english)
}

# Apply the ordering to df_nasal_breathing_therapy_long
df_nasal_breathing_therapy_long <- df_nasal_breathing_therapy_long %>%
  mutate(therapy_english = factor(therapy_english, levels = therapy_levels_ordered))

# Also apply the same ordering to count_df, then calculate xpos and ypos accordingly
count_df <- count_df %>%
  mutate(therapy_english = factor(therapy_english, levels = therapy_levels_ordered)) %>%
  arrange(therapy_english) %>%
  mutate(
    xpos = as.numeric(therapy_english),
    ypos = max(df_nasal_breathing_therapy_long$year_plot, na.rm = TRUE) + 3
  )

stripe_df <- data.frame(
  therapy_english = therapy_levels_ordered,
  xmin = seq_along(therapy_levels_ordered) - 0.5,
  xmax = seq_along(therapy_levels_ordered) + 0.5
)

set.seed(42)
nasal_breathing_therapy_plot <- ggplot() +
  geom_rect(data = stripe_df,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = rep(c("ivory2", "white"), length.out = nrow(stripe_df)),
            inherit.aes = FALSE) +
  geom_jitter(data = df_nasal_breathing_therapy_long,
              aes(x = therapy_english, y = year_plot, shape = factor(ongoing)),
              color = "black", # all points blackish
              width = 0.4, height = 0,
              alpha = 0.7,
              size = 2,
              show.legend = TRUE) +
  geom_text(data = count_df,
            aes(x = xpos, y = ypos, label = count),
            vjust = 0, size = 3, angle = 90) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(
    breaks = c(min(df_nasal_breathing_therapy_long$year_plot, na.rm = TRUE),
               seq(min(df_nasal_breathing_therapy_long$year_plot, na.rm = TRUE) + 1,
                   max(df_nasal_breathing_therapy_long$year_plot, na.rm = TRUE))),
    labels = c("No year",
               as.character(seq(min(df_nasal_breathing_therapy_long$year_plot, na.rm = TRUE) + 1,
                                max(df_nasal_breathing_therapy_long$year_plot, na.rm = TRUE))))
  ) +
  scale_shape_manual(
    name = "Ongoing therapy",
    values = c("TRUE" = 17, "FALSE" = 19), # different point shapes
    labels = c("Ongoing", "Not ongoing")
  ) +
  labs(x = "Therapy", y = "Year",
       title = "Therapies for nasal breathing problems over years") +
  theme_minimal(base_size = 8) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "black", size = 0.3, linetype = "dashed"),
    legend.position = c(0.2, 0.2),
    legend.background = element_rect(fill = ggplot2::alpha("white", 0.5), color = NA)
  ) +
  geom_hline(
    yintercept = seq(min(df_nasal_breathing_therapy_long$year_plot, na.rm = TRUE),
                     max(df_nasal_breathing_therapy_long$year_plot, na.rm = TRUE)),
    color = "grey80", size = 0.3, linetype = "dashed"
  )

print(nasal_breathing_therapy_plot)


ggsave(
  filename = "p_nasal_breathing_therapy_plot.svg",
  plot = nasal_breathing_therapy_plot,
  width = 12, height = 12
)

# ========================================================================== #
# 8. Create one-hot encoding for therapy ----------------------------------------
# ========================================================================== #

# Starting from df_nasal_breathing_therapy_long (one row per ID and therapy)
nasal_breating_therpies_one_hot <- df_nasal_breathing_therapy_long %>%
  distinct(ID, therapy_english) %>% # keep unique ID-therapy pairs
mutate(present = 1) %>% # flag presence
pivot_wider(
    names_from = therapy_english,
    values_from = present,
    values_fill = list(present = 0)
  )

nasal_breating_therpies_one_hot_prefixed <- nasal_breating_therpies_one_hot %>%
  rename_with(~paste0("Nasal_breathing_therapy_", .),
              .cols = -ID)

# ========================================================================== #
# 9. Save tables ----------------------------------------
# ========================================================================== #

nasal_breathing_df_clean <- nasal_breathing_df %>%
  dplyr::select(where(~!is.list(.)))

write.csv(as.data.frame(nasal_breathing_df_clean), "nasal_breating_prepocessed.csv", row.names = FALSE)
write.csv(nasal_breating_therpies_one_hot_prefixed, "nasal_breating_therpies_one_hot.csv", row.names = FALSE)



# ========================================================================== #
# B. ENT surgery analysis, with corrections form above
# ========================================================================== #

# ========================================================================== #
# 2. HELPER FUNCTIONS
# ========================================================================== #

#' Parse year-description pairs from ENT surgical procedure text strings
#'
#' @param text_vec Character vector of procedure descriptions (German text)
#' @param age Numeric vector of patient ages (same length as text_vec)
#' @param actual_year Numeric value or vector for reference year
#' @return Data frame with alternating year/desc columns (year1, desc1, year2, desc2, ...)
#' @details Handles multiple year formats:
#'   - Standard 4-digit year: "2015 Tonsillotomie"
#'   - Slash notation: "/74", "/2021"
#'   - MM/YY format: "05/21" (extracts year only)
#'   - Decade notation: "90er Jahre" (maps to 1995)
#'   - Childhood keywords: "Kindesalter", "Kindheit", "Kleinkind"
#'   - Multiple procedures with '+': "2015 Adenotomie + Tonsillotomie"
parse_procedure_descriptions <- function(text_vec, age, actual_year) {

  parse_one <- function(text, pat_age, pat_year) {
    items <- stringr::str_split(text, ",\\s*")[[1]]
    items <- items[items != ""]
    items <- stringr::str_trim(items)
    items <- items[items != "n"]

    extract_year <- function(s) {
      s <- stringr::str_trim(s)
      if (is.na(s) || s == "") return(NA_character_)

      y <- stringr::str_extract(s, "^\\d{4}")
      if (!is.na(y)) return(y)

      y <- stringr::str_extract(s, "/(\\d{2,4})")
      if (!is.na(y)) {
        y <- sub("/", "", y)
        if (nchar(y) == 2) {
          y_num <- as.numeric(y)
          y <- ifelse(y_num < 50, paste0("20", sprintf("%02d", y_num)), paste0("19", y))
        }
        return(y)
      }

      y <- stringr::str_extract(s, "\\d{2}/(\\d{2})")
      if (!is.na(y)) {
        yy <- sub("\\d{2}/", "", y)
        y_num <- as.numeric(yy)
        y <- ifelse(y_num < 50, paste0("20", sprintf("%02d", y_num)), paste0("19", yy))
        return(y)
      }

      if (stringr::str_detect(s, "90er Jahre")) return("1995")

      if (stringr::str_detect(s, "(Kindesalter|Kindheit|Kleinkind)")) {
        if (!is.na(pat_age) && !is.na(pat_year)) {
          est_year <- round(pat_year - pat_age)
          return(as.character(est_year))
        }
      }

      NA_character_
    }

    years <- sapply(items, extract_year)

    descs <- vector("list", length(items))
    for (i in seq_along(items)) {
      desc_raw <- items[i]

      # Safe check for childhood indicator (avoid NA)
      has_child_time <- isTRUE(stringr::str_detect(desc_raw, "(Kindesalter|Kindheit|Kleinkind)"))

      desc_raw <- stringr::str_remove(desc_raw, "^\\d{4}\\s*")
      desc_raw <- stringr::str_remove(desc_raw, "/\\d{2,4}")
      desc_raw <- stringr::str_remove(desc_raw, "\\d{2}/\\d{2}")
      desc_raw <- stringr::str_remove(desc_raw, "\\s*\\d{4}$")
      desc_raw <- stringr::str_remove(desc_raw, "90er Jahre")
      desc_raw <- stringr::str_trim(desc_raw)

      # Split on '+'
      desc_parts <- stringr::str_trim(stringr::str_split(desc_raw, "\\s*\\+\\s*")[[1]])

      # Apply childhood removal if pattern found
      if (has_child_time) {
        desc_parts <- stringr::str_remove(desc_parts, "(Kindesalter|Kindheit|Kleinkind)")
        desc_parts <- stringr::str_trim(desc_parts)
      }

      descs[[i]] <- desc_parts
    }

    years_expanded <- unlist(mapply(function(y, d) rep(y, length(d)), years, descs, SIMPLIFY = FALSE))
    descs_expanded <- unlist(descs)
    as.vector(rbind(years_expanded, descs_expanded))
  }

  parsed <- Map(parse_one, text_vec, age, actual_year)
  max_len <- max(sapply(parsed, length))
  parsed_padded <- lapply(parsed, function(x) { length(x) <- max_len; x })
  df <- as.data.frame(do.call(rbind, parsed_padded), stringsAsFactors = FALSE)
  row.names(df) <- NULL
  half <- max_len / 2
  colnames(df) <- paste0(rep(c("year", "desc"), half), rep(seq_len(half), each = 2))
  return(df)
}


# ========================================================================== #
# 3. DATA IMPORT
# ========================================================================== #

# Read the main corrected and translated dataset
trigeminale_daten_corrected_translated <- read.csv(
  "trigeminale_daten_corrected_translated.csv",
  check.names = FALSE
)

# ========================================================================== #
# 4. CLEAN AND PREPARE DATA
# ========================================================================== #

# Clean ENT surgery column: replace "n" (no), empty, NA, NaN with NA
clean_text_vec <- ifelse(
  trigeminale_daten_corrected_translated$`Surgery in ENT region` %in% c("n", "", "NA", "NaN"),
  NA_character_,
  trigeminale_daten_corrected_translated$`Surgery in ENT region`
)

# ========================================================================== #
# 5. PARSE SURGICAL PROCEDURE DESCRIPTIONS
# ========================================================================== #

# Parse procedure descriptions into year/description pairs
result_procedures <- parse_procedure_descriptions(clean_text_vec, age = trigeminale_daten_corrected_translated$Age, actual_year)

# ========================================================================== #
# 6. RESHAPE DATA TO LONG FORMAT
# ========================================================================== #

# Convert wide format to long format for analysis
# ID column is preserved as separate column (not pivoted)
ENT_surgeries_long_df <- result_procedures %>%
# Ensure ID is present and add row index
mutate(
    ID = trigeminale_daten_corrected_translated$ID,
    rowid = row_number()
  ) %>%
# Pivot all year/desc columns to long format, excluding ID and rowid
pivot_longer(
    cols = -c(ID, rowid),
    names_to = c(".value", "set"),
    names_pattern = "(year|desc)(\\d+)"
  ) %>%
# Remove empty or NA descriptions
filter(!is.na(desc) & desc != "") %>%
# Convert year to numeric (may result in NA for missing years)
mutate(year = as.numeric(year))

# ========================================================================== #
# 7. FIX YEAR EXTRACTION ISSUES
# ========================================================================== #
# Some years may remain in description text due to parsing edge cases
# This section catches and corrects these remaining year patterns

ENT_surgeries_long_df_fixed <- ENT_surgeries_long_df %>%
# Handle leading 4-digit year in description
mutate(
    year_new = ifelse(
      str_detect(desc, "^\\d{4}"),
      as.numeric(str_extract(desc, "^\\d{4}")),
      year
    ),
    desc = str_remove(desc, "^\\d{4}\\s*")
  ) %>%
# Handle two-digit leading year (e.g., "05 Adenotomie")
mutate(
    year_new = ifelse(
      is.na(year_new) & str_detect(desc, "^(\\d{2})\\s"),
      as.numeric(ifelse(
        as.numeric(str_extract(desc, "^(\\d{2})")) < 50,
        paste0("20", str_extract(desc, "^(\\d{2})")),
        paste0("19", str_extract(desc, "^(\\d{2})"))
      )),
      year_new
    ),
    desc = str_remove(desc, "^(\\d{2})\\s")
  ) %>%
# Special case: "er Tonsillotomie" (German decade notation) maps to 1995
mutate(
    year_new = ifelse(
      is.na(year_new) & str_detect(desc, "^er\\s*Tonsillotomie"),
      1995,
      year_new
    ),
    desc = str_remove(desc, "^er\\s*")
  ) %>%
# Clean up extra whitespace
mutate(desc = str_trim(desc)) %>%
# Replace original year with corrected year
mutate(year = year_new) %>%
  dplyr::select(-year_new) %>%
# Remove empty descriptions
dplyr::filter(!is.na(desc) & desc != "")

# ========================================================================== #
# 8. GERMAN TO ENGLISH TRANSLATION MAP
# ========================================================================== #

# Get condensed_ENT_dictionary from globals.R

# Translation dictionary for ENT surgical procedures
# Grouped by procedure type and normalized to consolidate variants
# breathing_and_ENT_surgery_therapy_dict <- c(
#   # Adenoid procedures - consolidated variants
#   "Adenotomie" = "Adenotomy",
#   "Adenotomie Kindesalter" = "Adenotomy",
#   "Kindesalter Adenotomie" = "Adenotomy",
#   "Adentomie" = "Adenotomy",  # Typo correction
#
#   # Nasal trauma/fractures
#   "Nasenfraktur" = "Nasal fracture",
#   "Nasen # OP" = "Nasal fracture surgery",
#
#   # Tonsil procedures - distinguished by type (full vs partial removal)
#   "Mandeln" = "Tonsillectomy",            # Full removal
#   "Mandeln OP" = "Tonsillectomy",
#   "Mandel OP" = "Tonsillectomy",
#   "Mandel Op" = "Tonsillectomy",
#   "Mandel-OP" = "Tonsillectomy",
#   "Kindesalter Mandeln" = "Tonsillectomy",
#   "Tonsillektomie" = "Tonsillectomy",     # Explicit full removal
#
#   "Tonsillotomie" = "Tonsillotomy",       # Partial removal
#   "Tonsillotomie Kindheit" = "Tonsillotomy",
#   "Kindesalter Tonsillotomie" = "Tonsillotomy",
#   "Tonsillotomie Kindesalter" = "Tonsillotomy",
#   "Tonsillotomie jugend" = "Tonsillotomy",
#   "er Tonsillotomie" = "Tonsillotomy",
#
#   # Nasal polyps - grouped and normalized
#   "Polypen" = "Nasal polyps surgery",
#   "Kleinkind Polypen" = "Nasal polyps surgery",
#   "Kindheit Polypen" = "Nasal polyps surgery",
#   "Kindesalter Polypen" = "Nasal polyps surgery",
#   "CRS Polypenentfernung 2x" = "Nasal polyps surgery",
#
#   # Nasal growth/sinus surgeries
#   "Wucherungen Nase" = "Nasal growth surgery",
#   "NNH" = "Paranasal sinus surgery",
#   "NNH Op" = "Paranasal sinus surgery",
#   "NNH links" = "Left paranasal sinus surgery",
#
#   # Nasal septum surgeries and corrections - grouped and normalized
#   "NSW" = "Nasal septum surgery",
#   "NSW OP" = "Nasal septum surgery",
#   "Nasenscheidewand" = "Nasal septum surgery",
#   "Nasenscheidewand OP" = "Nasal septum surgery",
#   "Begradigung Nasenscheidewand" = "Nasal septum correction",
#   "Septum OP" = "Nasal septum surgery",
#   "Septumkorrektur" = "Nasal septum correction",
#   "Septumdeviation" = "Septal deviation correction",
#
#   # Fenestration and Caldwell-Luc procedures
#   "Fensterung" = "Fenestration",
#   "Fensterung Kieferhöhle" = "Caldwell-Luc procedure",
#   "Fensterung Kieferhöhle li." = "Left Caldwell-Luc procedure",
#
#   # Ear surgeries
#   "Knorpelplastik Ohr" = "Cartilage plasty (ear)",
#   "Ohr li. Implantat" = "Left ear implant",
#   "Ohr anlegen" = "Otoplasty",
#   "Otoplastik" = "Otoplasty",
#   "Austausch Hämmerchen Ohr" = "Malleus replacement (middle ear)",
#
#   # Tympanic membrane surgeries
#   "Trommelfell" = "Tympanic membrane surgery",
#   "Trommelfell 2x" = "Tympanic membrane surgery",
#   "Trommelfell-OP" = "Tympanic membrane surgery",
#   "Trommelfellimplantat" = "Tympanic membrane implant",
#   "Trommelfellverschluss" = "Tympanic membrane closure",
#   "Trommelfellschnitt" = "Myringotomy",
#
#   "Knalltrauma" = "Acoustic trauma",
#
#   # Vocal cord surgeries
#   "Stimmband OP" = "Vocal cord surgery",
#   "Stimmknötchen" = "Vocal cord nodules surgery",
#
#   # Oncological procedures
#   "Mundbodenkarzinom" = "Mouth floor carcinoma surgery",
#   "Parotisentfernung" = "Parotidectomy",
#   "Karzinom Hals" = "Neck carcinoma surgery",
#   "Medulläres Schildrüsenkarzinom" = "Medullary thyroid carcinoma surgery",
#   "Karzinom Rachen li." = "Left pharyngeal carcinoma surgery",
#   "OP Halsbereich" = "Neck surgery",
#   "Abszess Rachen" = "Pharyngeal abscess surgery",
#   "Trauma Hals" = "Neck trauma surgery",
#
#   # Miscellaneous surgical procedures - grouped consistently
#   "Begradigung Nase" = "Nasal straightening surgery",
#   "Schilddrüse" = "Thyroid gland surgery",
#   "Thyreodektomie" = "Thyroidectomy",
#   "Kiefernzyste" = "Jaw cyst removal",
#   "OP Polypen Stirnhöhle" = "Frontal sinus polyp surgery",
#
#   # Turbinate reduction procedures - harmonized treatments
#   "Nasenmuscheln Verkleinerung" = "Inferior turbinate reduction",
#   "Reduktion Nasenmuscheln" = "Inferior turbinate reduction",
#   "Verkleinerung Nasenmuscheln" = "Inferior turbinate reduction",
#   "Laserung Nasenmuscheln" = "Inferior turbinate laser therapy",
#   "Verödung Nasenmuschel" = "Inferior turbinate cauterization",
#   "Verödung Nasenmuscheln" = "Inferior turbinate cauterization",
#   "Stromwellentherapie zum Abschwellen der Nasenmuscheln" = "Inferior turbinate radiofrequency therapy",
#
#   # Placeholder/invalid entries
#   "j" = "Unspecified"
# )

# ========================================================================== #
# 9. TRANSLATE DESCRIPTIONS TO ENGLISH
# ========================================================================== #

# Apply translation map to procedure descriptions
ENT_surgeries_long_df_fixed <- ENT_surgeries_long_df_fixed %>%
  mutate(desc_english = condensed_ENT_dictionary[desc]) %>%
  mutate(desc_english = ifelse(is.na(desc_english), "Unspecified", desc_english)) %>%
  mutate(desc_english = factor(desc_english, levels = sort(unique(desc_english))))


# ========================================================================== #
# 10. ADD MISSING SURGERIES FOR NASAL BREATHING PROBLEMS
# ========================================================================== #

# Your nasal breathing surgical procedures list
nasal_breathing_surgical_procedures <- c(
  "Adenotomy",
  "Nasal polyps surgery",
  "Frontal sinus polyp surgery",
  "Nasal fracture surgery",
  "Paranasal sinus surgery",
  "Left paranasal sinus surgery",
  "Nasal septum surgery",
  "Nasal septum surgery (septoplasty)",
  "Nasal septum correction",
  "Nasal septoplasty",
  "Left nasal septoplasty",
  "Septal deviation correction",
  "Septal perforation repair",
  "Inferior turbinate reduction",
  "Inferior turbinate cauterization",
  "Nasal straightening surgery",
  "Nasal septum fenestration",
  "Ethmoid sinus surgery",
  "Pansinusitis surgery",
  "Nasal mucosa excision",
  "Tympanostomy tubes insertion",
  "Tympanostomy tube surgery",
  "Middle ear bone implant",
  "Operation",
  "Unspecified surgery for nasal breathing problems"
)

# Prepare data frames (ensure character/numeric types)
df1 <- df_nasal_breathing_therapy_long %>%
  mutate(
    therapy_english = as.character(therapy_english),
    years_clean = as.numeric(years_clean)
  )

df2 <- ENT_surgeries_long_df_fixed %>%
  mutate(
    desc_english = as.character(desc_english),
    year = as.numeric(year)
  )

# Filter df1 to therapy_english in nasal breathing procedures list
filtered_df1 <- df1 %>%
  filter(therapy_english %in% nasal_breathing_surgical_procedures)

# Find missing combinations: those in filtered_df1 not in df2 by ID+procedure+year
missing_entries <- filtered_df1 %>%
  anti_join(
    df2,
    by = c("ID" = "ID", "therapy_english" = "desc_english", "years_clean" = "year")
  )

# Create new rows matching df2's structure
new_rows <- missing_entries %>%
  transmute(
    ID = ID,
    rowid = NA_integer_,
    set = NA_character_,
    year = years_clean,
    desc = therapy_english,
    desc_english = therapy_english,
    year_plot = years_clean
  )

# Append new rows to the existing ENT surgeries data frame
ENT_surgeries_long_df_fixed_completed <- bind_rows(df2, new_rows)

# Optionally reorder and assign new rowid
ENT_surgeries_long_df_fixed_completed <- ENT_surgeries_long_df_fixed_completed %>%
  arrange(ID, year) %>%
  mutate(rowid = row_number())

# ENT_surgeries_long_df_fixed_completed now contains all original and missing rows

# ========================================================================== #
# 10. STORE DATA FOR ONE-HOT ENCODING
# ========================================================================== #

# Store data frame before filtering for later one-hot encoding
ENT_surgeries_long_df_4_onehot <- ENT_surgeries_long_df_fixed_completed

# ========================================================================== #
# 11. PREPARE DATA FOR VISUALIZATION
# ========================================================================== #

# Determine minimum year for handling missing years
min_year <- min(ENT_surgeries_long_df_fixed_completed$year, na.rm = TRUE)

# Create plotting year: map NA years to pseudo-year below minimum
ENT_surgeries_long_df_fixed_completed <- ENT_surgeries_long_df_fixed_completed %>%
  mutate(
    year_plot = as.numeric(year),
    year_plot = ifelse(is.na(year_plot), min_year - 2, year_plot)
  )

# Remove unknown/untranslated procedures from visualization
ENT_surgeries_long_df_fixed_completed <- ENT_surgeries_long_df_fixed_completed %>%
  filter(desc_english != "Unknown")

# Reset factor levels after filtering
desc_levels <- sort(unique(ENT_surgeries_long_df_fixed_completed$desc_english))
ENT_surgeries_long_df_fixed_completed$desc_english <- factor(ENT_surgeries_long_df_fixed_completed$desc_english, levels = desc_levels)

# Create data frame for alternating background stripes in plot
stripe_df <- data.frame(
  desc_english = desc_levels,
  xmin = seq_along(desc_levels) - 0.5,
  xmax = seq_along(desc_levels) + 0.5
)

# Calculate procedure counts for annotation
count_df <- ENT_surgeries_long_df_fixed_completed %>%
  group_by(desc_english) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  mutate(
    xpos = as.numeric(desc_english),
    ypos = max(ENT_surgeries_long_df_fixed_completed$year_plot, na.rm = TRUE) + 3
  )

# ========================================================================== #
# 12. OPTIONAL: REORDER BY PROCEDURE FREQUENCY
# ========================================================================== #

if (order_plot) {
  # Order procedures by descending count
  desc_levels <- count_df %>%
    arrange(desc(count)) %>%
    pull(desc_english)

  # Apply new ordering to main data
  ENT_surgeries_long_df_fixed_completed$desc_english <- factor(ENT_surgeries_long_df_fixed_completed$desc_english, levels = desc_levels)

  # Reorder stripe and count data frames
  stripe_df <- data.frame(
    desc_english = desc_levels,
    xmin = seq_along(desc_levels) - 0.5,
    xmax = seq_along(desc_levels) + 0.5
  )

  count_df <- count_df %>%
    mutate(xpos = match(desc_english, desc_levels))
}

# ========================================================================== #
# 13. CREATE MAIN VISUALIZATION
# ========================================================================== #

# Set seed for reproducible jittering
set.seed(42)

# Generate plot
p_ENT_surgery <- ggplot() +
# Add alternating background stripes
geom_rect(
    data = stripe_df,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    fill = rep(c("ivory2", "white"), length.out = nrow(stripe_df)),
    inherit.aes = FALSE
  ) +
# Add jittered points for surgical procedures
geom_jitter(
    data = ENT_surgeries_long_df_fixed_completed,
    aes(x = desc_english, y = year_plot),
    width = 0.4,
    height = 0,
    alpha = 0.6,
    color = "black"
  ) +
# Add count labels above plot
geom_text(
    data = count_df,
    aes(x = xpos, y = ypos, label = count),
    vjust = 0,
    size = 3,
    angle = 90
  ) +
# Preserve all procedure categories on x-axis
scale_x_discrete(drop = FALSE) +
# Y-axis with "No year" label for missing years
scale_y_continuous(
    breaks = c(min_year - 2, seq(min_year, max(ENT_surgeries_long_df_fixed_completed$year_plot, na.rm = TRUE))),
    labels = c("No year", as.character(seq(min_year, max(ENT_surgeries_long_df_fixed_completed$year_plot, na.rm = TRUE))))
  ) +
# Axis labels and title
labs(
    x = "Disease",
    y = "Year",
    title = "ENT surgeries over years"
  ) +
# Minimal theme
theme_minimal(base_size = 8) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(
      color = "black",
      size = 0.3,
      linetype = "dashed"
    )
  ) +
# Add horizontal grid lines at each year
geom_hline(
    yintercept = seq(min_year, max(ENT_surgeries_long_df_fixed_completed$year_plot, na.rm = TRUE)),
    color = "grey80",
    size = 0.3,
    linetype = "dashed"
  )

# Display plot
print(p_ENT_surgery)

# Save plot as SVG
ggsave(
  paste0("p_ENT_surgery", ".svg"),
  p_ENT_surgery,
  width = 12,
  height = 12
)

# Print summary of disease counts
head(as.data.frame(count_df[order(-count_df$count),]), n = 20)

# ========================================================================== #
# 14. ONE-HOT ENCODING FOR STATISTICAL MODELING
# ========================================================================== #

# Convert English descriptions to character for pivoting
ENT_surgeries_long_df_clean <- ENT_surgeries_long_df_4_onehot %>%
  mutate(desc_english = as.character(desc_english))

# Get unique ID-procedure combinations
ENT_surgeries_long_df_unique <- ENT_surgeries_long_df_clean %>%
  distinct(ID, desc_english)

# Add presence indicator
ENT_surgeries_long_df_unique <- ENT_surgeries_long_df_unique %>%
  mutate(present = 1)

# Pivot to wide format: one column per procedure
one_hot_partial <- ENT_surgeries_long_df_unique %>%
  pivot_wider(
    id_cols = ID,
    names_from = desc_english,
    values_from = present,
    values_fill = list(present = 0),
    values_fn = list(present = max) # Handle duplicates by taking max
  )

# Ensure all participants are represented (fill missing with zeros)
all_ids <- tibble(ID = trigeminale_daten_corrected_translated$ID)

# Join and fill missing values with 0 (no procedure)
one_hot_encoded <- all_ids %>%
  left_join(one_hot_partial, by = "ID") %>%
  replace(is.na(.), 0)

# Remove "Unknown" column if present
if ("Unknown" %in% names(one_hot_encoded)) {
  one_hot_encoded <- one_hot_encoded %>%
    dplyr::select(-Unknown)
}

print(sum(rowSums(one_hot_encoded[, -1]) > 0))

# Add column for subjects with any ENT surgery
one_hot_encoded <- one_hot_encoded %>%
  mutate("has ENT surgery" = if_else(rowSums(dplyr::select(., - ID)) == 0, 0, 1))

# Export one-hot encoded data
write.csv(one_hot_encoded, "onehot_ENT_surgery.csv", row.names = FALSE)

# ========================================================================== #
# END OF SCRIPT
# ========================================================================== #
