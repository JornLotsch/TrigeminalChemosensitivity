################################################################################
# ENT Surgery History Analysis - Trigeminal Sensitivity Study
# Author: [Author Name]
# Date: [Date]
# Description: Comprehensive analysis of ear, nose, and throat surgical
#              procedures over time
#              - Parsing surgical procedure descriptions with years from free text
#              - German to English translation of medical procedures
#              - Temporal visualization of surgical interventions
#              - One-hot encoding for statistical modeling
################################################################################

# ========================================================================== #
# 1. LOAD REQUIRED LIBRARIES
# ========================================================================== #

library(stringr)   # String manipulation and regex operations
library(dplyr)     # Data manipulation and transformation
library(tidyr)     # Data tidying and reshaping (pivot operations)
library(ggplot2)   # Data visualization

# ========================================================================== #
# 2. CONFIGURATION SWITCHES
# ========================================================================== #

# Control whether to order plot by procedure frequency (descending)
# Set to TRUE to sort procedures by count, FALSE for alphabetical
order_plot <- TRUE

# ========================================================================== #
# 3. HELPER FUNCTIONS
# ========================================================================== #

#' Parse year-description pairs from ENT surgical procedure text strings
#'
#' @param text_vec Character vector of procedure descriptions (German text)
#' @return Data frame with alternating year/desc columns (year1, desc1, year2, desc2, ...)
#' @details Handles multiple year formats:
#'   - Standard 4-digit year: "2015 Tonsillotomie"
#'   - Slash notation: "/74", "/2021"
#'   - MM/YY format: "05/21" (extracts year only)
#'   - Decade notation: "90er Jahre" (maps to 1995)
#'   - Multiple procedures with '+': "2015 Adenotomie + Tonsillotomie"
parse_procedure_descriptions <- function(text_vec) {

  # Internal function to parse a single text entry
  parse_one <- function(text) {
    # Split by comma, trim whitespace, remove empty and "n" (no) entries
    items <- stringr::str_split(text, ",\\s*")[[1]]
    items <- items[items != ""]
    items <- stringr::str_trim(items)
    items <- items[items != "n"]

        #' Extract year from a procedure string
        #'
        #' @param s Single procedure string
        #' @return Year as character string or NA
    extract_year <- function(s) {
      s <- stringr::str_trim(s)
      if (is.na(s) || s == "") return(NA_character_)

      # Match standard 4-digit year at start (e.g., "2015 Adenotomie")
      y <- stringr::str_extract(s, "^\\d{4}")
      if (!is.na(y)) return(y)

      # Match year after slash: /74 or /2021
      y <- stringr::str_extract(s, "/(\\d{2,4})")
      if (!is.na(y)) {
        y <- sub("/", "", y)
        # Convert 2-digit year to 4-digit (< 50 = 20XX, >= 50 = 19XX)
        if (nchar(y) == 2) {
          y_num <- as.numeric(y)
          y <- ifelse(y_num < 50, paste0("20", sprintf("%02d", y_num)), paste0("19", y))
        }
        return(y)
      }

      # Match MM/YY format (e.g., "05/21"), extract only year part
      y <- stringr::str_extract(s, "\\d{2}/(\\d{2})")
      if (!is.na(y)) {
        yy <- sub("\\d{2}/", "", y)
        y_num <- as.numeric(yy)
        y <- ifelse(y_num < 50, paste0("20", sprintf("%02d", y_num)), paste0("19", yy))
        return(y)
      }

      # Match decade notation "90er Jahre" (map to mid-decade)
      if (stringr::str_detect(s, "90er Jahre")) return("1995")

      NA_character_
    }

    # Extract years for all items
    years <- sapply(items, extract_year)

    # Process descriptions: remove year patterns from text
    descs <- vector("list", length(items))
    for (i in seq_along(items)) {
      desc_raw <- items[i]
      # Remove all matched year patterns from string
      desc_raw <- stringr::str_remove(desc_raw, "^\\d{4}\\s*")     # Leading YYYY
      desc_raw <- stringr::str_remove(desc_raw, "/\\d{2,4}")       # /YY or /YYYY
      desc_raw <- stringr::str_remove(desc_raw, "\\d{2}/\\d{2}")   # MM/YY
      desc_raw <- stringr::str_remove(desc_raw, "\\s*\\d{4}$")     # Trailing YYYY
      desc_raw <- stringr::str_remove(desc_raw, "90er Jahre")      # Decade text
      desc_raw <- stringr::str_trim(desc_raw)

      # Split multiple procedures connected by '+' symbol
      desc_parts <- stringr::str_trim(stringr::str_split(desc_raw, "\\s*\\+\\s*")[[1]])
      descs[[i]] <- desc_parts
    }

    # Expand years to match number of descriptions
    # (repeat year for each procedure when connected by '+')
    years_expanded <- unlist(mapply(
      function(y, d) rep(y, length(d)),
      years, descs,
      SIMPLIFY = FALSE
    ))
    descs_expanded <- unlist(descs)

    # Interleave year and description: year1, desc1, year2, desc2, ...
    as.vector(rbind(years_expanded, descs_expanded))
  }

  # Apply parsing function to all text entries
  parsed <- lapply(text_vec, parse_one)
  max_len <- max(sapply(parsed, length))

  # Pad all entries to same length with NA for rectangular structure
  parsed_padded <- lapply(parsed, function(x) {
    length(x) <- max_len
    x
  })

  # Convert to data frame
  df <- as.data.frame(do.call(rbind, parsed_padded), stringsAsFactors = FALSE)

  # Create column names: year1, desc1, year2, desc2, ...
  half <- max_len / 2
  colnames(df) <- paste0(
    rep(c("year", "desc"), half),
    rep(seq_len(half), each = 2)
  )

  return(df)
}

# ========================================================================== #
# 4. DATA IMPORT
# ========================================================================== #

# Read the main corrected and translated dataset
trigeminale_daten_corrected_translated <- read.csv(
  "trigeminale_daten_corrected_translated.csv",
  check.names = FALSE
)

# ========================================================================== #
# 5. CLEAN AND PREPARE DATA
# ========================================================================== #

# Clean ENT surgery column: replace "n" (no), empty, NA, NaN with NA
clean_text_vec <- ifelse(
  trigeminale_daten_corrected_translated$`Surgery in ENT region` %in% c("n", "", "NA", "NaN"),
  NA_character_,
  trigeminale_daten_corrected_translated$`Surgery in ENT region`
)

# ========================================================================== #
# 6. PARSE SURGICAL PROCEDURE DESCRIPTIONS
# ========================================================================== #

# Parse procedure descriptions into year/description pairs
result_procedures <- parse_procedure_descriptions(clean_text_vec)

# ========================================================================== #
# 7. RESHAPE DATA TO LONG FORMAT
# ========================================================================== #

# Convert wide format to long format for analysis
# ID column is preserved as separate column (not pivoted)
long_df <- result_procedures %>%
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
# 8. FIX YEAR EXTRACTION ISSUES
# ========================================================================== #
# Some years may remain in description text due to parsing edge cases
# This section catches and corrects these remaining year patterns

long_df_fixed <- long_df %>%
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
# 9. GERMAN TO ENGLISH TRANSLATION MAP
# ========================================================================== #

# Comprehensive translation dictionary for ENT surgical procedures
# Grouped by procedure type and normalized to consolidate variants
translation_map <- c(
  # Adenoid procedures - consolidated variants
  "Adenotomie" = "Adenotomy",
  "Adenotomie Kindesalter" = "Adenotomy",
  "Kindesalter Adenotomie" = "Adenotomy",
  "Adentomie" = "Adenotomy",  # Typo correction

  # Nasal trauma/fractures
  "Nasenfraktur" = "Nasal fracture",
  "Nasen # OP" = "Nasal fracture surgery",

  # Tonsil procedures - distinguished by type (full vs partial removal)
  "Mandeln" = "Tonsillectomy",            # Full removal
  "Mandeln OP" = "Tonsillectomy",
  "Mandel OP" = "Tonsillectomy",
  "Mandel Op" = "Tonsillectomy",
  "Mandel-OP" = "Tonsillectomy",
  "Kindesalter Mandeln" = "Tonsillectomy",
  "Tonsillektomie" = "Tonsillectomy",     # Explicit full removal

  "Tonsillotomie" = "Tonsillotomy",       # Partial removal
  "Tonsillotomie Kindheit" = "Tonsillotomy",
  "Kindesalter Tonsillotomie" = "Tonsillotomy",
  "Tonsillotomie Kindesalter" = "Tonsillotomy",
  "Tonsillotomie jugend" = "Tonsillotomy",
  "er Tonsillotomie" = "Tonsillotomy",

  # Nasal polyps - grouped and normalized
  "Polypen" = "Nasal polyps surgery",
  "Kleinkind Polypen" = "Nasal polyps surgery",
  "Kindheit Polypen" = "Nasal polyps surgery",
  "Kindesalter Polypen" = "Nasal polyps surgery",
  "CRS Polypenentfernung 2x" = "Nasal polyps surgery",

  # Nasal growth/sinus surgeries
  "Wucherungen Nase" = "Nasal growth surgery",
  "NNH" = "Paranasal sinus surgery",
  "NNH Op" = "Paranasal sinus surgery",
  "NNH links" = "Left paranasal sinus surgery",

  # Nasal septum surgeries and corrections - grouped and normalized
  "NSW" = "Nasal septum surgery",
  "NSW OP" = "Nasal septum surgery",
  "Nasenscheidewand" = "Nasal septum surgery",
  "Nasenscheidewand OP" = "Nasal septum surgery",
  "Begradigung Nasenscheidewand" = "Nasal septum correction",
  "Septum OP" = "Nasal septum surgery",
  "Septumkorrektur" = "Nasal septum correction",
  "Septumdeviation" = "Septal deviation correction",

  # Fenestration and Caldwell-Luc procedures
  "Fensterung" = "Fenestration",
  "Fensterung Kieferhöhle" = "Caldwell-Luc procedure",
  "Fensterung Kieferhöhle li." = "Left Caldwell-Luc procedure",

  # Ear surgeries
  "Knorpelplastik Ohr" = "Cartilage plasty (ear)",
  "Ohr li. Implantat" = "Left ear implant",
  "Ohr anlegen" = "Otoplasty",
  "Otoplastik" = "Otoplasty",
  "Austausch Hämmerchen Ohr" = "Malleus replacement (middle ear)",

  # Tympanic membrane surgeries
  "Trommelfell" = "Tympanic membrane surgery",
  "Trommelfell 2x" = "Tympanic membrane surgery",
  "Trommelfell-OP" = "Tympanic membrane surgery",
  "Trommelfellimplantat" = "Tympanic membrane implant",
  "Trommelfellverschluss" = "Tympanic membrane closure",
  "Trommelfellschnitt" = "Myringotomy",

  "Knalltrauma" = "Acoustic trauma",

  # Vocal cord surgeries
  "Stimmband OP" = "Vocal cord surgery",
  "Stimmknötchen" = "Vocal cord nodules surgery",

  # Oncological procedures
  "Mundbodenkarzinom" = "Mouth floor carcinoma surgery",
  "Parotisentfernung" = "Parotidectomy",
  "Karzinom Hals" = "Neck carcinoma surgery",
  "Medulläres Schildrüsenkarzinom" = "Medullary thyroid carcinoma surgery",
  "Karzinom Rachen li." = "Left pharyngeal carcinoma surgery",
  "OP Halsbereich" = "Neck surgery",
  "Abszess Rachen" = "Pharyngeal abscess surgery",
  "Trauma Hals" = "Neck trauma surgery",

  # Miscellaneous surgical procedures - grouped consistently
  "Begradigung Nase" = "Nasal straightening surgery",
  "Schilddrüse" = "Thyroid gland surgery",
  "Thyreodektomie" = "Thyroidectomy",
  "Kiefernzyste" = "Jaw cyst removal",
  "OP Polypen Stirnhöhle" = "Frontal sinus polyp surgery",

  # Turbinate reduction procedures - harmonized treatments
  "Nasenmuscheln Verkleinerung" = "Inferior turbinate reduction",
  "Reduktion Nasenmuscheln" = "Inferior turbinate reduction",
  "Verkleinerung Nasenmuscheln" = "Inferior turbinate reduction",
  "Laserung Nasenmuscheln" = "Inferior turbinate laser therapy",
  "Verödung Nasenmuschel" = "Inferior turbinate cauterization",
  "Verödung Nasenmuscheln" = "Inferior turbinate cauterization",
  "Stromwellentherapie zum Abschwellen der Nasenmuscheln" = "Inferior turbinate radiofrequency therapy",

  # Placeholder/invalid entries
  "j" = "Unspecified"
)

# ========================================================================== #
# 10. TRANSLATE DESCRIPTIONS TO ENGLISH
# ========================================================================== #

# Apply translation map to procedure descriptions
long_df_fixed <- long_df_fixed %>%
  mutate(desc_english = translation_map[desc]) %>%
  mutate(desc_english = ifelse(is.na(desc_english), "Unknown", desc_english)) %>%
  mutate(desc_english = factor(desc_english, levels = sort(unique(desc_english))))

# ========================================================================== #
# 11. STORE DATA FOR ONE-HOT ENCODING
# ========================================================================== #

# Store data frame before filtering for later one-hot encoding
long_df_4_onehot <- long_df_fixed

# ========================================================================== #
# 12. PREPARE DATA FOR VISUALIZATION
# ========================================================================== #

# Determine minimum year for handling missing years
min_year <- min(long_df_fixed$year, na.rm = TRUE)

# Create plotting year: map NA years to pseudo-year below minimum
long_df_fixed <- long_df_fixed %>%
  mutate(
    year_plot = as.numeric(year),
    year_plot = ifelse(is.na(year_plot), min_year - 2, year_plot)
  )

# Remove unknown/untranslated procedures from visualization
long_df_fixed <- long_df_fixed %>%
  filter(desc_english != "Unknown")

# Reset factor levels after filtering
desc_levels <- sort(unique(long_df_fixed$desc_english))
long_df_fixed$desc_english <- factor(long_df_fixed$desc_english, levels = desc_levels)

# Create data frame for alternating background stripes in plot
stripe_df <- data.frame(
  desc_english = desc_levels,
  xmin = seq_along(desc_levels) - 0.5,
  xmax = seq_along(desc_levels) + 0.5
)

# Calculate procedure counts for annotation
count_df <- long_df_fixed %>%
  group_by(desc_english) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  mutate(
    xpos = as.numeric(desc_english),
    ypos = max(long_df_fixed$year_plot, na.rm = TRUE) + 3
  )

# ========================================================================== #
# 13. OPTIONAL: REORDER BY PROCEDURE FREQUENCY
# ========================================================================== #

if (order_plot) {
  # Order procedures by descending count
  desc_levels <- count_df %>%
    arrange(desc(count)) %>%
    pull(desc_english)

  # Apply new ordering to main data
  long_df_fixed$desc_english <- factor(long_df_fixed$desc_english, levels = desc_levels)

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
# 14. CREATE MAIN VISUALIZATION
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
    data = long_df_fixed,
    aes(x = desc_english, y = year_plot),
    width = 0.3,
    height = 0.3,
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
    breaks = c(min_year - 2, seq(min_year, max(long_df_fixed$year_plot, na.rm = TRUE))),
    labels = c("No year", as.character(seq(min_year, max(long_df_fixed$year_plot, na.rm = TRUE))))
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
    yintercept = seq(min_year, max(long_df_fixed$year_plot, na.rm = TRUE)),
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
# 15. ONE-HOT ENCODING FOR STATISTICAL MODELING
# ========================================================================== #

# Convert English descriptions to character for pivoting
long_df_clean <- long_df_4_onehot %>%
  mutate(desc_english = as.character(desc_english))

# Get unique ID-procedure combinations
long_df_unique <- long_df_clean %>%
  distinct(ID, desc_english)

# Add presence indicator
long_df_unique <- long_df_unique %>%
  mutate(present = 1)

# Pivot to wide format: one column per procedure
one_hot_partial <- long_df_unique %>%
  pivot_wider(
    id_cols = ID,
    names_from = desc_english,
    values_from = present,
    values_fill = list(present = 0),
    values_fn = list(present = max)  # Handle duplicates by taking max
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
    select(-Unknown)
}

print(sum(rowSums(one_hot_encoded[,-1]) > 0))

# Add column for subjects with no ENT surgeries
one_hot_encoded <- one_hot_encoded %>%
  mutate("no ENT surgery" = if_else(rowSums(select(., -ID)) == 0, 1, 0))

# Export one-hot encoded data
write.csv(one_hot_encoded, "onehot_ENT_surgery.csv", row.names = FALSE)

# ========================================================================== #
# END OF SCRIPT
# ========================================================================== #
