################################################################################
# Chronic and Neurological Disease Analysis - Trigeminal Sensitivity Study
# Author: Jorn Lotsch
# Description: Comprehensive analysis of chronic diseases and neurological
#              disorders over time
#              - Parsing disease descriptions with years from free text
#              - German to English translation
#              - Data correction for separately queried conditions
#              - Temporal visualization with color-coded neurological disorders
#              - One-hot encoding for statistical modeling
################################################################################

# ========================================================================== #
# 1. LOAD REQUIRED LIBRARIES
# ========================================================================== #

library(stringr)   # String manipulation and regex operations
library(dplyr)     # Data manipulation and transformation
library(tidyr)     # Data tidying and reshaping (pivot operations)
library(ggplot2)   # Data visualization
library(ggtext)    # Enhanced text rendering for plots
library(rlang)     # Tidy evaluation tools (for programming with dplyr)

# ========================================================================== #
# 2. CONFIGURATION SWITCHES
# ========================================================================== #

# Control whether to order plot by disease frequency (descending)
# Set to TRUE to sort diseases by count, FALSE for alphabetical
order_plot <- TRUE

# ========================================================================== #
# 3. HELPER FUNCTIONS
# ========================================================================== #

#' Parse year-description pairs from chronic disease text strings
#'
#' @param text_vec Character vector of disease descriptions (German text)
#' @return Data frame with alternating year/desc columns (year1, desc1, year2, desc2, ...)
#' @details Handles multiple formats:
#'   - Comma-separated entries: "2015 Asthma, 2020 Diabetes"
#'   - Multiple diseases with '+': "2015 Asthma + Diabetes"
#'   - "seit" (since) keyword: "seit 2015 Hypertonie"
#'   - Missing years (returns NA for year)
split_year_description <- function(text_vec) {

  # Internal function to parse a single text entry
  parse_one <- function(text) {
    # Split by comma, trim whitespace, remove "seit" (since) keyword
    items <- str_split(text, ",\\s*")[[1]]
    items <- str_trim(str_replace_all(items, "\\bseit\\b\\s*", ""))

    # Initialize storage vectors for years and descriptions
    years <- vector("character", length(items))
    descs <- vector("list", length(items))

    # Process each comma-separated item
    for (i in seq_along(items)) {
      # Extract ANY 4-digit year from the string (not just at start)
      year_match <- str_extract(items[i], "\\b\\d{4}\\b")
      years[i] <- year_match

      # Remove extracted year from description text
      desc_raw <- str_trim(str_replace(items[i], "\\b\\d{4}\\b", ""))

      # Split multiple diagnoses connected by '+' symbol
      desc_parts <- str_trim(str_split(desc_raw, "\\s*\\+\\s*")[[1]])
      descs[[i]] <- desc_parts
    }

    # Expand years to match number of descriptions
    # (repeat year for each disease when connected by '+')
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
# 5. PARSE DISEASE DESCRIPTIONS INTO STRUCTURED FORMAT
# ========================================================================== #

# Parse chronic disease descriptions into year/description pairs
result_chronic <- split_year_description(
  trigeminale_daten_corrected_translated$`Chronic disease`
)

# Parse neurological disease descriptions into year/description pairs
result_neurological <- split_year_description(
  trigeminale_daten_corrected_translated$`Neurological disorder`
)

# Rename neurological columns to continue numbering after chronic diseases
# This ensures unique column names when combining both datasets
names(result_neurological) <- paste0(
  rep(c("year", "desc"), dim(result_neurological)[2] / 2),
  rep(
    (dim(result_chronic)[2] / 2 + 1):(dim(result_chronic)[2] / 2 + dim(result_neurological)[2] / 2),
    each = 2
  )
)

# Combine chronic and neurological diseases with ID column
result_chronic_and_neurological <- cbind.data.frame(
  ID = trigeminale_daten_corrected_translated$ID,
  result_chronic,
  result_neurological
)

# Display first few rows for verification
head(result_chronic_and_neurological)

# ========================================================================== #
# 6. RESHAPE DATA TO LONG FORMAT
# ========================================================================== #

# Convert wide format to long format for analysis and plotting
# ID column is preserved as separate column (not pivoted)
long_df <- result_chronic_and_neurological %>%
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
# 7. GERMAN TO ENGLISH TRANSLATION MAP
# ========================================================================== #

# Comprehensive translation dictionary (including common typos)
translation_map <- c(
  # Gastrointestinal
  "Achalasie" = "Achalasia",

  # Cardiovascular
  "koromare Herzkrankheit" = "Coronary heart disease",
  "Hypercholesterinämie" = "Hypercholesterolemia",
  "Hypertonie" = "Hypertension",
  "essentielle Hypertonie" = "Essential hypertension",
  "Hyprtonie" = "Hypertension",  # Typo correction
  "arterielle Hypertonie" = "Arterial hypertension",
  "art. Hypertonie" = "Arterial hypertension",
  "Hypertonus" = "Hypertension (variant term)",
  "Arteriosklerose" = "Arteriosclerosis",
  "Herzinsuffizienz" = "Heart failure",
  "Herzklappeninsuffizienz l" = "Heart valve insufficiency (left)",
  "Herzrhythmusstörung" = "Cardiac arrhythmia",

  # Thyroid disorders
  "Thyreoektomie > behandelte Hypothyreose" = "Hypothyroidism (treated after thyroidectomy)",
  "Hypothyreose" = "Hypothyroidism",
  "Geburt Hypothyreose" = "Congenital hypothyroidism",
  "Mb. Basedow" = "Graves' disease",
  "Hashimoto" = "Hashimoto's thyroiditis",
  "Hasimodo" = "Hashimoto's thyroiditis",  # Typo correction
  "Hashimodo" = "Hashimoto's thyroiditis",  # Typo correction
  "Schilddrüsenprobleme" = "Thyroid problems",

  # Mental health
  "PTBS" = "Post-traumatic stress disorder",
  "Depression" = "Depression",
  "Depressionen" = "Depression",
  "chronische Schmerzstörung mit somatischen und psychischen Faktoren" = "Chronic pain disorder",
  "Dissoziation" = "Dissociation",
  "Zwangsstörung" = "Obsessive-compulsive disorder",

  # Respiratory - Asthma variants
  "allergisches Asthma" = "Allergic asthma",
  "Asthma" = "Asthma",
  "Asthma bronchiale" = "Bronchial asthma",
  "Asthma bronchiale Geburt" = "Congenital asthma",
  "Asthma Bronchiale Kindheit" = "Childhood bronchial asthma",
  "Asthma Kindesalter" = "Childhood asthma",
  "Asthma Geburt" = "Congenital asthma",
  "Asthma pollinose" = "Asthma with hay fever",
  "chron. Bronchitis" = "Chronic bronchitis",

  # Diabetes (including typos)
  "D.m." = "Diabetes mellitus (generic)",
  "D.m. Typ l" = "Diabetes mellitus type I",
  "D.m. Typ ll" = "Diabetes mellitus type II",
  "D.m. Typ I" = "Diabetes mellitus type I",
  "D.m. Typ II" = "Diabetes mellitus type II",
  "Diabethes" = "Diabetes mellitus (generic)",  # Typo
  "Diabethes mellitus" = "Diabetes mellitus (generic)",
  "Diabethes mellitus Typ ll" = "Diabetes mellitus type II",
  "Diabethes nellitus Typ 2" = "Diabetes mellitus type II",  # Typo
  "Diabethes Typ l" = "Diabetes mellitus type I",
  "Diabethes in Diagnostik" = "Diabetes under investigation",
  "Insulinresistenz (Diabethes)" = "Insulin resistance (diabetes type II)",

  # Inflammatory bowel disease
  "Mb. Chron" = "Crohn's disease",
  "Morbus Crohn" = "Crohn's disease",

  # Rheumatic conditions
  "Rheuma" = "Rheumatism",
  "rheumatische Atritis" = "Rheumatoid arthritis",
  "Arthrose" = "Osteoarthritis",
  "Weichteilrheuma" = "Soft tissue rheumatism",
  "Psoriasis-Arthritis" = "Psoriatic arthritis",
  "Fibromyalgie" = "Fibromyalgia",
  "Mb. Bechterew" = "Ankylosing spondylitis",
  "Mb. Scheuermann" = "Scheuermann's disease",

  # Skin conditions
  "Neurodermitis" = "Atopic dermatitis",
  "Neurodermitis Geburt" = "Congenital atopic dermatitis",
  "Schuppenflechte" = "Psoriasis",
  "3 LM atopisches Exzem" = "Atopic dermatitis",

  # Kidney disorders
  "z.N. Nierentransplantation" = "Kidney transplant (status post)",
  "Nieren-Op" = "Kidney surgery",
  "Niereninsuffizienz" = "Kidney insufficiency",
  "Nierenstauung" = "Kidney congestion",
  "Nierensteine" = "Kidney stones",
  "Nierensteine?" = "Kidney stones (suspected)",
  "Nierenreilresektion nach Tumor" = "Partial kidney resection after tumor",
  "schlechte Nierenwerte" = "Poor kidney function",

  # Liver
  "PBC" = "Primary biliary cirrhosis",

  # Neurological disorders
  "Migräne" = "Migraine",
  "Epilepsie" = "Epilepsy",
  "Tinitus" = "Tinnitus",
  "ADHS" = "ADHD",
  "MS" = "Multiple sclerosis",
  "Spinalkanalstenose" = "Spinal canal stenosis",
  "Polyneuropathie" = "Polyneuropathy",
  "Myastenia gravis" = "Myasthenia gravis",
  "cerebrale Ataxie" = "Cerebellar ataxia",
  "chron. Neuroboreliose" = "Chronic neuroborreliosis",
  "Mitochondriopathie" = "Mitochondriopathy",
  "einlaufende Demenz" = "Dementia (onset)",

  # Other conditions
  "Schlafapnoe" = "Sleep apnea",
  "Krebs" = "Cancer",
  "Hirntumor" = "Brain tumor",
  "Tumor kopf" = "Head tumor",
  "Asperger Autismus" = "Asperger's syndrome",
  "Raynaud Syndrom" = "Raynaud's syndrome",
  "Gürtelrose" = "Shingles",

  # Allergies
  "Pollinosis" = "Allergic rhinitis (hay fever)",
  "Heuschnupfen" = "Hay fever",
  "Heufieber" = "Hay fever",
  "Allergien" = "Allergic problems",
  "Allergie" = "Allergic problems",
  "Allergie Hausstaub" = "House dust allergy",
  "Hausstauballergie" = "House dust allergy",
  "Hausstaubmilbenallergie" = "House dust mite allergy",
  "Katze" = "Cat allergy",
  "Katzenhaarallergie" = "Cat hair allergy",
  "Lebensmittelallergien" = "Food allergies",
  "Histaminintoleranz" = "Histamine intolerance",

  # Miscellaneous
  "RR-Schwankungen" = "Blood pressure fluctuations",
  "Kreislaufprobleme" = "Circulatory problems",
  "Osteoporose" = "Osteoporosis",
  "Glaukom" = "Glaucoma",
  "Geruchsempfindungsminderung" = "Reduced sense of smell",
  "Faktor V Leiden" = "Factor V Leiden mutation",
  "MGUS" = "Monoclonal gammopathy of undetermined significance",
  "system. Lupus erythemadodes" = "Systemic lupus erythematosus",
  "Sarkoidose" = "Sarcoidosis",
  "Granulomatose mit Polyangiitis" = "Granulomatosis with polyangiitis",
  "PCOS" = "Polycystic ovary syndrome",
  "Thrombozytose" = "Thrombocytosis",
  "Zölliakie" = "Celiac disease",
  "chr. Entzündung der Schulter re." = "Chronic inflammation of right shoulder",
  "postvirales Fatiguesyndrom" = "Post-viral fatigue syndrome",
  "COPD" = "COPD",
  "Hodghin" = "M. Hodgkin",
  "j" = "Unspecified"
)

# ========================================================================== #
# 8. DEFINE DISEASE CATEGORIES FOR DATA CORRECTION
# ========================================================================== #

# Neurological disorders for separate handling and color-coding
neurological_disorders <- c(
  "ADHD",
  "Asperger's syndrome",
  "Brain tumor",
  "Cerebellar ataxia",
  "Chronic neuroborreliosis",
  "Dementia (onset)",
  "Dissociation",
  "Epilepsy",
  "Head tumor",
  "Migraine",
  "Multiple sclerosis",
  "Mitochondriopathy",
  "Myasthenia gravis",
  "Obsessive-compulsive disorder",
  "Polyneuropathy",
  "Post-traumatic stress disorder",
  "Spinal canal stenosis",
  "Tinnitus"
)

# Allergic diseases for validation against separately queried allergy data
allergic_diseases <- c(
  "Allergic asthma",
  "Allergic rhinitis (hay fever)",
  "Atopic dermatitis",
  "Bronchial asthma",
  "Cat allergy",
  "Food allergies",
  "Hay fever",
  "House dust allergy",
  "House dust mite allergy",
  "Jugendalter allergisches Asthma",
  "Childhood asthma",
  "Childhood bronchial asthma",
  "Congenital asthma",
  "Congenital atopic dermatitis",
  "Asthma with hay fever",
  "Asthma"
)

# ========================================================================== #
# 9. TRANSLATE DESCRIPTIONS TO ENGLISH
# ========================================================================== #

# Apply translation map to disease descriptions
long_df <- long_df %>%
  mutate(desc_english = translation_map[desc]) %>%
  mutate(desc_english = ifelse(is.na(desc_english), "Unknown", desc_english)) %>%
  mutate(desc_english = factor(desc_english, levels = sort(unique(desc_english))))

dim(long_df)

# ========================================================================== #
# 10. CORRECT FOR SEPARATELY QUERIED CONDITIONS
# ========================================================================== #
# Some conditions were asked about specifically in separate questions
# but may not have been mentioned in free text. Add placeholder entries.

#' Add flag for missing disease terms
#'
#' @param df Main data frame with disease information
#' @param check_column Column name to check for presence of terms
#' @param flag_text Text to add as flag when terms are missing
#' @param ids_to_check Restrict checking to specific IDs
#' @return Updated data frame with new rows for missing terms
add_flag_for_ids <- function(df, check_column, flag_text, ids_to_check) {
  check_col_sym <- sym(check_column)

  # IDs that already have the flag_text in check_column
  ids_with_flag <- df %>%
    filter(.data[[check_column]] == flag_text) %>%
    distinct(ID) %>%
    pull(ID)

  # IDs missing the flag_text row
  ids_missing_flag <- setdiff(ids_to_check, ids_with_flag)

  # Create new rows with NA except for ID and flag_text
  new_rows <- df %>%
    filter(ID %in% ids_missing_flag) %>%
    distinct(ID) %>%
    mutate(across(everything(), ~NA)) %>%
    mutate(ID = ids_missing_flag) %>%
    mutate(!!check_col_sym := flag_text)

  # Append new rows to original df
  df_updated <- bind_rows(df, new_rows)

  return(df_updated)
}

# --- Correct for separately queried allergies ---
# Get IDs with allergies from main questionnaire
has_allergy <- trigeminale_daten_corrected_translated %>%
  filter(!is.na(`Allergic problems`) & `Allergic problems` != "n") %>%
  pull(ID)

# Get IDs who mentioned specific allergies in free text
some_allergy_mentioned <- long_df %>%
  filter(desc_english %in% allergic_diseases) %>%
  distinct(ID) %>%
  pull(ID)

# Find IDs with allergies but no specific mention

has_allergy_all <- union(some_allergy_mentioned,has_allergy)

long_df <- add_flag_for_ids(
  df = long_df,
  check_column = "desc_english",
  flag_text = "Allergic problems",
  ids_to_check = has_allergy_all
)

# --- Correct for separately queried neurological disorders ---
has_neurological <- trigeminale_daten_corrected_translated %>%
  filter(!is.na(`Neurological disorder`) & `Neurological disorder` != "n") %>%
  pull(ID)

some_neurological_mentioned <- long_df %>%
  filter(desc_english %in% neurological_disorders) %>%
  distinct(ID) %>%
  pull(ID)

has_neurological_all <- union(has_neurological, some_neurological_mentioned)

long_df <- add_flag_for_ids(
  df = long_df,
  check_column = "desc_english",
  flag_text = "Neurological disorder",
  ids_to_check = has_neurological_all
)

# --- Correct for separately queried chronic sinusitis ---
has_chronic_sinusitis <- trigeminale_daten_corrected_translated %>%
  filter(!is.na(`Chronic sinusitis`) & `Chronic sinusitis` != "n") %>%
  pull(ID)

long_df <- add_flag_for_ids(
  df = long_df,
  check_column = "desc_english",
  flag_text = "Chronic sinusitis",
  ids_to_check = has_chronic_sinusitis
)

# --- Correct for separately queried migraines ---
has_migraine <- which(trigeminale_daten_corrected_translated$`How often do you have migraine per month` > 0)

long_df <- add_flag_for_ids(
  df = long_df,
  check_column = "desc_english",
  flag_text = "Migraine",
  ids_to_check = has_migraine
)


# Store corrected data frame for later one-hot encoding
long_df_4_onehot <- long_df

# ========================================================================== #
# 11. PREPARE DATA FOR VISUALIZATION
# ========================================================================== #

# Determine minimum year for handling missing years
min_year <- min(long_df$year, na.rm = TRUE)

# Create plotting year: map NA years to pseudo-year below minimum
long_df <- long_df %>%
  mutate(
    year_plot = as.numeric(year),
    year_plot = ifelse(is.na(year_plot), min_year - 2, year_plot)
  )

# Display diagnoses with years but unknown translation (for debugging)
long_df[!is.na(long_df$year) & long_df$desc_english == "Unknown", c("ID", "desc")]

# Remove unknown/untranslated diagnoses from visualization
long_df <- long_df %>%
  filter(desc_english != "Unknown")

# Reset factor levels after filtering
desc_levels <- sort(unique(long_df$desc_english))
long_df$desc_english <- factor(long_df$desc_english, levels = desc_levels)

# Create data frame for alternating background stripes in plot
stripe_df <- data.frame(
  desc_english = desc_levels,
  xmin = seq_along(desc_levels) - 0.5,
  xmax = seq_along(desc_levels) + 0.5
)

# Calculate disease counts for annotation
count_df <- long_df %>%
  group_by(desc_english) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  mutate(
    xpos = as.numeric(desc_english),
    ypos = max(long_df$year_plot, na.rm = TRUE) + 3
  )

# ========================================================================== #
# 12. OPTIONAL: REORDER BY DISEASE FREQUENCY
# ========================================================================== #

if (order_plot) {
  # Order diseases by descending count
  desc_levels <- count_df %>%
    arrange(desc(count)) %>%
    pull(desc_english)

  # Apply new ordering to main data
  long_df$desc_english <- factor(long_df$desc_english, levels = desc_levels)

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
# 13. FINALIZE PLOT ELEMENTS
# ========================================================================== #

# Color-code axis labels: red for neurological, black for others
label_colors <- ifelse(
  desc_levels %in% neurological_disorders,
  "red",
  "black"
)

# Add neurological indicator to plotting data
long_df <- long_df %>%
  mutate(is_neurological = desc_english %in% neurological_disorders)

# Set seed for reproducible jittering
set.seed(42)

# ========================================================================== #
# 14. CREATE MAIN VISUALIZATION
# ========================================================================== #

p_chronic_disseases <- ggplot() +
  # Add alternating background stripes
  geom_rect(
    data = stripe_df,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    fill = rep(c("ivory2", "white"), length.out = nrow(stripe_df)),
    inherit.aes = FALSE
  ) +
  # Add jittered points for diagnoses (color-coded by neurological status)
  geom_jitter(
    data = long_df,
    aes(x = desc_english, y = year_plot, color = is_neurological),
    width = 0.4,
    height = 0,
    alpha = 0.6,
    show.legend = FALSE
  ) +
  # Add count labels above plot
  geom_text(
    data = count_df,
    aes(x = xpos, y = ypos, label = count),
    vjust = 0,
    size = 3,
    angle = 90
  ) +
  # Preserve all disease categories on x-axis
  scale_x_discrete(drop = FALSE) +
  # Color neurological disorders red, others black
  scale_color_manual(
    values = c("TRUE" = "red", "FALSE" = "black"),
    labels = c("TRUE" = "Neurological", "FALSE" = "Other"),
    name = "Disease Type"
  ) +
  # Y-axis with "No year" label for missing years
  scale_y_continuous(
    breaks = c(min_year - 2, seq(min_year, max(long_df$year_plot, na.rm = TRUE))),
    labels = c("No year", as.character(seq(min_year, max(long_df$year_plot, na.rm = TRUE))))
  ) +
  # Axis labels and title
  labs(
    x = "Disease",
    y = "Year",
    title = "Diagnoses of chronic diseases over years"
  ) +
  # Minimal theme
  theme_minimal(base_size = 8) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      color = label_colors
    ),
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
    yintercept = seq(min_year, max(long_df$year_plot, na.rm = TRUE)),
    color = "grey80",
    size = 0.3,
    linetype = "dashed"
  )

# Display plot
print(p_chronic_disseases)

# Save plot as SVG
ggsave(
  paste0("p_chronic_disseases", ".svg"),
  p_chronic_disseases,
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

# Get unique ID-disease combinations
long_df_unique <- long_df_clean %>%
  distinct(ID, desc_english)

# Add presence indicator
long_df_unique <- long_df_unique %>%
  mutate(present = 1)

# Pivot to wide format: one column per disease
one_hot_encoded <- long_df_unique %>%
  pivot_wider(
    id_cols = ID,
    names_from = desc_english,
    values_from = present,
    values_fill = list(present = 0),
    values_fn = list(present = max)  # Handle duplicates by taking max
  )
# Remove "Unknown" column if present
one_hot_encoded <- one_hot_encoded %>%
  select(-Unknown)

dim(one_hot_encoded)
print(sum(rowSums(one_hot_encoded[,-1]) > 0))

# Add column for subjects with no diseases
one_hot_encoded <- one_hot_encoded %>%
  mutate("has chronic disease" = if_else(rowSums(select(., -ID)) == 0, 0, 1))


# Export one-hot encoded data
write.csv(one_hot_encoded, "onehot_chronic_diseases.csv", row.names = FALSE)

# ========================================================================== #
# END OF SCRIPT
# ========================================================================== #
