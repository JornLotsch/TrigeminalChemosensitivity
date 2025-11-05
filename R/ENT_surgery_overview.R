# ========================================================================== #
# Load required packages -----------------------------------------------------
# ========================================================================== #

library(readxl)    # Excel file reading
library(stringr)   # String manipulation utilities
library(dplyr)     # Data manipulation verbs
library(tidyr)     # Data tidying and reshaping


# ========================================================================== #
# Global switches -------------------------------------------------------------
# ========================================================================== #

order_plot <- FALSE  # Control ordering of plot elements (TRUE = order by frequency)


# ========================================================================== #
# Function: Parse procedure descriptions --------------------------------------
# ---------------------------------------------------------------------------
# Parses strings containing procedure descriptions with years,
# splitting into year-description pairs and handling various date formats
# ========================================================================== #

parse_procedure_descriptions <- function(text_vec) {

  parse_one <- function(text) {
    # Split by commas, trim spaces, and remove filler values
    items <- stringr::str_split(text, ",\\s*")[[1]]
    items <- items[items != ""]
    items <- stringr::str_trim(items)
    items <- items[items != "n"]  # Skip filler "n"

    # Helper function to extract year from string patterns
    extract_year <- function(s) {
      s <- stringr::str_trim(s)
      if (is.na(s) || s == "") return(NA_character_)

      # Match 4-digit year at start of string
      y <- stringr::str_extract(s, "^\\d{4}")
      if (!is.na(y)) return(y)

      # Match /YY or /YYYY patterns (convert 2-digit to 4-digit)
      y <- stringr::str_extract(s, "/(\\d{2,4})")
      if (!is.na(y)) {
        y <- sub("/", "", y)
        if (nchar(y) == 2) {
          y_num <- as.numeric(y)
          y <- ifelse(y_num < 50, paste0("20", sprintf("%02d", y_num)), paste0("19", y))
        }
        return(y)
      }

      # Match MM/YY format - take only year part
      y <- stringr::str_extract(s, "\\d{2}/(\\d{2})")
      if (!is.na(y)) {
        yy <- sub("\\d{2}/", "", y)
        y_num <- as.numeric(yy)
        y <- ifelse(y_num < 50, paste0("20", sprintf("%02d", y_num)), paste0("19", yy))
        return(y)
      }

      # Map "90er Jahre" to approximate mid-1990s
      if (stringr::str_detect(s, "90er Jahre")) return("1995")

      NA_character_
    }

    years <- sapply(items, extract_year)

    # Extract descriptions removing year patterns, split by '+'
    descs <- vector("list", length(items))
    for (i in seq_along(items)) {
      desc_raw <- items[i]

      desc_raw <- stringr::str_remove(desc_raw, "^\\d{4}\\s*")    # Leading YYYY
      desc_raw <- stringr::str_remove(desc_raw, "/\\d{2,4}")       # /YY or /YYYY
      desc_raw <- stringr::str_remove(desc_raw, "\\d{2}/\\d{2}")   # MM/YY
      desc_raw <- stringr::str_remove(desc_raw, "\\s*\\d{4}$")     # Trailing YYYY
      desc_raw <- stringr::str_remove(desc_raw, "90er Jahre")      # Remove text after mapping
      desc_raw <- stringr::str_trim(desc_raw)

      desc_parts <- stringr::str_trim(stringr::str_split(desc_raw, "\\s*\\+\\s*")[[1]])
      descs[[i]] <- desc_parts
    }

    # Repeat years for multiple descriptions per item
    years_expanded <- unlist(mapply(function(y, d) rep(y, length(d)), years, descs, SIMPLIFY = FALSE))
    descs_expanded <- unlist(descs)

    as.vector(rbind(years_expanded, descs_expanded))
  }

  parsed <- lapply(text_vec, parse_one)
  max_len <- max(sapply(parsed, length))

  # Pad vectors and create data frame with pairs of year and description columns
  parsed_padded <- lapply(parsed, function(x) { length(x) <- max_len; x })
  df <- as.data.frame(do.call(rbind, parsed_padded), stringsAsFactors = FALSE)

  half <- max_len / 2
  colnames(df) <- paste0(
    rep(c("year", "desc"), half),
    rep(seq_len(half), each = 2)
  )

  df
}


# ========================================================================== #
# Read data and initial processing -------------------------------------------
# ========================================================================== #

# Read Excel file with chronic disease data
trigeminale_daten_table1 <- read_excel(
  "/home/joern/Aktuell/TrigeminalSensitivity/09Originale/Bormann Trigeminale Studie Daten.xlsx",
  sheet = "Tabelle1"
)

# Identify character variables
character_vars <- names(trigeminale_daten_table1)[sapply(trigeminale_daten_table1, is.character)]
print(character_vars)

# Convert selected character variables to numeric
character_to_numeric_vars <- c(
  "Alter",
  "Körpergröße",
  "Riechvermögen unmittelbar nach Covid-19",
  "Wie oft Covid-19?",
  "Trinken Sie Alkohol?",
  "Wenn ich einen frischen Minzkaugummi gekaut habe, habe ich das Gefühl besser Luft durch die Nase zu bekommen"
)
trigeminale_daten_table1[character_to_numeric_vars] <-
  lapply(trigeminale_daten_table1[character_to_numeric_vars], as.numeric)


# ========================================================================== #
# Clean and parse ENT surgery descriptions -----------------------------------
# ========================================================================== #

clean_text_vec <- ifelse(
  trigeminale_daten_table1$`OP im HNO-Bereich` %in% c("n", "", "NA", "NaN"),
  NA_character_,
  trigeminale_daten_table1$`OP im HNO-Bereich`
)

result_procedures <- parse_procedure_descriptions(clean_text_vec)


# ========================================================================== #
# Reshape parsed data to long format -----------------------------------------
# ========================================================================== #

long_df <- result_procedures %>%
  mutate(rowid = row_number()) %>%
  pivot_longer(
    cols = everything(),
    names_to = c(".value", "set"),
    names_pattern = "(year|desc)(\\d+)"
  ) %>%
  # Convert year to numeric
  mutate(year = as.numeric(year))


# ========================================================================== #
# Fix years and descriptions -------------------------------------------------
# ========================================================================== #

long_df_fixed <- long_df %>%
  # Leading 4-digit year in desc column overrides year column
  mutate(
    year_new = ifelse(str_detect(desc, "^\\d{4}"),
                      as.numeric(str_extract(desc, "^\\d{4}")),
                      year),
    desc = str_remove(desc, "^\\d{4}\\s*")
  ) %>%
  # Handle two-digit leading year (e.g., "05 Something") in desc
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
  # Special case: "er Tonsillotomie" => 1995 if no year present
  mutate(
    year_new = ifelse(
      is.na(year_new) & str_detect(desc, "^er\\s*Tonsillotomie"),
      1995,
      year_new
    ),
    desc = str_remove(desc, "^er\\s*")
  ) %>%
  mutate(desc = str_trim(desc)) %>%
  mutate(year = year_new) %>%
  select(-year_new) %>%
  filter(!is.na(desc) & desc != "")


# ========================================================================== #
# Translation dictionary for surgery descriptions ----------------------------
# ========================================================================== #

translation_map <- c(
  "Pansinusitis" = "Pansinusitis",
  "Adenotomie" = "Adenotomy",
  "Adenotomie Kindesalter" = "Adenotomy",
  "Kindesalter Adenotomie" = "Adenotomy",
  "Adentomie" = "Adenotomy",
  "Nasenfraktur" = "Nasal fracture",
  "Mandel OP" = "Tonsil surgery",
  "Mandel Op" = "Tonsil surgery",
  "Mandeln" = "Tonsil surgery",
  "Mandeln OP" = "Tonsil surgery",
  "Mandel-OP" = "Tonsil surgery",
  "Tonsillotomie" = "Tonsillotomy",
  "Tonsillotomie Kindheit" = "Tonsillotomy",
  "Kindesalter Tonsillotomie" = "Tonsillotomy",
  "Tonsillotomie Kindesalter" = "Tonsillotomy",
  "Tonsillotomie jugend" = "Tonsillotomy",
  "er Tonsillotomie" = "Tonsillotomy",
  "Tonsillektomie" = "Tonsillectomy",
  "Polypen" = "Nasal polyps surgery",
  "Kleinkind Polypen" = "Nasal polyps surgery",
  "Kindheit Polypen" = "Nasal polyps surgery",
  "Kindesalter Polypen" = "Nasal polyps surgery",
  "Wucherungen Nase" = "Nasal growth surgery",
  "NNH" = "Paranasal sinus surgery",
  "NNH Op" = "Paranasal sinus surgery",
  "NNH links" = "Left paranasal sinus surgery",
  "NSW" = "Nasal septum surgery",
  "NSW OP" = "Nasal septum surgery",
  "NSW Fensterung" = "Nasal septum fenestration",
  "Nasenscheidewand" = "Nasal septum surgery",
  "Begradigung Nasenscheidewand" = "Nasal septum correction",
  "Nasenscheidewand OP" = "Nasal septum surgery",
  "Verödung Nasenscheidewand" = "Nasal septum cauterization",
  "Nasenmuscheln Verkleinerung" = "Inferior turbinate reduction",
  "Reduktion Nasenmuscheln" = "Inferior turbinate reduction",
  "Verkleinerung Nasenmuscheln" = "Inferior turbinate reduction",
  "Laserung Nasenmuscheln" = "Inferior turbinate laser therapy",
  "Verödung Nasenmuschel" = "Inferior turbinate cauterization",
  "Verödung Nasenmuscheln" = "Inferior turbinate cauterization",
  "Stromwellentherapie zum Abschwellen der Nasenmuscheln" = "Inferior turbinate radiofrequency therapy",
  "Nasen # OP" = "Nasal fracture surgery",
  "Septum OP" = "Nasal septum surgery",
  "Septumkorrektur" = "Nasal septum correction",
  "Septumdeviation" = "Septal deviation correction",
  "Fensterung Kieferhöhle" = "Caldwell-Luc procedure",
  "Fensterung Kieferhöhle li." = "Left Caldwell-Luc procedure",
  "Fensterung" = "Fenestration",
  "Knorpelplastik Ohr" = "Cartilage plasty (ear)",
  "Ohr li. Implantat" = "Left ear implant",
  "Ohr anlegen" = "Otoplasty",
  "Otoplastik" = "Otoplasty",
  "Austausch Hämmerchen Ohr" = "Malleus replacement (middle ear)",
  "Trommelfell" = "Tympanic membrane surgery",
  "Trommelfell 2x" = "Tympanic membrane surgery",
  "Trommelfell-OP" = "Tympanic membrane surgery",
  "Trommelfellimplantat" = "Tympanic membrane implant",
  "Trommelfellverschluss" = "Tympanic membrane closure",
  "Trommelfellschnitt" = "Myringotomy",
  "Knalltrauma" = "Acoustic trauma",
  "Mandeln OP" = "Tonsil surgery",
  "Stimmband OP" = "Vocal cord surgery",
  "Stimmknötchen" = "Vocal cord nodules surgery",
  "Polypen" = "Nasal polyps surgery",
  "Mandeln" = "Tonsil surgery",
  "Mundbodenkarzinom" = "Mouth floor carcinoma surgery",
  "Parotisentfernung" = "Parotidectomy",
  "Karzinom Hals" = "Neck carcinoma surgery",
  "Medulläres Schildrüsenkarzinom" = "Medullary thyroid carcinoma surgery",
  "Karzinom Rachen li." = "Left pharyngeal carcinoma surgery",
  "OP Halsbereich" = "Neck surgery",
  "Abszess Rachen" = "Pharyngeal abscess surgery",
  "Trauma Hals" = "Neck trauma surgery",
  "Fensterung" = "Fenestration",
  "Begradigung Nase" = "Nasal straightening surgery",
  "Schilddrüse" = "Thyroid gland surgery",
  "Thyreodektomie" = "Thyroidectomy",
  "Austausch Hämmerchen Ohr" = "Malleus replacement (middle ear)",
  "Kiefernzyste" = "Jaw cyst removal",
  "OP Polypen Stirnhöhle" = "Frontal sinus polyp surgery",
  "Kindesalter Mandeln" = "Tonsil surgery",
  "Kindheit Polypen" = "Nasal polyps surgery",
  "j" = NA_character_,
  "er Tonsillotomie" = "Tonsillotomy",
  "Erweiterung Durchgänge NNH" = "Widening of paranasal sinus passages",
  "Abtragung Nasenschleimhaut" = "Nasal mucosa removal",
  "Stromwellentherapie zum Abschwellen der Nasenmuscheln" = "Inferior turbinate radiofrequency therapy",
  "Adenotomie" = "Adenotomy",
  "Kindesalter Adenotomie" = "Adenotomy",
  "Adenotomie Kindesalter" = "Adenotomy",
  "Tonsillotomie jugend" = "Tonsillotomy",
  "CRS Polypenentfernung 2x" = "Chronic rhinosinusitis polyp removal",
  "Karzinom Rachen li." = "Left pharyngeal carcinoma surgery",
  "Mundbodenkarzinom" = "Mouth floor carcinoma surgery",
  "Kleinkind Polypen" = "Nasal polyps surgery",
  "Trommelfell-OP" = "Tympanic membrane surgery",
  "Siebbeinhöhle" = "Ethmoid sinus surgery",
  "NNH" = "Paranasal sinus surgery"
)


# ========================================================================== #
# Translate surgery descriptions ---------------------------------------------
# ========================================================================== #

long_df_fixed <- long_df_fixed %>%
  mutate(
    desc_english = translation_map[desc],
    desc_english = ifelse(is.na(desc_english), "Unknown", desc_english),
    desc_english = factor(desc_english, levels = sort(unique(desc_english)))
  )


# ========================================================================== #
# One-hot encode diagnoses for analysis --------------------------------------
# ========================================================================== #

# Create row_id assuming 'set' (from pivot_longer) reflects original rows
long_df_fixed <- long_df_fixed %>%
  mutate(row_id = as.integer(set))

# Filter out "Unknown" diagnoses
df_filtered <- long_df_fixed %>%
  filter(desc_english != "Unknown")

# One-hot encode per row_id and diagnosis
one_hot_partial <- df_filtered %>%
  distinct(row_id, desc_english) %>%
  mutate(present = 1) %>%
  pivot_wider(names_from = desc_english, values_from = present, values_fill = 0)

# Generate full set of row_ids and join to keep all rows, fill missing with zero
all_rows <- tibble(row_id = 1:max(long_df_fixed$row_id, na.rm = TRUE))

one_hot_full <- all_rows %>%
  left_join(one_hot_partial, by = "row_id") %>%
  replace(is.na(.), 0) %>%
  arrange(row_id) %>%
  select(-row_id)

print(dim(one_hot_full))
print(head(one_hot_full))


# ========================================================================== #
# Plot preparation and execution ---------------------------------------------
# ========================================================================== #

min_year <- min(long_df_fixed$year, na.rm = TRUE)

long_df_fixed <- long_df_fixed %>%
  mutate(year_plot = as.numeric(year),
         year_plot = ifelse(is.na(year_plot), min_year - 2, year_plot))

desc_levels <- sort(unique(long_df_fixed$desc_english))
long_df_fixed$desc_english <- factor(long_df_fixed$desc_english, levels = desc_levels)

stripe_df <- data.frame(
  desc_english = desc_levels,
  xmin = seq_along(desc_levels) - 0.5,
  xmax = seq_along(desc_levels) + 0.5
)

count_df <- long_df_fixed %>%
  group_by(desc_english) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  mutate(xpos = as.numeric(desc_english),
         ypos = max(long_df_fixed$year_plot, na.rm = TRUE) + 3)

if (order_plot) {
  desc_levels_ordered <- count_df %>%
    arrange(desc(count)) %>%
    pull(desc_english)

  long_df_fixed$desc_english <- factor(long_df_fixed$desc_english, levels = desc_levels_ordered)

  stripe_df <- data.frame(
    desc_english = desc_levels_ordered,
    xmin = seq_along(desc_levels_ordered) - 0.5,
    xmax = seq_along(desc_levels_ordered) + 0.5
  )

  count_df <- count_df %>%
    mutate(xpos = match(desc_english, desc_levels_ordered))
}

set.seed(42)  # For reproducible jitter

p_ENT_surgery <- ggplot() +
  geom_rect(data = stripe_df,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = rep(c("ivory2", "white"), length.out = nrow(stripe_df)),
            inherit.aes = FALSE) +
  geom_jitter(data = long_df_fixed,
              aes(x = desc_english, y = year_plot),
              width = 0.3, height = 0.3, alpha = 0.6, color = "black") +
  geom_text(data = count_df,
            aes(x = xpos, y = ypos, label = count),
            vjust = 0, size = 3, angle = 90) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(
    breaks = c(min_year - 2, seq(min_year, max(long_df_fixed$year_plot, na.rm = TRUE))),
    labels = c("No year", as.character(seq(min_year, max(long_df_fixed$year_plot, na.rm = TRUE))))
  ) +
  labs(x = "Disease", y = "Year", title = "ENT surgeries over years") +
  theme_minimal(base_size = 8) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "black", size = 0.3, linetype = "dashed")
  ) +
  geom_hline(yintercept = seq(min_year, max(long_df_fixed$year_plot, na.rm = TRUE)),
             color = "grey80", size = 0.3, linetype = "dashed")

print(p_ENT_surgery)

ggsave("p_ENT_surgery.svg", p_ENT_surgery, width = 12, height = 12)



