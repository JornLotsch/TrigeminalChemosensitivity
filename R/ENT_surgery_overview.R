library(stringr)
library(dplyr)

# Switches
order_plot <- FALSE

# Functions
parse_procedure_descriptions <- function(text_vec) {
  parse_one <- function(text) {
    # Split by comma, trim spaces, remove "seit"
    items <- stringr::str_split(text, ",\\s*")[[1]]
    items <- items[items != ""]
    items <- stringr::str_trim(items)
    items <- items[items != "n"] # also skip "n" here, if it denotes filler

    extract_year <- function(s) {
      s <- stringr::str_trim(s)
      if (is.na(s) || s == "") return(NA_character_)

      # Match normal 4-digit year at start
      y <- stringr::str_extract(s, "^\\d{4}")
      if (!is.na(y)) return(y)

      # Match year after slash, e.g. /74 or /2021
      y <- stringr::str_extract(s, "/(\\d{2,4})")
      if (!is.na(y)) {
        y <- sub("/", "", y)
        if (nchar(y) == 2) {
          y_num <- as.numeric(y)
          y <- ifelse(y_num < 50, paste0("20", sprintf("%02d", y_num)), paste0("19", y))
        }
        return(y)
      }

      # Match MM/YY format e.g. 05/21 (take only year)
      y <- stringr::str_extract(s, "\\d{2}/(\\d{2})")
      if (!is.na(y)) {
        yy <- sub("\\d{2}/", "", y)
        y_num <- as.numeric(yy)
        y <- ifelse(y_num < 50, paste0("20", sprintf("%02d", y_num)), paste0("19", yy))
        return(y)
      }

      # Match 90er Jahre style - roughly map to mid 1990s
      if (stringr::str_detect(s, "90er Jahre")) return("1995")

      NA_character_
    }

    years <- sapply(items, extract_year)

    descs <- vector("list", length(items))
    for (i in seq_along(items)) {
      desc_raw <- items[i]
      # Remove all matched year patterns from string
      desc_raw <- stringr::str_remove(desc_raw, "^\\d{4}\\s*") # leading 4 digit YYYY
      desc_raw <- stringr::str_remove(desc_raw, "/\\d{2,4}") # /YY or /YYYY
      desc_raw <- stringr::str_remove(desc_raw, "\\d{2}/\\d{2}") # MM/YY
      desc_raw <- stringr::str_remove(desc_raw, "\\s*\\d{4}$") # trailing YYYY
      desc_raw <- stringr::str_remove(desc_raw, "90er Jahre") # remove this text after mapping
      desc_raw <- stringr::str_trim(desc_raw)

      # Split multiple procs connected by '+'
      desc_parts <- stringr::str_trim(stringr::str_split(desc_raw, "\\s*\\+\\s*")[[1]])
      descs[[i]] <- desc_parts
    }

    years_expanded <- unlist(mapply(function(y, d) rep(y, length(d)), years, descs, SIMPLIFY = FALSE))
    descs_expanded <- unlist(descs)

    as.vector(rbind(years_expanded, descs_expanded))
  }

  parsed <- lapply(text_vec, parse_one)
  max_len <- max(sapply(parsed, length))

  parsed_padded <- lapply(parsed, function(x) { length(x) <- max_len; x })
  df <- as.data.frame(do.call(rbind, parsed_padded), stringsAsFactors = FALSE)

  half <- max_len / 2
  colnames(df) <- paste0(rep(c("year", "desc"), half), rep(seq_len(half), each = 2))

  df
}

# =============================== #
# Read data
# =============================== #

trigeminale_daten_corrected_translated <- read.csv("trigeminale_daten_corrected_translated.csv", check.names = FALSE)


# Clean-up the OP columns before extracting information
clean_text_vec <- ifelse(trigeminale_daten_corrected_translated$`Surgery in ENT region` %in% c("n", "", "NA", "NaN"), NA_character_,
                         trigeminale_daten_corrected_translated$`Surgery in ENT region` )


# Split ENT surgery descriptions into year/desc pairs
result_procedures <- parse_procedure_descriptions(clean_text_vec)

# Convert result to long format using tidyr pivot_longer
long_df <- result_procedures %>%
  mutate(rowid = row_number()) %>%
  pivot_longer(cols = everything(),
               names_to = c(".value", "set"),
               names_pattern = "(year|desc)(\\d+)") %>%
# filter(!is.na(desc) & desc != "") %>%       # Keep relevant desc entries
mutate(year = as.numeric(year))


long_df_fixed <- long_df %>%
# Handle leading year (4-digit)
mutate(
    year_new = ifelse(str_detect(desc, "^\\d{4}"),
                      as.numeric(str_extract(desc, "^\\d{4}")),
                      year),
    desc = str_remove(desc, "^\\d{4}\\s*")
  ) %>%
# Handle two-digit leading year (e.g., "05 ...")
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
# Special case: "er Tonsillotomie" means 1995 if no year is present
mutate(
    year_new = ifelse(
      is.na(year_new) & str_detect(desc, "^er\\s*Tonsillotomie"),
      1995,
      year_new
    ),
    desc = str_remove(desc, "^er\\s*") # Remove 'er' if it was used
  ) %>%
# Clean multiple spaces etc.
mutate(
    desc = str_trim(desc)
  ) %>%
# Replace year
mutate(year = year_new) %>%
  dplyr::select(-year_new) %>%
  dplyr::filter(!is.na(desc) & desc != "")


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

# Translate German descriptions to English
long_df_fixed <- long_df_fixed %>%
  mutate(desc_english = translation_map[desc]) %>%
  mutate(desc_english = ifelse(is.na(desc_english), "Unknown", desc_english)) %>%
  mutate(desc_english = factor(desc_english, levels = sort(unique(desc_english))))

## Create one hot data frame for later analysis ##
long_df_fixed_one_hot <- long_df_fixed

# Create unique row ids based on your original row indicator, assuming 'set' is row id (or create one if none)
df_long <- long_df_fixed_one_hot %>%
  mutate(row_id = as.integer(set))

# Filter out non-disease entries ("Unknown") from the one-hot process
df_filtered <- df_long %>%
  filter(desc_english != "Unknown")

# One-hot encode desc_english per row_id
one_hot_partial <- df_filtered %>%
  distinct(row_id, desc_english) %>% # unique diseases per row
mutate(present = 1) %>%
  pivot_wider(names_from = desc_english, values_from = present, values_fill = 0)

# Create a data frame with all original row_ids (1 to 1001)
all_rows <- tibble(row_id = 1:1001)

# Join to keep all rows, fill NAs with 0 meaning no presence
one_hot_full <- all_rows %>%
  left_join(one_hot_partial, by = "row_id") %>%
  replace(is.na(.), 0) %>%
  arrange(row_id)

# Optionally drop the row_id if unneeded or rename
one_hot_full <- one_hot_full %>% dplyr::select(-row_id)

print(dim(one_hot_full)) # Should be 1001 rows
print(head(one_hot_full))

## Continue with plot ##


# Determine the minimum observed year for later use
min_year <- min(long_df_fixed$year, na.rm = TRUE)

# Create plotting year column: use NA-year mapped to "pseudo-year" below the minimum
long_df_fixed <- long_df_fixed %>%
  mutate(year_plot = as.numeric(year),
         year_plot = ifelse(is.na(year_plot), min_year - 2, year_plot))

# Remove unknown diagnoses (e.g. those not in translation map)
long_df_fixed <- long_df_fixed %>%
  filter(desc_english != "Unknown")

# Reset factor levels to reflect actual data (post filtering)
desc_levels <- sort(unique(long_df_fixed$desc_english))
long_df_fixed$desc_english <- factor(long_df_fixed$desc_english, levels = desc_levels)

# Prepare data frame to add vertical stripes behind points
stripe_df <- data.frame(
  desc_english = desc_levels,
  xmin = seq_along(desc_levels) - 0.5,
  xmax = seq_along(desc_levels) + 0.5
)

# Calculate counts per disease for annotation above plot
count_df <- long_df_fixed %>%
  group_by(desc_english) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  mutate(xpos = as.numeric(desc_english),
         ypos = max(long_df_fixed$year_plot, na.rm = TRUE) + 3)

if (order_plot) {
  # Order by descending count
  desc_levels <- count_df %>%
    arrange(desc(count)) %>%
    pull(desc_english)

  # Apply new factor levels to long_df
  long_df$desc_english <- factor(long_df$desc_english, levels = desc_levels)

  # Also reorder stripe_df and count_df accordingly
  stripe_df <- data.frame(
    desc_english = desc_levels,
    xmin = seq_along(desc_levels) - 0.5,
    xmax = seq_along(desc_levels) + 0.5
  )
  count_df <- count_df %>%
    mutate(xpos = match(desc_english, desc_levels))
}

# Seed for reproducible jitter position
set.seed(42)

# Generate plot
p_ENT_surgery <- ggplot() +
# Vertical color stripes alternating ivory2 and white background
geom_rect(data = stripe_df,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = rep(c("ivory2", "white"), length.out = nrow(stripe_df)),
            inherit.aes = FALSE) +
# Jitter points for diagnosis occurrences by year and disease
geom_jitter(data = long_df_fixed,
              aes(x = desc_english, y = year_plot),
              width = 0.3, height = 0.3, alpha = 0.6, color = "black") +
# Vertical text counts above stripes (numbers of cases per disease)
geom_text(data = count_df,
            aes(x = xpos, y = ypos, label = count),
            vjust = 0, size = 3, angle = 90) +
# X axis categorical labels for disease, preserve all levels
scale_x_discrete(drop = FALSE) +
# Y axis with breaks including "No year" pseudo-year label below the min_year
scale_y_continuous(
    breaks = c(min_year - 2, seq(min_year, max(long_df_fixed$year_plot, na.rm = TRUE))),
    labels = c("No year", as.character(seq(min_year, max(long_df_fixed$year_plot, na.rm = TRUE))))
  ) +
# Axis labels and title
labs(x = "Disease", y = "Year", title = "ENT surgeries over years") +
# Minimal theme with smaller font size
theme_minimal(base_size = 8) +
# Customize axis text and grid
theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.grid.major.x = element_blank(), # remove vertical grid lines
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "black", size = 0.3, linetype = "dashed")
  ) +
# Add horizontal dashed lines at each year (behind points)
geom_hline(yintercept = seq(min_year, max(long_df_fixed$year_plot, na.rm = TRUE)),
             color = "grey80",
             size = 0.3,
             linetype = "dashed")

p_ENT_surgery

ggsave(paste0("p_ENT_surgery", ".svg"), p_ENT_surgery, width = 12, height = 12)
