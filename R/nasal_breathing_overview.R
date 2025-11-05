# Load necessary libraries
# ========================================================================== #
# Load required libraries ----------------------------------------------------
# ========================================================================== #

library(readxl)      # For reading Excel files
library(stringr)     # String manipulation functions
library(dplyr)       # Data manipulation verbs (filter, mutate, select, etc.)
library(tidyr)       # Data tidying and reshaping
library(ggplot2)     # Grammar of graphics plotting system
library(lubridate)   # Date and time manipulation
library(forcats)     # Factor level utilities
library(scales)      # Scaling functions for visualization
library(purrr)       # Functional programming utilities


# ========================================================================== #
# Read trigeminal sensitivity data from Excel --------------------------------
# ========================================================================== #

trigem_data <- read_excel(
  "/home/joern/Aktuell/TrigeminalSensitivity/09Originale/Bormann Trigeminale Studie Daten.xlsx",
  sheet = "Tabelle1"
)


# ========================================================================== #
# Convert selected character columns to numeric -----------------------------
# ========================================================================== #

# Specify character columns that should be numeric
character_to_numeric_vars <- c(
  "Alter",
  "Körpergröße",
  "Riechvermögen unmittelbar nach Covid-19",
  "Wie oft Covid-19?",
  "Trinken Sie Alkohol?",
  "Wenn ich einen frischen Minzkaugummi gekaut habe, habe ich das Gefühl besser Luft durch die Nase zu bekommen"
)

# Convert those columns to numeric, suppressing warnings if any
trigem_data[character_to_numeric_vars] <- lapply(
  trigem_data[character_to_numeric_vars],
  function(x) as.numeric(x)
)


# ========================================================================== #
# Subset and rename nasal breathing problem related variables ----------------
# ========================================================================== #

nasal_breathing_vars <- c(
  "R23", "R24", "R25",
  "Ärztliche Vorstellung wegen Problematik der Nasenatmung",
  "Wenn ja erfolgte eine Therapie? - Welche?"
)

nasal_breathing_df <- trigem_data[, nasal_breathing_vars] %>%
  rename(
    raw_text = `Ärztliche Vorstellung wegen Problematik der Nasenatmung`,
    therapy_text = `Wenn ja erfolgte eine Therapie? - Welche?`
  ) %>%
  mutate(age = as.integer(trigem_data$Alter))

current_year <- year(Sys.Date())


# ========================================================================== #
# Define function to extract years robustly from messy raw text --------------
# ========================================================================== #

extract_years_custom <- function(text, age, current_year = 2025) {
  if (is.na(text) || text == "n") {
    return(NA_character_)
  }

  text <- str_trim(text)

  # Remove leading "j " if present
  if (str_detect(text, "^j\\s+")) {
    text <- str_replace(text, "^j\\s+", "")
  }

  # Excel date serial number format
  if (str_detect(text, "^\\d{5,}$")) {
    return(format(as.Date(as.numeric(text), origin = "1899-12-30"), "%Y"))
  }

  # Year ranges (e.g., "2014-2016")
  if (str_detect(text, "\\d{4}-\\d{4}")) {
    years <- as.numeric(str_split(text, "-", simplify = TRUE))
    return(as.character(seq(years[1], years[2])))
  }

  # Slash-separated years (e.g., "2018/19")
  if (str_detect(text, "\\d{4}/\\d{2}")) {
    first_year <- as.numeric(str_sub(text, 1, 4))
    second_year <- as.numeric(paste0(str_sub(first_year, 1, 2), str_sub(text, 6, 7)))
    return(as.character(c(first_year, second_year)))
  }

  # Plus-separated years (e.g., "1997 + 2023")
  if (str_detect(text, "\\d{4}\\s*\\+\\s*\\d{4}")) {
    return(str_extract_all(text, "\\d{4}")[[1]])
  }

  # Single year, possibly with "ab " or "Winter"
  if (str_detect(text, "^(ab\\s*)?\\d{4}$")) {
    return(str_extract(text, "\\d{4}"))
  }
  if (str_detect(text, "^Winter\\s*\\d{4}$")) {
    return(str_extract(text, "\\d{4}"))
  }

  # Approximate years based on age
  if (text == "Grundschüler" && !is.na(age)) {
    return(as.character(current_year - age + 8))
  }
  if (text == "zur Geburt" && !is.na(age)) {
    return(as.character(current_year - age))
  }

  # Keywords indicating ongoing or periodic events
  if (str_detect(text, "regelmäßig|halbjährlich|Quartal|Kindesalter")) {
    return("ongoing")
  }

  # Plain 4-digit years
  if (str_detect(text, "^\\d{4}$")) {
    return(text)
  }

  # No match
  return(NA_character_)
}


# ========================================================================== #
# Apply year extraction and nasal breathing problem flag --------------------
# ========================================================================== #

nasal_breathing_df <- nasal_breathing_df %>%
  rowwise() %>%
  mutate(
    extracted_years = list(extract_years_custom(raw_text, age, current_year)),

    actual_nasal_breathing_problem = if_else(
      str_detect(raw_text, "regelmäßig|halbjährlich|Quartal|Kindesalter") |
        (raw_text == "j" & all(is.na(unlist(extracted_years)))),
      1, 0
    ),

    nasal_breathing_problem = if_else(
      !all(is.na(unlist(extracted_years))) |
        str_detect(raw_text, "regelmäßig|halbjährlich|Quartal|Kindesalter|Geburt|Grundschüler") |
        (raw_text == "j"),
      1, 0
    )
  ) %>%
  ungroup()

# Clean extracted years for easier readability
nasal_breathing_df$years_clean <- map_chr(nasal_breathing_df$extracted_years, ~paste(.x, collapse = ","))


# ========================================================================== #
# Translation dictionary for therapy text ------------------------------------
# ========================================================================== #

therapy_dict <- c(
  "n" = "No therapy",
  "j" = "Unknown",
  "OP Nasenscheidewand" = "Nasal septum surgery (septoplasty)",
  "OP NSW" = "Nasal septum surgery (septoplasty)",
  "OP" = "Operation",
  "Septumplastik" = "Septoplasty (nasal septum surgery)",
  "regelmäßige Untersuchungen" = "Regular checkups",
  "Desensibilisierung" = "Desensitization",
  "Cortisonspray" = "Cortisone spray",
  "Verödung" = "Cauterization",
  "Atemgerät Schlaf" = "Sleep breathing device",
  "Antibiotika" = "Antibiotics",
  "Nasenspray" = "Nasal spray",
  "Bedarfsspray Asthma" = "Asthma rescue spray",
  "Antiallergika" = "Antiallergics",
  "Polypen-Op" = "Polyp operation",
  "Asthmadiagnostik" = "Asthma diagnostics",
  "Op geplant" = "Operation planned",
  "Salbe" = "Ointment",
  "Doprident" = "Doprident",
  "Meersalzinhalation" = "Sea salt inhalation",
  "Cortison-NS" = "Cortisone nasal spray",
  "Creme" = "Cream",
  "Mometason" = "Mometasone",
  "AB" = "Antibiotics",
  "Rotlichthterapie" = "Red light therapy",
  "Stromwellentherapie" = "Radio wave therapy",
  "Trommelfellschnitt" = "Eardrum incision",
  "Akupunktur" = "Acupuncture",
  "Hyposensibilisierung" = "Hyposensitization",
  "Lasertherapie" = "Laser therapy",
  "Begradigung NSW" = "Nasal septum straightening",
  "Kortison-NS" = "Cortisone nasal spray"
)


# ========================================================================== #
# Translate therapy descriptions ---------------------------------------------
# ========================================================================== #

translate_therapy <- function(therapy_text) {
  if (is.na(therapy_text)) {
    return(NA_character_)
  }

  match_key <- names(therapy_dict)[
    map_lgl(names(therapy_dict), ~ str_detect(therapy_text, fixed(.x, ignore_case = TRUE)))
  ]

  if (length(match_key) > 0) {
    return(therapy_dict[match_key[1]])
  }

  # Default: Capitalize first letter of original text
  str_to_sentence(therapy_text)
}


# Translate all therapies in the dataset
nasal_breathing_df <- nasal_breathing_df %>%
  mutate(
    therapy_years = map(therapy_text, ~str_extract_all(.x, "\\b\\d{4}\\b")[[1]]),
    therapy_english = map_chr(therapy_text, translate_therapy)
  )


# ========================================================================== #
# Prepare data and plot nasal breathing problems timeline -------------------
# ========================================================================== #

library(tidyr) # needed for unnest

plot_nasal_breathing <- nasal_breathing_df %>%
  filter(nasal_breathing_problem == 1) %>%
  mutate(
    years_expanded = if_else(is.na(extracted_years), list(NA_character_), extracted_years)
  ) %>%
  unnest(years_expanded) %>%
  mutate(
    years_expanded = if_else(years_expanded == "ongoing", NA_character_, years_expanded),
    year_num = as.numeric(years_expanded),
    case_id = row_number()
  )

# Define x-axis positions for missing year points
min_year <- min(plot_nasal_breathing$year_num, na.rm = TRUE)
no_year_xpos <- min_year - 2

with_year_points <- plot_nasal_breathing %>%
  filter(!is.na(year_num)) %>%
  select(case_id, year_num)

without_year_points <- plot_nasal_breathing %>%
  filter(is.na(year_num) & actual_nasal_breathing_problem == 1) %>%
  distinct(case_id, .keep_all = TRUE) %>%
  mutate(year_num = no_year_xpos) %>%
  select(case_id, year_num)

all_points <- bind_rows(with_year_points, without_year_points)

unique_years <- sort(unique(with_year_points$year_num))
x_breaks <- c(no_year_xpos, unique_years)
x_labels <- c("No year", as.character(unique_years))


# Final nasal breathing problems plot
nasal_breathing_plot <- ggplot(all_points, aes(x = year_num, y = case_id)) +
  geom_point(size = 3) +
  scale_y_continuous(name = "Case number", breaks = unique(all_points$case_id)) +
  scale_x_continuous(
    name = "Year",
    breaks = x_breaks,
    labels = x_labels,
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  labs(title = "Nasal breathing problems by case and year") +
  theme_minimal(base_size = 8) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

print(nasal_breathing_plot)

# ========================================================================== #
# Prepare and plot therapy occurrence data -----------------------------------
# ========================================================================== #

df_therapy_long <- nasal_breathing_df %>%
  filter(!is.na(therapy_english)) %>%
  mutate(case_id = row_number()) %>%
  select(case_id, therapy_english, years_clean) %>%
  mutate(years_list = strsplit(as.character(years_clean), ",")) %>%
  unnest(years_list) %>%
  mutate(
    years_list = str_trim(years_list),
    year_num = suppressWarnings(as.numeric(years_list))
  )

pseudoyear_pos <- min_year - 2

df_therapy_long <- df_therapy_long %>%
  mutate(year_plot = if_else(is.na(year_num), pseudoyear_pos, year_num))

therapy_levels <- sort(unique(df_therapy_long$therapy_english))

df_therapy_long <- df_therapy_long %>%
  mutate(therapy_english = factor(therapy_english, levels = therapy_levels))

stripe_df <- data.frame(
  therapy_english = therapy_levels,
  xmin = seq_along(therapy_levels) - 0.5,
  xmax = seq_along(therapy_levels) + 0.5
)

count_df <- df_therapy_long %>%
  group_by(therapy_english) %>%
  summarise(count = n()) %>%
  mutate(
    xpos = as.numeric(therapy_english),
    ypos = max(df_therapy_long$year_plot, na.rm = TRUE) + 3
  )

nasal_breathing_therapy_plot <- ggplot() +
  geom_rect(data = stripe_df,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = rep(c("ivory2", "white"), length.out = nrow(stripe_df)),
            inherit.aes = FALSE) +
  geom_jitter(data = df_therapy_long,
              aes(x = therapy_english, y = year_plot),
              width = 0.3, height = 0.3, alpha = 0.6,
              color = "black", show.legend = FALSE) +
  geom_text(data = count_df,
            aes(x = xpos, y = ypos, label = count),
            vjust = 0, size = 3, angle = 90) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(
    breaks = c(pseudoyear_pos, seq(min_year, max(df_therapy_long$year_plot, na.rm = TRUE))),
    labels = c("No year", as.character(seq(min_year, max(df_therapy_long$year_plot, na.rm = TRUE))))
  ) +
  labs(x = "Therapy", y = "Year",
       title = "Therapies for nasal breathing problems over years") +
  theme_minimal(base_size = 8) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(
      color = "black", size = 0.3, linetype = "dashed"
    )
  ) +
  geom_hline(
    yintercept = seq(min_year, max(df_therapy_long$year_plot, na.rm = TRUE)),
    color = "grey80", size = 0.3, linetype = "dashed"
  )

print(nasal_breathing_therapy_plot)

ggsave(
  filename = "p_nasal_breathing_therapy_plot.svg",
  plot = nasal_breathing_therapy_plot,
  width = 12, height = 12
)


# ========================================================================== #
# Create therapy one-hot encoding matrix --------------------------------------
# ========================================================================== #

nasal_breathing_df$therapy_english <- factor(nasal_breathing_df$therapy_english)

therapy_one_hot_matrix <- model.matrix(~ 0 + therapy_english, data = nasal_breathing_df)

head(therapy_one_hot_matrix)

