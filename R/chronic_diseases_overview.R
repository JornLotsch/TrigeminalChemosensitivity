# Load necessary libraries
library(readxl)
library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggtext)

# Switches
order_plot <- FALSE

# Function to split chronic disease strings into alternating year/description pairs
split_year_description <- function(text_vec) {
  parse_one <- function(text) {
    # Split by comma, trim spaces, remove "seit"
    items <- str_split(text, ",\\s*")[[1]]
    items <- str_trim(str_replace_all(items, "\\bseit\\b\\s*", ""))

    years <- vector("character", length(items))
    descs <- vector("list", length(items))

    for (i in seq_along(items)) {
      # Extract ANY 4-digit year (not just at start)
      year_match <- str_extract(items[i], "\\b\\d{4}\\b")
      years[i] <- year_match

      # Remove year if present from description
      desc_raw <- str_trim(str_replace(items[i], "\\b\\d{4}\\b", ""))

      # Split multiple diagnoses connected by '+'
      desc_parts <- str_trim(str_split(desc_raw, "\\s*\\+\\s*")[[1]])
      descs[[i]] <- desc_parts
    }

    # Repeat year for each desc part, or NA if year missing
    years_expanded <- unlist(mapply(function(y, d) rep(y, length(d)), years, descs, SIMPLIFY = FALSE))
    descs_expanded <- unlist(descs)

    # Interleave year and desc as vector: year1, desc1, year2, desc2, ...
    as.vector(rbind(years_expanded, descs_expanded))
  }

  parsed <- lapply(text_vec, parse_one)
  max_len <- max(sapply(parsed, length))
  # Pad to rectangular structure
  parsed_padded <- lapply(parsed, function(x) {
    length(x) <- max_len
    x
  })
  df <- as.data.frame(do.call(rbind, parsed_padded), stringsAsFactors = FALSE)
  half <- max_len / 2
  colnames(df) <- paste0(rep(c("year", "desc"), half), rep(seq_len(half), each = 2))
  return(df)
}


# Read Excel file with chronic disease data
trigeminale_daten_table1 <- read_excel(
  "/home/joern/Aktuell/TrigeminalSensitivity/09Originale/Bormann Trigeminale Studie Daten.xlsx",
  sheet = "Tabelle1"
)

# Split chronic disease descriptions into year/desc pairs
result_chronic <- split_year_description(trigeminale_daten_table1$`chronische Erkrankungen`)

# Split neurological disease descriptions into year/desc pairs
result_neurological <- split_year_description(trigeminale_daten_table1$`neurologische Erkrankung`)
names(result_neurological) <- paste0(rep(c("year", "desc"), dim(result_neurological)[2] / 2),
                                     rep((dim(result_chronic)[2] / 2 + 1):(dim(result_chronic)[2] / 2 + dim(result_neurological)[2] / 2), each = 2))

result_chronic_and_neurological <- cbind.data.frame(result_chronic, result_neurological)

# Convert result to long format using tidyr pivot_longer
long_df <- result_chronic_and_neurological %>%
  mutate(rowid = row_number()) %>%
  pivot_longer(cols = everything(),
               names_to = c(".value", "set"),
               names_pattern = "(year|desc)(\\d+)") %>%
  filter(!is.na(desc) & desc != "") %>% # Keep relevant desc entries
mutate(year = as.numeric(year))



# Translation map from German to English (corrected typos included)
translation_map <- c(
  "Achalasie" = "Achalasia",
  "koromare Herzkrankheit" = "Coronary heart disease",
  "Hypercholesterinämie" = "Hypercholesterolemia",
  "Hypertonie" = "Hypertension",
  "essentielle Hypertonie" = "Essential hypertension",
  "Hyprtonie" = "Hypertension",
  "arterielle Hypertonie" = "Arterial hypertension",
  "art. Hypertonie" = "Arterial hypertension",
  "Hypertonus" = "Hypertension (variant term)",
  "Thyreoektomie > behandelte Hypothyreose" = "Hypothyroidism (treated after thyroidectomy)",
  "Hypothyreose" = "Hypothyroidism",
  "Geburt Hypothyreose" = "Congenital hypothyroidism",
  "PTBS" = "Post-traumatic stress disorder",
  "Depression" = "Depression",
  "Depressionen" = "Depression",
  "chronische Schmerzstörung mit somatischen und psychischen Faktoren" = "Chronic pain disorder",
  "allergisches Asthma" = "Allergic asthma",
  "Asthma" = "Asthma",
  "Asthma bronchiale" = "Bronchial asthma",
  "Asthma bronchiale Geburt" = "Congenital asthma",
  "Asthma Bronchiale Kindheit" = "Childhood bronchial asthma",
  "Asthma Kindesalter" = "Childhood asthma",
  "Asthma Geburt" = "Congenital asthma",
  "Asthma pollinose" = "Asthma with hay fever",
  "D.m." = "Diabetes mellitus (generic)",
  "D.m. Typ l" = "Diabetes mellitus type I",
  "D.m. Typ ll" = "Diabetes mellitus type II",
  "Diabethes" = "Diabetes mellitus (generic, typo)",
  "Diabethes mellitus" = "Diabetes mellitus (generic, typo)",
  "Diabethes mellitus Typ ll" = "Diabetes mellitus type II",
  "Diabethes nellitus Typ 2" = "Diabetes mellitus type II",
  "Diabethes Typ l" = "Diabetes mellitus type I",
  "Diabethes in Diagnostik" = "Diabetes under investigation",
  "Insulinresistenz (Diabethes)" = "Insulin resistance (diabetes type II)",
  "Mb. Chron" = "Crohn’s disease",
  "Morbus Crohn" = "Crohn’s disease",
  "Rheuma" = "Rheumatism",
  "rheumatische Atritis" = "Rheumatoid arthritis",
  "Arthrose" = "Osteoarthritis",
  "Weichteilrheuma" = "Soft tissue rheumatism",
  "Psoriasis-Arthritis" = "Psoriatic arthritis",
  "Neurodermitis" = "Atopic dermatitis",
  "Neurodermitis Geburt" = "Congenital atopic dermatitis",
  "Schuppenflechte" = "Psoriasis",
  "Mb. Basedow" = "Graves’ disease",
  "Hashimoto" = "Hashimoto’s thyroiditis",
  "Hasimodo" = "Hashimoto’s thyroiditis",
  "Hashimodo" = "Hashimoto’s thyroiditis",
  "Schilddrüsenprobleme" = "Thyroid problems",
  "z.N. Nierentransplantation" = "Kidney transplant (status post)",
  "Nieren-Op" = "Kidney surgery",
  "Niereninsuffizienz" = "Kidney insufficiency",
  "Nierenstauung" = "Kidney congestion",
  "Nierensteine" = "Kidney stones",
  "Nierensteine?" = "Kidney stones (suspected)",
  "Nierenreilresektion nach Tumor" = "Partial kidney resection after tumor",
  "schlechte Nierenwerte" = "Poor kidney function",
  "PBC" = "Primary biliary cirrhosis",
  "Arteriosklerose" = "Arteriosclerosis",
  "Herzinsuffizienz" = "Heart failure",
  "Herzklappeninsuffizienz l" = "Heart valve insufficiency (left)",
  "Herzrhythmusstörung" = "Cardiac arrhythmia",
  "koromare Herzkrankheit" = "Coronary heart disease",
  "Fibromyalgie" = "Fibromyalgia",
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
  "Dissoziation" = "Dissociation",
  "Zwangsstörung" = "Obsessive-compulsive disorder",
  "Schlafapnoe" = "Sleep apnea",
  "Krebs" = "Cancer",
  "Hirntumor" = "Brain tumor",
  "Tumor kopf" = "Head tumor",
  "Depressionen" = "Depression",
  "Asperger Autismus" = "Asperger’s syndrome",
  "Raynaud Syndrom" = "Raynaud’s syndrome",
  "Gürtelrose" = "Shingles",
  "Pollinosis" = "Allergic rhinitis (hay fever)",
  "Heuschnupfen" = "Hay fever",
  "Heufieber" = "Hay fever",
  "Allergien" = "Allergies",
  "Allergie" = "Allergy",
  "Allergie Hausstaub" = "House dust allergy",
  "Hausstauballergie" = "House dust allergy",
  "Hausstaubmilbenallergie" = "House dust mite allergy",
  "Katze" = "Cat allergy",
  "Katzenhaarallergie" = "Cat hair allergy",
  "Lebensmittelallergien" = "Food allergies",
  "Histaminintoleranz" = "Histamine intolerance",
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
  "Mb. Bechterew" = "Ankylosing spondylitis",
  "Mb. Scheuermann" = "Scheuermann’s disease",
  "PCOS" = "Polycystic ovary syndrome",
  "Thrombozytose" = "Thrombocytosis",
  "Zölliakie" = "Celiac disease",
  "einlaufende Demenz" = "Onset dementia",
  "chron. Bronchitis" = "Chronic bronchitis",
  "chr. Entzündung der Schulter re." = "Chronic inflammation of right shoulder"
)

## Create one hot data frame for later analysis ##
# Get names of the desc columns programmatically
desc_cols <- grep("^desc\\d+$", names(result_chronic_and_neurological), value = TRUE)

result_chronic_and_neurological_eng <- result_chronic_and_neurological[, desc_cols]

# Translation function to replace German terms with English, leave others unchanged
translate_vector <- function(vec, translation_map) {
  vec_char <- as.character(vec)
  translated <- translation_map[vec_char]
  translated[is.na(translated)] <- vec_char[is.na(translated)]
  return(translated)
}

# Apply translation function across all desc columns
result_chronic_and_neurological_eng <- result_chronic_and_neurological_eng %>%
  mutate(across(everything(), ~ translate_vector(.x, translation_map)))

library(dplyr)
library(tidyr)

result_chronic_and_neurological_eng <- result_chronic_and_neurological_eng %>%
  mutate(row_id = row_number())

long_df_one_hot <- result_chronic_and_neurological_eng %>%
  pivot_longer(cols = starts_with("desc"),
               names_to = "desc_col",
               values_to = "disease") %>%
  filter(!is.na(disease)) %>%
  filter(!(disease %in% c("n", "j")))

one_hot_partial <- long_df_one_hot %>%
  distinct(row_id, disease) %>%
  mutate(present = 1) %>%
  pivot_wider(names_from = disease,
              values_from = present,
              values_fill = list(present = 0))

all_rows <- tibble(row_id = 1:nrow(result_chronic_and_neurological_eng))

one_hot_full <- all_rows %>%
  left_join(one_hot_partial, by = "row_id") %>%
  replace(is.na(.), 0) %>%
  dplyr::select(-row_id)

# one_hot_full now has 1001 rows, all diseases as columns, zeros for rows with no diseases
print(dim(one_hot_full))
print(head(one_hot_full))

## Continue with plot ##

neurological_disorders <- c(
  "Migraine", "Epilepsy", "Tinnitus", "ADHD", "Multiple sclerosis",
  "Spinal canal stenosis", "Polyneuropathy", "Myasthenia gravis",
  "Cerebellar ataxia", "Chronic neuroborreliosis", "Mitochondriopathy",
  "Brain tumor", "Head tumor", "Asperger’s syndrome", "Onset dementia"
)

# Translate German descriptions to English
long_df <- long_df %>%
  mutate(desc_english = translation_map[desc]) %>%
  mutate(desc_english = ifelse(is.na(desc_english), "Unknown", desc_english)) %>%
  mutate(desc_english = factor(desc_english, levels = sort(unique(desc_english))))

dim(long_df)

# Determine the minimum observed year for later use
min_year <- min(long_df$year, na.rm = TRUE)

# Create plotting year column: use NA-year mapped to "pseudo-year" below the minimum
long_df <- long_df %>%
  mutate(year_plot = as.numeric(year),
         year_plot = ifelse(is.na(year_plot), min_year - 2, year_plot))

# Remove unknown diagnoses (e.g. those not in translation map)
long_df <- long_df %>%
  filter(desc_english != "Unknown")

# Reset factor levels to reflect actual data (post filtering)
desc_levels <- sort(unique(long_df$desc_english))
long_df$desc_english <- factor(long_df$desc_english, levels = desc_levels)

# Prepare data frame to add vertical stripes behind points
stripe_df <- data.frame(
  desc_english = desc_levels,
  xmin = seq_along(desc_levels) - 0.5,
  xmax = seq_along(desc_levels) + 0.5
)

# Calculate counts per disease for annotation above plot
count_df <- long_df %>%
  group_by(desc_english) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  mutate(xpos = as.numeric(desc_english),
         ypos = max(long_df$year_plot, na.rm = TRUE) + 3)

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

# Create colored labels: red for neurological disorders, black otherwise
label_colors <- ifelse(desc_levels %in% neurological_disorders, "red", "black")

# Add column to long_df indicating neurological disorder
long_df <- long_df %>%
  mutate(is_neurological = desc_english %in% neurological_disorders)

# Seed for reproducible jitter position
set.seed(42)

# Generate plot
p_chronic_disseases <- ggplot() +
  geom_rect(data = stripe_df,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = rep(c("ivory2", "white"), length.out = nrow(stripe_df)),
            inherit.aes = FALSE) +
  geom_jitter(data = long_df,
              aes(x = desc_english, y = year_plot, color = is_neurological),
              width = 0.3, height = 0.3, alpha = 0.6, show.legend = FALSE) +
  geom_text(data = count_df,
            aes(x = xpos, y = ypos, label = count),
            vjust = 0, size = 3, angle = 90) +
  scale_x_discrete(drop = FALSE) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"),
                     labels = c("TRUE" = "Neurological", "FALSE" = "Other"),
                     name = "Disease Type") +
  scale_y_continuous(
    breaks = c(min_year - 2, seq(min_year, max(long_df$year_plot, na.rm = TRUE))),
    labels = c("No year", as.character(seq(min_year, max(long_df$year_plot, na.rm = TRUE))))
  ) +
  labs(x = "Disease", y = "Year", title = "Diagnoses of chronic diseases over years") +
  theme_minimal(base_size = 8) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                               color = label_colors),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "black", size = 0.3, linetype = "dashed")
  ) +
  geom_hline(yintercept = seq(min_year, max(long_df$year_plot, na.rm = TRUE)),
             color = "grey80", size = 0.3, linetype = "dashed")

ggsave(paste0("p_chronic_disseases", ".svg"), p_chronic_disseases, width = 12, height = 12)
