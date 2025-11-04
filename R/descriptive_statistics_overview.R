################################################################################
# Trigeminal Sensitivity Analysis - Bormann Study
# Author: [YOUR NAME]
# Date: [DATE]
# Description: Analysis of trigeminal sensitivity study data, including
#              preprocessing, translation, categorization, descriptive stats,
#              and tabulations per variable category.
################################################################################

# ======================== #
# 1. Load Required Libraries
# ======================== #
library(readxl)
library(stringr)
library(dplyr)
library(tidyr)
library(psych)

# =============================== #
# 2. Data Import & Preparation
# =============================== #

# Import study data from Excel
trigeminale_daten_table1 <- read_excel(
  "/home/joern/Aktuell/TrigeminalSensitivity/09Originale/Bormann Trigeminale Studie Daten.xlsx",
  sheet = "Tabelle1"
)

# Import variable categorization and labels
variable_categories <- read_excel(
  "/home/joern/Aktuell/TrigeminalSensitivity/09Originale/Bormann Trigeminale Studie Daten_Categories.xlsx",
  sheet = "einteilung_4_R"
)

# Define category label mapping
category_labels <- c(
  "1" = "Demographics",
  "2" = "Disorders or health complaints",
  "3" = "COVID-19 history",
  "4" = "Smoking and alcohol use",
  "5" = "Subjective nasal chemosensory perception",
  "6" = "Rated olfactory function",
  "7" = "Ratings of nasal irritation and airflow",
  "8" = "Objective measurements"
)
variable_categories$Category_name <- category_labels[as.character(variable_categories$Category)]

# Variable name translations
translation_map <- c(
  # Category 1 - Demographics
  "Alter" = "Age",
  "Geschlecht" = "Gender",
  "Gewicht" = "Weight",
  "Körpergröße" = "Height",

  # Category 2 - Disorders or health complaints
  "Allergische Probleme" = "Allergic problems",
  "Infektion der oberen Atemwege" = "Upper respiratory infection",
  "Chronische Sinusitis" = "Chronic sinusitis",
  "OP im HNO-Bereich" = "Surgery in ENT region",
  "chronische Erkrankungen" = "Chronic disease",
  "neurologische Erkrankung" = "Neurological disorder",
  "Ärztliche Vorstellung wegen Problematik der Nasenatmung" = "Medical consultation for nasal breathing problems",
  "Wenn ja erfolgte eine Therapie? - Welche?" = "If yes, was therapy performed (and which one)",
  "Gesichtsschmerzen" = "Facial pain",
  "Wie oft haben Sie Gesichtsschmerzen?" = "How often do you have facial pain",
  "Wie sind die Gesichtsschmerzen?" = "What is the nature of facial pain",
  "Wie oft haben Sie im Monat Migräne?" = "How often do you have migraine per month",
  "Hat sich die Migräne in den letzten 10 Jahren verändert?" = "Has migraine changed over the last 10 years",

  # Category 3 - COVID-19 history
  "Waren Sie bereits an Covid erkrankt?" = "Have you had COVID-19?",
  "Besteht oder bestand eine Riechminderung nach Covid-19?" = "Is or was there smell reduction after COVID-19?",
  "Riechvermögen unmittelbar nach Covid-19" = "Smell ability immediately after COVID-19",
  "Wie oft Covid-19?" = "How many times have you had COVID-19?",
  "Zeitraum 1" = "Period 1",
  "Riechminderung?" = "Smell reduction 1",
  "Zeitraum 2" = "Period 2",
  "Riechminderung?2" = "Smell reduction 2",
  "Zeitraum 3" = "Period 3",
  "Riechminderung?3" = "Smell reduction 3",

  # Category 4 - Smoking and alcohol use
  "Rauchen Sie?" = "Do you smoke?",
  "Wenn Ja: Wie viele Zigaretten am Tag?" = "If yes: how many cigarettes per day?",
  "Wenn ja: Seit wann?" = "If yes: since when?",
  "Waren Sie je Raucher?" = "Have you ever smoked?",
  "Wenn ja: In welchem Zeitraum?" = "If yes: in what time period?",
  "Wenn ja: Wie viele Zigaretten am Tag?5" = "If yes: how many cigarettes per day?", # variant
  "Trinken Sie Alkohol?" = "Do you drink alcohol?",

  # Category 5 - Subjective nasal chemosensory perception (TriFunQ and related)
  "Stechende oder brennende Gerüche wie Rauch, Essig oder Nagellackentferner rufen starke Emotionen in mir hervor" =
    "Pungent or burning odors like smoke, vinegar, or nail polish remover elicit strong emotions in me",
  "Wenn ich einen frischen Minzkaugummi gekaut habe, habe ich das Gefühl besser Luft durch die Nase zu bekommen" =
    "After chewing a fresh mint gum, I feel I can breathe better through my nose",
  "Ich meide kohlensäurehaltige Getränke, weill sie in der Nase z.B. beim Aufstoßen zu sehr brennen" =
    "I avoid carbonated beverages because they burn my nose (e.g., on burping)",
  "Ich schätze meine Nasenatmung als sehr gut ein" =
    "I consider my nasal breathing to be very good",
  "Ich gehe nicht gerne in eine Sauna, weil ich heiße Luft in der Nase als brennend wahrnehme" =
    "I dislike going to saunas because I perceive hot air in my nose as burning",
  "Beim Zwiebelschneiden tränen mir die Augen stark" =
    "My eyes tear strongly when cutting onions",
  "Stechende oder brennende Gerüche lösen bei mir Husten oder Niesen aus" =
    "Pungent or burning odors cause me to cough or sneeze",
  "Ich meide brennende, stechende Gerüche (z.B. Ammoniak o. Chlor)" =
    "I avoid burning or pungent smells (e.g., ammonia or chlorine)",
  "Wenn ich etwas Beißendes oder Stechendes rieche, denke ich voller Panik an Situationen, in denen mir so etwas ähnliches passiert ist" =
    "When I smell something biting or pungent, I panic, remembering similar situations",
  "Im Winter ist mir die kalte Luft in der Nase ausgesprochen unangenehm" =
    "In winter, I find cold air in my nose extremely uncomfortable",
  "Wenn ich Meerrettich esse, finde ich das Beißen in der Nase besonders störend" =
    "When I eat horseradish, I find the burning in my nose especially bothersome",
  "Brennende oder stechende Gerüche können bei mir zu unangenehmen Empfindungen/ Schmerzen im Gesicht führen" =
    "Burning or pungent odors can cause unpleasant sensations or pain in my face",
  "Ich nehme die Kohlensäure in Getränken bewusst und intensiv wahr" =
    "I consciously and intensely perceive carbonation in drinks",
  "Wenn es um leicht kitzelnde oder stechende Gerüche geht, ist meine Nase viel empfindlicher für das Stechende und Beißende, als die Nase anderer Leute" =
    "When it comes to slightly tingling or pungent odors, my nose is much more sensitive than others'",
  "Ich benutze nur Zahnpasta mit sehr mildem Minzgeruch" =
    "I only use toothpaste with a very mild mint scent",
  "Wie oft schneiden Sie im Monat frische Zwiebeln?" =
    "How often do you cut fresh onions per month?",
  "Wenn Sie über das letzte halbe Jahr nachdenken, wie stark tränen Ihnen die Augen beim Zwiebelnschneiden?" =
    "Thinking back over the past six months, how much did your eyes tear when cutting onions?",
  "Hat sich das Tränen Ihrer Augen in den letzten 10 Jahren verändert?" =
    "Has your eye watering while cutting onions changed in the last 10 years?",

  # Category 6 - Rated olfactory function
  "Riechvermögen vor Covid" = "Smell ability before COVID-19",
  "Derzeitiges Riechvermögen" = "Current smell ability",
  "Riech- und Schmeckvermögen vermindert" = "Reduced smell and taste ability",
  "Wenn ja, wie begann das Problem" = "If yes, how did the problem start",
  "Wie hat sich das Problem verändert?" = "How has the problem changed",

  # Category 7 - Ratings of nasal irritation and airflow
  "Wie empfindlich ist Ihre Nase für stechendes/brennendes?" = "Sensitivity of the nose to stinging/burning stimuli",
  "Nasenatmung für beide Nasenlöcher" = "Nasal airflow (both nostrils)",
  "Nasenatmung für das rechte Nasenloch" = "Nasal airflow (right nostril)",
  "Nasenatmung für das linke Nasenloch" = "Nasal airflow (left nostril)",

  # Category 8 - Objective measurements
  "Lateralisierung (x/20)" = "Lateralization (x/20)",
  "AmmoLa - Intensität" = "AmmoLa intensity",
  "Identifikationstest (x/3)" = "Odor identification (x/3)",
  "CO2-Schwelle" = "CO2 threshold"
)


# Map German to English variable names, fallback to German if missing
variable_categories$Variable_english <- translation_map[variable_categories$Variable_german]
variable_categories$Variable_english[is.na(variable_categories$Variable_english)] <- variable_categories$Variable_german[is.na(variable_categories$Variable_english)]

# Create variable lists by category name
category_names <- make.names(unique(variable_categories$Category_name))
varlist <- split(variable_categories$Variable_english, variable_categories$Category_name)
names(varlist) <- category_names

# =============================== #
# 3. Descriptive statistics and tabulations per category
# =============================== #

category_vars <- split(variable_categories$Variable_german, variable_categories$Category_name)

# Descriptive statistics
desc_stats <- lapply(category_vars, function(vars) {
  existing_vars <- vars[vars %in% colnames(trigeminale_daten_table1)]
  if(length(existing_vars) > 0) {
    # Separate yes/no variables from others
    yes_no_vars <- c()
    numeric_vars <- c()

    for(var in existing_vars) {
      var_data <- trigeminale_daten_table1[[var]]
      unique_vals <- unique(na.omit(var_data))
      # Check if variable contains only j/n values OR only 0/1 values (ignoring NA)
      is_yes_no <- length(unique_vals) > 0 &&
        (all(unique_vals %in% c("j", "n")) || all(unique_vals %in% c(0, 1)))

      if(is_yes_no) {
        yes_no_vars <- c(yes_no_vars, var)
      } else {
        numeric_vars <- c(numeric_vars, var)
      }
    }

    # Get describe output for numeric variables
    if(length(numeric_vars) > 0) {
      desc_output <- describe(trigeminale_daten_table1[numeric_vars], na.rm = TRUE)
    } else {
      desc_output <- NULL
    }

    # Create describe-like output for yes/no variables
    if(length(yes_no_vars) > 0) {
      # First create the basic yes/no statistics
      yes_no_stats <- data.frame(
        n = sapply(yes_no_vars, function(v) sum(!is.na(trigeminale_daten_table1[[v]]))),
        count_yes = sapply(yes_no_vars, function(v) {
          var_data <- trigeminale_daten_table1[[v]]
          unique_vals <- unique(na.omit(var_data))
          # Count "j" for j/n variables, or 1 for 0/1 variables
          if(all(unique_vals %in% c("j", "n"))) {
            sum(var_data == "j", na.rm = TRUE)
          } else {
            sum(var_data == 1, na.rm = TRUE)
          }
        }),
        row.names = yes_no_vars,
        stringsAsFactors = FALSE
      )

      # If there's numeric output, match its structure
      if(!is.null(desc_output)) {
        # Get all column names from desc_output
        all_cols <- names(desc_output)

        # Create a complete data frame with all columns
        yes_no_output <- as.data.frame(matrix(NA, nrow = nrow(yes_no_stats), ncol = length(all_cols)))
        colnames(yes_no_output) <- all_cols
        rownames(yes_no_output) <- rownames(yes_no_stats)

        # Fill in the vars column (continuing numbering from numeric vars)
        yes_no_output$vars <- seq_along(yes_no_vars) + max(desc_output$vars)

        # Fill in the statistics we have
        yes_no_output$n <- yes_no_stats$n
        yes_no_output$max <- yes_no_stats$count_yes  # Count of "j" (or 1) in max column

        # Combine both outputs
        combined_output <- rbind(desc_output, yes_no_output)
        combined_output
      } else {
        # If no numeric vars, just add vars column to yes_no_stats
        yes_no_stats$vars <- seq_along(yes_no_vars)
        yes_no_stats <- yes_no_stats[, c("vars", "n", "count_yes")]
        yes_no_stats
      }
    } else {
      desc_output
    }
  } else {
    NULL
  }
})

# Add custom "Migraine yes" row to "Disorders or health complaints" category
if("Disorders or health complaints" %in% names(desc_stats)) {
  # Calculate migraine yes count
  migraine_yes_count <- sum(na.omit(trigeminale_daten_table1$`Wie oft haben Sie im Monat Migräne?` > 0))
  migraine_total_n <- sum(!is.na(trigeminale_daten_table1$`Wie oft haben Sie im Monat Migräne?`))

  # Create new row for Migraine yes
  migraine_row <- desc_stats$`Disorders or health complaints`[1, ]  # Copy structure
  migraine_row[] <- NA  # Set all to NA
  migraine_row$vars <- nrow(desc_stats$`Disorders or health complaints`) + 1
  migraine_row$n <- migraine_total_n
  migraine_row$max <- migraine_yes_count
  rownames(migraine_row) <- "Migraine yes"

  # Add to the data frame
  desc_stats$`Disorders or health complaints` <- rbind(
    desc_stats$`Disorders or health complaints`,
    migraine_row
  )
}

# Add custom "Alcohol yes" row to "Smoking and alcohol use" category
if("Smoking and alcohol use" %in% names(desc_stats)) {
  # Calculate alcohol yes count
  alcohol_yes_count <- sum(na.omit(trigeminale_daten_table1$`Trinken Sie Alkohol?` > 0))
  alcohol_total_n <- sum(!is.na(trigeminale_daten_table1$`Trinken Sie Alkohol?`))

  # Create new row for Alcohol yes
  alcohol_row <- desc_stats$`Smoking and alcohol use`[1, ]  # Copy structure
  alcohol_row[] <- NA  # Set all to NA
  alcohol_row$vars <- nrow(desc_stats$`Smoking and alcohol use`) + 1
  alcohol_row$n <- alcohol_total_n
  alcohol_row$max <- alcohol_yes_count
  rownames(alcohol_row) <- "Alcohol yes"

  # Add to the data frame
  desc_stats$`Smoking and alcohol use` <- rbind(
    desc_stats$`Smoking and alcohol use`,
    alcohol_row
  )
}

# Frequency tabulations
tabulations <- lapply(category_vars, function(vars) {
  existing_vars <- vars[vars %in% colnames(trigeminale_daten_table1)]
  lapply(existing_vars, function(var) table(trigeminale_daten_table1[[var]], useNA = "ifany"))
})

# =============================== #
# 4. Output
# =============================== #

# Print descriptive stats and tabulations
print(desc_stats$Demographics)
print(tabulations$Demographics)

# Print for Objective measurements as well
print(desc_stats$`Objective measurements`)
print(tabulations$`Objective measurements`)

# Etc
desc_stats$`Subjective nasal chemosensory perception`
tabulations$`Subjective nasal chemosensory perception`

desc_stats$`Disorders or health complaints`
tabulations$`Disorders or health complaints`

desc_stats$`Smoking and alcohol use`
tabulations$`Smoking and alcohol use`

desc_stats$`COVID-19 history`
tabulations$`COVID-19 history`

# Check which test categories have been applied in how many
lapply(lapply(desc_stats, "[[", "n"), max)

# Check the added questions applied only in a subset of participants
added_trigeminal_questions <- c("Wie oft schneiden Sie im Monat frische Zwiebeln?",
  "Wenn Sie über das letzte halbe Jahr nachdenken, wie stark tränen Ihnen die Augen beim Zwiebelnschneiden?",
  "Hat sich das Tränen Ihrer Augen in den letzten 10 Jahren verändert?")
max(desc_stats$`Subjective nasal chemosensory perception`$n[
  rownames(desc_stats$`Subjective nasal chemosensory perception`) %in% added_trigeminal_questions])

added_migraine_questions <- c("Wie oft haben Sie im Monat Migräne?",  "Hat sich die Migräne in den letzten 10 Jahren verändert?")
max(desc_stats$`Disorders or health complaints`$n[
  rownames(desc_stats$`Disorders or health complaints`) %in% added_migraine_questions])

