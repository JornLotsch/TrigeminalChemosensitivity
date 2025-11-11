################################################################################
# Trigeminal Sensitivity Analysis - Bormann Study
# Author: Jorn Lotsch
# Description: Analysis of trigeminal sensitivity study data, including
#              preprocessing, translation, categorization, descriptive stats,
#              and tabulations per variable category.
################################################################################


# ======================== #
# 1. Load Required Libraries
# ======================== #
library(stringr)
library(dplyr)
library(readr)
library(tidyr)
library(psych)

# =============================== #
# 2. Read data
# =============================== #

trigeminale_daten_corrected_translated <- read.csv("trigeminale_daten_corrected_translated.csv", check.names = FALSE)

character_vars <- names(trigeminale_daten_corrected_translated)[sapply(trigeminale_daten_corrected_translated, is.character)]
print(character_vars)

variable_categories <- read.csv("trigeminale_daten_variable_categories.csv")

# Remove duplicated columns, keeping the first occurrence
trigeminale_daten_corrected_translated <- trigeminale_daten_corrected_translated[ , !duplicated(names(trigeminale_daten_corrected_translated))]


# =============================== #
# 3. Descriptive statistics and tabulations per category
# =============================== #

# Helper function to classify variable type
get_var_type <- function(x) {
  if (is.numeric(x)) {
    if (any(x %% 1 != 0, na.rm = TRUE)) return("interval")
    else return("ordinal")
  }
  else return("non-numeric")
}

# Helper to check if numeric vector is consecutive 0..max
is_consecutive_sequence <- function(x_num) {
  x_sorted <- sort(unique(x_num))
  identical(x_sorted, seq(from = 0, to = max(x_sorted)))
}

summarize_var <- function(x, var_type, var_name = NULL) {
  x_char <- tolower(as.character(x))
  x_char[x_char == "n.b."] <- NA_character_
  n_valid <- sum(!is.na(x_char))
  x_no_na <- x_char[!is.na(x_char)]
  unique_items <- unique(x_no_na)
  n_unique <- length(unique_items)

  var_name_lower <- tolower(var_name %||% "")

  special_separately_vars <- c(
    "period 1", "period 2", "period 3",
    "if yes: since when?",
    "if yes: in what time period?",
    "if yes: how many cigarettes per day?",
    "if yes: how many cigarettes then per day?"
  )

  special_count_vars <- c(
    "how often do you have facial pain",
    "what is the nature of facial pain",
    "odor identification (x/3)"
  )

  special_numeric_vars <- c(
  )

  if (var_name_lower %in% special_separately_vars) {
    return(tibble(
      Summary = "Summarized separately",
      var_type = "separately summarized"
    ))
  }

  # Strict yes/no detection
  values_no_na <- unique(x_no_na)
  is_strict_yes_no <- length(values_no_na) == 2 && all(values_no_na %in% c("j", "n"))
  if (is_strict_yes_no) {
    no_count <- sum(x_no_na == "n")
    yes_count <- sum(x_no_na == "j")
    return(tibble(
      Summary = sprintf("no(%d), yes(%d), n=%d", no_count, yes_count, n_valid),
      var_type = "yes/no"
    ))
  }

  if (var_name_lower %in% special_count_vars) {
    counts <- as.data.frame(table(x_no_na), stringsAsFactors = FALSE)
    items <- paste0(counts$x_no_na, "(", counts$Freq, ")")
    return(tibble(
      Summary = sprintf("%s, n=%d", paste(items, collapse = ", "), n_valid),
      var_type = "count"
    ))
  }

  if (var_name_lower %in% special_numeric_vars) {
    x_num_no_na <- as.numeric(x[!is.na(x)])
    mean_x <- mean(x_num_no_na)
    sd_x <- sd(x_num_no_na)
    median_x <- median(x_num_no_na)
    q <- quantile(x_num_no_na, probs = c(0.25, 0.75), na.rm = TRUE)
    rng <- range(x_num_no_na)
    return(tibble(
      Summary = sprintf(
        "mean=%.2f, sd=%.2f, median=%.2f, IQR=%.2f-%.2f, range=%s, n=%d",
        mean_x, sd_x, median_x, q[1], q[2], paste(rng, collapse = "-"), length(x_num_no_na)
      ),
      var_type = "interval/ordinal"
    ))
  }

  # if (var_name_lower %in% special_numeric_vars) {
  #   x_num_no_na <- as.numeric(x[!is.na(x)])
  #   mean_x <- mean(x_num_no_na)
  #   sd_x <- sd(x_num_no_na)
  #   median_x <- median(x_num_no_na)
  #   iqr_x <- IQR(x_num_no_na)
  #   rng <- range(x_num_no_na)
  #   return(tibble(
  #     Summary = sprintf("mean=%.2f, sd=%.2f, median=%.2f, iqr=%.2f, range=%s, n=%d",
  #                       mean_x, sd_x, median_x, iqr_x, paste(rng, collapse = "-"), length(x_num_no_na)),
  #     var_type = "interval/ordinal"
  #   ))
  # }

  if (is.numeric(x)) {
    x_no_na_num <- x[!is.na(x)]
    n_valid_num <- length(x_no_na_num)
    unique_nums <- unique(x_no_na_num)
    n_unique_num <- length(unique_nums)
    rng <- range(x_no_na_num)

    if (var_type == "ordinal" && n_unique_num <= 10 && is_consecutive_sequence(unique_nums)) {
      counts <- as.data.frame(table(x_no_na_num), stringsAsFactors = FALSE)
      items <- paste0(counts$x_no_na_num, "(", counts$Freq, ")")
      return(tibble(
        Summary = sprintf("%s, range=%s, n=%d",
                          paste(items, collapse = ", "), paste(rng, collapse = "-"), n_valid_num),
        var_type = "count"
      ))
    }

    if (var_type == "interval" || (var_type == "ordinal" && n_unique_num > 5)) {
      mean_x <- mean(x_no_na_num)
      sd_x <- sd(x_no_na_num)
      median_x <- median(x_no_na_num)
      q <- quantile(x_no_na_num, probs = c(0.25, 0.75), na.rm = TRUE)
      rng <- range(x_no_na_num)
      return(tibble(
        Summary = sprintf(
          "mean=%.2f, sd=%.2f, median=%.2f, IQR=%.2f-%.2f, range=%s, n=%d",
          mean_x, sd_x, median_x, q[1], q[2], paste(rng, collapse = "-"), n_valid_num
        ),
        var_type = "interval/ordinal"
      ))
    }

    # if (var_type == "interval" || (var_type == "ordinal" && n_unique_num > 5)) {
    #   mean_x <- mean(x_no_na_num)
    #   sd_x <- sd(x_no_na_num)
    #   median_x <- median(x_no_na_num)
    #   iqr_x <- IQR(x_no_na_num)
    #   return(tibble(
    #     Summary = sprintf("mean=%.2f, sd=%.2f, median=%.2f, iqr=%.2f, range=%s, n=%d",
    #                       mean_x, sd_x, median_x, iqr_x, paste(rng, collapse = "-"), n_valid_num),
    #     var_type = "interval/ordinal"
    #   ))
    # }

    if (var_type == "ordinal" && n_unique_num <= 5) {
      counts <- as.data.frame(table(x_no_na_num), stringsAsFactors = FALSE)
      items <- paste0(counts$x_no_na_num, "(", counts$Freq, ")")
      return(tibble(
        Summary = sprintf("%s, range=%s, n=%d",
                          paste(items, collapse = ", "), paste(rng, collapse = "-"), n_valid_num),
        var_type = "count"
      ))
    }
  }

  if (all(grepl("^[a-z]$", unique_items))) {
    counts <- as.data.frame(table(x_no_na), stringsAsFactors = FALSE)
    items <- paste0(counts$x_no_na, "(", counts$Freq, ")")
    return(tibble(
      Summary = sprintf("%s, n=%d", paste(items, collapse = ", "), n_valid),
      var_type = "count"
    ))
  }

  # Relaxed yes/no detection fallback
  no_values <- c("n", "nein", "0")
  yes_values <- c("j", "je", "ja", "1")
  is_yes_no_relaxed <- n_valid > 0 && all(x_no_na == "n" | x_no_na != "n")

  if (is_yes_no_relaxed) {
    no_count <- sum(x_no_na %in% no_values)
    yes_count <- sum(!(x_no_na %in% no_values))
    return(tibble(
      Summary = sprintf("no(%d), yes(%d), n=%d", no_count, yes_count, n_valid),
      var_type = "yes/no"
    ))
  }

  if (n_unique <= 5) {
    counts <- as.data.frame(table(x_no_na), stringsAsFactors = FALSE)
    items <- paste0(counts$x_no_na, "(", counts$Freq, ")")
    return(tibble(
      Summary = sprintf("%s, n=%d", paste(items, collapse = ", "), n_valid),
      var_type = "count"
    ))
  }

  tibble(Summary = NA_character_, var_type = var_type)
}

# Main routine

desc_stats <- tibble()
for (v in names(trigeminale_daten_corrected_translated[ , -1])) {  # skip ID column if any
  x <- trigeminale_daten_corrected_translated[[v]]
  vt <- get_var_type(x)
  summ_and_type <- summarize_var(x, vt, v)
  desc_stats <- bind_rows(desc_stats, tibble(
    Variable_english = v,
    Summary = summ_and_type$Summary,
    Variable_type_1 = vt,
    Variable_type = summ_and_type$var_type
  ))
}

# Join with category names
desc_stats <- desc_stats %>%
  left_join(variable_categories %>% select(Variable_english, Category_name), by = "Variable_english") %>%
  select(Category_name, Variable_english, Variable_type, Summary)

desc_stats <- desc_stats[!duplicated(desc_stats), ]

correct_facial_pain_summary <- function(desc_stats) {
  translations <- c(
    "beim passivrauchen(1)" = "during passive smoking (1)",
    "et(5)" = "once a day (5)",
    "i(2)" = "always (2)",
    "m(23)" = "monthly (23)",
    "mt(5)" = "several times a day (5)",
    "nur beim apfelessen(1)" = "only when eating apples (1)",
    "w(14)" = "weekly (14)"
  )

  desc_stats %>%
    mutate(Summary = if_else(
      Variable_english == "How often do you have facial pain",
      {
        summary_text <- Summary
        for (pattern in names(translations)) {
          summary_text <- gsub(pattern, translations[[pattern]], summary_text, fixed = TRUE)
        }
        summary_text
      },
      Summary
    ))
}

correct_facial_pain_nature_summary <- function(desc_stats) {
  translations <- c(
    "z" = "pulling",
    "s" = "stabbing",
    "sz" = "stabbing,pulling",
    "d" = "pressing",
    "sd" = "stabbing,pressing",
    "b" = "burning",
    "db" = "pressing,burning",
    "dz" = "pressing,pulling",
    "zb" = "pulling,burning",
    "sdzb" = "stabbing,pressing,pulling,burning"
  )

  desc_stats %>%
    mutate(Summary = if_else(
      Variable_english == "What is the nature of facial pain",
      {
        summary_text <- Summary
        for (pattern in names(translations)) {
          # Replace exact matches and within parentheses, case-insensitive
          summary_text <- gsub(paste0("\\b", pattern, "\\b"), translations[[pattern]], summary_text, ignore.case = TRUE)
        }
        summary_text
      },
      Summary
    ))
}

correct_types <- tibble::tribble(
  ~Variable_english, ~Correct_type,
  "Age", "interval/ordinal",
  "Gender", "nominal",
  "Weight", "interval/ordinal",
  "Height", "interval/ordinal",
  "Allergic problems", "yes/no",
  "Upper respiratory infection", "yes/no",
  "Chronic sinusitis", "yes/no",
  "Surgery in ENT region", "yes/no",
  "Chronic disease", "yes/no",
  "Neurological disorder", "yes/no",
  "Medical consultation for nasal breathing problems", "yes/no",
  "If yes, was therapy performed (and which one)", "yes/no",
  "Facial pain", "yes/no",
  "How often do you have facial pain", "nominal",
  "What is the nature of facial pain", "nominal",
  "How often do you have migraine per month", "interval/ordinal",
  "Has migraine changed over the last 10 years", "interval/ordinal",
  "Have you had COVID-19?", "yes/no",
  "Is or was there smell reduction after COVID-19?", "yes/no",
  "Smell ability before COVID-19", "interval/ordinal",
  "Smell ability immediately after COVID-19", "interval/ordinal",
  "How many times have you had COVID-19?", "count",
  "Period 1", "Year(s), separately summarized",
  "Smell reduction 1", "yes/no",
  "Period 2", "Year(s), separately summarized",
  "Smell reduction 2", "yes/no",
  "Period 3", "Year(s), separately summarized",
  "Smell reduction 3", "yes/no",
  "Current smell ability", "interval/ordinal",
  "Reduced smell and taste ability", "yes/no",
  "If yes, how did the problem start", "interval/ordinal",
  "How has the problem changed", "interval/ordinal",
  "Sensitivity of the nose to stinging/burning stimuli", "interval/ordinal",
  "Nasal airflow (both nostrils)", "interval/ordinal",
  "Nasal airflow (right nostril)", "interval/ordinal",
  "Nasal airflow (left nostril)", "interval/ordinal",
  "Do you smoke?", "yes/no",
  "If yes: how many cigarettes per day?", "Count or range, separately summarized",
  "If yes: since when?", "separately summarized",
  "Have you ever smoked?", "yes/no",
  "If yes: in what time period?", "Year(s), separately summarized",
  "If yes: how many cigarettes then per day?", "Count or range, separately summarized",
  "Do you drink alcohol?", "interval/ordinal",
  "Pungent or burning odors like smoke, vinegar, or nail polish remover elicit strong emotions in me", "interval/ordinal",
  "After chewing a fresh mint gum, I feel I can breathe better through my nose", "interval/ordinal",
  "I avoid carbonated beverages because they burn my nose (e.g., on burping)", "interval/ordinal",
  "I consider my nasal breathing to be very good", "interval/ordinal",
  "I dislike going to saunas because I perceive hot air in my nose as burning", "interval/ordinal",
  "My eyes tear strongly when cutting onions", "interval/ordinal",
  "Pungent or burning odors cause me to cough or sneeze", "interval/ordinal",
  "I avoid burning or pungent smells (e.g., ammonia or chlorine)", "interval/ordinal",
  "When I smell something biting or pungent, I panic, remembering similar situations", "interval/ordinal",
  "In winter, I find cold air in my nose extremely uncomfortable", "interval/ordinal",
  "When I eat horseradish, I find the burning in my nose especially bothersome", "interval/ordinal",
  "Burning or pungent odors can cause unpleasant sensations or pain in my face", "interval/ordinal",
  "I consciously and intensely perceive carbonation in drinks", "interval/ordinal",
  "When it comes to slightly tingling or pungent odors, my nose is much more sensitive than others'", "interval/ordinal",
  "I only use toothpaste with a very mild mint scent", "interval/ordinal",
  "How often do you cut fresh onions per month?", "interval/ordinal",
  "Thinking back over the past six months, how much did your eyes tear when cutting onions?", "interval/ordinal",
  "Has your eye watering while cutting onions changed in the last 10 years?", "interval/ordinal",
  "Lateralization (x/20)", "interval/ordinal",
  "AmmoLa intensity", "interval/ordinal",
  "Odor identification (x/3)", "count",
  "CO2 threshold", "interval/ordinal"
)

# Replace Variable_type
desc_stats <- desc_stats %>%
  left_join(correct_types, by = "Variable_english") %>%
  mutate(Variable_type = ifelse(!is.na(Correct_type), Correct_type, Variable_type)) %>%
  select(-Correct_type)


desc_stats <- correct_facial_pain_summary(desc_stats)
desc_stats <- correct_facial_pain_nature_summary(desc_stats)
desc_stats <- desc_stats %>%
  mutate(Category_name = factor(Category_name, levels = unique(variable_categories$Category_name))) %>%
  arrange(Category_name)
table(desc_stats$Category_name)

# Output as CSV
write_csv(desc_stats, "descriptive_statistics_by_category.csv")


