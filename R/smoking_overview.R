################################################################################
# Smoking Overview Analysis
# Description: Analysis and visualization of smoking patterns including
#              current and former smokers with temporal information
################################################################################

# Load necessary libraries
library(readxl)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(psych)

# ========================================================================== #
# 1. DATA IMPORT
# ========================================================================== #

# Read data from Excel file
trigeminale_daten_table1 <- read_excel(
  "/home/joern/Aktuell/TrigeminalSensitivity/09Originale/Bormann Trigeminale Studie Daten.xlsx",
  sheet = "Tabelle1"
)

# ========================================================================== #
# Convert selected character columns to numeric -----------------------------
# ========================================================================== #

character_vars <- names(trigeminale_daten_table1)[sapply(trigeminale_daten_table1, is.character)]
print(character_vars)
character_to_numeric_vars <- c(
  "Alter",
  "Körpergröße",
  "Riechvermögen unmittelbar nach Covid-19",
  "Wie oft Covid-19?",
  "Trinken Sie Alkohol?",
  "Wenn ich einen frischen Minzkaugummi gekaut habe, habe ich das Gefühl besser Luft durch die Nase zu bekommen"
)

trigeminale_daten_table1[character_to_numeric_vars] <- lapply(trigeminale_daten_table1[character_to_numeric_vars], function(x) as.numeric(x))


# Define smoking-related column names
col_current_smoker <- "Rauchen Sie?"
col_cigs_current <- "Wenn Ja: Wie viele Zigaretten am Tag?"
col_since_current <- "Wenn ja: Seit wann?"
col_former_smoker <- "Waren Sie je Raucher?"
col_period_former <- "Wenn ja: In welchem Zeitraum?"
col_cigs_former <- "Wenn ja: Wie viele Zigaretten am Tag?5"

# ========================================================================== #
# 2. FIX SWAPPED COLUMNS (cigarette count vs. year)
# ========================================================================== #

trigeminale_daten_table1_fixed <- trigeminale_daten_table1 %>%
  mutate(
    # Check if columns are swapped
    cigs_is_year = suppressWarnings(as.numeric(.data[[col_cigs_current]])) > 1900,
    since_is_year = suppressWarnings(as.numeric(.data[[col_since_current]])) > 1900,
    swapped = cigs_is_year & !since_is_year
  ) %>%
  mutate(
    temp = ifelse(swapped, .data[[col_cigs_current]], .data[[col_since_current]]),
    !!col_cigs_current := ifelse(swapped, .data[[col_since_current]], .data[[col_cigs_current]]),
    !!col_since_current := ifelse(swapped, temp, .data[[col_since_current]])
  ) %>%
  dplyr::select(-c(cigs_is_year, since_is_year, swapped, temp))

# ========================================================================== #
# 3. HELPER FUNCTIONS
# ========================================================================== #

# Function to convert cigarette count strings to numeric min/max daily values
convert_to_daily <- function(x) {
  min_day <- numeric(length(x))
  max_day <- numeric(length(x))

  for (i in seq_along(x)) {
    val <- str_trim(as.character(x[i]))

    if (val == "" || is.na(val)) {
      min_day[i] <- NA
      max_day[i] <- NA
      next
    }

    # Handle ranges like '4 bis 8'
    if (str_detect(val, "\\d+\\s*bis\\s*\\d+")) {
      nums <- as.numeric(str_extract_all(val, "\\d+")[[1]])
      min_day[i] <- nums[1]
      max_day[i] <- nums[2]
      next
    }

    # Handle per week
    if (str_detect(val, "/Woche")) {
      num <- as.numeric(str_extract(val, "\\d+"))
      min_day[i] <- num / 7
      max_day[i] <- num / 7
      next
    }

    # Handle per month
    if (str_detect(val, "/Monat")) {
      num <- as.numeric(str_extract(val, "\\d+"))
      min_day[i] <- num / 30
      max_day[i] <- num / 30
      next
    }

    # Single numbers
    if (str_detect(val, "^\\d+$")) {
      num <- as.numeric(val)
      min_day[i] <- num
      max_day[i] <- num
      next
    }

    # Fallback
    min_day[i] <- NA
    max_day[i] <- NA
  }

  tibble(min_per_day = min_day, max_per_day = max_day)
}

# Function to parse year ranges for former smokers
parse_former_period <- function(zeitraum_str, fallback_year) {
  zeitraum_str <- as.character(zeitraum_str)

  # Missing or "n.b." -> return fallback
  if (is.na(zeitraum_str) || zeitraum_str == "" || zeitraum_str == "n.b.") {
    return(data.frame(start = fallback_year, end = fallback_year, period = 1))
  }

  # Duration only (e.g., "10 Jahre", "16 Jahre") -> fallback
  if (grepl("^\\d+\\s+Jahre?$", zeitraum_str) ||
    (grepl("Jahre", zeitraum_str) && !grepl("-", zeitraum_str))) {
    return(data.frame(start = fallback_year, end = fallback_year, period = 1))
  }

  # "bis YYYY" format
  if (grepl("^bis\\s+\\d{4}$", zeitraum_str)) {
    end_year <- as.numeric(sub("bis\\s+", "", zeitraum_str))
    return(data.frame(start = fallback_year, end = end_year, period = 1))
  }

  # Interrupted periods like "2000-2013 + 2016-2023"
  if (grepl("\\+", zeitraum_str)) {
    periods <- strsplit(zeitraum_str, "\\s*\\+\\s*")[[1]]
    result <- data.frame()

    for (p in seq_along(periods)) {
      period_str <- trimws(periods[p])
      years <- as.numeric(unlist(strsplit(period_str, "-")))

      if (length(years) == 2 && !any(is.na(years))) {
        result <- rbind(result, data.frame(start = years[1], end = years[2], period = p))
      }
    }

    if (nrow(result) > 0) {
      return(result)
    }
  }

  # Standard "YYYY-YYYY" format
  if (grepl("\\d{4}\\s*-\\s*\\d{4}", zeitraum_str)) {
    years <- as.numeric(unlist(strsplit(zeitraum_str, "\\s*-\\s*")))

    if (length(years) == 2 && !any(is.na(years))) {
      return(data.frame(start = years[1], end = years[2], period = 1))
    }
  }

  # Fallback if nothing worked
  return(data.frame(start = fallback_year, end = fallback_year, period = 1))
}

# ========================================================================== #
# 4. PROCESS CURRENT SMOKERS
# ========================================================================== #

current_year <- as.numeric(format(Sys.Date(), "%Y"))

# Convert cigarette counts
zig_current <- convert_to_daily(trigeminale_daten_table1_fixed[[col_cigs_current]])

# Clean smoking_since: convert 2-digit years to 4-digit
smoking_since_cleaned <- as.numeric(trigeminale_daten_table1_fixed[[col_since_current]])
smoking_since_cleaned <- ifelse(
  !is.na(smoking_since_cleaned) & smoking_since_cleaned < 100,
  ifelse(smoking_since_cleaned < 50,
         2000 + smoking_since_cleaned,
         1900 + smoking_since_cleaned),
  smoking_since_cleaned
)

# Build current smokers dataframe
current_smokers <- data.frame(
  rowid = seq_len(nrow(trigeminale_daten_table1_fixed)),
  is_current = trigeminale_daten_table1_fixed[[col_current_smoker]] == "j",
  smoking_since = smoking_since_cleaned,
  min_cigs = zig_current$min_per_day,
  max_cigs = zig_current$max_per_day,
  stringsAsFactors = FALSE
)

# Filter to actual current smokers
current_smokers <- current_smokers %>%
  filter(is_current) %>%
  mutate(
    mean_cigs = rowMeans(cbind(min_cigs, max_cigs), na.rm = TRUE),
    smoker_type = "current"
  )

# ========================================================================== #
# 5. PROCESS FORMER SMOKERS
# ========================================================================== #

# Convert cigarette counts for former smokers
zig_former <- convert_to_daily(trigeminale_daten_table1_fixed[[col_cigs_former]])

# Build former smokers dataframe
former_smokers_raw <- data.frame(
  rowid = seq_len(nrow(trigeminale_daten_table1_fixed)),
  is_former = trigeminale_daten_table1_fixed[[col_former_smoker]] == "j",
  zeitraum = trigeminale_daten_table1_fixed[[col_period_former]],
  min_cigs = zig_former$min_per_day,
  max_cigs = zig_former$max_per_day,
  stringsAsFactors = FALSE
)

# Filter to actual former smokers
former_smokers_raw <- former_smokers_raw %>%
  filter(is_former) %>%
  mutate(mean_cigs = rowMeans(cbind(min_cigs, max_cigs), na.rm = TRUE))

# Determine fallback year for unspecified periods
# Use minimum year from current smokers minus 3
valid_current_years <- current_smokers$smoking_since[!is.na(current_smokers$smoking_since)]
if (length(valid_current_years) > 0) {
  fallback_year <- min(valid_current_years) - 3
} else {
  fallback_year <- current_year - 50  # Default fallback
}

# Parse former smoker periods
former_smokers_list <- list()
for (i in seq_len(nrow(former_smokers_raw))) {
  periods <- parse_former_period(former_smokers_raw$zeitraum[i], fallback_year)
  periods$rowid <- former_smokers_raw$rowid[i]
  periods$mean_cigs <- former_smokers_raw$mean_cigs[i]
  periods$smoker_type <- "former"
  former_smokers_list[[i]] <- periods
}

former_smokers <- bind_rows(former_smokers_list)

# ========================================================================== #
# 6. COMBINE AND PREPARE PLOT DATA
# ========================================================================== #

# Prepare current smokers for plotting
current_plot <- current_smokers %>%
  mutate(
    bar_start = ifelse(is.na(smoking_since), fallback_year, smoking_since),
    bar_end = current_year,
    period = 1,
    has_period = !is.na(smoking_since)
  ) %>%
  dplyr::select(rowid, bar_start, bar_end, mean_cigs, smoker_type, period, has_period)

# Prepare former smokers for plotting
former_plot <- former_smokers %>%
  mutate(
    bar_start = start,
    bar_end = end,
    has_period = (start != fallback_year | end != fallback_year)
  ) %>%
  dplyr::select(rowid, bar_start, bar_end, mean_cigs, smoker_type, period, has_period)

# Combine all smokers
all_smokers <- bind_rows(current_plot, former_plot)

# Validate years (must be reasonable)
all_smokers <- all_smokers %>%
  mutate(
    bar_start = ifelse(bar_start < 1900 | bar_start > current_year + 1, fallback_year, bar_start),
    bar_end = ifelse(bar_end < 1900 | bar_end > current_year + 1, fallback_year, bar_end)
  )

# Create unique group ID for plotting (to handle interrupted periods)
all_smokers <- all_smokers %>%
  mutate(group_id = paste(rowid, period, sep = "_"))

# ========================================================================== #
# 7. CREATE VISUALIZATION
# ========================================================================== #

# Define x-axis breaks and labels
x_breaks <- c(fallback_year, seq(1950, current_year, by = 10))
x_labels <- c("Period Not\nspecified", as.character(seq(1950, current_year, by = 10)))

# Create plot
p_smoking <- ggplot(all_smokers, aes(
  y = reorder(interaction(factor(rowid), - rowid), mean_cigs),
  color = mean_cigs
)) +
  # Bars for those with specified periods
  geom_errorbarh(
    data = all_smokers %>% filter(has_period),
    aes(xmin = bar_start, xmax = bar_end),
    height = 0.3,
    linewidth = 1
  ) +
  # Points for those without specified periods
  geom_point(
    data = all_smokers %>% filter(!has_period),
    aes(x = bar_start, y = group_id),
    size = 2,
    shape = 16
  ) +
  scale_x_continuous(
    breaks = x_breaks,
    labels = x_labels,
    limits = c(fallback_year - 1, current_year + 1)
  ) +
  scale_y_discrete(labels = function(x) {
    # Extract the first part before the dot/dash from interaction labels
    sapply(strsplit(as.character(x), "\\.|\\-"), function(y) y[1])
  }) +
  scale_color_gradient(
    low = "gold",
    high = "red",
    na.value = "grey50",
    name = "Cigarettes/Day"
  )  +
  labs(
    x = "Year",
    y = "Participant",
    title = "Smoking duration: Current and former smokers"
  ) +
  theme_minimal(base_size = 8) +
  theme(
    legend.position = "right",
    axis.text.y = element_text(size = 3),
   # axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank()
  )

# Display plot
print(p_smoking)

# Save plot
ggsave("p_smoking.svg", p_smoking, width = 12, height = 12)

# ========================================================================== #
# 8. SUMMARY STATISTICS
# ========================================================================== #


cat("\n=== Smoking Analysis Summary ===\n")

# Summary by smoker type
cat("\n--- Smoker Counts by Type ---\n")
smoker_counts <- all_smokers %>%
  group_by(smoker_type) %>%
  summarise(
    total = n_distinct(rowid),
    with_period = sum(has_period),
    without_period = sum(!has_period)
  )
print(smoker_counts)

# Overall summary
cat("\n--- Overall Summary ---\n")
cat("Total unique smokers in plot:", n_distinct(all_smokers$rowid), "\n")
cat("Total smoking periods/bars plotted:", nrow(all_smokers), "\n")
cat("Fallback year for unspecified periods:", fallback_year, "\n")

# Descriptive statistics for cigarette consumption by smoker type
cat("\n--- Cigarette Consumption Statistics ---\n")
cat("\nCURRENT SMOKERS:\n")
if (nrow(current_smokers) > 0) {
  current_cigs <- current_smokers %>%
    dplyr::select(min_cigs, max_cigs, mean_cigs) %>%
    filter(!is.na(mean_cigs))

  if (nrow(current_cigs) > 0) {
    print(psych::describe(current_cigs, na.rm = TRUE))
  } else {
    cat("No cigarette consumption data available for current smokers.\n")
  }
} else {
  cat("No current smokers in dataset.\n")
}

cat("\nALL SMOKERS COMBINED:\n")
all_cigs <- all_smokers %>%
  group_by(rowid, smoker_type) %>%
  slice(1) %>%  # Take one row per person
  ungroup() %>%
  dplyr::select(mean_cigs) %>%
  filter(!is.na(mean_cigs))

if (nrow(all_cigs) > 0) {
  print(psych::describe(all_cigs, na.rm = TRUE))
} else {
  cat("No cigarette consumption data available.\n")
}

# Smoking period length statistics
cat("\n--- Smoking Period Length (Years) ---\n")

# Calculate period lengths (excluding fallback year entries)
period_lengths <- all_smokers %>%
  filter(bar_start != fallback_year | bar_end != fallback_year) %>%
  mutate(period_length = bar_end - bar_start) %>%
  dplyr::select(rowid, smoker_type, period_length, has_period)

cat("\nCURRENT SMOKERS (period length):\n")
current_periods <- period_lengths %>%
  filter(smoker_type == "current") %>%
  dplyr::select(period_length)

if (nrow(current_periods) > 0) {
  print(psych::describe(current_periods, na.rm = TRUE))
} else {
  cat("No period length data available for current smokers.\n")
}

cat("\nFORMER SMOKERS (period length):\n")
former_periods <- period_lengths %>%
  filter(smoker_type == "former") %>%
  dplyr::select(period_length)

if (nrow(former_periods) > 0) {
  print(psych::describe(former_periods, na.rm = TRUE))
} else {
  cat("No period length data available for former smokers.\n")
}

cat("\nALL SMOKERS COMBINED (period length):\n")
if (nrow(period_lengths) > 0) {
  all_periods <- period_lengths %>%
    dplyr::select(period_length)
  print(psych::describe(all_periods, na.rm = TRUE))
} else {
  cat("No period length data available.\n")
}

cat("\n=== End of Smoking Analysis ===\n")
