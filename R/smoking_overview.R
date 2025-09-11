# Load necessary libraries
library(readxl)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)

# Read data from Excel file
trigeminale_daten_table1_fixed <- read_excel(
  "/home/joern/Aktuell/TrigeminalSensitivity/09Originale/Bormann Trigeminale Studie Daten.xlsx",
  sheet = "Tabelle1"
)

# Define smoking-related column names
col_cigs <- "Wenn Ja: Wie viele Zigaretten am Tag?"
col_since <- "Wenn ja: Seit wann?"

# Detect and swap values of columns if needed based on year presence
trigeminale_daten_table1_fixed <- trigeminale_daten_table1_fixed %>%
  mutate(
# Check if the smoking count column contains a year (likely swapped with date column)
    cigs_is_year = suppressWarnings(as.numeric(.data[[col_cigs]])) > 1900,
    since_is_year = suppressWarnings(as.numeric(.data[[col_since]])) > 1900,
    swapped = cigs_is_year & !since_is_year
  ) %>%
  mutate(
    temp = ifelse(swapped, .data[[col_cigs]], .data[[col_since]]),
# Swap the values if swap condition is true
    !!col_cigs := ifelse(swapped, .data[[col_since]], .data[[col_cigs]]),
    !!col_since := ifelse(swapped, temp, .data[[col_since]])
  ) %>%
# Remove intermediate helper columns
dplyr::select(-c(cigs_is_year, since_is_year, swapped, temp))

# Function to extract all 4-digit years from strings in a vector
extract_years <- function(x) {
  lapply(x, function(s) {
    years <- stringr::str_extract_all(s, "\\d{4}")[[1]]
    as.numeric(years)
  })
}

# Extract years from 'Wenn ja: In welchem Zeitraum?'
years2 <- lapply(trigeminale_daten_table1_fixed$`Wenn ja: In welchem Zeitraum?`, extract_years)

# Determine length of year entries per record and max length for matrix
lengths_vec <- sapply(years2, function(x) length(x[[1]]))
max_length <- max(lengths_vec, na.rm = TRUE)

# Initialize matrix to hold years and fill from extracted years list
year_matrix <- matrix(NA, nrow = length(years2), ncol = max_length)
for (i in seq_along(years2)) {
  y <- years2[[i]][[1]]
  if (!is.null(y)) {
    year_matrix[i, seq_len(length(y))] <- y
  }
}

# Convert matrix to data frame and assign column names
years_df <- as.data.frame(year_matrix)
colnames(years_df) <- paste0("year", seq_len(max_length))

# Add smoking-related info from original data frame to years_df
years_df$actual_somker <- trigeminale_daten_table1_fixed$`Rauchen Sie?`
years_df$smoking_since <- trigeminale_daten_table1_fixed$`Wenn ja: Seit wann?`

# Convert smoking count strings to numeric min/max daily values
convert_to_daily <- function(x) {
  min_day <- numeric(length(x))
  max_day <- numeric(length(x))

  for (i in seq_along(x)) {
    val <- str_trim(x[i])

    # Handle empty or NA
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

    # Handle frequencies per week/month
    if (str_detect(val, "/Woche")) {
      num <- as.numeric(str_extract(val, "\\d+"))
      min_day[i] <- num / 7
      max_day[i] <- num / 7
      next
    }
    if (str_detect(val, "/Monat")) {
      num <- as.numeric(str_extract(val, "\\d+"))
      min_day[i] <- num / 30
      max_day[i] <- num / 30
      next
    }

    # Handle single numbers
    if (str_detect(val, "^\\d+$")) {
      num <- as.numeric(val)
      min_day[i] <- num
      max_day[i] <- num
      next
    }

    # Unknown format fallback
    min_day[i] <- NA
    max_day[i] <- NA
  }

  tibble(min_per_day = min_day, max_per_day = max_day)
}

# Apply conversion function to smoking count column
zig <- convert_to_daily(trigeminale_daten_table1_fixed[[col_cigs]])
years_df$min_cigarettes_per_day <- zig$min_per_day
years_df$max_cigarettes_per_day <- zig$max_per_day

# Add former smoker info and row id
years_df$former_somker <- trigeminale_daten_table1_fixed$`Waren Sie je Raucher?`
years_df$rowid <- seq_len(nrow(years_df))

# Clean 'smoking_since' converting 2-digit years to 4-digit years
years_df_cleaned <- years_df %>%
  mutate(
    smoking_since = as.numeric(smoking_since),
    smoking_since = ifelse(
      !is.na(smoking_since) & smoking_since < 100,
      ifelse(smoking_since < 50, 2000 + smoking_since, 1900 + smoking_since),
      smoking_since
    )
  )

# Set current year for plotting
current_year <- as.numeric(format(Sys.Date(), "%Y"))

# Filter active smokers, compute bar plotting values and mean cigarettes per day
smoker_df <- years_df_cleaned %>%
  filter(actual_somker == "j", !is.na(smoking_since)) %>%
  mutate(
    bar_start = smoking_since,
    bar_end = current_year,
    mean_cigs = rowMeans(
      cbind(as.numeric(min_cigarettes_per_day), as.numeric(max_cigarettes_per_day)),
      na.rm = TRUE
    ),
    person = row_number()
  )

rownames(smoker_df) <- smoker_df$rowid

# Plot smoking duration bars colored by mean daily cigarettes
p_smoking <- ggplot(smoker_df, aes(
  y = reorder(interaction(factor(rowid), - rowid), mean_cigs),
  xmin = bar_start,
  xmax = bar_end,
  color = mean_cigs
)) +
  geom_errorbarh(height = 0.3, linewidth = 2) +
#scale_color_viridis_c(option = "C", end = 0.9, na.value = "grey80") +
labs(
    x = "Year",
    y = "Person",
    color = "Cigarettes/Day",
    title = "Active Smokers' Smoking Duration Colored by Mean Cigarettes/Day"
  ) +
  theme_minimal(base_size = 8) +
  theme(legend.position.inside = TRUE, legend.position = c(.1, .5)) +
  scale_y_discrete(labels = smoker_df$rowid) +

# 1. Use other options of viridis with scale_color_viridis_c()
#  scale_color_viridis_c(option = "magma")  # Warm perceptually-uniform palette
#  scale_color_viridis_c(option = "inferno") # Dark, high-contrast palette
# scale_color_viridis_c(option = "plasma")  # Purple-yellow palette
# scale_color_viridis_c(option = "cividis") # Colorblind-friendly blue-yellow palette
#
# # 2. Use ColorBrewer palettes with scale_color_distiller() (note these have less smooth gradients)
# scale_color_distiller(palette = "YlGnBu")  # Yellow-Green-Blue palette
# scale_color_distiller(palette = "RdYlBu")  # Diverging Red-Yellow-Blue palette
# scale_color_distiller(palette = "PuRd")    # Purple-Red gradient
#
# # 4. Use base R gradients (less recommended but available)
scale_color_gradient(low = "yellow", high = "red")
# scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = median_value)
#

p_smoking

ggsave(paste0("p_smoking", ".svg"), p_smoking, width = 12, height = 12)

