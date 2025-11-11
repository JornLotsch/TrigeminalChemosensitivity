# Load required libraries
library(dplyr)
library(ggplot2)
library(opGMMassessment)

# Define zero-invariant log transform function
slog <- function(x, base = 10) {
  absX <- abs(x)
  s <- sign(x)
  if (base == 0) {
    return(s * log1p(absX))
  } else if (base == 2) {
    return(s * log2(absX + 1))
  } else if (base == 10) {
    return(s * log10(absX + 1))
  } else {
    factor <- log(base)
    return(s * log1p(absX) * factor)
  }
}

# Original skewed data vector
ammo_intensity <- trigeminale_daten_corrected_translated$R28

# Reflect the data (for left skewness) to positive domain
ammo_intensity_reflect <- max(ammo_intensity, na.rm = TRUE) + 1 - ammo_intensity

# Apply zero-invariant log transform to reflected data
ammo_intensity_slog <- slog(ammo_intensity_reflect, base = 10)

# Shapiro-Wilk normality tests on original and transformed data
shapiro_original <- shapiro.test(ammo_intensity)
shapiro_slog <- shapiro.test(ammo_intensity_slog)

print(shapiro_original)
print(shapiro_slog)

# Density plots for visualization
par(mfrow = c(1, 2))
plot(density(ammo_intensity, na.rm = TRUE), main = "Original Data Density")
plot(density(ammo_intensity_slog, na.rm = TRUE), main = "Reflected & Slog Transformed Density")
par(mfrow = c(1, 1))

# Gaussian Mixture Model assessment on transformed data
gmm_result <- opGMMassessment::opGMMassessment(
  ammo_intensity_slog,
  MaxModes = 4,
  MaxCores = parallel::detectCores() - 1,
  FitAlg = "DO",
  Seed = 42
)

print(gmm_result)

# Retransform GMM boundaries, means, and SDs back to original data scale
# Function to inverse slog and reflection
inv_slog <- function(y, base = 10, max_x = max(ammo_intensity, na.rm = TRUE)) {
  # Inverse of slog for base 10
  # y = s * log10(|x| + 1)
  # => |x| = 10^{|y|} - 1
  abs_y <- abs(y)
  s <- sign(y)
  x_unreflect <- s * (10 ^ abs_y - 1)
  # Undo reflection: original x = max_x + 1 - x_unreflect
  x_orig <- max_x + 1 - x_unreflect
  return(x_orig)
}

# Apply inverse transform to GMM boundaries, means, SDs
gmm_original_boundaries <- inv_slog(gmm_result$Boundaries)
gmm_original_means <- inv_slog(gmm_result$Means)
gmm_original_sds <- gmm_result$SDs # SD needs careful handling, approximate here

# Print retransformed GMM parameters in original scale
cat("Retransformed GMM Boundaries (original data scale):\n")
print(gmm_original_boundaries)

cat("Retransformed GMM Means (original data scale):\n")
print(gmm_original_means)

cat("GMM SDs (on transformed scale, approximate on original scale):\n")
print(gmm_original_sds)
