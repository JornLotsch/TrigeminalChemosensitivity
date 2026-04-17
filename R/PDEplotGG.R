################################################################################
# Pareto Density Estimation Plot - ggplot2
# Description: Creates a ggplot of Pareto density estimation curves for each
#              column in a data object. Missing values are removed before
#              estimation, and each variable is plotted as a separate color.
################################################################################

# ======================== #
# Functions
# ======================== #

#' Plot Pareto Density Estimation curves
#'
#' This function converts the input to a matrix, computes Pareto density
#' estimation for each column, and returns a ggplot object with one curve per
#' variable.
#'
#' @param Data A data frame or matrix containing numeric variables.
#' @return A ggplot object showing Pareto density estimation curves.
#' @export
PDEplotGG <- function(Data) {
  Data = as.matrix(Data)
  m <- matrix(NA, nrow = 0, ncol = 3)
  require(DataVisualizations)
  for (i in 1:ncol(Data)) {
    PDE <- DataVisualizations::ParetoDensityEstimation(as.vector(na.omit(Data[,i])))
    m2 <- PDE$kernels
    m3 <- PDE$paretoDensity
    m1 <- rep(i, length(m2))
    m <- rbind(m,cbind(m1,m2,m3))
  }
  mdf <- data.frame(m)
  require(ggplot2)
  p <- ggplot(data = mdf, aes(x = m2, y = m3, colour = factor(m1))) +
    geom_line(aes(linewidth=factor(1))) + guides(linewidth = FALSE)
  return(p)
}