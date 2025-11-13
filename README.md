# Trigeminal sensitivity and its modulation in a general population sample

A research project analyzing trigeminal sensitivity and related factors using machine learning, statistical analysis and data visualization in R.

<img src="./wip-figure.png " width="30%">

## Project Overview

This project analyzes medical data related to trigeminal sensitivity, focusing on various risk factors and patient characteristics including:

- **Smoking history and patterns** - Duration, intensity, and temporal patterns
- **ENT surgical history** - Types, timing, and frequency of ear, nose, and throat procedures
- **Chronic disease patterns** - Comorbidities and their relationships
- **Facial pain characteristics** - Pain distribution and intensity analysis
- **COVID-19 impact** - Effects on trigeminal sensitivity
- **Statistical distributions** - Advanced statistical modeling of patient data

## Repository Structure
```
📂 .
├── 📂 R
│   ├── 📄 smoking_overview.R                   # Smoking history analysis and visualization
│   ├── 📄 ENT_surgery_overview.R               # ENT surgical history analysis
│   ├── 📄 facial_pain_overview.R               # Facial pain pattern analysis
│   ├── 📄 chronic_diseases_overview.R          # Chronic disease analysis
│   ├── 📄 nasal_breathing_overview.R           # Nasal breathing problems analysis
│   ├── 📄 covid_overview.R                     # COVID-19 impact analysis
│   ├── 📄 trigeminal_measures_distributions_correlations.R   # Distribution and correlation analysis of trigeminal measures
│   └── 📄 ammo_distribution_analysis.R         # Statistical distribution analysis
├── 📂 Python                                   # Python analysis scripts (future)
│   ├── 📄 trigeminal_measures_exploreTukey.py  # Trigeminal distrib. and transformation
├── 📂 Matlab                                   # MATLAB analysis scripts (option)
└── 📄 README.md                                # This file
```

## Key Features
### 🚬 Smoking Analysis (`smoking_overview.R`)
- Comprehensive smoking history parsing and validation
- Timeline visualization of smoking duration
- Daily cigarette consumption analysis with range handling
- Data cleaning for inconsistent date formats
- Color-coded visualization by smoking intensity

### 🏥 ENT Surgery Analysis (`ENT_surgery_overview.R`)
- Multi-language procedure description parsing (German to English)
- Temporal analysis of surgical interventions
- One-hot encoding for statistical modeling
- Comprehensive procedure categorization and translation
- Year-wise surgery distribution visualization

### 📊 Statistical Modeling (`ammo_distribution_analysis.R`)
- Advanced distribution analysis using Gaussian Mixture Models
- Zero-invariant logarithmic transformations
- Normality testing and data transformation pipelines
- Multi-modal distribution detection
- Data reflection techniques for skewed distributions

### 🎯 Data Processing Features
- Robust date parsing and validation
- Multi-format data cleaning (Excel integration)
- Missing data handling strategies
- Flexible visualization with multiple color schemes
- Reproducible analysis with set seeds

### 📈 Distribution and Correlation Analysis (`trigeminal_measures_distributions_correlations.R`)
- Transformation of trigeminal sensitivity measures
- Correlation analysis of trigeminal sensitivity measures
- Statistical relationship exploration between different assessment methods
- Multi-variable dependency analysis
- Clinical measure validation and cross-validation

### 🐍 Python Analysis (`trigeminal_measures_exploreTukey.py`)
- Tukey ladder of powers exploratory data analysis
- Advanced data transformation techniques
- Systematic evaluation of multiple trigeminal sensitivity measures
- Automated handling of data subsets (CO2 threshold analysis by time periods)
- Integration with R-generated CSV data exports



## Requirements

### R Dependencies
```r
# Data importing and wrangling
library(readxl)        # Read Excel files (.xls and .xlsx)
library(dplyr)         # Data manipulation (filter, mutate, summarise, select, etc.)
library(tidyr)         # Data tidying (pivot, unnest, etc.)
library(stringr)       # String processing (regex, pattern matching, etc.)
library(forcats)       # Manipulating factor levels
library(purrr)         # Functional programming (map, walk, etc.)
library(lubridate)     # Date and time manipulation

# Data visualization
library(ggplot2)       # Grammar of graphics for data visualization
library(scales)        # Axis labeling and scaling utilities for ggplot2
library(ggtext)        # Rich text formatting for ggplot2 axis and titles
library(viridis)       # Color palettes for publication-ready plots

# Advanced analysis (only if used in your scripts; include/exclude as needed)
library(MASS)              # Robust regression and statistical modeling
library(ComplexHeatmap)    # Advanced heatmap visualization
library(circlize)          # Circular visualizations
library(ggthemes)          # Extra themes for ggplot2
library(grid)              # Grid graphics for layout/customization
library(vcd)               # Categorical data visualization
library(cvms)              # Cross-validation tools for model selection
library(reshape2)          # Legacy data reshaping (use tidyr where possible)
library(opGMMassessment)   # Gaussian Mixture Model assessment and analysis
```

## Data Source
The project analyzes data from:
- **Source**: Bormann Trigeminale Studie Daten.xlsx
- **Primary Sheet**: Tabelle1
- **Sample Size**: 1001+ patient records
- **Data Types**: Medical histories, surgical records, smoking patterns, pain assessments

## Key Visualizations
1. **Smoking Timeline**: Horizontal bar charts showing smoking duration colored by daily cigarette consumption
2. **ENT Surgery Distribution**: Scatter plots with jittered points showing surgical procedures over time
3. **Statistical Distributions**: Density plots comparing original vs. transformed data distributions
4. **Pain Pattern Analysis**: Comprehensive facial pain mapping and intensity visualization

## Usage
Each R script can be run independently:
``` r
# Example: Run smoking analysis
source("R/smoking_overview.R")

```
## Data Processing Pipeline
1. **Data Import**: Excel file reading with sheet specification
2. **Data Validation**: Automatic detection and correction of data entry errors
3. **Feature Engineering**: Creation of derived variables and time-based features
4. **Statistical Analysis**: Advanced modeling and distribution fitting
5. **Visualization**: Multiple chart types with customizable aesthetics

## Research Applications
This codebase supports research in:
- **Medical Epidemiology**: Risk factor analysis for trigeminal sensitivity
- **Clinical Research**: Patient characteristic profiling
- **Statistical Methodology**: Advanced distribution modeling techniques
- **Data Visualization**: Medical data presentation and interpretation

## Future Development
- **Python Integration**: Planned Python scripts for machine learning applications
- **MATLAB Analysis**: Advanced signal processing and statistical modeling
- **Web Dashboard**: Interactive visualization platform
- **Database Integration**: Automated data pipeline development


## License
This project is part of ongoing research into trigeminal sensitivity and related medical factors.
