# Suboptimal Response to Statins (SoRS) Meta-analysis

This repository contains the R Markdown workflow and accompanying resources for the study:  
**"Prevalence and Predictors of Suboptimal Response to Statins in Primary Prevention: A Systematic Review and Meta-analysis"** by Elias Yeshitila.

## Contents
- `SoRS_meta_analysis.Rmd`
- `StatinMeta-analysis.csv`: Input dataset
- `renv.lock`: Captures exact R + package versions for reproducibility

## Instructions
1. Clone this repository.
2. Open `SoRS_meta_analysis.Rmd` in RStudio.
3. Run the renv setup:
   ```r
   if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")
   renv::restore()
   ```
4. Knit the `.Rmd` file to reproduce the analysis.

## Requirements
- R (≥ 4.2)
- All package versions are automatically managed via `renv.lock`.

## Contact
For any questions, contact **Elias Yeshitila** at ey270@cam.ac.uk.
