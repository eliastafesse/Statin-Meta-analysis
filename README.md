# Suboptimal Response to Statins (SoRS) Meta-analysis

This repository contains the R Markdown workflow and accompanying resources for the study:  
**"Magnitude, Predictors, and Consequences of Suboptimal Response to Statins in Primary Prevention of Cardiovascular Disease: A Systematic Review and Meta-Analysis "** by Elias Yeshitila, MD, MPhil.

## Contents
- `R_code_for_meta-analysis.R`
- `StatinMeta-analysis.csv`: Input dataset
- `renv.lock`: Captures exact R + package versions for reproducibility

## Instructions
1. Clone this repository.
2. Open `R_code_for_meta-analysis.R` in RStudio.
3. Run the renv setup:
   ```r
   if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")
   renv::restore()
   ```
4. Run the `R_code_for_meta-analysis.R` file to reproduce the analysis.

## Requirements
- All package versions are automatically managed via `renv.lock`.

## Contact
For any questions, contact **Elias Yeshitila** at ey270@cam.ac.uk.
