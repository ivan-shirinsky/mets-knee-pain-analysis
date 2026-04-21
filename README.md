# Metabolic Dysregulation and Knee Pain Study

This repository contains the analytical code used for cross-sectional and longitudinal analyses of metabolic dysregulation and knee pain in a population-based cohort.

## Study overview

This project investigates whether metabolic dysregulation, assessed using a continuous metabolic syndrome score (`cMetS_resid`), is associated with knee pain in adults from a population-based cohort.

The analysis includes both cross-sectional and longitudinal models, with a focus on:
- sex-specific associations
- interaction effects (cMetS × sex)
- potential bias due to loss to follow-up

## Analytical scope

The repository includes code for:

- Logistic regression models  
- Cross-sectional analyses  
- Longitudinal analyses (baseline → follow-up)  
- Sex-stratified analyses  
- Interaction testing (cMetS × sex)  
- Inverse probability weighting (IPW) to address non-response bias  
- Sensitivity analyses:
  - BMI-adjusted models  
  - Truncated IPW weights  

## Repository structure

- `analysis/` — R scripts for statistical analysis  
- `docs/` — study notes and protocol (optional)  

## Data availability

The individual-level data used in this study are not publicly available due to institutional and ethical restrictions.

No real data or study results are included in this repository.

## Expected data structure

The main analysis dataset is expected to include variables such as:

- `idn` — participant identifier  
- `cMetS_resid` — continuous metabolic syndrome score  
- `knee_pain_broad` — knee pain status (`0 = no`, `1 = yes`)  
- `sex` — participant sex (`"male"`, `"female"`)  
- `age_f1` — age at baseline  

Additional variables are required for extended and sensitivity analyses (e.g., cardiometabolic components, follow-up indicators, baseline covariates).

## Example (dummy) dataset

```r
set.seed(123)

df_analysis <- data.frame(
  idn = 1:10,
  cMetS_resid = rnorm(10),
  knee_pain_broad = sample(0:1, 10, replace = TRUE),
  sex = sample(c("male", "female"), 10, replace = TRUE),
  age_f1 = runif(10, 25, 45)
)

This example is provided for illustration only and does not reproduce the study results.

## Reproducibility

The code reflects the original analytical workflow used in the study. Full reproducibility requires access to the original dataset.

## How to run

1. Prepare an input dataset with the required structure  
2. Place the dataset in the working directory  
3. Adjust file paths if needed  
4. Run the main analysis script:

```r
source("analysis/main_analysis.R")


## Purpose
This repository is intended to provide transparency of the analytical workflow and to document the computational component of the study.



