# Metabolic dysregulation and knee pain study

This repository contains the analytical code used for cross-sectional and longitudinal analyses of metabolic dysregulation and knee pain in a population-based cohort.

## Contents

* R scripts for statistical analysis

## Data availability

The data used in this study are not publicly available due to institutional regulations governing access to individual-level data.

## Expected data structure

The main analysis dataset is expected to include variables such as:

* `idn` — participant identifier
* `cMetS_resid` — continuous metabolic syndrome score
* `knee_pain_broad` — knee pain status (0 = no, 1 = yes)
* `sex` — sex ("male", "female")
* `age_f1` — age at baseline

Additional variables are required for extended and sensitivity analyses (e.g., cardiometabolic components, follow-up measurements).

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
```

This example is provided for illustration only and does not reproduce the study results.

## Reproducibility

The code reflects the original analytical workflow used in the study.
Full reproducibility requires access to the original dataset.

## How to run

1. Prepare input data with the required structure
2. Update file paths if needed
3. Run the analysis script

---
