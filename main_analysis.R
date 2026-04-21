# =========================================================
# Metabolic Dysregulation and Knee Pain Study
# Full analytical pipeline 
#
# Note: data are not included in this repository.
# =========================================================



suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(haven)
  library(lubridate)
  library(openxlsx)
  library(broom)
})

resid_z <- function(x, age, sex) {
  ok <- complete.cases(x, age, sex)
  out <- rep(NA_real_, length(x))
  fit <- lm(x[ok] ~ age[ok] + sex[ok])
  out[ok] <- as.numeric(scale(residuals(fit)))
  out
}

fmt_p <- function(p) {
  ifelse(is.na(p), NA, ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

fmt_or_ci <- function(or, lcl, ucl) {
  sprintf("%.2f (%.2f to %.2f)", or, lcl, ucl)
}

mean_sd <- function(x, digits = 1) {
  sprintf("%.*f ± %.*f",
          digits, mean(x, na.rm = TRUE),
          digits, sd(x, na.rm = TRUE))
}

n_perc <- function(x, level) {
  x_chr <- as.character(x)
  n <- sum(x_chr == as.character(level), na.rm = TRUE)
  d <- sum(!is.na(x_chr))
  sprintf("%d (%.1f%%)", n, 100 * n / d)
}

missing_n_perc <- function(x) {
  n <- sum(is.na(x))
  d <- length(x)
  sprintf("%d (%.1f%%)", n, 100 * n / d)
}

extract_or <- function(model, var) {
  est <- coef(model)[var]
  se  <- sqrt(diag(vcov(model)))[var]
  c(
    OR  = exp(est),
    LCL = exp(est - 1.96 * se),
    UCL = exp(est + 1.96 * se),
    p   = summary(model)$coefficients[var, "Pr(>|z|)"]
  )
}

## --- Load raw data --------------------------------------------------

df_new <- read_sav(
  "data/joints_03.04.2026.sav"
)

df_full <- read_sav(
  "data/baseline_full_cohort.sav"
)

additional_var <- read_sav(
  "data/additional_baseline_ipw.sav"
)

## --- Prepare baseline full cohort with cMetS ------------------------

df_full <- df_full %>%
  mutate(
    sex = factor(sex, levels = c(1, 2), labels = c("male", "female"))
  ) %>%
  mutate(
    z_waist_res = resid_z(waist, age_f, sex),
    z_tg_res    = resid_z(tg, age_f, sex),
    z_hdl_res   = resid_z(hdl, age_f, sex),
    z_gluc_res  = resid_z(gluc, age_f, sex),
    z_syst_res  = resid_z(syst, age_f, sex),
    cMetS_resid = z_waist_res + z_tg_res - z_hdl_res + z_gluc_res + z_syst_res
  )

## --- Lookup table ---------------------------------------------------

df_cmets <- df_full %>%
  select(
    idn,
    cMetS_resid,
    z_waist_res,
    z_tg_res,
    z_hdl_res,
    z_gluc_res,
    z_syst_res
  )

## --- Main analysis dataset ------------------------------------------

df_analysis <- df_new %>%
  rename(
    exam_date_2 = `Examination date`
  ) %>%
  mutate(
    age_f1 = age_f,
    sex = factor(sex, levels = c(1, 2), labels = c("male", "female")),
    knee_pain_broad = if_else(
      v_2 == 1 | v_1.10 == 1 | v_1.14 == 1,
      1L, 0L
    ),
    follow_up_years = as.numeric(difftime(dexam_2, dexam_1, units = "days")) / 365.25
  ) %>%
  left_join(df_cmets, by = "idn")

## --- Visit 2 cMetS dataset ------------------------------------------

df_v2 <- df_analysis %>%
  mutate(
    z_waist_res_2 = resid_z(waist_2, age_f2, sex),
    z_tg_res_2    = resid_z(tg_2, age_f2, sex),
    z_hdl_res_2   = resid_z(hdl_2, age_f2, sex),
    z_gluc_res_2  = resid_z(gluc_2, age_f2, sex),
    z_syst_res_2  = resid_z(syst_2, age_f2, sex),
    cMetS_resid_2 = z_waist_res_2 + z_tg_res_2 - z_hdl_res_2 + z_gluc_res_2 + z_syst_res_2
  )

analysis_initialized <- TRUE

cat("Session initialized successfully\n")
cat("df_new:", nrow(df_new), "rows\n")
cat("df_full:", nrow(df_full), "rows\n")
cat("df_analysis:", nrow(df_analysis), "rows\n")
cat("df_v2:", nrow(df_v2), "rows\n")
cat("Knee pain cases:", sum(df_analysis$knee_pain_broad == 1, na.rm = TRUE), "\n")

# Determining the dates of visits


library(dplyr)
library(lubridate)

## --- 1. Dates range -----------------------------------------------

cat("=== Date range ===\n")

cat("\nVisit 1 (baseline):\n")
print(range(df_analysis$dexam_1, na.rm = TRUE))

cat("\nVisit 2 (follow-up):\n")
print(range(df_analysis$dexam_2, na.rm = TRUE))

## --- 2. Year range -----------------------------------------------

cat("\n=== Year range ===\n")

cat("\nVisit 1 (baseline):\n")
print(range(year(df_analysis$dexam_1), na.rm = TRUE))

cat("\nVisit 2 (follow-up):\n")
print(range(year(df_analysis$dexam_2), na.rm = TRUE))

## --- 3. Months distribution -------------------------------------

cat("\n=== Month distribution ===\n")

cat("\nVisit 1:\n")
print(table(month(df_analysis$dexam_1)))

cat("\nVisit 2:\n")
print(table(month(df_analysis$dexam_2)))

## --- 4. Follow-up duration  ------------------------------------

df_analysis <- df_analysis %>%
  mutate(
    follow_up_years = as.numeric(difftime(dexam_2, dexam_1, units = "days")) / 365.25
  )

cat("\n=== Follow-up duration (years) ===\n")
print(summary(df_analysis$follow_up_years))


# Calculating broad vs narrow knee pain definition

stopifnot(exists("df_analysis"))

library(dplyr)

## =========================================================
## 1. Create components of knee pain definition
## =========================================================

df_knee <- df_analysis %>%
  mutate(
    knee_pain_direct  = if_else(v_2 == 1, 1L, 0L),
    knee_pain_manikin = if_else(v_1.10 == 1 | v_1.14 == 1, 1L, 0L),
    knee_pain_group = case_when(
      knee_pain_manikin == 0 & knee_pain_direct == 0 ~ "Neither",
      knee_pain_manikin == 1 & knee_pain_direct == 0 ~ "Manikin only",
      knee_pain_manikin == 0 & knee_pain_direct == 1 ~ "Direct question only",
      knee_pain_manikin == 1 & knee_pain_direct == 1 ~ "Both"
    )
  )

## =========================================================
## 2. Overall distribution (n = 249)
## =========================================================

cat("\n=== Overall distribution ===\n")

overall_tab <- df_knee %>%
  count(knee_pain_group) %>%
  mutate(percent = round(100 * n / sum(n), 1))

print(overall_tab)

## =========================================================
## 3. Distribution among cases only (knee pain = 1)
## =========================================================

cat("\n=== Among participants with knee pain ===\n")

cases_tab <- df_knee %>%
  filter(knee_pain_broad == 1) %>%
  count(knee_pain_group) %>%
  mutate(percent = round(100 * n / sum(n), 1))

print(cases_tab)

## =========================================================
## 4. Quick cross-tab (for sanity check)
## =========================================================

cat("\n=== Cross-tab (direct vs manikin) ===\n")

print(
  with(df_knee, table(knee_pain_direct, knee_pain_manikin))
)

## =========================================================
## 5. Extract numbers for manuscript
## =========================================================

cases_summary <- cases_tab %>%
  select(knee_pain_group, n) %>%
  tidyr::pivot_wider(
    names_from = knee_pain_group,
    values_from = n,
    values_fill = 0
  )

cat("\n=== Numbers for manuscript ===\n")
print(cases_summary)

# Table 1 baseline characteristics 

library(dplyr)
library(tibble)
library(openxlsx)

mean_sd <- function(x, digits = 1) {
  sprintf(
    "%.*f ± %.*f",
    digits, mean(x, na.rm = TRUE),
    digits, sd(x, na.rm = TRUE)
  )
}

n_perc <- function(x, level) {
  x_chr <- as.character(x)
  n <- sum(x_chr == as.character(level), na.rm = TRUE)
  d <- sum(!is.na(x_chr))
  sprintf("%d (%.1f%%)", n, 100 * n / d)
}

p_cont <- function(x, g) {
  keep <- !is.na(x) & !is.na(g)
  x <- x[keep]
  g <- g[keep]
  if (length(unique(g)) < 2) return("")
  tryCatch(sprintf("%.3f", t.test(x ~ g)$p.value), error = function(e) "")
}

p_cat <- function(x, g) {
  keep <- !is.na(x) & !is.na(g)
  tbl <- table(x[keep], g[keep])
  if (all(dim(tbl) > 1)) {
    sprintf("%.3f", chisq.test(tbl)$p.value)
  } else {
    ""
  }
}

df_tab1 <- df_analysis %>%
  mutate(
    sex = case_when(
      as.character(sex) %in% c("1", "male") ~ "male",
      as.character(sex) %in% c("2", "female") ~ "female",
      TRUE ~ NA_character_
    ),
    sex = factor(sex, levels = c("male", "female")),
    tg_1_mmol  = tg_1 / 88.57,
    hdl_1_mmol = hdl_1 * 0.02586,
    mets_waist_1 = case_when(
      sex == "male"   & waist_1 >= 94 ~ 1,
      sex == "female" & waist_1 >= 80 ~ 1,
      TRUE ~ 0
    ),
    mets_tg_1 = if_else(tg_1_mmol >= 1.7, 1, 0, missing = 0),
    mets_hdl_1 = case_when(
      sex == "male"   & hdl_1_mmol < 1.03 ~ 1,
      sex == "female" & hdl_1_mmol < 1.29 ~ 1,
      TRUE ~ 0
    ),
    mets_bp_1   = if_else(syst_1 >= 130 | diast_1 >= 85, 1, 0, missing = 0),
    mets_gluc_1 = if_else(gluc_1 >= 5.6, 1, 0, missing = 0),
    mets_components_1 = mets_tg_1 + mets_hdl_1 + mets_bp_1 + mets_gluc_1,
    mets_idf_1 = if_else(mets_waist_1 == 1 & mets_components_1 >= 2, 1, 0),
    knee_pain = factor(knee_pain_broad, levels = c(0, 1), labels = c("No", "Yes"))
  ) %>%
  filter(!is.na(knee_pain))

df_no  <- df_tab1 %>% filter(knee_pain == "No")
df_yes <- df_tab1 %>% filter(knee_pain == "Yes")

table1 <- tibble(
  Variable = c(
    "Participants, n",
    "Age, years",
    "Female, n (%)",
    "BMI, kg/m²",
    "Waist circumference, cm",
    "Triglycerides, mmol/L",
    "HDL-C, mmol/L",
    "Glucose, mmol/L",
    "Systolic BP, mmHg",
    "MetS by IDF, n (%)",
    "cMetS"
  ),

  `Overall` = c(
    nrow(df_tab1),
    mean_sd(df_tab1$age_f1, 1),
    n_perc(df_tab1$sex, "female"),
    mean_sd(df_tab1$bmi_1, 1),
    mean_sd(df_tab1$waist_1, 1),
    mean_sd(df_tab1$tg_1_mmol, 2),
    mean_sd(df_tab1$hdl_1_mmol, 2),
    mean_sd(df_tab1$gluc_1, 2),
    mean_sd(df_tab1$syst_1, 1),
    n_perc(df_tab1$mets_idf_1, 1),
    mean_sd(df_tab1$cMetS_resid, 2)
  ),

  `No knee pain` = c(
ppp nrow(df_no),
    mean_sd(df_no$age_f1, 1),
    n_perc(df_no$sex, "female"),
    mean_sd(df_no$bmi_1, 1),
    mean_sd(df_no$waist_1, 1),
    mean_sd(df_no$tg_1_mmol, 2),
    mean_sd(df_no$hdl_1_mmol, 2),
    mean_sd(df_no$gluc_1, 2),
    mean_sd(df_no$syst_1, 1),
    n_perc(df_no$mets_idf_1, 1),
    mean_sd(df_no$cMetS_resid, 2)
  ),

  `Knee pain` = c(
    nrow(df_yes),
    mean_sd(df_yes$age_f1, 1),
    n_perc(df_yes$sex, "female"),
    mean_sd(df_yes$bmi_1, 1),
    mean_sd(df_yes$waist_1, 1),
    mean_sd(df_yes$tg_1_mmol, 2),
    mean_sd(df_yes$hdl_1_mmol, 2),
    mean_sd(df_yes$gluc_1, 2),
    mean_sd(df_yes$syst_1, 1),
    n_perc(df_yes$mets_idf_1, 1),
    mean_sd(df_yes$cMetS_resid, 2)
  )
)

print(table1, n = Inf)
write.csv(table1, "table1.csv", row.names = FALSE, na = "")

## --- Export Table 1 to Excel ---------------------------------------

wb1 <- createWorkbook()
addWorksheet(wb1, "Table 1")

title_text_1 <- "Table 1. Baseline characteristics of participants according to knee pain status at follow-up"
note_text_1 <- paste(
"Note: Continuous variables are presented as mean ± standard deviation and categorical variables as number (%).",
"Knee pain status was missing for one participant at follow-up; this individual was excluded from the analysis.",
"MetS by IDF was defined using clinical and laboratory criteria; information on antihypertensive treatment was not included.",
"Abbreviations: BMI = body mass index; HDL-C = high-density lipoprotein cholesterol; BP = blood pressure; MetS = metabolic syndrome; cMetS = continuous metabolic syndrome risk score; IDF = International Diabetes Federation."
)

title_style <- createStyle(textDecoration = "bold", fontSize = 12, wrapText = TRUE)
header_style <- createStyle(textDecoration = "bold", border = "Bottom", wrapText = TRUE)
text_style <- createStyle(wrapText = TRUE, valign = "top")

writeData(wb1, "Table 1", title_text_1, startRow = 1, startCol = 1)
addStyle(wb1, "Table 1", title_style, rows = 1, cols = 1)

writeData(wb1, "Table 1", table1, startRow = 3, startCol = 1, headerStyle = header_style)

note_row_1 <- 3 + nrow(table1) + 3

writeData(wb1, "Table 1", note_text_1, startRow = note_row_1, startCol = 1)
addStyle(wb1, "Table 1", text_style, rows = note_row_1, cols = 1, gridExpand = TRUE)

mergeCells(wb1, "Table 1", cols = 1:4, rows = 1)
mergeCells(wb1, "Table 1", cols = 1:4, rows = note_row_1)

setColWidths(wb1, "Table 1", cols = 1:4, widths = c(28, 18, 18, 18))
freezePane(wb1, "Table 1", firstActiveRow = 4)

saveWorkbook(wb1, "Table_1.xlsx", overwrite = TRUE)

# Table 2 pain longitudinal: main analysis

stopifnot(exists("analysis_initialized"))
stopifnot(isTRUE(analysis_initialized))

## --- Packages --------------------------------------------------------
library(dplyr)
library(tidyr)
library(openxlsx)

## --- Complete-case dataset ------------------------------------------

df_main <- df_analysis %>%
  select(knee_pain_broad, cMetS_resid, sex, age_f1) %>%
  filter(
    !is.na(knee_pain_broad),
    !is.na(cMetS_resid),
    !is.na(sex),
    !is.na(age_f1)
  )

cat("N in main complete-case analysis:\n")
print(nrow(df_main))

## --- Models ----------------------------------------------------------

model_crude <- glm(
  knee_pain_broad ~ cMetS_resid,
  data   = df_main,
  family = binomial
)

model_overall <- glm(
  knee_pain_broad ~ cMetS_resid + sex + age_f1,
  data   = df_main,
  family = binomial
)

model_adj <- glm(
  knee_pain_broad ~ cMetS_resid * sex + age_f1,
  data   = df_main,
  family = binomial
)

coef_crude   <- summary(model_crude)$coefficients
coef_overall <- summary(model_overall)$coefficients
coef_adj     <- summary(model_adj)$coefficients
vc_adj       <- vcov(model_adj)

## --- Crude effect ----------------------------------------------------

beta_crude <- coef_crude["cMetS_resid", "Estimate"]
se_crude   <- coef_crude["cMetS_resid", "Std. Error"]
p_crude    <- coef_crude["cMetS_resid", "Pr(>|z|)"]

or_crude  <- exp(beta_crude)
lcl_crude <- exp(beta_crude - 1.96 * se_crude)
ucl_crude <- exp(beta_crude + 1.96 * se_crude)

## --- Overall adjusted effect ----------------------------------------

beta_overall <- coef_overall["cMetS_resid", "Estimate"]
se_overall   <- coef_overall["cMetS_resid", "Std. Error"]
p_overall    <- coef_overall["cMetS_resid", "Pr(>|z|)"]

or_overall  <- exp(beta_overall)
lcl_overall <- exp(beta_overall - 1.96 * se_overall)
ucl_overall <- exp(beta_overall + 1.96 * se_overall)

## --- Sex-specific estimates from interaction model ------------------

beta_male <- coef_adj["cMetS_resid", "Estimate"]
se_male   <- coef_adj["cMetS_resid", "Std. Error"]
p_male    <- coef_adj["cMetS_resid", "Pr(>|z|)"]

or_male  <- exp(beta_male)
lcl_male <- exp(beta_male - 1.96 * se_male)
ucl_male <- exp(beta_male + 1.96 * se_male)

interaction_name <- "cMetS_resid:sexfemale"

beta_female <- coef_adj["cMetS_resid", "Estimate"] +
  coef_adj[interaction_name, "Estimate"]

var_female <- vc_adj["cMetS_resid", "cMetS_resid"] +
  vc_adj[interaction_name, interaction_name] +
  2 * vc_adj["cMetS_resid", interaction_name]

se_female <- sqrt(var_female)
z_female  <- beta_female / se_female
p_female  <- 2 * pnorm(-abs(z_female))

or_female  <- exp(beta_female)
lcl_female <- exp(beta_female - 1.96 * se_female)
ucl_female <- exp(beta_female + 1.96 * se_female)

## --- Interaction -----------------------------------------------------

beta_interaction <- coef_adj[interaction_name, "Estimate"]
se_interaction   <- coef_adj[interaction_name, "Std. Error"]
p_interaction    <- coef_adj[interaction_name, "Pr(>|z|)"]

or_interaction  <- exp(beta_interaction)
lcl_interaction <- exp(beta_interaction - 1.96 * se_interaction)
ucl_interaction <- exp(beta_interaction + 1.96 * se_interaction)

interaction_measure_text <- paste0(
  "Measure of interaction on multiplicative scale: ratio of ORs (95% CI) = ",
  fmt_or_ci(or_interaction, lcl_interaction, ucl_interaction),
  "; P = ",
  fmt_p(p_interaction)
)

## --- Table for manuscript -------------------------------------------

tab_n <- df_analysis %>%
  filter(!is.na(sex), !is.na(knee_pain_broad)) %>%
  group_by(sex, knee_pain_broad) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(
    names_from  = knee_pain_broad,
    values_from = n,
    values_fill = 0
  ) %>%
  rename(
    no_pain = `0`,
    pain    = `1`
  )

tab_final <- tab_n %>%
  mutate(
    Sex = case_when(
      sex == "male"   ~ "Male",
      sex == "female" ~ "Female",
      TRUE            ~ as.character(sex)
    ),
    `Knee pain, n`    = pain,
    `No knee pain, n` = no_pain,
    `OR per 1-unit increase in cMetS (95% CI)` = case_when(
      sex == "male"   ~ fmt_or_ci(or_male, lcl_male, ucl_male),
      sex == "female" ~ fmt_or_ci(or_female, lcl_female, ucl_female),
      TRUE            ~ NA_character_
    ),
    `P-value` = case_when(
      sex == "male"   ~ fmt_p(p_male),
      sex == "female" ~ fmt_p(p_female),
      TRUE            ~ NA_character_
    )
  ) %>%
  select(
    Sex,
    `Knee pain, n`,
    `No knee pain, n`,
    `OR per 1-unit increase in cMetS (95% CI)`,
    `P-value`
  ) %>%
  arrange(factor(Sex, levels = c("Male", "Female")))

cat("\nTable 2. Sex-specific association between baseline cMetS and knee pain at 9-year follow-up\n\n")
print(tab_final, row.names = FALSE)

cat("\n", interaction_measure_text, "\n", sep = "")
cat("Adjusted for: age\n")
cat("Note: Odds ratios represent the association between cMetS and knee pain per 1-unit increase in cMetS. In the model without interaction, there was no clear evidence of an overall association between cMetS and knee pain.\n", "Abbreviations: cMetS = continuous metabolic syndrome score; OR = odds ratio.\n")

## --- Excel export ----------------------------------------------------

wb <- createWorkbook()
addWorksheet(wb, "Table 2")

title_text <- "Table 2. Sex-specific association between baseline cMetS and knee pain at 9-year follow-up"
adjusted_text <- "Adjusted for: age"
note_text <- "Note: Odds ratios represent the association between cMetS and knee pain per 1-unit increase in cMetS. In the model without interaction, there was no clear evidence of an overall association between cMetS and knee pain. Abbreviations: cMetS = continuous metabolic syndrome risk score; OR = odds ratio; CI = confidence interval."

title_style  <- createStyle(textDecoration = "bold", fontSize = 12, wrapText = TRUE)
header_style <- createStyle(textDecoration = "bold", border = "Bottom", wrapText = TRUE)
text_style   <- createStyle(wrapText = TRUE, valign = "top")

writeData(wb, "Table 2", title_text, startRow = 1, startCol = 1)
addStyle(wb, "Table 2", title_style, rows = 1, cols = 1)

writeData(wb, "Table 2", tab_final, startRow = 3, startCol = 1, headerStyle = header_style)

note_row <- 3 + nrow(tab_final) + 3

writeData(wb, "Table 2", interaction_measure_text, startRow = note_row,     startCol = 1)
writeData(wb, "Table 2", adjusted_text,            startRow = note_row + 1, startCol = 1)
writeData(wb, "Table 2", note_text,                startRow = note_row + 2, startCol = 1)

addStyle(wb, "Table 2", text_style, rows = note_row:(note_row + 2), cols = 1, gridExpand = TRUE)

mergeCells(wb, "Table 2", cols = 1:5, rows = 1)
mergeCells(wb, "Table 2", cols = 1:5, rows = note_row)
mergeCells(wb, "Table 2", cols = 1:5, rows = note_row + 1)
mergeCells(wb, "Table 2", cols = 1:5, rows = note_row + 2)

setColWidths(wb, "Table 2", cols = 1:5, widths = c(12, 14, 16, 35, 12))
freezePane(wb, "Table 2", firstActiveRow = 4)

saveWorkbook(wb, "Table_2.xlsx", overwrite = TRUE)

cat("\nExcel file saved: Table_2.xlsx\n")

# Table 3 cross-sectional: main analysis

stopifnot(exists("analysis_initialized"))
stopifnot(isTRUE(analysis_initialized))
stopifnot(exists("fmt_or_ci"))
stopifnot(exists("fmt_p"))

library(dplyr)
library(tidyr)
library(openxlsx)

cat("\n=== Cross-sectional analysis at 9-year follow-up ===\n")

model_cross <- glm(
  knee_pain_broad ~ cMetS_resid_2 * sex + age_f2,
  data   = df_v2,
  family = binomial
)

coef_cross <- summary(model_cross)$coefficients
vc_cross   <- vcov(model_cross)

interaction_name_cross <- "cMetS_resid_2:sexfemale"

## Male effect
beta_male_cross <- coef_cross["cMetS_resid_2", "Estimate"]
se_male_cross   <- coef_cross["cMetS_resid_2", "Std. Error"]
p_male_cross    <- coef_cross["cMetS_resid_2", "Pr(>|z|)"]

or_male_cross  <- exp(beta_male_cross)
lcl_male_cross <- exp(beta_male_cross - 1.96 * se_male_cross)
ucl_male_cross <- exp(beta_male_cross + 1.96 * se_male_cross)

## Female effect
beta_female_cross <- coef_cross["cMetS_resid_2", "Estimate"] +
  coef_cross[interaction_name_cross, "Estimate"]

var_female_cross <- vc_cross["cMetS_resid_2", "cMetS_resid_2"] +
  vc_cross[interaction_name_cross, interaction_name_cross] +
  2 * vc_cross["cMetS_resid_2", interaction_name_cross]

se_female_cross <- sqrt(var_female_cross)
z_female_cross  <- beta_female_cross / se_female_cross
p_female_cross  <- 2 * pnorm(-abs(z_female_cross))

or_female_cross  <- exp(beta_female_cross)
lcl_female_cross <- exp(beta_female_cross - 1.96 * se_female_cross)
ucl_female_cross <- exp(beta_female_cross + 1.96 * se_female_cross)

## Interaction
beta_interaction_cross <- coef_cross[interaction_name_cross, "Estimate"]
se_interaction_cross   <- coef_cross[interaction_name_cross, "Std. Error"]
p_interaction_cross    <- coef_cross[interaction_name_cross, "Pr(>|z|)"]

or_interaction_cross  <- exp(beta_interaction_cross)
lcl_interaction_cross <- exp(beta_interaction_cross - 1.96 * se_interaction_cross)
ucl_interaction_cross <- exp(beta_interaction_cross + 1.96 * se_interaction_cross)

interaction_measure_text_cross <- paste0(
  "Measure of interaction on multiplicative scale: ratio of ORs (95% CI) = ",
  fmt_or_ci(or_interaction_cross, lcl_interaction_cross, ucl_interaction_cross),
  "; P = ",
  fmt_p(p_interaction_cross)
)

## Counts for table
tab_n_cross <- df_v2 %>%
  filter(!is.na(sex), !is.na(knee_pain_broad)) %>%
  group_by(sex, knee_pain_broad) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(
    names_from  = knee_pain_broad,
    values_from = n,
    values_fill = 0
  ) %>%
  rename(
    no_pain = `0`,
    pain    = `1`
  )

tab_final_cross <- tab_n_cross %>%
  mutate(
    Sex = case_when(
      sex == "male"   ~ "Male",
      sex == "female" ~ "Female",
      TRUE            ~ as.character(sex)
    ),
    `Knee pain, n`    = pain,
    `No knee pain, n` = no_pain,
    `OR per 1-unit increase in cMetS (95% CI)` = case_when(
      sex == "male"   ~ fmt_or_ci(or_male_cross, lcl_male_cross, ucl_male_cross),
      sex == "female" ~ fmt_or_ci(or_female_cross, lcl_female_cross, ucl_female_cross),
      TRUE            ~ NA_character_
    ),
    `P-value` = case_when(
      sex == "male"   ~ fmt_p(p_male_cross),
      sex == "female" ~ fmt_p(p_female_cross),
      TRUE            ~ NA_character_
    )
  ) %>%
  select(
    Sex,
    `Knee pain, n`,
    `No knee pain, n`,
    `OR per 1-unit increase in cMetS (95% CI)`,
    `P-value`
  ) %>%
  arrange(factor(Sex, levels = c("Male", "Female")))

cat("\nTable 3. Sex-specific cross-sectional association between cMetS and knee pain at 9-year follow-up\n\n")
print(tab_final_cross, row.names = FALSE)

cat("\n", interaction_measure_text_cross, "\n", sep = "")
cat("Adjusted for: age\n")
cat("Note: Odds ratios represent the cross-sectional association between cMetS and knee pain per 1-unit increase in cMetS at 9-year follow-up. Abbreviations: cMetS = continuous metabolic syndrome score; OR = odds ratio; CI = confidence interval.\n")

## Excel export
wb3 <- createWorkbook()
addWorksheet(wb3, "Table 3")

title_text_cross <- "Table 3. Sex-specific cross-sectional association between cMetS and knee pain at 9-year follow-up"
adjusted_text_cross <- "Adjusted for: age"
note_text_cross <- "Note: Odds ratios represent the cross-sectional association between cMetS and knee pain per 1-unit increase in cMetS at 9-year follow-up. Abbreviations: cMetS = continuous metabolic syndrome score; OR = odds ratio; CI = confidence interval."

title_style  <- createStyle(textDecoration = "bold", fontSize = 12, wrapText = TRUE)
header_style <- createStyle(textDecoration = "bold", border = "Bottom", wrapText = TRUE)
text_style   <- createStyle(wrapText = TRUE, valign = "top")

writeData(wb3, "Table 3", title_text_cross, startRow = 1, startCol = 1)
addStyle(wb3, "Table 3", title_style, rows = 1, cols = 1)

writeData(wb3, "Table 3", tab_final_cross, startRow = 3, startCol = 1, headerStyle = header_style)

note_row_cross <- 3 + nrow(tab_final_cross) + 3

writeData(wb3, "Table 3", interaction_measure_text_cross, startRow = note_row_cross,     startCol = 1)
writeData(wb3, "Table 3", adjusted_text_cross,            startRow = note_row_cross + 1, startCol = 1)
writeData(wb3, "Table 3", note_text_cross,                startRow = note_row_cross + 2, startCol = 1)

addStyle(wb3, "Table 3", text_style,
         rows = note_row_cross:(note_row_cross + 2),
         cols = 1, gridExpand = TRUE)

mergeCells(wb3, "Table 3", cols = 1:5, rows = 1)
mergeCells(wb3, "Table 3", cols = 1:5, rows = note_row_cross)
mergeCells(wb3, "Table 3", cols = 1:5, rows = note_row_cross + 1)
mergeCells(wb3, "Table 3", cols = 1:5, rows = note_row_cross + 2)

setColWidths(wb3, "Table 3", cols = 1:5, widths = c(12, 14, 16, 35, 12))
freezePane(wb3, "Table 3", firstActiveRow = 4)

saveWorkbook(wb3, "Table_3.xlsx", overwrite = TRUE)

cat("\nExcel file saved: Table_3.xlsx\n")

# Supplementary Table S1: longitudinal models full set

stopifnot(exists("analysis_initialized"))
stopifnot(isTRUE(analysis_initialized))

library(dplyr)
library(openxlsx)

## --- Complete-case dataset ------------------------------------------

df_main <- df_analysis %>%
  select(knee_pain_broad, cMetS_resid, sex, age_f1) %>%
  filter(
    !is.na(knee_pain_broad),
    !is.na(cMetS_resid),
    !is.na(sex),
    !is.na(age_f1)
  )

cat("N in main complete-case analysis:\n")
print(nrow(df_main))

## --- Models ----------------------------------------------------------

# Model 1: no interaction
m1_crude <- glm(
  knee_pain_broad ~ cMetS_resid,
  data = df_main,
  family = binomial
)

m1_adj <- glm(
  knee_pain_broad ~ cMetS_resid + sex + age_f1,
  data = df_main,
  family = binomial
)

# Model 2: interaction
m2_crude <- glm(
  knee_pain_broad ~ cMetS_resid * sex,
  data = df_main,
  family = binomial
)

m2_adj <- glm(
  knee_pain_broad ~ cMetS_resid * sex + age_f1,
  data = df_main,
  family = binomial
)

## --- Helper extraction ----------------------------------------------

extract_main <- function(model, var) {
  cf <- summary(model)$coefficients
  est <- cf[var, "Estimate"]
  se  <- cf[var, "Std. Error"]
  p   <- cf[var, "Pr(>|z|)"]
  c(
    OR  = exp(est),
    LCL = exp(est - 1.96 * se),
    UCL = exp(est + 1.96 * se),
    p   = p
  )
}

extract_int <- function(model) {
  cf <- summary(model)$coefficients
  vc <- vcov(model)

  main <- "cMetS_resid"
  int  <- "cMetS_resid:sexfemale"

  beta_m <- cf[main, "Estimate"]
  se_m   <- cf[main, "Std. Error"]
  p_m    <- cf[main, "Pr(>|z|)"]

  beta_f <- beta_m + cf[int, "Estimate"]
  var_f  <- vc[main, main] + vc[int, int] + 2 * vc[main, int]
  se_f   <- sqrt(var_f)
  z_f    <- beta_f / se_f
  p_f    <- 2 * pnorm(-abs(z_f))

  beta_i <- cf[int, "Estimate"]
  se_i   <- cf[int, "Std. Error"]
  p_i    <- cf[int, "Pr(>|z|)"]

  list(
    male = c(
      OR  = exp(beta_m),
      LCL = exp(beta_m - 1.96 * se_m),
      UCL = exp(beta_m + 1.96 * se_m),
      p   = p_m
    ),
    female = c(
      OR  = exp(beta_f),
      LCL = exp(beta_f - 1.96 * se_f),
      UCL = exp(beta_f + 1.96 * se_f),
      p   = p_f
    ),
    interaction = c(
      OR  = exp(beta_i),
      LCL = exp(beta_i - 1.96 * se_i),
      UCL = exp(beta_i + 1.96 * se_i),
      p   = p_i
    )
  )
}

## --- Extract estimates ----------------------------------------------

m1c <- extract_main(m1_crude, "cMetS_resid")
m1a <- extract_main(m1_adj,   "cMetS_resid")

m2c <- extract_int(m2_crude)
m2a <- extract_int(m2_adj)

## --- Table -----------------------------------------------------------

tab_s2 <- bind_rows(
  data.frame(
    Model = "Model 1: no interaction",
    Adjustment = "Crude",
    Effect = "cMetS",
    `OR (95% CI)` = fmt_or_ci(m1c["OR"], m1c["LCL"], m1c["UCL"]),
    `P-value` = fmt_p(m1c["p"]),
    check.names = FALSE
  ),
  data.frame(
    Model = "Model 1: no interaction",
    Adjustment = "Adjusted",
    Effect = "cMetS",
    `OR (95% CI)` = fmt_or_ci(m1a["OR"], m1a["LCL"], m1a["UCL"]),
    `P-value` = fmt_p(m1a["p"]),
    check.names = FALSE
  ),
  data.frame(
    Model = rep("Model 2: interaction", 3),
    Adjustment = rep("Crude", 3),
    Effect = c("cMetS (males)", "cMetS (females)", "cMetS × sex"),
    `OR (95% CI)` = c(
      fmt_or_ci(m2c$male["OR"],        m2c$male["LCL"],        m2c$male["UCL"]),
      fmt_or_ci(m2c$female["OR"],      m2c$female["LCL"],      m2c$female["UCL"]),
      fmt_or_ci(m2c$interaction["OR"], m2c$interaction["LCL"], m2c$interaction["UCL"])
    ),
    `P-value` = c(
      fmt_p(m2c$male["p"]),
      fmt_p(m2c$female["p"]),
      fmt_p(m2c$interaction["p"])
    ),
    check.names = FALSE
  ),
  data.frame(
    Model = rep("Model 2: interaction", 3),
    Adjustment = rep("Adjusted", 3),
    Effect = c("cMetS (males)", "cMetS (females)", "cMetS × sex"),
    `OR (95% CI)` = c(
      fmt_or_ci(m2a$male["OR"],        m2a$male["LCL"],        m2a$male["UCL"]),
      fmt_or_ci(m2a$female["OR"],      m2a$female["LCL"],      m2a$female["UCL"]),
      fmt_or_ci(m2a$interaction["OR"], m2a$interaction["LCL"], m2a$interaction["UCL"])
    ),
    `P-value` = c(
      fmt_p(m2a$male["p"]),
      fmt_p(m2a$female["p"]),
      fmt_p(m2a$interaction["p"])
    ),
    check.names = FALSE
  )
)

cat("\nSupplementary Table S1\n\n")
print(tab_s2, row.names = FALSE)

## --- Interaction text -----------------------------------------------

interaction_text_crude <- paste0(
  "Model 2 crude: Measure of interaction on multiplicative scale: ratio of ORs (95% CI) = ",
  fmt_or_ci(m2c$interaction["OR"], m2c$interaction["LCL"], m2c$interaction["UCL"]),
  "; P = ", fmt_p(m2c$interaction["p"])
)

interaction_text_adj <- paste0(
  "Model 2 adjusted: Measure of interaction on multiplicative scale: ratio of ORs (95% CI) = ",
  fmt_or_ci(m2a$interaction["OR"], m2a$interaction["LCL"], m2a$interaction["UCL"]),
  "; P = ", fmt_p(m2a$interaction["p"])
)

cat("\n", interaction_text_crude, "\n", sep = "")
cat(interaction_text_adj, "\n", sep = "")

## --- Excel export ----------------------------------------------------

wb <- createWorkbook()
addWorksheet(wb, "Table S1")

title_text <- "Supplementary Table S1. Longitudinal associations between baseline cMetS and knee pain at 9-year follow-up: full set of crude and adjusted logistic regression models with and without cMetS-by-sex interaction"
adjusted_text <- "Model 1: without interaction. Model 2: with a cMetS-by-sex interaction term."
note_text <- paste(
  "Note: Crude models are unadjusted.",
  "Adjusted models include age and sex where appropriate.",
  "For interaction models, sex-specific effects are shown for males and females together with the multiplicative interaction term.",
  "Abbreviations: cMetS = continuous metabolic syndrome score; OR = odds ratio; CI = confidence interval."
)

title_style  <- createStyle(textDecoration = "bold", fontSize = 12, wrapText = TRUE)
header_style <- createStyle(textDecoration = "bold", border = "Bottom", wrapText = TRUE)
text_style   <- createStyle(wrapText = TRUE, valign = "top")

writeData(wb, "Table S1", title_text, startRow = 1, startCol = 1)
addStyle(wb, "Table S1", title_style, rows = 1, cols = 1)

writeData(wb, "Table S1", tab_s2, startRow = 3, startCol = 1, headerStyle = header_style)

note_row <- 3 + nrow(tab_s2) + 3

writeData(wb, "Table S1", adjusted_text,         startRow = note_row,     startCol = 1)
writeData(wb, "Table S1", interaction_text_crude,startRow = note_row + 1, startCol = 1)
writeData(wb, "Table S1", interaction_text_adj,  startRow = note_row + 2, startCol = 1)
writeData(wb, "Table S1", note_text,             startRow = note_row + 3, startCol = 1)

addStyle(wb, "Table S1", text_style, rows = note_row:(note_row + 3), cols = 1, gridExpand = TRUE)

mergeCells(wb, "Table S1", cols = 1:5, rows = 1)
mergeCells(wb, "Table S1", cols = 1:5, rows = note_row)
mergeCells(wb, "Table S1", cols = 1:5, rows = note_row + 1)
mergeCells(wb, "Table S1", cols = 1:5, rows = note_row + 2)
mergeCells(wb, "Table S1", cols = 1:5, rows = note_row + 3)

setColWidths(wb, "Table S1", cols = 1:5, widths = c(24, 12, 18, 22, 12))
freezePane(wb, "Table S1", firstActiveRow = 4)

saveWorkbook(wb, "Table_S1_longitudinal_associations_full_set.xlsx", overwrite = TRUE)

cat("\nExcel file saved: Table_S1_longitudinal_models.xlsx\n")


# Supplementary Table S2: longitudinal individual components

stopifnot(exists("analysis_initialized"))
stopifnot(isTRUE(analysis_initialized))

library(dplyr)
library(openxlsx)

components1 <- c("z_waist_res", "z_tg_res", "z_hdl_res", "z_gluc_res", "z_syst_res")

labels1 <- c(
  z_waist_res = "Waist circumference (baseline residual z)",
  z_tg_res    = "Triglycerides (baseline residual z)",
  z_hdl_res   = "HDL-C (baseline residual z)",
  z_gluc_res  = "Glucose (baseline residual z)",
  z_syst_res  = "Systolic BP (baseline residual z)"
)

results_long_comp <- list()

for (comp in components1) {
  df_comp <- df_analysis %>%
    select(knee_pain_broad, all_of(comp), sex, age_f1) %>%
    filter(
      !is.na(knee_pain_broad),
      !is.na(.data[[comp]]),
      !is.na(sex),
      !is.na(age_f1)
    )

  form <- as.formula(paste("knee_pain_broad ~", comp, "* sex + age_f1"))
  fit  <- glm(form, data = df_comp, family = binomial)

  coef_tab <- summary(fit)$coefficients
  vc       <- vcov(fit)

  main_name <- comp
  int_name  <- paste0(comp, ":sexfemale")
  if (!int_name %in% rownames(coef_tab)) {
    int_name <- paste0("sexfemale:", comp)
  }

  ## Male
  est_m <- coef_tab[main_name, "Estimate"]
  se_m  <- coef_tab[main_name, "Std. Error"]
  p_m   <- coef_tab[main_name, "Pr(>|z|)"]

  or_m  <- exp(est_m)
  lcl_m <- exp(est_m - 1.96 * se_m)
  ucl_m <- exp(est_m + 1.96 * se_m)

  ## Female
  est_f <- coef_tab[main_name, "Estimate"] + coef_tab[int_name, "Estimate"]
  var_f <- vc[main_name, main_name] +
    vc[int_name, int_name] +
    2 * vc[main_name, int_name]
  se_f <- sqrt(var_f)
  z_f  <- est_f / se_f
  p_f  <- 2 * pnorm(-abs(z_f))

  or_f  <- exp(est_f)
  lcl_f <- exp(est_f - 1.96 * se_f)
  ucl_f <- exp(est_f + 1.96 * se_f)

  ## Interaction
  est_i <- coef_tab[int_name, "Estimate"]
  se_i  <- coef_tab[int_name, "Std. Error"]
  p_i   <- coef_tab[int_name, "Pr(>|z|)"]

  or_i  <- exp(est_i)
  lcl_i <- exp(est_i - 1.96 * se_i)
  ucl_i <- exp(est_i + 1.96 * se_i)

  results_long_comp[[comp]] <- data.frame(
    Exposure = labels1[[comp]],
    N = nrow(df_comp),
    `Male OR (95% CI)` = fmt_or_ci(or_m, lcl_m, ucl_m),
    `Male P-value` = fmt_p(p_m),
    `Female OR (95% CI)` = fmt_or_ci(or_f, lcl_f, ucl_f),
    `Female P-value` = fmt_p(p_f),
    `Interaction ratio of ORs (95% CI)` = fmt_or_ci(or_i, lcl_i, ucl_i),
    `Interaction P-value` = fmt_p(p_i),
    check.names = FALSE
  )
}

tab_components_clean <- bind_rows(results_long_comp)

cat("\nSupplementary Table S2. Sex-specific longitudinal associations of individual baseline components of cMetS with knee pain at 9-year follow-up\n\n")
print(tab_components_clean, row.names = FALSE)

## --- Excel export ----------------------------------------------------

wb <- createWorkbook()
addWorksheet(wb, "Table S2")

title_text <- "Supplementary Table S2. Sex-specific longitudinal associations of individual baseline components of cMetS with knee pain at 9-year follow-up"
adjusted_text <- "Adjusted for: age; sex-specific estimates derived from models including a component-by-sex interaction term"
note_text <- paste(
  "Note: Odds ratios are reported per 1 SD increase in each baseline residualized metabolic component.",
  "The interaction term is presented as the ratio of odds ratios on the multiplicative scale.",
  "Abbreviations: OR = odds ratio; CI = confidence interval; HDL-C = high-density lipoprotein cholesterol; BP = blood pressure."
)

title_style  <- createStyle(textDecoration = "bold", fontSize = 12, wrapText = TRUE)
header_style <- createStyle(textDecoration = "bold", border = "Bottom", wrapText = TRUE)
text_style   <- createStyle(wrapText = TRUE, valign = "top")

writeData(wb, "Table S2", title_text, startRow = 1, startCol = 1)
addStyle(wb, "Table S2", title_style, rows = 1, cols = 1)

writeData(
  wb,
  "Table S2",
  tab_components_clean,
  startRow = 3,
  startCol = 1,
  headerStyle = header_style
)

note_row <- 3 + nrow(tab_components_clean) + 3

writeData(wb, "Table S2", adjusted_text, startRow = note_row,     startCol = 1)
writeData(wb, "Table S2", note_text,     startRow = note_row + 1, startCol = 1)

addStyle(wb, "Table S2", text_style, rows = note_row:(note_row + 1), cols = 1, gridExpand = TRUE)

mergeCells(wb, "Table S2", cols = 1:8, rows = 1)
mergeCells(wb, "Table S2", cols = 1:8, rows = note_row)
mergeCells(wb, "Table S2", cols = 1:8, rows = note_row + 1)

setColWidths(wb, "Table S2", cols = 1:8, widths = c(38, 8, 18, 10, 18, 10, 24, 12))
setRowHeights(wb, "Table S2", rows = c(1, note_row:(note_row + 1)), heights = 70)
freezePane(wb, "Table S2", firstActiveRow = 4)

saveWorkbook(wb, "Table_S2_longitudinal_components.xlsx", overwrite = TRUE)

cat("\nExcel file saved: Table_S2_longitudinal_components.xlsx\n")

* Supplementary Table S3: responders vs non-responders

stopifnot(exists("analysis_initialized"))
stopifnot(isTRUE(analysis_initialized))
stopifnot(exists("df_new"))
stopifnot(exists("df_full"))
stopifnot(exists("additional_var"))

library(dplyr)
library(tibble)
library(openxlsx)

## --- 0. Helper functions ---------------------------------------------

mean_sd <- function(x, digits = 1) {
  x <- as.numeric(x)
  if (all(is.na(x))) return("NA")
  sprintf(
    "%.*f ± %.*f",
    digits, mean(x, na.rm = TRUE),
    digits, sd(x, na.rm = TRUE)
  )
}

missing_n_perc <- function(x) {
  n <- sum(is.na(x))
  d <- length(x)
  if (d == 0) return("0 (NA)")
  sprintf("%d (%.1f%%)", n, 100 * n / d)
}

n_perc <- function(x, level) {
  x_chr <- as.character(x)
  n <- sum(x_chr == as.character(level), na.rm = TRUE)
  d <- sum(!is.na(x_chr))
  if (d == 0) return("0 (NA)")
  sprintf("%d (%.1f%%)", n, 100 * n / d)
}

smd_cont <- function(x1, x0) {
  x1 <- as.numeric(x1)
  x0 <- as.numeric(x0)
  x1 <- x1[!is.na(x1)]
  x0 <- x0[!is.na(x0)]
  if (length(x1) < 2 || length(x0) < 2) return("")
  s1 <- sd(x1)
  s0 <- sd(x0)
  sp <- sqrt((s1^2 + s0^2) / 2)
  if (is.na(sp) || sp == 0) return("")
  sprintf("%.3f", abs((mean(x1) - mean(x0)) / sp))
}

smd_binary <- function(x1, x0, level) {
  p1 <- mean(as.character(x1) == as.character(level), na.rm = TRUE)
  p0 <- mean(as.character(x0) == as.character(level), na.rm = TRUE)
  p  <- (p1 + p0) / 2
  denom <- sqrt(p * (1 - p))
  if (is.na(denom) || denom == 0) return("")
  sprintf("%.3f", abs((p1 - p0) / denom))
}

smd_multicat <- function(x1, x0) {
  x1 <- x1[!is.na(x1)]
  x0 <- x0[!is.na(x0)]
  if (length(x1) == 0 || length(x0) == 0) return("")
  levs <- union(levels(factor(x1)), levels(factor(x0)))
  p1 <- prop.table(table(factor(x1, levels = levs)))
  p0 <- prop.table(table(factor(x0, levels = levs)))
  smd <- sqrt(sum((as.numeric(p1) - as.numeric(p0))^2))
  sprintf("%.3f", smd)
}

make_cont_rows <- function(var, label, digits = 1, data_resp, data_non) {
  tibble(
    Variable = c(label, "  Missing"),
    `Completed follow-up` = c(
      mean_sd(data_resp[[var]], digits),
      missing_n_perc(data_resp[[var]])
    ),
    `Not assessed at follow-up` = c(
      mean_sd(data_non[[var]], digits),
      missing_n_perc(data_non[[var]])
    ),
    SMD = c(
      smd_cont(data_resp[[var]], data_non[[var]]),
      ""
    )
  )
}

make_binary_rows <- function(var, label, level, data_resp, data_non) {
  tibble(
    Variable = c(label, "  Missing"),
    `Completed follow-up` = c(
      n_perc(data_resp[[var]], level),
      missing_n_perc(data_resp[[var]])
    ),
    `Not assessed at follow-up` = c(
      n_perc(data_non[[var]], level),
      missing_n_perc(data_non[[var]])
    ),
    SMD = c(
      smd_binary(data_resp[[var]], data_non[[var]], level),
      ""
    )
  )
}

make_cat_rows <- function(varname, label, levels_order, data_resp, data_non) {
  tibble(
    Variable = c(label, paste0("  ", levels_order), "  Missing"),
    `Completed follow-up` = c(
      "",
      sapply(levels_order, function(l) n_perc(data_resp[[varname]], l)),
      missing_n_perc(data_resp[[varname]])
    ),
    `Not assessed at follow-up` = c(
      "",
      sapply(levels_order, function(l) n_perc(data_non[[varname]], l)),
      missing_n_perc(data_non[[varname]])
    ),
    SMD = c(
      smd_multicat(data_resp[[varname]], data_non[[varname]]),
      rep("", length(levels_order)),
      ""
    )
  )
}

## --- 1. Additional baseline variables --------------------------------

additional_var_s1 <- additional_var %>%
  select(idn, v12, v16, v18, v19, v27) %>%
  mutate(
    idn = as.character(idn),

    education = factor(
      if_else(as.numeric(v12) %in% c(1, 2, 3), as.numeric(v12), NA_real_),
      levels = c(1, 2, 3),
      labels = c(
        "incomplete secondary",
        "secondary or vocational",
        "higher education"
      )
    ),

    employment_status = factor(
      if_else(as.numeric(v16) %in% c(1, 2, 3, 4, 5), as.numeric(v16), NA_real_),
      levels = c(1, 2, 3, 4, 5),
      labels = c(
        "employed",
        "student",
        "homemaker",
        "unemployed",
        "retired/disabled"
      )
    ),

    self_rated_health = factor(
      if_else(as.numeric(v18) %in% c(1, 2, 3, 4, 5), as.numeric(v18), NA_real_),
      levels = c(1, 2, 3, 4, 5),
      labels = c(
        "very good",
        "good",
        "fair",
        "poor",
        "very poor"
      )
    ),

    physical_activity = factor(
      if_else(as.numeric(v19) %in% c(1, 2, 3), as.numeric(v19), NA_real_),
      levels = c(1, 2, 3),
      labels = c("yes", "sometimes", "no")
    ),

    smoking = case_when(
      as.numeric(v27) == 1 ~ "never smoked",
      as.numeric(v27) == 2 ~ "former smoker",
      as.numeric(v27) %in% c(3, 4, 5) ~ "current smoker",
      TRUE ~ NA_character_
    ),
    smoking = factor(
      smoking,
      levels = c("never smoked", "former smoker", "current smoker")
    )
  ) %>%
  select(
    idn,
    education,
    employment_status,
    self_rated_health,
    physical_activity,
    smoking
  )

## --- 2. Construct responder indicator --------------------------------

followup_ids <- df_new %>%
  transmute(idn = as.character(idn)) %>%
  distinct()

df_s1 <- df_full %>%
  mutate(
    idn = as.character(idn),
    response = if_else(idn %in% followup_ids$idn, 1L, 0L)
  ) %>%
  left_join(additional_var_s1, by = "idn") %>%
  mutate(
    sex = case_when(
      as.character(sex) %in% c("1", "male") ~ "male",
      as.character(sex) %in% c("2", "female") ~ "female",
      TRUE ~ NA_character_
    ),
    sex = factor(sex, levels = c("male", "female")),

    bmi_baseline = as.numeric(imt),
    tg_mmol = as.numeric(tg_mm),
    hdl_mmol = as.numeric(hdl_mm),

    mets_waist = case_when(
      sex == "male"   & waist >= 94 ~ 1,
      sex == "female" & waist >= 80 ~ 1,
      TRUE ~ 0
    ),
    mets_tg = if_else(tg_mmol >= 1.7, 1, 0, missing = 0),
    mets_hdl = case_when(
      sex == "male"   & hdl_mmol < 1.03 ~ 1,
      sex == "female" & hdl_mmol < 1.29 ~ 1,
      TRUE ~ 0
    ),
    mets_bp   = if_else(syst >= 130 | diast >= 85, 1, 0, missing = 0),
    mets_gluc = if_else(gluc >= 5.6, 1, 0, missing = 0),
    mets_components = mets_tg + mets_hdl + mets_bp + mets_gluc,
    mets_idf = if_else(mets_waist == 1 & mets_components >= 2, 1, 0)
  )

df_resp <- df_s1 %>% filter(response == 1)
df_non  <- df_s1 %>% filter(response == 0)

cat("Responders:", nrow(df_resp), "\n")
cat("Non-responders:", nrow(df_non), "\n")

## --- 3. Build table --------------------------------------------------

table_s1_main <- bind_rows(
  tibble(
    Variable = "Participants, n",
    `Completed follow-up` = as.character(nrow(df_resp)),
    `Not assessed at follow-up` = as.character(nrow(df_non)),
    SMD = ""
  ),
  make_cont_rows("age_f", "Age, years", 1, df_resp, df_non),
  make_binary_rows("sex", "Female, n (%)", "female", df_resp, df_non),
  make_cont_rows("bmi_baseline", "BMI, kg/m²", 1, df_resp, df_non),
  make_cont_rows("waist", "Waist circumference, cm", 1, df_resp, df_non),
  make_cont_rows("tg_mmol", "Triglycerides, mmol/L", 2, df_resp, df_non),
  make_cont_rows("hdl_mmol", "HDL-C, mmol/L", 2, df_resp, df_non),
  make_cont_rows("gluc", "Glucose, mmol/L", 2, df_resp, df_non),
  make_cont_rows("syst", "Systolic BP, mmHg", 1, df_resp, df_non),
  make_binary_rows("mets_idf", "MetS by IDF, n (%)", 1, df_resp, df_non),
  make_cont_rows("cMetS_resid", "cMetS", 2, df_resp, df_non)
)

table_s1_education <- make_cat_rows(
  "education",
  "Education level, n (%)",
  c("incomplete secondary", "secondary or vocational", "higher education"),
  df_resp, df_non
)

table_s1_employment <- make_cat_rows(
  "employment_status",
  "Employment status, n (%)",
  c("employed", "student", "homemaker", "unemployed", "retired/disabled"),
  df_resp, df_non
)

table_s1_health <- make_cat_rows(
  "self_rated_health",
  "Self-rated health, n (%)",
  c("very good", "good", "fair", "poor", "very poor"),
  df_resp, df_non
)

table_s1_physical <- make_cat_rows(
  "physical_activity",
  "Physical activity, n (%)",
  c("yes", "sometimes", "no"),
  df_resp, df_non
)

table_s1_smoking <- make_cat_rows(
  "smoking",
  "Smoking status, n (%)",
  c("never smoked", "former smoker", "current smoker"),
  df_resp, df_non
)

table_s1 <- bind_rows(
  table_s1_main,
  table_s1_education,
  table_s1_employment,
  table_s1_health,
  table_s1_physical,
  table_s1_smoking
)

## --- 4. Print --------------------------------------------------------

table_title <- "Supplementary Table S3. Baseline characteristics of participants according to follow-up assessment status"

cat("\n", table_title, "\n\n", sep = "")
print(table_s1, n = Inf)

## --- 5. Excel export -------------------------------------------------

stopifnot(exists("table_s1"))
stopifnot(nrow(table_s1) > 0)

note_text <- paste(
  "Note: Continuous variables are presented as mean ± standard deviation and categorical variables as number (%).",
  "Missing values are shown for all variables.",
  "Standardized mean differences (SMDs) were used to compare baseline characteristics between participants who completed the follow-up assessment and those not assessed at follow-up.",
  "This table is intended to evaluate potential selection bias.",
  "MetS by IDF was defined using clinical and laboratory criteria; information on antihypertensive treatment was not included.",
  "Abbreviations: BMI = body mass index; HDL-C = high-density lipoprotein cholesterol; BP = blood pressure; MetS = metabolic syndrome; cMetS = continuous metabolic syndrome risk score; IDF = International Diabetes Federation; SMD = standardized mean difference."
)

wb <- createWorkbook()
addWorksheet(wb, "Table S3")

writeData(wb, "Table S3", table_title, startRow = 1, startCol = 1)
writeData(wb, "Table S3", table_s1, startRow = 3, startCol = 1, colNames = TRUE)

note_row <- nrow(table_s1) + 5
writeData(wb, "Table S3", note_text, startRow = note_row, startCol = 1)

title_style    <- createStyle(textDecoration = "bold", fontSize = 12)
header_style   <- createStyle(textDecoration = "bold", border = "bottom")
body_style     <- createStyle(halign = "center")
firstcol_style <- createStyle(halign = "left")
note_style     <- createStyle(wrapText = TRUE, valign = "top")

addStyle(wb, "Table S3", title_style, rows = 1, cols = 1, gridExpand = TRUE)
addStyle(wb, "Table S3", header_style, rows = 3, cols = 1:ncol(table_s1), gridExpand = TRUE)

addStyle(
  wb, "Table S3", body_style,
  rows = 4:(nrow(table_s1) + 3),
  cols = 2:ncol(table_s1),
  gridExpand = TRUE
)

addStyle(
  wb, "Table S3", firstcol_style,
  rows = 4:(nrow(table_s1) + 3),
  cols = 1,
  gridExpand = TRUE
)

setColWidths(wb, "Table S3", cols = 1, widths = 38)
setColWidths(wb, "Table S3", cols = 2:4, widths = "auto")

mergeCells(wb, "Table S3", cols = 1:ncol(table_s1), rows = note_row)
addStyle(wb, "Table S3", note_style, rows = note_row, cols = 1, gridExpand = TRUE)
setRowHeights(wb, "Table S3", rows = note_row, heights = 110)

outfile <- "Table_S3_responders_vs_nonresponders.xlsx"
saveWorkbook(wb, outfile, overwrite = TRUE)

cat("\nExcel file saved:", outfile, "\n")



stopifnot(exists("analysis_initialized"))
stopifnot(isTRUE(analysis_initialized))
stopifnot(exists("df_full"))
stopifnot(exists("df_new"))
stopifnot(exists("df_analysis"))
stopifnot(exists("additional_var"))

library(dplyr)
library(haven)
library(broom)

## =========================================================
## 1. Helper functions
## =========================================================

extract_or_table <- function(model) {
  est <- broom::tidy(model, conf.int = TRUE, exponentiate = TRUE)
  est %>%
    mutate(
      OR_CI = sprintf("%.2f (%.2f to %.2f)", estimate, conf.low, conf.high),
      p = ifelse(p.value < 0.001, "<0.001", sprintf("%.3f", p.value))
    ) %>%
    select(term, OR_CI, p)
}

summarize_weights <- function(x) {
  c(
    min = min(x, na.rm = TRUE),
    p25 = quantile(x, 0.25, na.rm = TRUE),
    median = median(x, na.rm = TRUE),
    mean = mean(x, na.rm = TRUE),
    p75 = quantile(x, 0.75, na.rm = TRUE),
    max = max(x, na.rm = TRUE)
  )
}

## =========================================================
## 2. Construct baseline response dataset
## =========================================================

followup_ids <- df_new %>%
  transmute(idn = as.character(idn)) %>%
  distinct()

additional_var_ipw <- additional_var %>%
  select(idn, v12, v16, v18, v19, v27) %>%
  mutate(
    idn = as.character(idn),

    education = factor(
      if_else(as.numeric(v12) %in% c(1, 2, 3), as.numeric(v12), NA_real_),
      levels = c(1, 2, 3),
      labels = c(
        "incomplete secondary",
        "secondary or vocational",
        "higher education"
      )
    ),

    employment_status = factor(
      if_else(as.numeric(v16) %in% c(1, 2, 3, 4, 5), as.numeric(v16), NA_real_),
      levels = c(1, 2, 3, 4, 5),
      labels = c(
        "employed",
        "student",
        "homemaker",
        "unemployed",
        "retired/disabled"
      )
    ),

    self_rated_health = factor(
      if_else(as.numeric(v18) %in% c(1, 2, 3, 4, 5), as.numeric(v18), NA_real_),
      levels = c(1, 2, 3, 4, 5),
      labels = c(
        "very good",
        "good",
        "fair",
        "poor",
        "very poor"
      )
    ),

    physical_activity = factor(
      if_else(as.numeric(v19) %in% c(1, 2, 3), as.numeric(v19), NA_real_),
      levels = c(1, 2, 3),
      labels = c("yes", "sometimes", "no")
    ),

    smoking = case_when(
      as.numeric(v27) == 1 ~ "never smoked",
      as.numeric(v27) == 2 ~ "former smoker",
      as.numeric(v27) %in% c(3, 4, 5) ~ "current smoker",
      TRUE ~ NA_character_
    ),
    smoking = factor(
      smoking,
      levels = c("never smoked", "former smoker", "current smoker")
    )
  ) %>%
  mutate(
    ## More stable recodings for IPW
    education_ipw = case_when(
      education == "higher education" ~ "higher education",
      education %in% c("incomplete secondary", "secondary or vocational") ~ "no higher education",
      TRUE ~ NA_character_
    ),
    education_ipw = factor(
      education_ipw,
      levels = c("no higher education", "higher education")
    ),

    employment_ipw = case_when(
      employment_status == "employed" ~ "employed",
      employment_status %in% c("student", "homemaker", "unemployed", "retired/disabled") ~ "not employed/other",
      TRUE ~ NA_character_
    ),
    employment_ipw = factor(
      employment_ipw,
      levels = c("employed", "not employed/other")
    ),

    self_rated_health_ipw = case_when(
      self_rated_health %in% c("very good", "good") ~ "good/very good",
      self_rated_health == "fair" ~ "fair",
      self_rated_health %in% c("poor", "very poor") ~ "poor/very poor",
      TRUE ~ NA_character_
    ),
    self_rated_health_ipw = factor(
      self_rated_health_ipw,
      levels = c("good/very good", "fair", "poor/very poor")
    ),

    physical_activity_ipw = case_when(
      physical_activity == "yes" ~ "yes",
      physical_activity %in% c("sometimes", "no") ~ "sometimes/no",
      TRUE ~ NA_character_
    ),
    physical_activity_ipw = factor(
      physical_activity_ipw,
      levels = c("yes", "sometimes/no")
    ),

    smoking_ipw = case_when(
      smoking == "never smoked" ~ "never",
      smoking %in% c("former smoker", "current smoker") ~ "ever",
      TRUE ~ NA_character_
    ),
    smoking_ipw = factor(
      smoking_ipw,
      levels = c("never", "ever")
    )
  ) %>%
  select(
    idn,
    education, employment_status, self_rated_health, physical_activity, smoking,
    education_ipw, employment_ipw, self_rated_health_ipw, physical_activity_ipw, smoking_ipw
  )

df_ipw_base <- df_full %>%
  mutate(
    idn = as.character(idn),
    response = if_else(idn %in% followup_ids$idn, 1L, 0L),
    sex = factor(sex, levels = c("male", "female"))
  ) %>%
  left_join(additional_var_ipw, by = "idn")

cat("Total baseline sample:", nrow(df_ipw_base), "\n")
cat("Responders:", sum(df_ipw_base$response == 1, na.rm = TRUE), "\n")
cat("Non-responders:", sum(df_ipw_base$response == 0, na.rm = TRUE), "\n")

## =========================================================
## 3. Response model datasets
## =========================================================

df_resp_model_main <- df_ipw_base %>%
  select(idn, response, age_f, sex, waist, tg, hdl, gluc, syst) %>%
  na.omit()

cat("\nN in main response model:", nrow(df_resp_model_main), "\n")
cat("Responders in main response model:", sum(df_resp_model_main$response == 1), "\n")
cat("Non-responders in main response model:", sum(df_resp_model_main$response == 0), "\n")

df_resp_model_soc <- df_ipw_base %>%
  select(
    idn, response, age_f, sex, waist, tg, hdl, gluc, syst,
    education_ipw, employment_ipw, self_rated_health_ipw,
    physical_activity_ipw, smoking_ipw
  ) %>%
  na.omit()

cat("\nN in socio-clinical response model:", nrow(df_resp_model_soc), "\n")
cat("Responders in socio-clinical response model:", sum(df_resp_model_soc$response == 1), "\n")
cat("Non-responders in socio-clinical response model:", sum(df_resp_model_soc$response == 0), "\n")

## Optional: inspect recoded distributions
cat("\n=== employment_ipw by response ===\n")
print(table(df_resp_model_soc$employment_ipw, df_resp_model_soc$response, useNA = "ifany"))

cat("\n=== self_rated_health_ipw by response ===\n")
print(table(df_resp_model_soc$self_rated_health_ipw, df_resp_model_soc$response, useNA = "ifany"))

cat("\n=== physical_activity_ipw by response ===\n")
print(table(df_resp_model_soc$physical_activity_ipw, df_resp_model_soc$response, useNA = "ifany"))

cat("\n=== smoking_ipw by response ===\n")
print(table(df_resp_model_soc$smoking_ipw, df_resp_model_soc$response, useNA = "ifany"))

## =========================================================
## 4. Fit response models
## =========================================================

model_response_main <- glm(
  response ~ age_f + sex + waist + tg + hdl + gluc + syst,
  family = binomial,
  data = df_resp_model_main
)

cat("\n=== Main response model ===\n")
print(summary(model_response_main))

cat("\n=== ORs: Main response model ===\n")
print(extract_or_table(model_response_main))

model_response_soc <- glm(
  response ~ age_f + sex + waist + tg + hdl + gluc + syst +
    education_ipw + employment_ipw + self_rated_health_ipw +
    physical_activity_ipw + smoking_ipw,
  family = binomial,
  data = df_resp_model_soc
)

cat("\n=== Socio-clinical response model ===\n")
print(summary(model_response_soc))

cat("\n=== ORs: Socio-clinical response model ===\n")
print(extract_or_table(model_response_soc))

## =========================================================
## 5. Predicted probabilities and IPW
## =========================================================

df_resp_model_main <- df_resp_model_main %>%
  mutate(
    prob_main = predict(model_response_main, type = "response"),
    ipw_main  = 1 / prob_main
  )

cat("\n=== Main IPW summary ===\n")
print(summary(df_resp_model_main$ipw_main))
cat("\n")
print(summarize_weights(df_resp_model_main$ipw_main))

df_resp_model_soc <- df_resp_model_soc %>%
  mutate(
    prob_soc = predict(model_response_soc, type = "response"),
    ipw_soc  = 1 / prob_soc
  )

cat("\n=== Socio-clinical IPW summary ===\n")
print(summary(df_resp_model_soc$ipw_soc))
cat("\n")
print(summarize_weights(df_resp_model_soc$ipw_soc))

## =========================================================
## 6. Attach weights back to responder analysis dataset
## =========================================================

df_weights_main <- df_resp_model_main %>%
  filter(response == 1) %>%
  select(idn, prob_main, ipw_main)

df_weights_soc <- df_resp_model_soc %>%
  filter(response == 1) %>%
  select(idn, prob_soc, ipw_soc)

cat("\nRows in main weights table:", nrow(df_weights_main), "\n")
cat("Rows in socio-clinical weights table:", nrow(df_weights_soc), "\n")

df_analysis_ipw <- df_analysis %>%
  mutate(idn = as.character(idn)) %>%
  left_join(df_weights_main, by = "idn") %>%
  left_join(df_weights_soc, by = "idn")

cat("\nRows in df_analysis_ipw:", nrow(df_analysis_ipw), "\n")
cat("Responders with main IPW:", sum(!is.na(df_analysis_ipw$ipw_main)), "\n")
cat("Responders with socio-clinical IPW:", sum(!is.na(df_analysis_ipw$ipw_soc)), "\n")

## =========================================================
## 7. Outcome models
## =========================================================

model_unweighted <- glm(
  knee_pain_broad ~ cMetS_resid * sex + age_f,
  family = binomial,
  data = df_analysis
)

cat("\n=== Unweighted outcome model ===\n")
print(summary(model_unweighted))

model_ipw_main <- glm(
  knee_pain_broad ~ cMetS_resid * sex + age_f,
  family = binomial,
  data = df_analysis_ipw,
  weights = ipw_main,
  subset = !is.na(ipw_main)
)

cat("\n=== Weighted outcome model (main IPW) ===\n")
print(summary(model_ipw_main))

model_ipw_soc <- glm(
  knee_pain_broad ~ cMetS_resid * sex + age_f,
  family = binomial,
  data = df_analysis_ipw,
  weights = ipw_soc,
  subset = !is.na(ipw_soc)
)

cat("\n=== Weighted outcome model (socio-clinical IPW) ===\n")
print(summary(model_ipw_soc))

## =========================================================
## 8. Compact OR tables
## =========================================================

cat("\n=== ORs: Unweighted outcome model ===\n")
print(extract_or_table(model_unweighted))

cat("\n=== ORs: Weighted outcome model (main IPW) ===\n")
print(extract_or_table(model_ipw_main))

cat("\n=== ORs: Weighted outcome model (socio-clinical IPW) ===\n")
print(extract_or_table(model_ipw_soc))


# IPW sensitivity analysis

stopifnot(exists("analysis_initialized"))
stopifnot(isTRUE(analysis_initialized))
stopifnot(exists("df_full"))
stopifnot(exists("df_new"))
stopifnot(exists("df_analysis"))
stopifnot(exists("additional_var"))

library(dplyr)
library(haven)
library(broom)

extract_or_table <- function(model) {
  est <- broom::tidy(model, conf.int = TRUE, exponentiate = TRUE)
  est %>%
    mutate(
      OR_CI = sprintf("%.2f (%.2f to %.2f)", estimate, conf.low, conf.high),
      p = ifelse(p.value < 0.001, "<0.001", sprintf("%.3f", p.value))
    ) %>%
    select(term, OR_CI, p)
}

summarize_weights <- function(x) {
  c(
    min = min(x, na.rm = TRUE),
    p25 = quantile(x, 0.25, na.rm = TRUE),
    median = median(x, na.rm = TRUE),
    mean = mean(x, na.rm = TRUE),
    p75 = quantile(x, 0.75, na.rm = TRUE),
    max = max(x, na.rm = TRUE)
  )
}

followup_ids <- df_new %>%
  transmute(idn = as.character(idn)) %>%
  distinct()

additional_var_ipw <- additional_var %>%
  select(idn, v12, v16, v18, v19, v27) %>%
  mutate(
    idn = as.character(idn),

    education = factor(
      if_else(as.numeric(v12) %in% c(1, 2, 3), as.numeric(v12), NA_real_),
      levels = c(1, 2, 3),
      labels = c(
        "incomplete secondary",
        "secondary or vocational",
        "higher education"
      )
    ),

    employment_status = factor(
      if_else(as.numeric(v16) %in% c(1, 2, 3, 4, 5), as.numeric(v16), NA_real_),
      levels = c(1, 2, 3, 4, 5),
      labels = c(
        "employed",
        "student",
        "homemaker",
        "unemployed",
        "retired/disabled"
      )
    ),

    self_rated_health = factor(
      if_else(as.numeric(v18) %in% c(1, 2, 3, 4, 5), as.numeric(v18), NA_real_),
      levels = c(1, 2, 3, 4, 5),
      labels = c(
        "very good",
        "good",
        "fair",
        "poor",
        "very poor"
      )
    ),

    physical_activity = factor(
      if_else(as.numeric(v19) %in% c(1, 2, 3), as.numeric(v19), NA_real_),
      levels = c(1, 2, 3),
      labels = c("yes", "sometimes", "no")
    ),

    smoking = case_when(
      as.numeric(v27) == 1 ~ "never smoked",
      as.numeric(v27) == 2 ~ "former smoker",
      as.numeric(v27) %in% c(3, 4, 5) ~ "current smoker",
      TRUE ~ NA_character_
    ),
    smoking = factor(
      smoking,
      levels = c("never smoked", "former smoker", "current smoker")
    )
  ) %>%
  mutate(
    education_ipw = case_when(
      education == "higher education" ~ "higher education",
      education %in% c("incomplete secondary", "secondary or vocational") ~ "no higher education",
      TRUE ~ NA_character_
    ),
    education_ipw = factor(education_ipw, levels = c("no higher education", "higher education")),

    employment_ipw = case_when(
      employment_status == "employed" ~ "employed",
      employment_status %in% c("student", "homemaker", "unemployed", "retired/disabled") ~ "not employed/other",
      TRUE ~ NA_character_
    ),
    employment_ipw = factor(employment_ipw, levels = c("employed", "not employed/other")),

    self_rated_health_ipw = case_when(
      self_rated_health %in% c("very good", "good") ~ "good/very good",
      self_rated_health == "fair" ~ "fair",
      self_rated_health %in% c("poor", "very poor") ~ "poor/very poor",
      TRUE ~ NA_character_
    ),
    self_rated_health_ipw = factor(self_rated_health_ipw, levels = c("good/very good", "fair", "poor/very poor")),

    physical_activity_ipw = case_when(
      physical_activity == "yes" ~ "yes",
      physical_activity %in% c("sometimes", "no") ~ "sometimes/no",
      TRUE ~ NA_character_
    ),
    physical_activity_ipw = factor(physical_activity_ipw, levels = c("yes", "sometimes/no")),

    smoking_ipw = case_when(
      smoking == "never smoked" ~ "never",
      smoking %in% c("former smoker", "current smoker") ~ "ever",
      TRUE ~ NA_character_
    ),
    smoking_ipw = factor(smoking_ipw, levels = c("never", "ever"))
  ) %>%
  select(
    idn,
    education_ipw, employment_ipw, self_rated_health_ipw,
    physical_activity_ipw, smoking_ipw
  )

df_ipw_base <- df_full %>%
  mutate(
    idn = as.character(idn),
    response = if_else(idn %in% followup_ids$idn, 1L, 0L),
    sex = factor(sex, levels = c("male", "female"))
  ) %>%
  left_join(additional_var_ipw, by = "idn")

df_resp_model_main <- df_ipw_base %>%
  select(idn, response, age_f, sex, waist, tg, hdl, gluc, syst) %>%
  na.omit()

df_resp_model_soc <- df_ipw_base %>%
  select(
    idn, response, age_f, sex, waist, tg, hdl, gluc, syst,
    education_ipw, employment_ipw, self_rated_health_ipw,
    physical_activity_ipw, smoking_ipw
  ) %>%
  na.omit()

model_response_main <- glm(
  response ~ age_f + sex + waist + tg + hdl + gluc + syst,
  family = binomial,
  data = df_resp_model_main
)

model_response_soc <- glm(
  response ~ age_f + sex + waist + tg + hdl + gluc + syst +
    education_ipw + employment_ipw + self_rated_health_ipw +
    physical_activity_ipw + smoking_ipw,
  family = binomial,
  data = df_resp_model_soc
)

df_resp_model_main <- df_resp_model_main %>%
  mutate(
    prob_main = predict(model_response_main, type = "response"),
    ipw_main = 1 / prob_main
  )

df_resp_model_soc <- df_resp_model_soc %>%
  mutate(
    prob_soc = predict(model_response_soc, type = "response"),
    ipw_soc = 1 / prob_soc
  )

df_weights_main <- df_resp_model_main %>%
  filter(response == 1) %>%
  select(idn, prob_main, ipw_main)

df_weights_soc <- df_resp_model_soc %>%
  filter(response == 1) %>%
  select(idn, prob_soc, ipw_soc)

df_analysis_ipw <- df_analysis %>%
  mutate(idn = as.character(idn)) %>%
  left_join(df_weights_main, by = "idn") %>%
  left_join(df_weights_soc, by = "idn")

cat("Objects created:\n")
cat("model_response_main:", exists("model_response_main"), "\n")
cat("model_response_soc:", exists("model_response_soc"), "\n")
cat("df_resp_model_main:", exists("df_resp_model_main"), "\n")
cat("df_resp_model_soc:", exists("df_resp_model_soc"), "\n")

# Supplementary Table S4. Predictors of follow-up participation and inverse probability weights


stopifnot(exists("analysis_initialized"))
stopifnot(isTRUE(analysis_initialized))
stopifnot(exists("df_full"))
stopifnot(exists("df_new"))
stopifnot(exists("df_analysis"))
stopifnot(exists("additional_var"))

library(dplyr)
library(haven)
library(broom)

## =========================================================
## 1. Construct baseline response dataset
## =========================================================

followup_ids <- df_new %>%
  transmute(idn = as.character(idn)) %>%
  distinct()

## Additional baseline variables for optional extended model
additional_var_ipw <- additional_var %>%
  select(idn, v12, v16, v18, v19, v27) %>%
  mutate(
    idn = as.character(idn),

    education = factor(
      if_else(as.numeric(v12) %in% c(1, 2, 3), as.numeric(v12), NA_real_),
      levels = c(1, 2, 3),
      labels = c(
        "incomplete secondary",
        "secondary or vocational",
        "higher education"
      )
    ),

    employment_status = factor(
      if_else(as.numeric(v16) %in% c(1, 2, 3, 4, 5), as.numeric(v16), NA_real_),
      levels = c(1, 2, 3, 4, 5),
      labels = c(
        "employed",
        "student",
        "homemaker",
        "unemployed",
        "retired/disabled"
      )
    ),

    self_rated_health = factor(
      if_else(as.numeric(v18) %in% c(1, 2, 3, 4, 5), as.numeric(v18), NA_real_),
      levels = c(1, 2, 3, 4, 5),
      labels = c(
        "very good",
        "good",
        "fair",
        "poor",
        "very poor"
      )
    ),

    physical_activity = factor(
      if_else(as.numeric(v19) %in% c(1, 2, 3), as.numeric(v19), NA_real_),
      levels = c(1, 2, 3),
      labels = c("yes", "sometimes", "no")
    ),

    smoking = case_when(
      as.numeric(v27) == 1 ~ "never smoked",
      as.numeric(v27) == 2 ~ "former smoker",
      as.numeric(v27) %in% c(3, 4, 5) ~ "current smoker",
      TRUE ~ NA_character_
    ),
    smoking = factor(
      smoking,
      levels = c("never smoked", "former smoker", "current smoker")
    )
  ) %>%
  select(
    idn,
    education,
    employment_status,
    self_rated_health,
    physical_activity,
    smoking
  )

## Build baseline cohort with response indicator
df_ipw_base <- df_full %>%
  mutate(
    idn = as.character(idn),
    response = if_else(idn %in% followup_ids$idn, 1L, 0L),
    sex = factor(sex, levels = c("male", "female"))
  ) %>%
  left_join(additional_var_ipw, by = "idn")

cat("Total baseline sample:", nrow(df_ipw_base), "\n")
cat("Responders:", sum(df_ipw_base$response == 1, na.rm = TRUE), "\n")
cat("Non-responders:", sum(df_ipw_base$response == 0, na.rm = TRUE), "\n")

## =========================================================
## 2. Response model dataset
## =========================================================

## Main response model: core cardiometabolic variables only
df_resp_model_main <- df_ipw_base %>%
  select(idn, response, age_f, sex, waist, tg, hdl, gluc, syst) %>%
  na.omit()

cat("\nN in main response model:", nrow(df_resp_model_main), "\n")
cat("Responders in main response model:", sum(df_resp_model_main$response == 1), "\n")
cat("Non-responders in main response model:", sum(df_resp_model_main$response == 0), "\n")

## Optional extended response model with extra baseline variables
df_resp_model_ext <- df_ipw_base %>%
  select(
    idn, response, age_f, sex, waist, tg, hdl, gluc, syst,
    education, employment_status, self_rated_health, physical_activity, smoking
  ) %>%
  na.omit()

cat("\nN in extended response model:", nrow(df_resp_model_ext), "\n")
cat("Responders in extended response model:", sum(df_resp_model_ext$response == 1), "\n")
cat("Non-responders in extended response model:", sum(df_resp_model_ext$response == 0), "\n")

## =========================================================
## 3. Fit response models
## =========================================================

model_response_main <- glm(
  response ~ age_f + sex + waist + tg + hdl + gluc + syst,
  family = binomial,
  data = df_resp_model_main
)

cat("\n=== Main response model ===\n")
print(summary(model_response_main))

cat("\nOdds ratios (main response model):\n")
print(exp(cbind(OR = coef(model_response_main), confint(model_response_main))))

model_response_ext <- glm(
  response ~ age_f + sex + waist + tg + hdl + gluc + syst +
    education + employment_status + self_rated_health + physical_activity + smoking,
  family = binomial,
  data = df_resp_model_ext
)

cat("\n=== Extended response model ===\n")
print(summary(model_response_ext))

cat("\nOdds ratios (extended response model):\n")
print(exp(cbind(OR = coef(model_response_ext), confint(model_response_ext))))

## =========================================================
## 4. Predicted probabilities and IPW
## =========================================================

df_resp_model_main <- df_resp_model_main %>%
  mutate(
    prob_main = predict(model_response_main, type = "response"),
    ipw_main  = 1 / prob_main
  )

cat("\n=== Main IPW summary ===\n")
print(summary(df_resp_model_main$ipw_main))

df_resp_model_ext <- df_resp_model_ext %>%
  mutate(
    prob_ext = predict(model_response_ext, type = "response"),
    ipw_ext  = 1 / prob_ext
  )

cat("\n=== Extended IPW summary ===\n")
print(summary(df_resp_model_ext$ipw_ext))

## =========================================================
## 5. Attach weights back to responder analysis dataset
## =========================================================

df_weights_main <- df_resp_model_main %>%
  filter(response == 1) %>%
  select(idn, prob_main, ipw_main)

df_weights_ext <- df_resp_model_ext %>%
  filter(response == 1) %>%
  select(idn, prob_ext, ipw_ext)

cat("\nRows in main weights table:", nrow(df_weights_main), "\n")
cat("Rows in extended weights table:", nrow(df_weights_ext), "\n")

df_analysis_ipw <- df_analysis %>%
  mutate(idn = as.character(idn)) %>%
  left_join(df_weights_main, by = "idn") %>%
  left_join(df_weights_ext, by = "idn")

cat("\nRows in df_analysis_ipw:", nrow(df_analysis_ipw), "\n")
cat("Responders with main IPW:", sum(!is.na(df_analysis_ipw$ipw_main)), "\n")
cat("Responders with extended IPW:", sum(!is.na(df_analysis_ipw$ipw_ext)), "\n")

## =========================================================
## 6. Example weighted outcome models
## =========================================================

## Unweighted model
model_unweighted <- glm(
  knee_pain_broad ~ cMetS_resid * sex + age_f,
  family = binomial,
  data = df_analysis
)

cat("\n=== Unweighted outcome model ===\n")
print(summary(model_unweighted))

## Weighted model using main IPW
model_ipw_main <- glm(
  knee_pain_broad ~ cMetS_resid * sex + age_f,
  family = binomial,
  data = df_analysis_ipw,
  weights = ipw_main,
  subset = !is.na(ipw_main)
)

cat("\n=== Weighted outcome model (main IPW) ===\n")
print(summary(model_ipw_main))

## Weighted model using extended IPW
model_ipw_ext <- glm(
  knee_pain_broad ~ cMetS_resid * sex + age_f,
  family = binomial,
  data = df_analysis_ipw,
  weights = ipw_ext,
  subset = !is.na(ipw_ext)
)

cat("\n=== Weighted outcome model (extended IPW) ===\n")
print(summary(model_ipw_ext))

## =========================================================
## 7. Compact coefficient tables
## =========================================================

extract_or_table <- function(model) {
  est <- broom::tidy(model, conf.int = TRUE, exponentiate = TRUE)
  est %>%
    mutate(
      OR_CI = sprintf("%.2f (%.2f to %.2f)", estimate, conf.low, conf.high),
      p = ifelse(p.value < 0.001, "<0.001", sprintf("%.3f", p.value))
    ) %>%
    select(term, OR_CI, p)
}

cat("\n=== ORs: Unweighted model ===\n")
print(extract_or_table(model_unweighted))

cat("\n=== ORs: Weighted model (main IPW) ===\n")
print(extract_or_table(model_ipw_main))

cat("\n=== ORs: Weighted model (extended IPW) ===\n")
print(extract_or_table(model_ipw_ext))


stopifnot(exists("analysis_initialized"))
stopifnot(isTRUE(analysis_initialized))
stopifnot(exists("model_response_main"))
stopifnot(exists("model_response_soc"))
stopifnot(exists("df_resp_model_main"))
stopifnot(exists("df_resp_model_soc"))

library(dplyr)
library(broom)
library(openxlsx)

## =========================================================
## Helper
## =========================================================

fmt_or <- function(est, lcl, ucl) {
  ifelse(
    is.na(est) | is.na(lcl) | is.na(ucl),
    "",
    sprintf("%.2f (%.2f to %.2f)", est, lcl, ucl)
  )
}

fmt_p <- function(p) {
  ifelse(
    is.na(p),
    "",
    ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
  )
}

tidy_model <- function(model) {
  broom::tidy(model, conf.int = TRUE, exponentiate = TRUE)
}

## =========================================================
## Tidy models
## =========================================================

soc_tidy  <- tidy_model(model_response_soc)
main_tidy <- tidy_model(model_response_main)

tab_models <- full_join(
  soc_tidy %>%
    select(term, estimate, conf.low, conf.high, p.value),
  main_tidy %>%
    select(term, estimate, conf.low, conf.high, p.value),
  by = "term",
  suffix = c("_soc", "_main")
)

## =========================================================
## Format model results
## =========================================================

tab_models <- tab_models %>%
  mutate(
    `Primary model OR (95% CI)` =
      fmt_or(estimate_soc, conf.low_soc, conf.high_soc),
    `Primary model P-value` =
      fmt_p(p.value_soc),

    `Additional model OR (95% CI)` =
      fmt_or(estimate_main, conf.low_main, conf.high_main),
    `Additional model P-value` =
      fmt_p(p.value_main)
  ) %>%
  mutate(
    Variable = case_when(
      term == "(Intercept)" ~ "Intercept",
      term == "age_f" ~ "Age (per year)",
      term == "sexfemale" ~ "Female (vs male)",
      term == "waist" ~ "Waist circumference, cm",
      term == "tg" ~ "Triglycerides, mg/dL",
      term == "hdl" ~ "HDL-C, mg/dL",
      term == "gluc" ~ "Glucose, mmol/L",
      term == "syst" ~ "Systolic BP, mmHg",
      term == "education_ipwhigher education" ~
        "Higher education (vs no higher education)",
      term == "employment_ipwnot employed/other" ~
        "Not employed/other (vs employed)",
      term == "self_rated_health_ipwfair" ~
        "Self-rated health: fair (vs good/very good)",
      term == "self_rated_health_ipwpoor/very poor" ~
        "Self-rated health: poor/very poor (vs good/very good)",
      term == "physical_activity_ipwsometimes/no" ~
        "Physical activity: sometimes/no (vs yes)",
      term == "smoking_ipwever" ~
        "Ever smoker (vs never)",
      TRUE ~ term
    )
  ) %>%
  select(
    Variable,
    `Primary model OR (95% CI)`,
    `Primary model P-value`,
    `Additional model OR (95% CI)`,
    `Additional model P-value`
  )

## Optional ordering
variable_order <- c(
  "Intercept",
  "Age (per year)",
  "Female (vs male)",
  "Waist circumference, cm",
  "Triglycerides, mg/dL",
  "HDL-C, mg/dL",
  "Glucose, mmol/L",
  "Systolic BP, mmHg",
  "Higher education (vs no higher education)",
  "Not employed/other (vs employed)",
  "Self-rated health: fair (vs good/very good)",
  "Self-rated health: poor/very poor (vs good/very good)",
  "Physical activity: sometimes/no (vs yes)",
  "Ever smoker (vs never)"
)

tab_models <- tab_models %>%
  mutate(Variable = factor(Variable, levels = variable_order)) %>%
  arrange(Variable) %>%
  mutate(Variable = as.character(Variable))

## =========================================================
## Weight summaries
## =========================================================

weights_summary <- tibble(
  Variable = c(
    "N in response model",
    "Responders",
    "Non-responders",
    "IPW: min",
    "IPW: 25th percentile",
    "IPW: median",
    "IPW: mean",
    "IPW: 75th percentile",
    "IPW: max"
  ),

  `Primary model OR (95% CI)` = c(
    nrow(df_resp_model_soc),
    sum(df_resp_model_soc$response == 1),
    sum(df_resp_model_soc$response == 0),
    sprintf("%.2f", min(df_resp_model_soc$ipw_soc, na.rm = TRUE)),
    sprintf("%.2f", quantile(df_resp_model_soc$ipw_soc, 0.25, na.rm = TRUE)),
    sprintf("%.2f", median(df_resp_model_soc$ipw_soc, na.rm = TRUE)),
    sprintf("%.2f", mean(df_resp_model_soc$ipw_soc, na.rm = TRUE)),
    sprintf("%.2f", quantile(df_resp_model_soc$ipw_soc, 0.75, na.rm = TRUE)),
    sprintf("%.2f", max(df_resp_model_soc$ipw_soc, na.rm = TRUE))
  ),

  `Primary model P-value` = "",

  `Additional model OR (95% CI)` = c(
    nrow(df_resp_model_main),
    sum(df_resp_model_main$response == 1),
    sum(df_resp_model_main$response == 0),
    sprintf("%.2f", min(df_resp_model_main$ipw_main, na.rm = TRUE)),
    sprintf("%.2f", quantile(df_resp_model_main$ipw_main, 0.25, na.rm = TRUE)),
    sprintf("%.2f", median(df_resp_model_main$ipw_main, na.rm = TRUE)),
    sprintf("%.2f", mean(df_resp_model_main$ipw_main, na.rm = TRUE)),
    sprintf("%.2f", quantile(df_resp_model_main$ipw_main, 0.75, na.rm = TRUE)),
    sprintf("%.2f", max(df_resp_model_main$ipw_main, na.rm = TRUE))
  ),

  `Additional model P-value` = ""
)

## =========================================================
## Combine
## =========================================================

blank_row <- tibble(
  Variable = "",
  `Primary model OR (95% CI)` = "",
  `Primary model P-value` = "",
  `Additional model OR (95% CI)` = "",
  `Additional model P-value` = ""
)

tab_final <- bind_rows(
  tab_models,
  blank_row,
  weights_summary
)

cat(
  "\nSupplementary Table S4. Predictors of follow-up participation and inverse probability weights\n\n"
)
print(tab_final, n = Inf)

## =========================================================
## Excel export
## =========================================================

wb <- createWorkbook()
addWorksheet(wb, "Table S4")

title <- "Supplementary Table S4. Predictors of follow-up participation and inverse probability weights"

writeData(wb, "Table S4", title, startRow = 1, startCol = 1)
writeData(wb, "Table S4", tab_final, startRow = 3, startCol = 1)

note <- paste(
  "Note: Odds ratios (ORs) are derived from logistic regression models of follow-up participation.\n",
  "The primary response model included age, sex, waist circumference, triglycerides, HDL-C, glucose, systolic blood pressure, as well as education, employment status, self-rated health, physical activity, and smoking (recoded categories).\n",
  "A simpler cardiometabolic model including age, sex, waist circumference, triglycerides, HDL-C, glucose, and systolic blood pressure was used as an additional sensitivity analysis.\n",
  "Inverse probability weights (IPW) were calculated as the inverse of the predicted probability of follow-up participation."
)

note_row <- nrow(tab_final) + 5
writeData(wb, "Table S4", note, startRow = note_row, startCol = 1)

## Styles
title_style <- createStyle(
  textDecoration = "bold",
  fontSize = 12,
  wrapText = TRUE
)

header_style <- createStyle(
  textDecoration = "bold",
  border = "bottom",
  wrapText = TRUE,
 halign = "center"
)

body_style <- createStyle(
  wrapText = TRUE,
  valign = "top",
  halign = "center"
)

firstcol_style <- createStyle(
  wrapText = TRUE,
  valign = "top",
  halign = "left"
)

note_style <- createStyle(
  wrapText = TRUE,
  valign = "top"
)

addStyle(wb, "Table S4", title_style, rows = 1, cols = 1, gridExpand = TRUE)
addStyle(wb, "Table S4", header_style, rows = 3, cols = 1:5, gridExpand = TRUE)

addStyle(
  wb, "Table S4", body_style,
  rows = 4:(nrow(tab_final) + 3),
  cols = 2:5,
  gridExpand = TRUE
)

addStyle(
  wb, "Table S4", firstcol_style,
  rows = 4:(nrow(tab_final) + 3),
  cols = 1,
  gridExpand = TRUE
)

addStyle(
  wb, "Table S4", note_style,
  rows = note_row,
  cols = 1,
  gridExpand = TRUE
)

## Layout
mergeCells(wb, "Table S4", cols = 1:5, rows = 1)
mergeCells(wb, "Table S4", cols = 1:5, rows = note_row)

setColWidths(wb, "Table S4", cols = 1, widths = 48)
setColWidths(wb, "Table S4", cols = 2, widths = 24)
setColWidths(wb, "Table S4", cols = 3, widths = 12)
setColWidths(wb, "Table S4", cols = 4, widths = 24)
setColWidths(wb, "Table S4", cols = 5, widths = 12)

setRowHeights(wb, "Table S4", rows = 1, heights = 28)
setRowHeights(wb, "Table S4", rows = 3, heights = 45)
setRowHeights(wb, "Table S4", rows = 4:(nrow(tab_final) + 3), heights = 18)
setRowHeights(wb, "Table S4", rows = note_row, heights = 95)

freezePane(wb, "Table S4", firstActiveRow = 4)

saveWorkbook(wb, "Table_S4_IPW_models.xlsx", overwrite = TRUE)

cat("\nExcel file saved: Table_S4_IPW_models.xlsx\n")

# Supplementary Table S5. Sensitivity analyses IPW

stopifnot(
  exists("analysis_initialized"),
  isTRUE(analysis_initialized),
  exists("df_analysis"),
  exists("model_response_main"),
  exists("model_response_soc"),
  exists("df_resp_model_main"),
  exists("df_resp_model_soc")
)

library(dplyr)
library(tidyr)
library(openxlsx)
library(tibble)

## =========================================================
## Rebuild predicted probabilities, IPW, and df_analysis_ipw
## =========================================================

if (!"prob_main" %in% names(df_resp_model_main)) {
  df_resp_model_main <- df_resp_model_main %>%
    mutate(
      prob_main = predict(model_response_main, newdata = ., type = "response")
    )
}

if (!"ipw_main" %in% names(df_resp_model_main)) {
  df_resp_model_main <- df_resp_model_main %>%
    mutate(
      ipw_main = 1 / prob_main
    )
}

if (!"prob_soc" %in% names(df_resp_model_soc)) {
  df_resp_model_soc <- df_resp_model_soc %>%
    mutate(
      prob_soc = predict(model_response_soc, newdata = ., type = "response")
    )
}

if (!"ipw_soc" %in% names(df_resp_model_soc)) {
  df_resp_model_soc <- df_resp_model_soc %>%
    mutate(
      ipw_soc = 1 / prob_soc
    )
}

df_weights_main <- df_resp_model_main %>%
  filter(response == 1) %>%
  select(idn, prob_main, ipw_main)

df_weights_soc <- df_resp_model_soc %>%
  filter(response == 1) %>%
  select(idn, prob_soc, ipw_soc)

df_analysis_ipw <- df_analysis %>%
  mutate(idn = as.character(idn)) %>%
  left_join(df_weights_main, by = "idn") %>%
  left_join(df_weights_soc, by = "idn")

cat("df_analysis rows:", nrow(df_analysis), "\n")
cat("df_analysis_ipw rows:", nrow(df_analysis_ipw), "\n")
cat("Non-missing ipw_soc:", sum(!is.na(df_analysis_ipw$ipw_soc)), "\n")
cat("Non-missing ipw_main:", sum(!is.na(df_analysis_ipw$ipw_main)), "\n")

## =========================================================
## Helpers
## =========================================================

fmt_p <- function(p) {
  ifelse(
    is.na(p),
    "",
    ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
  )
}

fmt_or_ci <- function(or, lcl, ucl) {
  ifelse(
    is.na(or) | is.na(lcl) | is.na(ucl),
    "",
    sprintf("%.2f (%.2f to %.2f)", or, lcl, ucl)
  )
}

extract_sex_specific <- function(model, sex_var = "sexfemale", exposure = "cMetS_resid") {
  co <- summary(model)$coefficients
  vc <- vcov(model)

  interaction_name <- paste0(exposure, ":", sex_var)
  if (!interaction_name %in% rownames(co)) {
    interaction_name <- paste0(sex_var, ":", exposure)
  }

  beta_male <- co[exposure, "Estimate"]
  se_male   <- co[exposure, "Std. Error"]
  p_male    <- co[exposure, "Pr(>|z|)"]

  or_male  <- exp(beta_male)
  lcl_male <- exp(beta_male - 1.96 * se_male)
  ucl_male <- exp(beta_male + 1.96 * se_male)

  beta_female <- co[exposure, "Estimate"] + co[interaction_name, "Estimate"]

  var_female <- vc[exposure, exposure] +
    vc[interaction_name, interaction_name] +
    2 * vc[exposure, interaction_name]

  se_female <- sqrt(var_female)
  z_female  <- beta_female / se_female
  p_female  <- 2 * pnorm(-abs(z_female))

  or_female  <- exp(beta_female)
  lcl_female <- exp(beta_female - 1.96 * se_female)
  ucl_female <- exp(beta_female + 1.96 * se_female)

  beta_interaction <- co[interaction_name, "Estimate"]
  se_interaction   <- co[interaction_name, "Std. Error"]
  p_interaction    <- co[interaction_name, "Pr(>|z|)"]

  or_interaction  <- exp(beta_interaction)
  lcl_interaction <- exp(beta_interaction - 1.96 * se_interaction)
  ucl_interaction <- exp(beta_interaction + 1.96 * se_interaction)

  tibble(
    `Male OR (95% CI)` = fmt_or_ci(or_male, lcl_male, ucl_male),
    `Male P-value` = fmt_p(p_male),
    `Female OR (95% CI)` = fmt_or_ci(or_female, lcl_female, ucl_female),
    `Female P-value` = fmt_p(p_female),
    `Interaction ratio of ORs (95% CI)` = fmt_or_ci(or_interaction, lcl_interaction, ucl_interaction),
    `Interaction P-value` = fmt_p(p_interaction)
  )
}

make_counts <- function(data) {
  data %>%
    filter(!is.na(sex), !is.na(knee_pain_broad)) %>%
    group_by(sex, knee_pain_broad) %>%
    summarise(n = n(), .groups = "drop") %>%
    pivot_wider(
      names_from = knee_pain_broad,
      values_from = n,
      values_fill = 0
    ) %>%
    rename(
      no_pain = `0`,
      pain = `1`
    ) %>%
    mutate(
      sex = case_when(
        sex == "male" ~ "Male",
        sex == "female" ~ "Female",
        TRUE ~ as.character(sex)
      )
    ) %>%
    arrange(factor(sex, levels = c("Male", "Female")))
}

counts_to_text <- function(tab) {
  male_row <- tab %>% filter(sex == "Male")
  female_row <- tab %>% filter(sex == "Female")

  male_text <- if (nrow(male_row) == 1) paste0(male_row$pain, "/", male_row$no_pain) else "NA/NA"
  female_text <- if (nrow(female_row) == 1) paste0(female_row$pain, "/", female_row$no_pain) else "NA/NA"

  paste0("Male: ", male_text, "; Female: ", female_text)
}

## =========================================================
## Datasets
## =========================================================

df_unweighted <- df_analysis %>%
  select(knee_pain_broad, cMetS_resid, sex, age_f) %>%
  filter(
    !is.na(knee_pain_broad),
    !is.na(cMetS_resid),
    !is.na(sex),
    !is.na(age_f)
  )

df_ipw_primary <- df_analysis_ipw %>%
  select(knee_pain_broad, cMetS_resid, sex, age_f, ipw_soc) %>%
  filter(
    !is.na(knee_pain_broad),
    !is.na(cMetS_resid),
    !is.na(sex),
    !is.na(age_f),
    !is.na(ipw_soc)
  )

df_ipw_additional <- df_analysis_ipw %>%
  select(knee_pain_broad, cMetS_resid, sex, age_f, ipw_main) %>%
  filter(
    !is.na(knee_pain_broad),
    !is.na(cMetS_resid),
    !is.na(sex),
    !is.na(age_f),
    !is.na(ipw_main)
  )

cat("N in unweighted analysis:", nrow(df_unweighted), "\n")
cat("N in primary IPW analysis:", nrow(df_ipw_primary), "\n")
cat("N in additional IPW analysis:", nrow(df_ipw_additional), "\n")

## =========================================================
## Models
## =========================================================

model_unweighted_s5 <- glm(
  knee_pain_broad ~ cMetS_resid * sex + age_f,
  family = binomial,
  data = df_unweighted
)

model_ipw_primary_s5 <- glm(
  knee_pain_broad ~ cMetS_resid * sex + age_f,
  family = binomial,
  data = df_ipw_primary,
  weights = ipw_soc
)

model_ipw_additional_s5 <- glm(
  knee_pain_broad ~ cMetS_resid * sex + age_f,
  family = binomial,
  data = df_ipw_additional,
  weights = ipw_main
)

## =========================================================
## Counts
## =========================================================

counts_unweighted <- make_counts(df_unweighted)
counts_primary    <- make_counts(df_ipw_primary)
counts_additional <- make_counts(df_ipw_additional)

## =========================================================
## Final table
## =========================================================

tab_final_s5 <- bind_rows(
  tibble(
    Analysis = "Main analysis (unweighted)",
    N = nrow(df_unweighted),
    `Knee pain / no knee pain by sex, n` = counts_to_text(counts_unweighted)
  ) %>% bind_cols(extract_sex_specific(model_unweighted_s5)),

  tibble(
    Analysis = "IPW-weighted analysis, extended socio-clinical response model",
    N = nrow(df_ipw_primary),
    `Knee pain / no knee pain by sex, n` = counts_to_text(counts_primary)
  ) %>% bind_cols(extract_sex_specific(model_ipw_primary_s5)),

  tibble(
    Analysis = "IPW-weighted analysis, cardiometabolic response model",
    N = nrow(df_ipw_additional),
    `Knee pain / no knee pain by sex, n` = counts_to_text(counts_additional)
  ) %>% bind_cols(extract_sex_specific(model_ipw_additional_s5))
)

cat(
  "\nSupplementary Table S5. Sensitivity analyses of the association between baseline cMetS and knee pain using inverse probability weighting\n\n"
)
print(tab_final_s5, n = Inf)

## =========================================================
## Excel export
## =========================================================

tab_export_s5 <- tab_final_s5
colnames(tab_export_s5) <- c(
  "Analysis",
  "N",
  "Knee pain / no knee pain\nby sex, n",
  "Male\nOR (95% CI)",
  "Male\nP-value",
  "Female\nOR (95% CI)",
  "Female\nP-value",
  "Interaction ratio of ORs\n(95% CI)",
  "Interaction\nP-value"
)

wb <- createWorkbook()
addWorksheet(wb, "Table S5")

title_text <- "Supplementary Table S5. Sensitivity analyses of the association between baseline cMetS and knee pain using inverse probability weighting"
writeData(wb, "Table S5", title_text, startRow = 1, startCol = 1)
writeData(wb, "Table S5", tab_export_s5, startRow = 3, startCol = 1)

note_text <- paste(
  "Note: Odds ratios represent the association between baseline cMetS and knee pain per 1-unit increase in cMetS.",
  "Sex-specific estimates were derived from logistic regression models including a cMetS × sex interaction term and adjusted for age.",
  "IPW analyses were performed as sensitivity analyses using two alternative response models: a socio-clinical response model and a simpler cardiometabolic response model.",
  "The interaction term is presented as the ratio of odds ratios on the multiplicative scale."
)

note_row <- nrow(tab_export_s5) + 5
writeData(wb, "Table S5", note_text, startRow = note_row, startCol = 1)

title_style <- createStyle(textDecoration = "bold", fontSize = 12, wrapText = TRUE)
header_style <- createStyle(textDecoration = "bold", border = "bottom", wrapText = TRUE, halign = "center", valign = "center")
body_style <- createStyle(wrapText = TRUE, valign = "top", halign = "center")
firstcol_style <- createStyle(wrapText = TRUE, valign = "top", halign = "left")
note_style <- createStyle(wrapText = TRUE, valign = "top")

addStyle(wb, "Table S5", title_style, rows = 1, cols = 1, gridExpand = TRUE)
addStyle(wb, "Table S5", header_style, rows = 3, cols = 1:ncol(tab_export_s5), gridExpand = TRUE)
addStyle(wb, "Table S5", body_style, rows = 4:(nrow(tab_export_s5) + 3), cols = 2:ncol(tab_export_s5), gridExpand = TRUE)
addStyle(wb, "Table S5", firstcol_style, rows = 4:(nrow(tab_export_s5) + 3), cols = 1, gridExpand = TRUE)
addStyle(wb, "Table S5", note_style, rows = note_row, cols = 1, gridExpand = TRUE)

mergeCells(wb, "Table S5", cols = 1:ncol(tab_export_s5), rows = 1)
mergeCells(wb, "Table S5", cols = 1:ncol(tab_export_s5), rows = note_row)

setColWidths(wb, "Table S5", cols = 1, widths = 38)
setColWidths(wb, "Table S5", cols = 2, widths = 8)
setColWidths(wb, "Table S5", cols = 3, widths = 28)
setColWidths(wb, "Table S5", cols = 4, widths = 18)
setColWidths(wb, "Table S5", cols = 5, widths = 10)
setColWidths(wb, "Table S5", cols = 6, widths = 18)
setColWidths(wb, "Table S5", cols = 7, widths = 10)
setColWidths(wb, "Table S5", cols = 8, widths = 22)
setColWidths(wb, "Table S5", cols = 9, widths = 12)

setRowHeights(wb, "Table S5", rows = 1, heights = 40)
setRowHeights(wb, "Table S5", rows = 3, heights = 50)
setRowHeights(wb, "Table S5", rows = 4:(nrow(tab_export_s5) + 3), heights = 35)
setRowHeights(wb, "Table S5", rows = note_row, heights = 90)

freezePane(wb, "Table S5", firstActiveRow = 4)

saveWorkbook(wb, "Table_S5_IPW_sensitivity.xlsx", overwrite = TRUE)

cat("\nExcel file saved: Table_S5_IPW_sensitivity.xlsx\n")

stopifnot(exists("analysis_initialized"))
stopifnot(isTRUE(analysis_initialized))
stopifnot(exists("df_analysis"))

library(dplyr)
library(tidyr)
library(openxlsx)

## =========================================================
## Helpers (local)
## =========================================================

fmt_p_s6 <- function(p) {
  ifelse(is.na(p), "", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

fmt_or_ci_s6 <- function(or, lcl, ucl) {
  ifelse(
    is.na(or) | is.na(lcl) | is.na(ucl),
    "",
    sprintf("%.2f (%.2f to %.2f)", or, lcl, ucl)
  )
}

## =========================================================
## Dataset (new object)
## =========================================================

df_s6 <- df_analysis %>%
  select(knee_pain_broad, cMetS_resid, sex, age_f1, bmi_1) %>%
  filter(
    !is.na(knee_pain_broad),
    !is.na(cMetS_resid),
    !is.na(sex),
    !is.na(age_f1),
    !is.na(bmi_1)
  )

cat("N in BMI-adjusted complete-case analysis:\n")
print(nrow(df_s6))

## =========================================================
## Model
## =========================================================

model_s6 <- glm(
  knee_pain_broad ~ cMetS_resid * sex + age_f1 + bmi_1,
  data   = df_s6,
  family = binomial
)

coef_s6 <- summary(model_s6)$coefficients
vc_s6   <- vcov(model_s6)

## =========================================================
## Sex-specific estimates
## =========================================================

beta_male_s6 <- coef_s6["cMetS_resid", "Estimate"]
se_male_s6   <- coef_s6["cMetS_resid", "Std. Error"]
p_male_s6    <- coef_s6["cMetS_resid", "Pr(>|z|)"]

or_male_s6  <- exp(beta_male_s6)
lcl_male_s6 <- exp(beta_male_s6 - 1.96 * se_male_s6)
ucl_male_s6 <- exp(beta_male_s6 + 1.96 * se_male_s6)

interaction_name_s6 <- "cMetS_resid:sexfemale"

beta_female_s6 <- coef_s6["cMetS_resid", "Estimate"] +
  coef_s6[interaction_name_s6, "Estimate"]

var_female_s6 <- vc_s6["cMetS_resid", "cMetS_resid"] +
  vc_s6[interaction_name_s6, interaction_name_s6] +
  2 * vc_s6["cMetS_resid", interaction_name_s6]

se_female_s6 <- sqrt(var_female_s6)
z_female_s6  <- beta_female_s6 / se_female_s6
p_female_s6  <- 2 * pnorm(-abs(z_female_s6))

or_female_s6  <- exp(beta_female_s6)
lcl_female_s6 <- exp(beta_female_s6 - 1.96 * se_female_s6)
ucl_female_s6 <- exp(beta_female_s6 + 1.96 * se_female_s6)

## =========================================================
## Interaction
## =========================================================

beta_int_s6 <- coef_s6[interaction_name_s6, "Estimate"]
se_int_s6   <- coef_s6[interaction_name_s6, "Std. Error"]
p_int_s6    <- coef_s6[interaction_name_s6, "Pr(>|z|)"]

or_int_s6  <- exp(beta_int_s6)
lcl_int_s6 <- exp(beta_int_s6 - 1.96 * se_int_s6)
ucl_int_s6 <- exp(beta_int_s6 + 1.96 * se_int_s6)

interaction_text_s6 <- paste0(
  "Measure of interaction on multiplicative scale: ratio of ORs (95% CI) = ",
  fmt_or_ci_s6(or_int_s6, lcl_int_s6, ucl_int_s6),
  "; P = ",
  fmt_p_s6(p_int_s6)
)

## =========================================================
## Table
## =========================================================

tab_s6 <- df_s6 %>%
  group_by(sex, knee_pain_broad) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(
    names_from  = knee_pain_broad,
    values_from = n,
    values_fill = 0
  ) %>%
  rename(no_pain = `0`, pain = `1`) %>%
  mutate(
    Sex = ifelse(sex == "male", "Male", "Female"),
    `Knee pain, n` = pain,
    `No knee pain, n` = no_pain,
    `OR per 1-unit increase in cMetS (95% CI)` = ifelse(
      sex == "male",
      fmt_or_ci_s6(or_male_s6, lcl_male_s6, ucl_male_s6),
      fmt_or_ci_s6(or_female_s6, lcl_female_s6, ucl_female_s6)
    ),
    `P-value` = ifelse(
      sex == "male",
      fmt_p_s6(p_male_s6),
      fmt_p_s6(p_female_s6)
    )
  ) %>%
  select(
    Sex,
    `Knee pain, n`,
    `No knee pain, n`,
    `OR per 1-unit increase in cMetS (95% CI)`,
    `P-value`
  ) %>%
  arrange(factor(Sex, levels = c("Male", "Female")))

cat("\nSupplementary Table S6. BMI-adjusted sensitivity analysis\n\n")
print(tab_s6)

cat("\n", interaction_text_s6, "\n")
cat("Adjusted for: age, BMI\n")

## =========================================================
## Excel
## =========================================================

wb <- createWorkbook()
addWorksheet(wb, "Table S6")

writeData(wb, "Table S6",
          "Supplementary Table S6. BMI-adjusted sensitivity analysis",
          startRow = 1, startCol = 1)

writeData(wb, "Table S6", tab_s6, startRow = 3, startCol = 1)

saveWorkbook(wb, "Table_S6_BMI_sensitivity.xlsx", overwrite = TRUE)

cat("\nExcel file saved: Table_S6_BMI_sensitivity.xlsx\n")


# Sensitivity analysis: trimming of IPW (1–99%)


stopifnot(exists("df_ipw_primary"))
stopifnot(exists("df_ipw_additional"))

library(dplyr)

## --- 1. Create trimmed weights -------------------------------------

df_ipw_primary_trim <- df_ipw_primary |>
  mutate(
    ipw_soc_trim = pmin(
      pmax(ipw_soc, quantile(ipw_soc, 0.01, na.rm = TRUE)),
      quantile(ipw_soc, 0.99, na.rm = TRUE)
    )
  )

df_ipw_additional_trim <- df_ipw_additional |>
  mutate(
    ipw_main_trim = pmin(
      pmax(ipw_main, quantile(ipw_main, 0.01, na.rm = TRUE)),
      quantile(ipw_main, 0.99, na.rm = TRUE)
    )
  )

## --- 2. Fit models --------------------------------------------------

model_ipw_primary_trim <- glm(
  knee_pain_broad ~ cMetS_resid * sex + age_f,
  family = binomial,
  data = df_ipw_primary_trim,
  weights = ipw_soc_trim
)

model_ipw_additional_trim <- glm(
  knee_pain_broad ~ cMetS_resid * sex + age_f,
  family = binomial,
  data = df_ipw_additional_trim,
  weights = ipw_main_trim
)

## --- 3. Extract results --------------------------------------------

cat("Socio-clinical model (trimmed weights):\n")
print(extract_sex_specific(model_ipw_primary_trim))

cat("\nCardiometabolic model (trimmed weights):\n")
print(extract_sex_specific(model_ipw_additional_trim))



