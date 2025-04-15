# =============================================================================
# Title:    IV Analysis, Sensitivity Analysis, and Summary Tables for SSC Study
# Author:   Sangyun Lee
# Date:     2025-04-12
# Purpose:  Data preprocessing, OLS regression and 2SLS regression (IV), 
#           sensitivity analysis, and creation of summary tables and IV diagnostics.
# =============================================================================

# ---------------------------#
# 1. Required Libraries      #
# ---------------------------#
library(dplyr)
library(ivreg)         # For instrumental variable regression
library(sandwich)      # For robust covariance estimation
library(iv.sensemakr)  # For sensitivity analysis for IV
library(tidyr)         # For pivot_wider
library(knitr)         # For table printing

# ---------------------------#
# 2. Reproducibility         #
# ---------------------------#
set.seed(123)

# ---------------------------#
# 3. Data Preprocessing      #
# ---------------------------#
## Load datasets and select required variables
data1 <- read.csv("baseline survey.csv") %>%
  select(HID, MMAR, MWORK01, ARA01, HCUL, HSECU, MPGN04, MPGN05, MPGN06, MPGN08,
         CBIRTH_O, TWIN)

data2 <- read.csv("wave 1 main survey.csv") %>%
  select(HID, MPGN12, MBRF03, MBRF06_01, MBRF06_02, MBRF06_03)

## Merge datasets by HID, select and rename variables, and filter out cesarean responses
data <- merge(data1, data2, by = "HID") %>%
  select(MPGN12, MBRF03, MBRF06_01, MBRF06_02, MBRF06_03,
         MMAR, MWORK01, ARA01, HCUL, HSECU, MPGN04, MPGN05, MPGN06, MPGN08, CBIRTH_O, TWIN) %>%
  rename(
    Z = MPGN12, A = MBRF03, 
    Y1 = MBRF06_01, Y2 = MBRF06_02, Y3 = MBRF06_03,
    C1 = MMAR, C2 = MWORK01, C3 = ARA01, 
    C4 = HCUL, C5 = HSECU, C6 = MPGN04, 
    C7 = MPGN05, C8 = MPGN06, C9 = MPGN08, 
    C10 = CBIRTH_O, C11 = TWIN
  ) %>%
  filter(Z != 2)  # Exclude cesarean delivery responses

# ---------------------------#
# 4. Recode Variables        #
# ---------------------------#
## Binary recoding for IV, treatment, and outcomes (0/1)
data <- data %>%
  mutate(
    Z = ifelse(Z == 1, 1, 0),
    A = ifelse(A == 1, 1, 0),
    Y1 = ifelse(Y1 == 1, 1, 0),
    Y2 = ifelse(Y2 == 1, 1, 0),
    Y3 = ifelse(Y3 == 1, 1, 0)
  )

## Recode covariates and convert to factors
data <- data %>%
  mutate(
    C1 = ifelse(C1 %in% c(1), 1, 0),
    C2 = ifelse(C2 %in% c(1, 2), 1, 0),
    C7 = ifelse(C7 %in% c(1), 1, 0),
    C8 = ifelse(C8 %in% c(1), 1, 0),
    C9 = ifelse(C9 %in% c(1), 1, 0),
    C10 = ifelse(C10 %in% c(4, 5), 4, C10),
    across(c(C1, C2, C3, C4, C5, C6, C7, C8, C9, C10, C11), as.factor)
  )

# ---------------------------#
# 5. Table 1: Summary Statistics (Baseline Characteristics & Outcomes) #
# ---------------------------#

## Define variable labels
var_label <- c(
  "C1"  = "Maternal marital status",
  "C2"  = "Maternal employment or education status",
  "C3"  = "Urbanicity",
  "C4"  = "Multicultural household status",
  "C5"  = "Low-income household status",
  "C6"  = "ART use",
  "C7"  = "Pre-pregnancy vaccination status",
  "C8"  = "Pre-pregnancy prenatal exam status",
  "C9"  = "Prenatal education experience status",
  "C10" = "Birth order of the child",
  "C11" = "Multiple birth status",
  "Y1"  = "1 month Exclusive breastfeeding",
  "Y2"  = "2 month Exclusive breastfeeding",
  "Y3"  = "3 month Exclusive breastfeeding"
)

vars <- names(var_label)

## Calculate group sizes for the SSC (treatment) variable A
n0 <- sum(data$A == 0, na.rm = TRUE)
n1 <- sum(data$A == 1, na.rm = TRUE)

## Loop over each variable to compute frequency, percentage, and chi-square p-value
summary_list <- list()

for (v in vars) {
  tbl <- table(data[[v]], data$A)
  p_val <- suppressWarnings(chisq.test(tbl)$p.value)
  
  df <- as.data.frame(tbl)
  colnames(df) <- c("Level", "SSC", "Freq")
  
  df <- df %>%
    group_by(SSC) %>%
    mutate(Percent = round(Freq / sum(Freq) * 100, 2)) %>%
    ungroup()
  
  # For binary variables (0/1), convert to "Yes"/"No" with order: Yes then No
  unique_vals <- sort(unique(df$Level))
  if (all(unique_vals %in% c("0", "1")) && length(unique_vals) == 2) {
    df <- df %>%
      mutate(Level = ifelse(Level == "1", "Yes", "No"),
             Level = factor(Level, levels = c("Yes", "No")))
  } else {
    df$Level <- as.character(df$Level)
  }
  
  df <- df %>%
    mutate(FreqPct = paste0(Freq, " (", Percent, "%)"))
  
  df_wide <- df %>%
    select(Level, SSC, FreqPct) %>%
    pivot_wider(names_from = SSC, values_from = FreqPct) %>%
    mutate(Variable = var_label[[v]],
           p_value = round(p_val, 3))
  
  df_wide <- df_wide %>%
    arrange(Level) %>%
    mutate(Variable = ifelse(row_number() == 1, Variable, ""),
           p_value = ifelse(row_number() == 1, as.character(p_value), ""))
  
  summary_list[[v]] <- df_wide
}

final_table <- bind_rows(summary_list) %>%
  rename(
    "Characteristics and Outcomes" = Variable,
    "No skin-to-skin (n = 0)" = `0`,
    "Skin-to-skin (n = 1)" = `1`,
    "p-value" = p_value
  )

colnames(final_table) <- gsub("\\(n = 0\\)",
                              paste0("(n = ", n0, ")"), colnames(final_table))
colnames(final_table) <- gsub("\\(n = 1\\)",
                              paste0("(n = ", n1, ")"), colnames(final_table))

table1 <- final_table %>%
  select(`Characteristics and Outcomes`, everything())

print(table1, row.names = FALSE)

# ---------------------------#
# 6. Table 2: OLS and IV (2SLS) Estimates (Excluding Anderson-Rubin CI) #
# ---------------------------#

## Helper function for OLS models
run_ols <- function(outcome, covars = NULL) {
  formula <- as.formula(paste(outcome, "~ A", if (!is.null(covars)) paste("+", covars)))
  model <- lm(formula, data = data)
  list(
    coef = round(summary(model)$coef[2, c(1, 4)], 3),
    conf = round(confint(model)[2, ], 3)
  )
}

## Helper function for IV (2SLS) models
run_iv <- function(outcome, covars = NULL) {
  if (is.null(covars)) {
    formula <- as.formula(paste(outcome, "~ A | Z"))
  } else {
    formula <- as.formula(paste(outcome, "~ A +", covars, "| Z +", covars))
  }
  model <- ivreg(formula, data = data)
  list(
    coef = round(summary(model, vcov = sandwich, df = Inf)$coefficients[2, c(1, 4)], 3),
    conf = round(confint(model)[2, ], 3)
  )
}

# Specify covariates for adjusted models
covariates <- "C1+C2+C3+C4+C5+C6+C7+C8+C9+C10+C11"

# Run OLS models for outcomes Y1, Y2, Y3
ols_results <- lapply(c("Y1", "Y2", "Y3"), function(y) {
  list(
    unadjusted = run_ols(y),
    adjusted   = run_ols(y, covariates)
  )
})

# Run IV (2SLS) models for outcomes Y1, Y2, Y3
iv_results <- lapply(c("Y1", "Y2", "Y3"), function(y) {
  list(
    unadjusted = run_iv(y),
    adjusted   = run_iv(y, covariates)
  )
})

## Helper function to format estimates and 95% confidence intervals as a string
format_est_ci <- function(res) {
  paste0(res$coef[1], " [", res$conf[1], ", ", res$conf[2], "]")
}

table2 <- data.frame(
  Outcome = c("1 month", "2 month", "3 month"),
  `OLS (unadjusted) [95% CI]` = sapply(ols_results, function(x) format_est_ci(x$unadjusted)),
  `OLS (adjusted) [95% CI]`   = sapply(ols_results, function(x) format_est_ci(x$adjusted)),
  `2SLS (unadjusted) [95% CI]` = sapply(iv_results, function(x) format_est_ci(x$unadjusted)),
  `2SLS (adjusted) [95% CI]`   = sapply(iv_results, function(x) format_est_ci(x$adjusted)),
  stringsAsFactors = FALSE
)

knitr::kable(
  table2,
  align   = c("l", "c", "c", "c", "c"),
  caption = "Effect of Immediate Skin-to-Skin Contact on Exclusive Breastfeeding: OLS and 2SLS Estimates"
)

# =============================================================================
# Anderson-Rubin Confidence Intervals for IV Models (Outcomes: Y1, Y2, Y3)
# =============================================================================

# Helper function to compute Anderson–Rubin confidence intervals for a given outcome
compute_AR_CI <- function(outcome, covars = NULL) {
  # Extract the outcome, treatment, IV, and covariate matrix
  y <- data[[outcome]]
  d <- data$A
  z <- as.numeric(data$Z)
  x <- model.matrix(~ C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C10 + C11, data = data)
  
  # Fit the IV model using iv_fit (from iv.sensemakr)
  mod.fit <- iv_fit(y, d, z, x)
  
  # Extract the Anderson-Rubin confidence interval for the coefficient on A
  # Using confint() with method = "anderson-rubin" returns a vector with names "lwr" and "upr"
  ar_ci <- confint(mod.fit)
  
  # Convert the vector to numeric and round to three decimal places
  ar_ci_num <- round(as.numeric(ar_ci), 3)
  return(ar_ci_num)
}

# Define outcomes for which to compute AR CIs
outcomes <- c("Y1", "Y2", "Y3")

# Compute AR CIs for each outcome and store the results in a list
ar_results <- lapply(outcomes, compute_AR_CI)
names(ar_results) <- outcomes

# Create a data frame with the AR CI results
df_AR <- data.frame(
  Outcome = outcomes,
  `AR Confidence Interval` = sapply(ar_results, function(ci) paste0("[", ci[1], ", ", ci[2], "]")),
  stringsAsFactors = FALSE
)

knitr::kable(
  df_AR,
  align = c("l", "c"),
  caption = "Anderson-Rubin Confidence Intervals for IV Models",
  row.names = FALSE
)



# ---------------------------#
# 8. F-test Diagnostics for IV Models ---------------------------
# ---------------------------#

## Unadjusted IV model F-test
iv_mod_unadj <- ivreg(Y1 ~ A | Z, data = data)
iv_sum_unadj <- summary(iv_mod_unadj, diagnostics = TRUE)
weak_inst_unadj <- iv_sum_unadj$diagnostics["Weak instruments", ]
cat("Unadjusted IV Model - Weak Instruments:\n")
print(weak_inst_unadj)

## Adjusted IV model F-test
iv_mod_adj <- ivreg(Y1 ~ A + C1+C2+C3+C4+C5+C6+C7+C8+C9+C10+C11 |
                      Z + C1+C2+C3+C4+C5+C6+C7+C8+C9+C10+C11, data = data)
iv_sum_adj <- summary(iv_mod_adj, diagnostics = TRUE)
weak_inst_adj <- iv_sum_adj$diagnostics["Weak instruments", ]
cat("Adjusted IV Model - Weak Instruments:\n")
print(weak_inst_adj)

# ---------------------------#
# 9. Table 3: Sensitivity Analysis ---------------------------
# ---------------------------#

# Helper function for Sensitivity Analysis using iv.sensemakr
run_sensitivity <- function(outcome) {
  y <- data[[outcome]]
  d <- data$A
  z <- as.numeric(data$Z)
  x <- model.matrix(~ C1+C2+C3+C4+C5+C6+C7+C8+C9+C10+C11, data = data)
  mod.fit <- iv_fit(y, d, z, x)  # Assumes iv_fit is defined (e.g., from iv.sensemakr)
  IV.sens <- sensemakr(mod.fit, benchmark_covariates = c("C111"), kz = 1:3)
  round(summary(IV.sens)$bounds$iv[, c(4, 5)], 3)
}

# Run sensitivity analysis for outcomes Y1, Y2, Y3
sensitivity_results <- lapply(c("Y1", "Y2", "Y3"), run_sensitivity)

# Outcome labels
outcome_names <- c("1 month", "2 month", "3 month")

# Create a data frame to store the sensitivity analysis results
table3 <- data.frame(
  Outcome = character(),
  `1x adjusted 95% CI` = character(),
  `2x adjusted 95% CI` = character(),
  `3x adjusted 95% CI` = character(),
  stringsAsFactors = FALSE
)

for (i in seq_along(sensitivity_results)) {
  df_sens <- sensitivity_results[[i]]
  ci_1x <- paste0("(", df_sens$lwr[1], ", ", df_sens$upr[1], ")")
  ci_2x <- paste0("(", df_sens$lwr[2], ", ", df_sens$upr[2], ")")
  ci_3x <- paste0("(", df_sens$lwr[3], ", ", df_sens$upr[3], ")")
  
  table3 <- rbind(table3, data.frame(
    Outcome = outcome_names[i],
    `1x adjusted 95% CI` = ci_1x,
    `2x adjusted 95% CI` = ci_2x,
    `3x adjusted 95% CI` = ci_3x,
    stringsAsFactors = FALSE
  ))
}

knitr::kable(
  table3,
  align = c("l", "c", "c", "c"),
  caption = "Sensitivity Analysis for Violation of the Exclusion Restriction Assumption"
)
