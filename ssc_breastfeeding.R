# Required Libraries ------------------------------------------------------
library(dplyr)
library(ivreg)
library(sandwich)
library(iv.sensemakr)

# Reproducibility ---------------------------------------------------------
set.seed(123)

# Data Preprocessing ------------------------------------------------------

## Load datasets
data1 <- read.csv("baseline survey.csv") %>%
  select(HID, MAGE, CBIRTH_O, MPGN04, MPGN06, MPGN08, ARA01, MEDU, TWIN)

data2 <- read.csv("wave 1 main survey.csv") %>%
  select(HID, MPGN12, MBRF03, MBRF06_01, MBRF06_02, MBRF06_03)

## Merge datasets by HID and rename variables
data <- merge(data1, data2, by = "HID") %>%
  select(MPGN12, MBRF03, MBRF06_01, MBRF06_02, MBRF06_03,
         MAGE, CBIRTH_O, MPGN04, MPGN06, MPGN08, ARA01, MEDU, TWIN) %>%
  rename(Z = MPGN12, A = MBRF03, 
         Y1 = MBRF06_01, Y2 = MBRF06_02, Y3 = MBRF06_03,
         C1 = MAGE, C2 = CBIRTH_O, C3 = MPGN04, C4 = MPGN06, 
         C5 = MPGN08, C6 = ARA01, C7 = MEDU, C8 = TWIN) %>%
  filter(Z != 2) ### Filter out cesarean delivery responses

# Recode Variables --------------------------------------------------------

## Binary recoding for IV, treatment, outcomes
data <- data %>%
  mutate(
    Z = ifelse(Z == 1, 1, 0),
    A = ifelse(A == 1, 1, 0),
    Y1 = ifelse(Y1 == 1, 1, 0),
    Y2 = ifelse(Y2 == 1, 1, 0),
    Y3 = ifelse(Y3 == 1, 1, 0)
  )

## Recode covariates
data <- data %>%
  mutate(
    C2 = ifelse(C2 %in% c(4, 5), 4, C2),
    C7 = case_when(
      C7 %in% c(2, 3) ~ 1,
      C7 == 4 ~ 2,
      C7 %in% c(5, 6) ~ 3,
      TRUE ~ 4
    ),
    across(c(C2, C3, C4, C5, C6, C7, C8), as.factor)
  )

# Table 1: Summary Statistics ---------------------------------------------

table_data <- data %>%
  group_by(A) %>%
  summarise(
    Maternal_age_mean = mean(C1, na.rm = TRUE),
    First_birth_order = sum(C2 == 1, na.rm = TRUE),
    Second_birth_order = sum(C2 == 2, na.rm = TRUE),
    Third_birth_order = sum(C2 == 3, na.rm = TRUE),
    Fourth_or_more = sum(as.numeric(C2) >= 4, na.rm = TRUE),
    Natural_conception = sum(C3 == 0, na.rm = TRUE),
    ART_IUI_IVF = sum(C3 == 1, na.rm = TRUE),
    Childbirth_education = sum(C4 == 1, na.rm = TRUE),
    Prenatal_testing = sum(C5 == 1, na.rm = TRUE),
    Metropolitan = sum(C6 == 1, na.rm = TRUE),
    Small_medium_city = sum(C6 == 2, na.rm = TRUE),
    Town_rural = sum(C6 == 3, na.rm = TRUE),
    Less_than_high_school = sum(C7 == 1, na.rm = TRUE),
    High_school = sum(C7 == 2, na.rm = TRUE),
    College = sum(C7 == 3, na.rm = TRUE),
    Graduate_school = sum(C7 == 4, na.rm = TRUE),
    Twin_pregnancy = sum(C8 == 1, na.rm = TRUE)
  ) %>%
  mutate(A = ifelse(A == 0, "No skin-to-skin", "Skin-to-skin"))

## Calculate p-values for covariates
p_values <- list(
  Age = t.test(C1 ~ A, data = data)$p.value,
  Birth_order = chisq.test(data$C2, data$A)$p.value,
  ART = chisq.test(data$C3, data$A)$p.value,
  Childbirth_education = chisq.test(data$C4, data$A)$p.value,
  Prenatal_testing = chisq.test(data$C5, data$A)$p.value,
  Urbanization = chisq.test(data$C6, data$A)$p.value,
  Education = chisq.test(data$C7, data$A)$p.value,
  Twin = chisq.test(data$C8, data$A)$p.value
)

# Main Analysis: Table 2 --------------------------------------------------

## Helper function for OLS models
run_ols <- function(outcome, covars = NULL) {
  formula <- as.formula(paste(outcome, "~ A", if (!is.null(covars)) paste("+", covars)))
  model <- lm(formula, data = data)
  list(
    coef = round(summary(model)$coef[2, c(1, 4)], 3),
    conf = round(confint(model)[2, ], 3)
  )
}

## Run OLS Models
covariates <- "C1+C2+C3+C4+C5+C6+C7+C8"
ols_results <- lapply(c("Y1", "Y2", "Y3"), function(y) {
  list(
    unadjusted = run_ols(y),
    adjusted = run_ols(y, covariates)
  )
})

## Helper function for IV models
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


## Run IV Models
iv_results <- lapply(c("Y1", "Y2", "Y3"), function(y) {
  list(
    unadjusted = run_iv(y),
    adjusted = run_iv(y, covariates)
  )
})

# Sensitivity Analysis: Table 3 -------------------------------------------

run_sensitivity <- function(outcome) {
  y <- data[[outcome]]
  d <- data$A
  z <- as.numeric(data$Z)
  x <- model.matrix(~ C1+C2+C3+C4+C5+C6+C7+C8, data = data)
  mod.fit <- iv_fit(y, d, z, x)
  IV.sens <- sensemakr(mod.fit, benchmark_covariates = c("C81"), kz = 1:3)
  round(summary(IV.sens)$bounds$iv[, c(4, 5)], 3)
}

sensitivity_results <- lapply(c("Y1", "Y2", "Y3"), run_sensitivity)
