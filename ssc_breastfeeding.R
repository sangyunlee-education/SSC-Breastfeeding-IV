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
  select(HID, MMAR, MWORK01, ARA01, HCUL, HSECU, MPGN04, MPGN05, MPGN06, MPGN08,
         CBIRTH_O, TWIN)

data2 <- read.csv("wave 1 main survey.csv") %>%
  select(HID, MPGN12, MBRF03, MBRF06_01, MBRF06_02, MBRF06_03)

## Merge datasets by HID and rename variables
data <- merge(data1, data2, by = "HID") %>%
  select(MPGN12, MBRF03, MBRF06_01, MBRF06_02, MBRF06_03,
         MMAR, MWORK01, ARA01, HCUL, HSECU, MPGN04, 
         MPGN05, MPGN06, MPGN08, CBIRTH_O, TWIN) %>%
  rename(Z = MPGN12, A = MBRF03, 
         Y1 = MBRF06_01, Y2 = MBRF06_02, Y3 = MBRF06_03,
         C1 = MMAR, C2 = MWORK01, C3 = ARA01, 
         C4 = HCUL, C5 = HSECU, C6 = MPGN04, C7 = MPGN05, C8 = MPGN06, 
         C9 = MPGN08, C10 = CBIRTH_O, C11 = TWIN) %>%
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
    C1 = ifelse(C1 %in% c(1), 1, 0),
    C2 = ifelse(C2 %in% c(1, 2), 1, 0),
    C7 = ifelse(C7 %in% c(1), 1, 0),
    C8 = ifelse(C8 %in% c(1), 1, 0),
    C9 = ifelse(C9 %in% c(1), 1, 0),
    C10 = ifelse(C10 %in% c(4, 5), 4, C10),
    across(c(C1, C2, C3, C4, C5, C6, C7, C8, C9, C10, C11), as.factor)
  )


# Table 1: Summary Statistics ---------------------------------------------

table_data <- data %>%
  group_by(A) %>%
  summarise(
    Total = n(),
    Marital = sum(as.numeric(as.character(C1)), na.rm = TRUE),
    Marital_Ratio = Marital / Total,
    Employment_study = sum(as.numeric(as.character(C2)), na.rm = TRUE),
    Employment_study_Ratio = Employment_study / Total,
    Large_city = sum(C3 == 1, na.rm = TRUE),
    Large_city_Ratio = Large_city / Total,
    Mediumsmall_city = sum(C3 == 2, na.rm = TRUE),
    Mediumsmall_city_Ratio = Mediumsmall_city / Total,
    Rural_area = sum(C3 == 3, na.rm = TRUE),
    Rural_area_Ratio = Rural_area / Total,
    Multicultural = sum(as.numeric(as.character(C4)), na.rm = TRUE),
    Multicultural_Ratio = Multicultural / Total,
    Lowincome = sum(as.numeric(as.character(C5)), na.rm = TRUE),
    Lowincome_Ratio = Lowincome / Total,
    Natural_conception = sum(C6 == 1, na.rm = TRUE),
    Natural_conception_Ratio = Natural_conception / Total,
    ART_Os_OI = sum(C6 == 2, na.rm = TRUE),
    ART_Os_OI_Ratio = ART_Os_OI / Total,
    ART_IUI_IVF = sum(C6 == 3, na.rm = TRUE),
    ART_IUI_IVF_Ratio = ART_IUI_IVF / Total,
    Preconception_vaccination = sum(as.numeric(as.character(C7)), na.rm = TRUE),
    Preconception_vaccination_Ratio = Preconception_vaccination / Total,
    Preconception_check = sum(as.numeric(as.character(C8)), na.rm = TRUE),
    Preconception_check_Ratio = Preconception_check / Total,
    Participation_education = sum(as.numeric(as.character(C9)), na.rm = TRUE),
    Participation_education_Ratio = Participation_education / Total,
    Firstborn = sum(C10 == 1, na.rm = TRUE),
    Firstborn_Ratio = Firstborn / Total,
    Second_order = sum(C10 == 2, na.rm = TRUE),
    Second_order_Ratio = Second_order / Total,
    Third_order = sum(C10 == 3, na.rm = TRUE),
    Third_order_Ratio = Third_order / Total,
    Fourth_or_more = sum(C10 == 4, na.rm = TRUE),
    Fourth_or_more_Ratio = Fourth_or_more / Total,
    Twins_multiples = sum(as.numeric(as.character(C11)), na.rm = TRUE),
    Twins_multiples_Ratio = Twins_multiples / Total
  ) %>%
  mutate(A = ifelse(A == 0, "No skin-to-skin", "Skin-to-skin"))


## Calculate p-values for covariates
p_values <- list(
  Marital = chisq.test(data$C1, data$A)$p.value,
  Employment_study = chisq.test(data$C2, data$A)$p.value,
  Urbanicity = chisq.test(data$C3, data$A)$p.value,
  Multicultural = chisq.test(data$C4, data$A)$p.value,
  Low_income = chisq.test(data$C5, data$A)$p.value,
  ART = chisq.test(data$C6, data$A)$p.value,
  Preconception_vaccination = chisq.test(data$C7, data$A)$p.value,
  Preconception_check = chisq.test(data$C8, data$A)$p.value,
  participation_education = chisq.test(data$C9, data$A)$p.value,
  Birth_order = chisq.test(data$C10, data$A)$p.value,
  Twins = chisq.test(data$C11, data$A)$p.value
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
covariates <- "C1+C2+C3+C4+C5+C6+C7+C8+C9+C10+C11"
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
  x <- model.matrix(~ C1+C2+C3+C4+C5+C6+C7+C8+C9+C10+C11, data = data)
  mod.fit <- iv_fit(y, d, z, x)
  IV.sens <- sensemakr(mod.fit, benchmark_covariates = c("C111"), kz = 1:3)
  round(summary(IV.sens)$bounds$iv[, c(4, 5)], 3)
}

sensitivity_results <- lapply(c("Y1", "Y2", "Y3"), run_sensitivity)
