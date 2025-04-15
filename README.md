# The Effect of Immediate Skin-to-Skin Contact on Exclusive Breastfeeding: An Instrumental Variable Approach

This repository contains R code to analyze the causal effect of immediate skin-to-skin contact (SSC) at birth on exclusive breastfeeding outcomes using Instrumental Variable (IV) methods. The analysis also includes a sensitivity analysis to assess the robustness of the IV estimates under hypothetical violations of the exclusion restriction.

---

## **Overview**  
The code performs the following steps:  
1. **Data Preprocessing**: Merges, cleans, and recodes the data.  
2. **Descriptive Analysis**: Generates summary statistics for baseline characteristics and outcomes by treatment group.  
3. **Main Analysis**: Estimates the effect of SSC on exclusive breastfeeding using:  
   - Ordinary Least Squares (OLS) regression  
   - Two-Stage Least Squares (2SLS) IV models  
4. **Sensitivity Analysis**: Evaluates the robustness of the IV estimates under potential violations of the exclusion restriction.

The IV approach is used to address possible unobserved confounding, and sensitivity analysis further supports the robustness of the findings.

---

## **Data Availability**  

The data that support the findings of this study are openly available at:  
[https://kicce.re.kr/kececp](https://kicce.re.kr/kececp)

---

## **File Structure**  

- **ssc_breastfeeding.R**: The complete R code for data preprocessing, regression analyses, sensitivity analysis, and table generation.

---

## **Requirements**  

To execute this code, ensure the following R packages are installed:

```r
install.packages(c("dplyr", "ivreg", "sandwich", "iv.sensemakr", "tidyr", "knitr"))
