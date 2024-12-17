# The Effect of Immediate Skin-to-Skin Contact on Exclusive Breastfeeding: An Instrumental Variable Approach  

This repository contains R code to analyze the causal effect of immediate skin-to-skin contact at birth on exclusive breastfeeding using Instrumental Variable (IV) methods. The analysis also includes a sensitivity analysis to assess the robustness of the IV estimates.  

---

## **Overview**  
The code performs the following steps:  
1. **Data Preprocessing**: Merges, cleans, and recodes the data.  
2. **Descriptive Analysis**: Generates summary statistics for covariates by treatment group.  
3. **Main Analysis**: Estimates the effect of skin-to-skin contact using:  
   - Ordinary Least Squares (OLS)  
   - Two-Stage Least Squares (2SLS) Instrumental Variable models  
4. **Sensitivity Analysis**: Evaluates the robustness of the IV results.

---

## **Data Availability**  

The data that support the findings of this study are openly available at:  
[https://kicce.re.kr/kececp](https://kicce.re.kr/kececp)

---

## **File Structure**  

- **main_analysis.R**: The complete R code for the analysis.  

---

## **Requirements**  
To execute this code, ensure the following R packages are installed:  

```r
install.packages(c("dplyr", "ivreg", "sandwich", "iv.sensemakr"))
