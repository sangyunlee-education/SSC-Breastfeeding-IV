# Immediate Skin-to-Skin Contact and Exclusive Breastfeeding: An Instrumental Variable Analysis

This repository contains the R code used for the data analysis in the study:
**"The Effect of Immediate Skin-to-Skin Contact on Exclusive Breastfeeding: An Instrumental Variable Approach"**

## Overview
This study investigates the causal effect of immediate skin-to-skin contact (SSC) after birth on exclusive breastfeeding rates using instrumental variable (IV) analysis. The repository provides the full analytical workflow, ensuring transparency and reproducibility of the results.

---

## Files and Directory Structure

- `data/`
  *(Not included)*: Placeholder for the dataset used in the study. The dataset is not shared due to privacy concerns.
- `scripts/`
  - **`data_preprocessing.R`**: Code for cleaning and preparing the raw data for analysis.
  - **`descriptive_statistics.R`**: Code for calculating summary statistics and visualizations.
  - **`iv_analysis.R`**: Implementation of the instrumental variable analysis.
  - **`sensitivity_analysis.R`**: Additional analyses to assess robustness of the results.
- `results/`
  *(Optional)*: Folder to save output files like regression tables and plots.

---

## How to Use

### Prerequisites
1. **R**: Version X.X.X or later
2. R Packages:
   - `tidyverse`
   - `AER`
   - `ggplot2`
   - Additional packages listed in `scripts/iv_analysis.R`

### Steps to Reproduce the Analysis
1. Clone this repository:
   ```bash
   git clone https://github.com/yourusername/SSC_IV_Breastfeeding.git
   ```
2. Set your working directory to the repository:
   ```R
   setwd("path_to_cloned_repo")
   ```
3. Run the scripts in the following order:
   - `data_preprocessing.R`
   - `descriptive_statistics.R`
   - `iv_analysis.R`
   - `sensitivity_analysis.R`

4. Output tables and plots will be saved in the `results/` folder.

---

## Citation
If you use this repository in your research or projects, please cite the associated study:

> **Author(s)**: [Sangyun Lee, Yongnam Kim]  
> **Title**: *The Effect of Immediate Skin-to-Skin Contact on Exclusive Breastfeeding: An Instrumental Variable Approach*  
> **Year**: 2025  
> **DOI**: [Add if available]

---

## Notes and Limitations

- **Data Availability**: The original dataset is proprietary and cannot be shared. You may use a similar dataset to replicate the workflow.
- **Ethical Considerations**: All analyses were conducted in compliance with ethical guidelines.

---

## Contact
For questions or collaborations, please contact:
[Sangyun Lee]  
[lee10260@snu.ac.kr]
