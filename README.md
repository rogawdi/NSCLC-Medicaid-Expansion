# NSCLC-Medicaid-Expansion
Generalized Difference-in-Differences Analysis with Propensity Score Matching from the SEER database

This repository contains R code and related materials for the study:
"Effects of State-Wide Medicaid Expansion on Resectable Non–Small Cell Lung Cancer Survival"

This study used SEER data (2006–2019) and a quasi-experimental, generalized difference-in-differences approach to evaluate whether state-level Medicaid expansion was associated with improved 2-year survival and early-stage diagnosis among patients with resectable non-small cell lung cancer (NSCLC).

Content:
  1. Prepares SEER data for analysis (state-based grouping, filtering, formatting, variable creation)
  2. Descriptive Statistics of cleaned data
  3. Kaplan-Meier survival analysis of pre-matched expansion groups with log rank and Benjamini-Hochberg adjusted pairwise comparisons
  4. Generalized difference-in-differences analysis of pre-matched populations for early post-implementation, late post-implementation, and total time for each expansion against control
  5. Propensity score modeling and matching by expansion group and era with love plots and pre-/post-matching propensity disribution graphs
  6.  Cox proportional hazards models of matched groups with generalized difference-in-differences models for survival
  7.  Placebo Falsification Modeling for Survival across entire cohort - to verify parallel trends assumption
  8.  Kaplan-Meier survival analysis of post-matching expansion groups log rank and Benjamini-Hochberg adjusted pairwise comparisons
  9.  Stage at diagnosis of matched cohorts with placebo falsification modeling - to verify parallel trends assumption
  10.  State map plot
	•	Functions: "survival", "survminer" for Kaplan-Meier analyses. "Matching", "MatchIt" for propensity score matching, "lmtest" for linear model assumptions, "usmap" for map generation. Other custom R functions for matching, covariate balance assessment, and plotting. 
	•	README.md: This file.
	•	LICENSE: MIT License (or other license of your choice).
This study used the SEER (Surveillance, Epidemiology, and End Results) public-use dataset. Access to SEER data is available upon request and approval from the National Cancer Institute: https://seer.cancer.gov/data/

Note: No patient-level data are included in this repository. Only code and summary-level outputs are shared.

Reproducibility

To reproduce the results:
	1.	Register and export SEER Research Plus data (https://seer.cancer.gov/data/) as a csv file with included variables from 
	2.	Import csv file and format  data to match the necessary structure using step 1. 
	3.	Run each step of the code sequentially, or modify paths and parameters as needed. Create necessary export directories as they appear. 
	4.	Use the provided plotting scripts to replicate survival curves, covariate balance plots, and tables.

Citation
If you use this code or build upon this work, please cite our paper (once published) and acknowledge this repository.

Rohin Gawdi. Effects of State-Wide Medicaid Expansion on Resectable Non–Small Cell Lung Cancer Survival. [Journal once available]. 2025. [DOI once available]
