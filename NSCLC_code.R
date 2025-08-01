# Load necessary libraries
library(survival)
library(survminer)
library(dplyr)
library(broom)
library(Matching)
library(MatchIt)
library(nnet)
library(cobalt)
library(data.table)
library(rgenoud)
library(tidyr)
library(gridExtra)
library(grid)
library(ggplot2)
library(broom)
library(scales)
library(sandwich)
library(lmtest)
library(purrr)
library(ggpubr)
library(usmap)

### 1. Variable Coding ###

# Read in the dataset
seer_final <- read.csv("inputs/nsclc_stage3a_2.csv")


# Define Medicaid expansion categories
seer_final$Medicaid.Expansion.Status <- tolower(seer_final$Medicaid.Expansion.Status)

seer_final$Medicaid_Expansion_Status <- as.factor(seer_final$Medicaid.Expansion.Status)

seer_final <- seer_final %>%
  mutate(dx.yr_grouped = case_when(
    dx.yr >= 2006 & dx.yr <= 2010 ~ "2006-2010",
    dx.yr >= 2011 & dx.yr <= 2013 ~ "2011-2013",
    dx.yr >= 2014 & dx.yr <= 2016 ~ "2014-2016",
    dx.yr >= 2017 & dx.yr <= 2019 ~ "2017-2019",
    TRUE ~ NA_character_
  ))
seer_final$dx.yr_grouped <- factor(seer_final$dx.yr_grouped)

seer_final <- seer_final %>%
  mutate(
    Medicaid_Expansion_Status = case_when(
      # 2006-2010: Assign pre-expansion groups based on future expansion
      dx.yr_grouped == "2006-2010" & Medicaid.Expansion.Status %in% c("california", "connecticut", "new jersey", "seattle (puget sound)") ~ 1,  # future EARLY expander
      dx.yr_grouped == "2006-2010" & Medicaid.Expansion.Status %in% c("iowa", "hawaii", "kentucky", "new Mexico") ~ 2,  # future ON TIME expander
      dx.yr_grouped == "2006-2010" & Medicaid.Expansion.Status %in% c("louisiana", "alaska natives") ~ 3,  # future LATE expander
      dx.yr_grouped == "2006-2010" ~ 0,  # Never expanders stay 0
      
      # 2010-2013: First expansion
      dx.yr_grouped == "2011-2013" & Medicaid.Expansion.Status %in% c("california", "connecticut", "new jersey", "seattle (puget sound)") ~ 1,  # NEW expander
      dx.yr_grouped == "2011-2013" & Medicaid.Expansion.Status %in% c("iowa", "hawaii", "kentucky", "new Mexico") ~ 2,  # future expander
      dx.yr_grouped == "2011-2013" & Medicaid.Expansion.Status %in% c("louisiana", "alaska natives") ~ 3,  # future expander
      dx.yr_grouped == "2011-2013" ~ 0,  
      
      # 2014-2016: Second expansion
      dx.yr_grouped == "2014-2016" & Medicaid.Expansion.Status %in% c("iowa", "hawaii", "kentucky", "new Mexico") ~ 2,  #new expander
      dx.yr_grouped == "2014-2016" & Medicaid.Expansion.Status %in% c("california", "connecticut", "new jersey", "seattle (puget sound)") ~ 1,  #already expanded
      dx.yr_grouped == "2014-2016" & Medicaid.Expansion.Status %in% c("louisiana", "alaska natives") ~ 3,  # future expander
      dx.yr_grouped == "2014-2016" ~ 0,  
      
      # 2017-2019: Final expansion
      dx.yr_grouped == "2017-2019" & Medicaid.Expansion.Status %in% c("louisiana", "alaska natives") ~ 3,  # new expander
      dx.yr_grouped == "2017-2019" & Medicaid.Expansion.Status %in% c("iowa", "hawaii", "kentucky", "new Mexico") ~ 2,  #already expanded
      dx.yr_grouped == "2017-2019" & Medicaid.Expansion.Status %in% c("california", "connecticut", "new jersey", "seattle (puget sound)") ~ 1,  #already expanded
      dx.yr_grouped == "2017-2019" ~ 0,  
      
      TRUE ~ 0  # Default: No expansion
    )
  )

# seer_final$Medicaid.Expansion.Status <- relevel(seer_final$Medicaid.Expansion.Status, ref = "0")



# Check if mapping was successful:
table(seer_final$Medicaid_Expansion_Status)

# Convert categorical variables into factors
seer_final$Age <- factor(seer_final$Age)  # Age grouped by 5 years
seer_final$Combined.Stage <- factor(seer_final$Combined.Stage)
stage_map <- data.frame(
  stage_range = c("Stage I", "Stage IIA", "Stage IIB", "Stage IIIA"),
  stage_group = c("I", "II+III", "II+III", "II+III")
) 

# Merge SEER dataset with age categories
seer_final <- seer_final %>%
  left_join(stage_map, by = c("Combined.Stage" = "stage_range"))

# # View new age categories
table(seer_final$stage_group, useNA = "ifany")

seer_final$Sex <- factor(seer_final$Sex, levels = c("Male", "Female"))
seer_final$Median.household.income.inflation.adj.to.2022 <- factor(seer_final$Median.household.income.inflation.adj.to.2022)  # Median household 
seer_final$Rural.Urban.Continuum.Code <- factor(seer_final$Rural.Urban.Continuum.Code)  # Rural/urban continuum
#seer_final <- seer_final %>%
#  filter(!grepl("Unknown/missing/no match", trimws(Rural.Urban.Continuum.Code)))
seer_final$Marital.status.at.diagnosis <- factor(seer_final$Marital.status.at.diagnosis)  # Marital status
# seer_final$Marital.status.at.diagnosis <- relevel(seer_final$Marital.status.at.diagnosis, ref = "Married (including common law)") # make Married the reference

seer_final$Race.and.origin.recode..NHW..NHB..NHAIAN..NHAPI..Hispanic. <- factor(seer_final$Race.and.origin.recode..NHW..NHB..NHAIAN..NHAPI..Hispanic.)  # Race/ethnicity
# seer_final$Race.and.origin.recode..NHW..NHB..NHAIAN..NHAPI..Hispanic. <- relevel(seer_final$Race.and.origin.recode..NHW..NHB..NHAIAN..NHAPI..Hispanic., ref = "Non-Hispanic White") # make NH White the reference


# Group Radiation recodes
seer_final$Radiation.recode <- factor(seer_final$Radiation.recode)  # Radiation
rad_map <- data.frame(
  rad_range = c("Beam radiation", 
                "Combination of beam with implants or isotopes", 
                "None/Unknown",
                "Radiation, NOS  method or source not specified",
                "Radioactive implants (includes brachytherapy) (1988+)",
                "Radioisotopes (1988+)","Recommended, unknown if administered",
                "Refused (1988+)"),
  rad_group = c("Radiation recieved", "Radiation recieved","none/unknown/refused",
                "Radiation recieved",
                "Radiation recieved", "Radiation recieved",
                "none/unknown/refused","none/unknown/refused")
)
# Merge SEER dataset with radiation categories
seer_final <- seer_final %>%
  left_join(rad_map, by = c("Radiation.recode" = "rad_range"))

# # View new rad categories
table(seer_final$rad_group, useNA = "ifany")

# seer_final$rad_group <- relevel(seer_final$rad_group, ref = "none/unknown/refused") # should radiation performed be baseline? mmm actually no i think should not be is baseline since binary. 


# Group Surgery recodes
seer_final$Surgery.Performed <- factor(seer_final$Surgery.Performed)  # Surgery Status

seer_final$Surgery.Performed <- factor(seer_final$Surgery.Performed)  # Surgery
surg_map <- data.frame(
  surg_range = c("Surgery performed", 
                 "Not performed, patient died prior to recommended surgery", 
                 "Not recommended",
                 "Not recommended, contraindicated due to other cond; autopsy only (1973-2002)",
                 "Recommended but not performed, patient refused",
                 "Recommended but not performed, unknown reason",
                 "Recommended, unknown if performed",
                 "Unknown; death certificate; or autopsy only (2003+)"),
  surgery_status = c("Recommended + performed", 
                     "Recommended, not performed",
                     "Not recommended/unknown",
                     "Not recommended/unknown",
                     "Recommended, not performed",
                     "Recommended, not performed",
                     "Recommended, not performed",
                     "Not recommended/unknown")
)

# Merge SEER dataset with radiation categories
seer_final <- seer_final %>%
  left_join(surg_map, by = c("Surgery.Performed" = "surg_range"))

# Convert surgery_status to a factor if it's not already
seer_final$surgery_status <- factor(seer_final$surgery_status)

# Now, relevel the factor
# seer_final$surgery_status <- relevel(seer_final$surgery_status, ref = "Recommended + performed")

# # View new age categories
table(seer_final$surgery_status, useNA = "ifany")



# Create a binary event indicator (1 = death, 0 = alive)
seer_final$Event <- ifelse(seer_final$Vital.status.recode..study.cutoff.used. == "Dead", 1, 0)

# Convert Time into numeric variable
seer_final$Survival.months <- as.numeric(as.character(seer_final$Survival.months))
seer_final <- seer_final[!is.na(seer_final$Survival.months), ]

# INCOME QUINTILE
# Define income ranges and assign midpoints
income_map <- data.frame(
  income_range = c("< $40,000", "$40,000 - $44,999", "$45,000 - $49,999",
                   "$50,000 - $54,999", "$55,000 - $59,999", "$60,000 - $64,999",
                   "$65,000 - $69,999", "$70,000 - $74,999", "$75,000 - $79,999",
                   "$80,000 - $84,999", "$85,000 - $89,999", "$90,000 - $94,999",
                   "$95,000 - $99,999", "$100,000 - $109,999", "$110,000 - $119,999",
                   "$120,000+"),
  midpoint = c(35000, 42500, 47500, 52500, 57500, 62500, 67500, 72500, 77500, 
               82500, 87500, 92500, 97500, 105000, 115000, 125000)
)

seer_final <- seer_final %>%
  left_join(income_map, by = c("Median.household.income.inflation.adj.to.2022" = "income_range"))

quintile_cutoffs <- quantile(seer_final$midpoint, probs = seq(0, 1, 0.2), na.rm = TRUE)

seer_final <- seer_final %>%
  mutate(income_quintile = cut(midpoint, 
                               breaks = c(-Inf, quintile_cutoffs[-1]), 
                               labels = c("Q1", "Q2", "Q3", "Q4", "Q5"),
                               include.lowest = TRUE))

# seer_final$income_quintile <- relevel(seer_final$income_quintile, ref = "Q3")

seer_final <- seer_final %>%
  filter(!grepl("<NA>", trimws(income_quintile)))

# # View result
table(seer_final$income_quintile, useNA = "ifany")


### State by State population income quintiles ### 

# 1. Create a new variable for state-specific income quintiles
seer_final <- seer_final %>%
  group_by(Medicaid.Expansion.Status) %>%
  mutate(
    state_income_quintile = ntile(midpoint, 5)  # Split income into quintiles within each state
  ) %>%
  ungroup()

# 2. Create a 'low_income' flag for the lowest 2 quintiles within each state
seer_final <- seer_final %>%
  mutate(
    low_income_flag = ifelse(state_income_quintile %in% c(1, 2), 1, 0)
  )

# 3. Optional: If you want to see the quintile ranges by state
state_income_ranges <- seer_final %>%
  filter(!is.na(midpoint)) %>%  # Exclude NA values
  group_by(Medicaid.Expansion.Status, state_income_quintile) %>%
  summarize(
    min_income = min(midpoint, na.rm = TRUE),
    max_income = max(midpoint, na.rm = TRUE),
    .groups = "drop"
  )

# 5. View or export the quintile cutoffs table
print(state_quintile_cutoffs)
# write.csv(state_quintile_cutoffs, "state_income_quintile_cutoffs.csv", row.names = FALSE)


# AGE GROUPS
# Create a mapping of age groups
age_map <- data.frame(
  age_range = c("20-24 years", "25-29 years", "30-34 years", "35-39 years",
                "40-44 years", "45-49 years", "50-54 years", "55-59 years", "60-64 years"),
  age_category = c("Young Adult", "Young Adult", "Young Adult", "Young Adult",
                   "Young Adult", "Young Adult", "Young Adult", "Middle-Aged",
                   "Older Adult")
)

# Merge SEER dataset with age categories
seer_final <- seer_final %>%
  left_join(age_map, by = c("Age" = "age_range"))

# # View new age categories
table(seer_final$age_category, useNA = "ifany")

#Truncate survival at 24 months
seer_final <- seer_final %>%
  mutate(capped_time = pmin(Survival.months, 24))


seer_final$Marital_Grouped <- factor(ifelse(
  seer_final$Marital.status.at.diagnosis %in% c("Married (including common law)", "Unmarried or Domestic Partner"), "Married",
  ifelse(seer_final$Marital.status.at.diagnosis %in% c("Divorced", "Separated", "Widowed"), "Previously Married", "Single/Unknown")
))

table(seer_final$Marital_Grouped)

seer_final$Marital_Grouped <- relevel(seer_final$Marital_Grouped, ref = "Single/Unknown")

seer_final$capped_time <- pmin(seer_final$Survival.months, 24)
seer_final$capped_event <- ifelse(seer_final$Survival.months > 24, 0, seer_final$Event)



write.csv(state_quintile_cutoffs, "state_income_quintile_cutoffs.csv", row.names = FALSE)



#### 2. Descriptive Statistics #####


# Variables to include in the analysis
# Clean the data by removing rows with missing values in the necessary columns
cleaned_data <- seer_final[complete.cases(seer_final[, c("age_category", "Sex", "Marital_Grouped", 
                                                         "surgery_status", "rad_group", "stage_group", 
                                                         "Rural.Urban.Continuum.Code", "Medicaid_Expansion_Status", "income_quintile",
                                                         "Race.and.origin.recode..NHW..NHB..NHAIAN..NHAPI..Hispanic.",
                                                         "Medicaid.Expansion.Status")]), ]


# Trim whitespace and apply the filter again
cleaned_data$Rural.Urban.Continuum.Code <- trimws(cleaned_data$Rural.Urban.Continuum.Code)

# Apply the filter again
cleaned_data <- cleaned_data %>%
  filter(!grepl("Unknown/missing/no match|Unknown/missing", cleaned_data$Rural.Urban.Continuum.Code, ignore.case = TRUE))


# Check the unique values again
table(cleaned_data$Rural.Urban.Continuum.Code)



# Convert necessary variables to factors
cleaned_data$Race.and.origin.recode..NHW..NHB..NHAIAN..NHAPI..Hispanic. <- as.factor(cleaned_data$Race.and.origin.recode..NHW..NHB..NHAIAN..NHAPI..Hispanic.)
cleaned_data$age_category <- as.factor(cleaned_data$age_category)
cleaned_data$Sex <- as.factor(cleaned_data$Sex)
cleaned_data$income_quintile <- as.factor(cleaned_data$income_quintile)
cleaned_data$Marital_Grouped <- as.factor(cleaned_data$Marital_Grouped)
cleaned_data$surgery_status <- as.factor(cleaned_data$surgery_status)
cleaned_data$rad_group <- as.factor(cleaned_data$rad_group)
cleaned_data$stage_group <- as.factor(cleaned_data$stage_group)
cleaned_data$Rural.Urban.Continuum.Code <- as.factor(cleaned_data$Rural.Urban.Continuum.Code)
cleaned_data$Medicaid_Expansion_Status <- as.factor(cleaned_data$Medicaid_Expansion_Status)
cleaned_data$Medicaid.Expansion.Status <- as.factor(cleaned_data$Medicaid.Expansion.Status)

# List of variables to include in the analysis
vars <- c("age_category", "Sex", "Marital_Grouped", "surgery_status", 
          "rad_group", "stage_group", "Rural.Urban.Continuum.Code", "Race.and.origin.recode..NHW..NHB..NHAIAN..NHAPI..Hispanic.", "income_quintile")


# Convert 'Rural.Urban.Continuum.Code' to a factor (if it's not already)
cleaned_data$Rural.Urban.Continuum.Code <- as.factor(cleaned_data$Rural.Urban.Continuum.Code)
# Convert 'income_quintile' to a factor (if it's not already)
cleaned_data$income_quintile <- as.factor(cleaned_data$income_quintile)
# Convert 'Race.and.origin.recode..NHW..NHB..NHAIAN..NHAPI..Hispanic.' to a factor
cleaned_data$Race.and.origin.recode..NHW..NHB..NHAIAN..NHAPI..Hispanic. <- as.factor(cleaned_data$Race.and.origin.recode..NHW..NHB..NHAIAN..NHAPI..Hispanic.)



# Create the contingency table for Rural.Urban.Continuum.Code and Medicaid_Expansion_Status
table_data <- table(cleaned_data$Rural.Urban.Continuum.Code, cleaned_data$Medicaid_Expansion_Status)

# Create the contingency table for income_quintile and Medicaid_Expansion_Status
table_data_income <- table(cleaned_data$income_quintile, cleaned_data$Medicaid_Expansion_Status)

# Create the contingency table for Rural.Urban.Continuum.Code and Medicaid_Expansion_Status
table_data_race <- table(cleaned_data$Race.and.origin.recode..NHW..NHB..NHAIAN..NHAPI..Hispanic., cleaned_data$Medicaid_Expansion_Status)

# Check the table
print(table_data_income)
print(table_data_race)
print(table_data)

# You can adjust the columns based on your analysis

# Initialize an empty list to store results
final_results <- list()


# Check that other relevant columns are factors (e.g., Medicaid Expansion Status, Stage Group, etc.)
cleaned_data$Medicaid_Expansion_Status <- as.factor(cleaned_data$Medicaid_Expansion_Status)
cleaned_data$stage_group <- as.factor(cleaned_data$stage_group)
table_data_income <- table(cleaned_data$income_quintile, cleaned_data$Medicaid_Expansion_Status)

# Run the chi-squared test
chi_squared_test <- chisq.test(table_data)

# Run the chi-squared test
chi_squared_income <- chisq.test(table_data_income)

chi_squared_race <- chisq.test(table_data_race)

# Display the results
print(chi_squared_test)
print(chi_squared_income)
print(chi_squared_race)



# Loop through each variable
for (var in vars) {
  # Create a contingency table for the variable by Medicaid Expansion Status
  temp_table <- table(cleaned_data[[var]], cleaned_data$Medicaid_Expansion_Status)
  
  # Print counts for each level of the variable
  counts <- as.data.frame(temp_table)
  
  # Initialize an empty list to store percentages
  percentage_list <- list()
  
  # Loop through each Medicaid Expansion Status group (0: Non-expansion, 1: Late, 2: 2014, 3: Early)
  for (Medicaid_Expansion_Status in 0:3) {
    # Get the total count for the current Medicaid Expansion Status
    total_count <- sum(temp_table[, Medicaid_Expansion_Status + 1])  # +1 because R indexes start from 1
    
    # Calculate percentages for the current expansion status group
    percentage_df <- temp_table[, Medicaid_Expansion_Status + 1] / total_count * 100
    
    # Store the percentages in the list
    percentage_list[[as.character(Medicaid_Expansion_Status)]] <- percentage_df
  }
  
  # Combine the counts and percentages into one data frame
  combined_table <- cbind(counts, Percentage = unlist(percentage_list))
  
  # Run Chi-Squared test (for categorical variables)
  if (is.factor(cleaned_data[[var]])) {
    chi_sq_test <- chisq.test(temp_table)
    p_value <- chi_sq_test$p.value
  } else {
    # ANOVA for continuous variables
    anova_test <- aov(cleaned_data[[var]] ~ cleaned_data$Medicaid_Expansion_Status)
    p_value <- summary(anova_test)[[1]]["Pr(>F)"][1]
  }
  
  # Add the p-value to the combined table
  combined_table$Chi_Sq_ANOVA_p_value <- p_value
  
  # Add the result to the final list
  final_results[[var]] <- combined_table
}

# Combine all results into a single data frame
final_results_df <- do.call(rbind, final_results)

# Pivot the results to widen the data frame
pivoted_results <- final_results_df %>%
  pivot_wider(names_from = Var2, values_from = c(Freq, Percentage, Chi_Sq_ANOVA_p_value))

# Fix Freq and Percentage columns to merge them properly
pivoted_results <- pivoted_results %>%
  mutate(
    across(starts_with("Freq"), 
           ~ paste0(.x, " (", 
                    get(paste0("Percentage_", gsub("^Freq_", "", cur_column()))), 
                    ")", sep = ""))
  )

# Export the results to a CSV file
write.csv(pivoted_results, "coxdata/descriptive_stats_expansion.csv", row.names = FALSE)



# Create the contingency table for Rural.Urban.Continuum.Code and Medicaid.Expansion.Status
table_state <- table(cleaned_data$Rural.Urban.Continuum.Code, cleaned_data$Medicaid.Expansion.Status)

write.csv(table_state, "supplement_stats/state_rurality.csv", row.names = FALSE)

# Create the contingency table for income_quintile and Medicaid.Expansion.Status
table_state_income <- table(cleaned_data$income_quintile, cleaned_data$Medicaid.Expansion.Status)
write.csv(table_state_income, "supplement_stats/incomes.csv", row.names = FALSE)

# Create the contingency table for Rural.Urban.Continuum.Code and Medicaid.Expansion.Status
table_state_race <- table(cleaned_data$Race.and.origin.recode..NHW..NHB..NHAIAN..NHAPI..Hispanic., cleaned_data$Medicaid.Expansion.Status)
write.csv(table_state_race, "supplement_stats/state_race.csv", row.names = FALSE)
# Check the table
print(table_state_income)
print(table_state_race)
print(table_state)


# Initialize an empty list to store results
final_results <- list()


# Check that other relevant columns are factors (e.g., Medicaid Expansion Status, Stage Group, etc.)
cleaned_data$Medicaid.Expansion.Status <- as.factor(cleaned_data$Medicaid.Expansion.Status)
cleaned_data$stage_group <- as.factor(cleaned_data$stage_group)
table_state_income <- table(cleaned_data$income_quintile, cleaned_data$Medicaid.Expansion.Status)

# Run the chi-squared test
chi_squared_test <- chisq.test(table_state)

# Run the chi-squared test
chi_squared_income <- chisq.test(table_state_income)

chi_squared_race <- chisq.test(table_state_race)

# Display the results
print(chi_squared_test)
print(chi_squared_income)
print(chi_squared_race)



# Check the table
print(table_state)


# Run the chi-squared test
chi_squared_test <- chisq.test(table_data)

# Run the chi-squared test
chi_squared_income <- chisq.test(table_data_income)

chi_squared_race <- chisq.test(table_data_race)

# Display the results
print(chi_squared_test)
print(chi_squared_income)
print(chi_squared_race)
# Initialize an empty list to store results
state_results <- list()


# Run the chi-squared test
chi_squared_test <- chisq.test(table_data)

# Display the results
print(chi_squared_test)

# Loop through each Medicaid Expansion Status group (0: Non-expansion, 1: Late, 2: 2014, 3: Early)
for (Medicaid_Expansion_Status in 0:3) {
  # Get the total count for the current Medicaid Expansion Status
  total_count <- sum(temp_table[, Medicaid_Expansion_Status + 1])  # +1 because R indexes start from 1
  
  # Calculate percentages for the current expansion status group
  percentage_df <- temp_table[, Medicaid_Expansion_Status + 1] / total_count * 100
  
  # Store the percentages in the list
  percentage_list[[as.character(Medicaid_Expansion_Status)]] <- percentage_df
}

# Combine the counts and percentages into one data frame
combined_table <- cbind(counts, Percentage = unlist(percentage_list))

# Run Chi-Squared test (for categorical variables)
if (is.factor(cleaned_data[[var]])) {
  chi_sq_test <- chisq.test(temp_table)
  p_value <- chi_sq_test$p.value
} else {
  # ANOVA for continuous variables
  anova_test <- aov(cleaned_data[[var]] ~ cleaned_data$Medicaid_Expansion_Status)
  p_value <- summary(anova_test)[[1]]["Pr(>F)"][1]
}

# Add the p-value to the combined table
combined_table$Chi_Sq_ANOVA_p_value <- p_value

# Add the result to the final list
state_results[[var]] <- combined_table




# Combine all results into a single data frame
state_results_df <- do.call(rbind, state_results)

pivoted_results <- state_results_df %>%
  pivot_wider(names_from = Var2, values_from = c(Freq, Percentage, Chi_Sq_ANOVA_p_value))

# Fix Freq and Percentage columns to merge them properly
pivoted_results <- pivoted_results %>%
  mutate(
    across(starts_with("Freq"), 
           ~ paste0(.x, " (", 
                    get(paste0("Percentage_", gsub("^Freq_", "", cur_column()))), 
                    ")", sep = ""))
  )

# Export the results to a CSV file
write.csv(pivoted_results, "coxdata/descriptive_stats_state.csv", row.names = FALSE)

print("Results have been exported to 'descriptive_stats_state.csv'")

# Transpose the data frame
transposed_data <- t(pivoted_results)

# If you want to convert it back to a data frame after transposing
transposed_data <- as.data.frame(transposed_data)


# # View the first few rows of the transposed data
head(transposed_data)

# Export the results to a CSV file
write.csv(transposed_data, "coxdata/descriptive_stats_state_transposed.csv", row.names = TRUE)



#### 3. Un-adjusted KM ##### 

# Combine all the plots using `arrange_ggsurvplots` (supports risk tables!)
combined_plot <- arrange_ggsurvplots(plots, print = FALSE, ncol = 2, nrow = ceiling(length(plots)/2))

# # Save to vector-quality PDF
# # ggsave("figures/Figure1_KM_Survival_by_Medicaid_Status.pdf", 
#        plot = combined_plot,
#        width = 18, height = 12, device = cairo_pdf)

#Step 1: Cap survival at 24 months
cleaned_data <- cleaned_data %>%
   mutate(capped_survival = pmin(Survival.months, 24))

# Step 2: Group and summarize
survival_summary_noadj <- cleaned_data %>%
  group_by(dx.yr_grouped, Medicaid_Expansion_Status) %>%
  summarise(
    mean_survival_months = mean(capped_survival, na.rm = TRUE),
    sd_survival_months = sd(capped_survival, na.rm = TRUE),
    n = n()
  ) %>%
  ungroup()

# Step 3: Print the table
print(survival_summary_noadj)

# Save the table as a CSV file
write.csv(survival_summary_noadj, "figures/NONADJ_months_survived_capped.csv", row.names = FALSE)


# Cap survival time at 24 months
cleaned_data$capped_time <- pmin(cleaned_data$Survival.months, 24)
cleaned_data$capped_event <- ifelse(cleaned_data$Survival.months > 24, 0, cleaned_data$Event)


cleaned_data$dx.yr_grouped <- factor(cleaned_data$dx.yr_grouped, 
                                   levels = c("2006-2010", "2011-2013", "2014-2016", "2017-2019"))

cleaned_data$dx.yr_grouped <- as.factor(as.character(cleaned_data$dx.yr_grouped))

# Fit KM model
km_fit <- survfit(Surv(capped_time, capped_event) ~ Medicaid_Expansion_Status +
                    dx.yr_grouped, data = cleaned_data)



# Custom colors for JAMA aesthetic
jama_colors <- c("black", "#003f5c", "#bc5090", "#ffa600")  # black, navy, muted magenta, soft orange

# Split data by time period
seer_split <- split(cleaned_data, cleaned_data$dx.yr_grouped)


# Ensure factor levels are ordered for clean legend and facets
cleaned_data$Medicaid_Expansion_Status <- factor(cleaned_data$Medicaid_Expansion_Status,
                                               levels = c(0, 1, 2, 3),
                                               labels = c("Non-Expansion", "Early Expansion", "2014 Expansion", "Late Expansion"))


cleaned_data$capped_time <- pmin(cleaned_data$Survival.months, 24)
cleaned_data$capped_event <- ifelse(cleaned_data$Survival.months > 24, 0, cleaned_data$Event)


data_by_era <- split(cleaned_data, cleaned_data$dx.yr_grouped)

plots <- map(names(data_by_era), function(era) {
  data_subset <- data_by_era[[era]]
  
  data_subset <- data_subset %>%
    mutate(
      capped_time = pmin(Survival.months, 24),
      capped_event = ifelse(Survival.months > 24, 0, Event)
    )
  
  km <- survfit(Surv(capped_time, capped_event) ~ Medicaid_Expansion_Status, data = data_subset)
  
  ggsurvplot(
    km,
    data = data_subset,
    title = era,
    conf.int = TRUE,
    pval = TRUE,
    risk.table = TRUE,
    risk.table.height = 0.25,
    tables.theme = theme_cleantable() + theme(text = element_text(size = 2, family = "Times")),
    palette = jama_colors,
    xlim = c(0, 24),
    ylim = c(0.5, 1.0),
    break.time.by = 6,
    xlab = "Months Since Diagnosis",
    ylab = "Overall Survival Probability",
    legend = "none"
  )
})


print(plots)

# Log-rank test for comparing survival between groups (Medicaid expansion status)
log_rank_test <- survdiff(Surv(capped_time, capped_event) ~ Medicaid_Expansion_Status + dx.yr_grouped, data = seer_final)

log_rank_results <- data.frame(
  dx.yr_grouped = unique(seer_final$dx.yr_grouped),
  p_value = sapply(unique(seer_final$dx.yr_grouped), function(period) {
    period_data <- seer_final %>% filter(dx.yr_grouped == period)
    log_rank_test <- survdiff(Surv(capped_time, capped_event) ~ Medicaid_Expansion_Status, data = period_data)
    log_rank_test$pvalue
  })
)


# Save the results to a CSV file
write.csv(log_rank_results, "stats_results/log_rank_test_results.csv", row.names = FALSE)

# # View the results
# # # view(log_rank_results)
# Create an empty data frame to store pairwise log-rank test results
pairwise_log_rank_results <- data.frame(dx.yr_grouped = character(),
                                        group_comparison = character(),
                                        p_value = numeric(),
                                        p_value_holm = numeric(),  # Store adjusted p-value
                                        p_value_bh = numeric(),    # Store Benjamini-Hochberg adjusted p-value
                                        p_value_bonferroni = numeric(),  # Store Bonferroni adjusted p-value
                                        stringsAsFactors = FALSE)

# Iterate over each dx.yr_grouped period
for (period in unique(seer_final$dx.yr_grouped)) {
  
  # Filter data for the specific dx.yr_grouped period
  period_data <- seer_final %>% filter(dx.yr_grouped == period)
  
  # Check if there are more than one Medicaid expansion group in the period
  expansion_groups_in_period <- unique(period_data$Medicaid_Expansion_Status)
  
  # Skip the period if there are fewer than two groups
  if (length(expansion_groups_in_period) < 2) {
    next  # Skip to the next iteration
  }

  # Create the Kaplan-Meier survival fit
  km_fit <- survfit(Surv(capped_time, capped_event) ~ Medicaid_Expansion_Status, data = period_data)
  
  # Create an empty vector to store p-values
  p_values <- c()
  group_comparisons <- c()
  
  # Perform pairwise log-rank tests between each pair of Medicaid_Expansion_Status groups
  for (i in 1:(length(expansion_groups_in_period) - 1)) {
    for (j in (i + 1):length(expansion_groups_in_period)) {
      group1 <- expansion_groups_in_period[i]
      group2 <- expansion_groups_in_period[j]
      
      # Subset data for each pair of groups
      pairwise_data <- period_data %>%
        filter(Medicaid_Expansion_Status %in% c(group1, group2))
      
      # Run log-rank test
      log_rank_test <- survdiff(Surv(capped_time, capped_event) ~ Medicaid_Expansion_Status, data = pairwise_data)
      
      # Extract p-value
      p_value <- log_rank_test$pvalue
      
      # Store p-value and the group comparison
      p_values <- c(p_values, p_value)
      group_comparisons <- c(group_comparisons, paste("Group", group1, "vs Group", group2))
    }
  }

  
  # Apply Holm-Bonferroni correction to the p-values
  p_values_holm <- p.adjust(p_values, method = "holm")
  
  # Apply Benjamini-Hochberg correction to the p-values
  p_values_bh <- p.adjust(p_values, method = "BH")
  
  # Apply Bonferroni correction to the p-values
  p_values_bonferroni <- p.adjust(p_values, method = "bonferroni")
  
  # Append results to the data frame
  for (i in 1:length(p_values)) {
    pairwise_log_rank_results <- rbind(pairwise_log_rank_results,
                                       data.frame(dx.yr_grouped = period,
                                                  group_comparison = group_comparisons[i],
                                                  p_value = p_values[i],
                                                  p_value_holm = p_values_holm[i],
                                                  p_value_bh = p_values_bh[i],
                                                  p_value_bonferroni = p_values_bonferroni[i]))
  }
  
}

#### 4. Unadjusted DiD ####
seer_final <- seer_final %>%
  mutate(expansion1_time = case_when(
    dx.yr_grouped %in% c("2006-2010") ~ "Pre",
    dx.yr_grouped %in% c("2011-2013", "2014-2016", "2017-2019") ~ "Post",
    TRUE ~ NA_character_  # Assign NA for years outside the study period
  ))

# Create expansion1_group variable (1 = expansion1 states, 0 = non-expansion states)
seer_final <- seer_final %>%
  mutate(expansion1_group = case_when(
    Medicaid_Expansion_Status == 1 ~ 1,  # Expansion 1 (early)
    Medicaid_Expansion_Status == 0 ~ 0,  # Non-expansion (control)
    TRUE ~ NA_real_  # Assign NA for other groups
  ))

# EARLY - consider early months (post-implementation periods)
# Step 1: Subset data to focus on 2014-2019
seer_final_did_unadj_1 <- seer_final %>%
  filter(dx.yr_grouped %in% c("2006-2010", "2011-2013", "2014-2016", "2017-2019"))

# Step 2: Update the time variable for 2014-2016 and 2017-2019 periods
#seer_final_did_unadj_1 <- seer_final_did_early_1 %>%
seer_final_did_unadj_1 <- seer_final_did_unadj_1 %>% 
 mutate(did_unadj_time = ifelse(dx.yr_grouped %in% c("2006-2010"), "Pre", "Post"))

# Create treatment variable (same as before)
seer_final_did_unadj_1 <- seer_final_did_unadj_1 %>%
  mutate(did_unadj_treat = ifelse(Medicaid_Expansion_Status == 1, 1, 0))

# Step 3: Ensure 'did_unadj_time' is a factor with 'Pre' as reference
seer_final_did_unadj_1 <- seer_final_did_unadj_1 %>%
  mutate(
    did_unadj_time = factor(did_unadj_time, levels = c("Pre", "Post")),
    did_unadj_treat = factor(did_unadj_treat, levels = c(0, 1))
    
  )
# Check reference levels
contrasts(seer_final_did_unadj_1$did_unadj_time)
contrasts(seer_final_did_unadj_1$did_unadj_treat)

seer_final_did_unadj_1$capped_time <- pmin(seer_final_did_unadj_1$Survival.months, 24)
seer_final_did_unadj_1$capped_event <- ifelse(seer_final_did_unadj_1$Survival.months > 24, 0, seer_final_did_unadj_1$Event)


capped_cox_final_did_unadj_1 <- coxph(Surv(capped_time, capped_event) ~ did_unadj_time * did_unadj_treat + cluster(Medicaid.Expansion.Status),
                                data = seer_final_did_unadj_1)

summary(capped_cox_final_did_unadj_1)

capped_cox_final_did_unadj_1 <- tidy(capped_cox_final_did_unadj_1, conf.int = TRUE) %>%
  mutate(HR = exp(estimate),  # Convert coefficients to hazard ratios
         CI_lower = exp(conf.low),
         CI_upper = exp(conf.high)) %>%
  dplyr::select(term, HR, CI_lower, CI_upper, p.value) 

# view(capped_cox_final_did_unadj_1)

write.csv(capped_cox_final_did_unadj_1, "coxdata/unadj_0v1.csv")




# EARLY - consider early months (post-implementation periods)
# Step 1: Subset data to focus on 2014-2019
seer_final_did_unadj_early_1 <- seer_final %>%
  filter(dx.yr_grouped %in% c("2006-2010", "2011-2013"))

# Step 2: Update the time variable for 2014-2016 and 2017-2019 periods
seer_final_did_unadj_early_1 <- seer_final_did_unadj_early_1 %>%
  mutate(did_unadj_time = ifelse(dx.yr_grouped %in% c("2011-2013"), "Post", "Pre"))

# Create treatment variable (same as before)
seer_final_did_unadj_early_1 <- seer_final_did_unadj_early_1 %>%
  mutate(did_unadj_treat = ifelse(Medicaid_Expansion_Status == 1, 1, 0))

# Step 3: Ensure 'did_unadj_time' is a factor with 'Pre' as reference
seer_final_did_unadj_early_1 <- seer_final_did_unadj_early_1 %>%
  mutate(
    did_unadj_time = factor(did_unadj_time, levels = c("Pre", "Post")),
    did_unadj_treat = factor(did_unadj_treat, levels = c(0, 1))
    
  )
# Check reference levels
contrasts(seer_final_did_unadj_early_1$did_unadj_time)
contrasts(seer_final_did_unadj_early_1$did_unadj_treat)

# Step 4: Run Cox regression model
capped_cox_final_early_1 <- coxph(Surv(capped_time, capped_event) ~ did_unadj_time * did_unadj_treat + cluster(Medicaid.Expansion.Status),
                                  data = seer_final_did_unadj_early_1)


# Step 5: Summary of the model
summary(capped_cox_final_early_1)

capped_cox_final_did_unadj_early_1 <- tidy(capped_cox_final_early_1, conf.int = TRUE) %>%
  mutate(HR = exp(estimate),  # Convert coefficients to hazard ratios
         CI_lower = exp(conf.low),
         CI_upper = exp(conf.high)) %>%
  dplyr::select(term, HR, CI_lower, CI_upper, p.value)

# # # view(capped_cox_final_did_unadj_early_1)

write.csv(capped_cox_final_did_unadj_early_1, "coxdata/unadj_0v1_early.csv")



#LATER
# Step 1: Subset data to focus on 2014-2019
seer_final_did_unadj_later_1 <- seer_final %>%
  filter(dx.yr_grouped %in% c("2006-2010", "2014-2016", "2017-2019"))

# Step 2: Update the time variable for 2014-2016 and 2017-2019 periods
seer_final_did_unadj_later_1 <- seer_final_did_unadj_later_1 %>%
  mutate(did_unadj_time = ifelse(dx.yr_grouped %in% c("2014-2016", "2017-2019"), "Post", "Pre"))

# Create treatment variable (same as before)
seer_final_did_unadj_later_1 <- seer_final_did_unadj_later_1 %>%
  mutate(did_unadj_treat = ifelse(Medicaid_Expansion_Status == 1, 1, 0))

# Step 3: Ensure 'did_unadj_time' is a factor with 'Pre' as reference
seer_final_did_unadj_later_1 <- seer_final_did_unadj_later_1 %>%
  mutate(
    did_unadj_time = factor(did_unadj_time, levels = c("Pre", "Post")),
    did_unadj_treat = factor(did_unadj_treat, levels = c(0, 1))
    
  )
# Check reference levels
contrasts(seer_final_did_unadj_later_1$did_unadj_time)
contrasts(seer_final_did_unadj_later_1$did_unadj_treat)

# Step 4: Run Cox regression model
capped_cox_final_later_1 <- coxph(Surv(capped_time, capped_event) ~ did_unadj_time * did_unadj_treat + cluster(Medicaid.Expansion.Status),
                                  data = seer_final_did_unadj_later_1)


# Step 5: Summary of the model
summary(capped_cox_final_later_1)

capped_cox_final_did_unadj_later_1 <- tidy(capped_cox_final_later_1, conf.int = TRUE) %>%
  mutate(HR = exp(estimate),  # Convert coefficients to hazard ratios
         CI_lower = exp(conf.low),
         CI_upper = exp(conf.high)) %>%
  dplyr::select(term, HR, CI_lower, CI_upper, p.value)

# # # view(capped_cox_final_did_unadj_later_1)

write.csv(capped_cox_final_did_unadj_later_1, "coxdata/unadj_0v1_later.csv")


# 0 vs 2
# Create expansion1_time variable (Pre = 2006-2013, Post = 2014-2019)
seer_final <- seer_final %>%
  mutate(expansion2_time = case_when(
    dx.yr_grouped %in% c("2006-2010", "2011-2013") ~ "Pre",
    dx.yr_grouped %in% c("2014-2016", "2017-2019") ~ "Post",
    TRUE ~ NA_character_  # Assign NA for years outside the study period
  ))

# Create expansion1_group variable (1 = expansion2 states, 0 = non-expansion states)
seer_final <- seer_final %>%
  mutate(expansion2_group = case_when(
    Medicaid_Expansion_Status == 2 ~ 1,  # Expansion 2 (on-time)
    Medicaid_Expansion_Status == 0 ~ 0,  # Non-expansion (control)
    TRUE ~ NA_real_  # Assign NA for other groups
  ))



# Subset to only include subjects from Expansion Group 0 and 1
seer_final_did <- seer_final %>%
  filter(Medicaid_Expansion_Status %in% c(0, 2)) %>%
  # Assign treatment variable: 1 for expansion (2), 0 for non-expansion (0)
  mutate(did_unadj_treat_2 = ifelse(Medicaid_Expansion_Status == 2, 1, 0),
         # Define Pre (2006-2013) and Post (2014+)
         did_unadj_time_2 = ifelse(dx.yr_grouped %in% c("2006-2010", "2011-2013"), "Pre", "Post"))

# Convert did_unadj_time and did_unadj_treat to factors with explicit ordering:
seer_final_did <- seer_final_did %>%
  mutate(
    did_unadj_time_2 = factor(did_unadj_time_2, levels = c("Pre", "Post")),
    did_unadj_treat_2 = factor(did_unadj_treat_2, levels = c(0, 1))
  )

seer_final_did <- seer_final_did %>%
  filter(!grepl("Unknown Race", trimws(Race.and.origin.recode..NHW..NHB..NHAIAN..NHAPI..Hispanic.)))


capped_cox_final_did_unadj_2 <- coxph(Surv(capped_time, capped_event) ~ did_unadj_time_2 * did_unadj_treat_2 + cluster(Medicaid.Expansion.Status),
                                data = seer_final_did)

summary(capped_cox_final_did_unadj_2)

capped_cox_final_did_unadj_2 <- tidy(capped_cox_final_did_unadj_2, conf.int = TRUE) %>%
  mutate(HR = exp(estimate),  # Convert coefficients to hazard ratios
         CI_lower = exp(conf.low),
         CI_upper = exp(conf.high)) %>%
  dplyr::select(term, HR, CI_lower, CI_upper, p.value)

# # # view(capped_cox_final_did_unadj_2)

write.csv(capped_cox_final_did_unadj_2, "coxdata/unadj_0v2.csv")


# EARLY - consider early months (post-implementation periods)
# Step 1: Subset data to focus on 2014-2019
seer_final_did_unadj_early_2 <- seer_final %>%
  filter(dx.yr_grouped %in% c("2006-2010", "2011-2013", "2014-2016"))

# Step 2: Update the time variable for 2014-2016 and 2017-2019 periods
seer_final_did_unadj_early_2 <- seer_final_did_unadj_early_2 %>%
  mutate(did_unadj_time = ifelse(dx.yr_grouped %in% c("2014-2016"), "Post", "Pre"))

# Create treatment variable (same as before)
seer_final_did_unadj_early_2 <- seer_final_did_unadj_early_2 %>%
  mutate(did_unadj_treat = ifelse(Medicaid_Expansion_Status == 2, 1, 0))

# Step 3: Ensure 'did_unadj_time' is a factor with 'Pre' as reference
seer_final_did_unadj_early_2 <- seer_final_did_unadj_early_2 %>%
  mutate(
    did_unadj_time = factor(did_unadj_time, levels = c("Pre", "Post")),
    did_unadj_treat = factor(did_unadj_treat, levels = c(0, 1))
    
  )
# Check reference levels
contrasts(seer_final_did_unadj_early_2$did_unadj_time)
contrasts(seer_final_did_unadj_early_2$did_unadj_treat)

# Step 4: Run Cox regression model
capped_cox_final_early_2 <- coxph(Surv(capped_time, capped_event) ~ did_unadj_time * did_unadj_treat + cluster(Medicaid.Expansion.Status),
                                  data = seer_final_did_unadj_early_2)


# Step 5: Summary of the model
summary(capped_cox_final_early_2)

capped_cox_final_did_unadj_early_2 <- tidy(capped_cox_final_early_2, conf.int = TRUE) %>%
  mutate(HR = exp(estimate),  # Convert coefficients to hazard ratios
         CI_lower = exp(conf.low),
         CI_upper = exp(conf.high)) %>%
  dplyr::select(term, HR, CI_lower, CI_upper, p.value)

# # # view(capped_cox_final_did_unadj_early_2)

write.csv(capped_cox_final_did_unadj_early_2, "coxdata/unadj_0v2_early.csv")


# LATE - consider later months (post-implementation periods)
# Step 1: Subset data to focus on 2017-2019
seer_final_did_unadj_later_2 <- seer_final %>%
  filter(dx.yr_grouped %in% c("2006-2010", "2011-2013", "2017-2019"))

# Step 2: Update the time variable for 2014-2016 and 2017-2019 periods
seer_final_did_unadj_later_2 <- seer_final_did_unadj_later_2 %>%
  mutate(did_unadj_time_2 = ifelse(dx.yr_grouped %in% c("2017-2019"), "Post", "Pre"))

# Create treatment variable (same as before)
seer_final_did_unadj_later_2 <- seer_final_did_unadj_later_2 %>%
  mutate(did_unadj_treat_2 = ifelse(Medicaid_Expansion_Status == 2, 1, 0))

# Step 3: Ensure 'did_unadj_time' is a factor with 'Pre' as reference
seer_final_did_unadj_later_2 <- seer_final_did_unadj_later_2 %>%
  mutate(
    did_unadj_time_2 = factor(did_unadj_time_2, levels = c("Pre", "Post")),
    did_unadj_treat_2 = factor(did_unadj_treat_2, levels = c(0, 1))
  )

# Check reference levels
contrasts(seer_final_did_unadj_later_2$did_unadj_time_2)
contrasts(seer_final_did_unadj_later_2$did_unadj_treat_2)

# Run Cox regression model
capped_cox_final_later_2 <- coxph(Surv(capped_time, capped_event) ~ did_unadj_time_2 * did_unadj_treat_2 + cluster(Medicaid.Expansion.Status),
                                  data = seer_final_did_unadj_later_2)


# Step 5: Summary of the model
summary(capped_cox_final_later_2)

capped_cox_final_did_unadj_later_2 <- tidy(capped_cox_final_later_2, conf.int = TRUE) %>%
  mutate(HR = exp(estimate),  # Convert coefficients to hazard ratios
         CI_lower = exp(conf.low),
         CI_upper = exp(conf.high)) %>%
  dplyr::select(term, HR, CI_lower, CI_upper, p.value)

# # # view(capped_cox_final_did_unadj_later_2)

write.csv(capped_cox_final_did_unadj_later_2, "coxdata/unadj_0v2_later.csv")

# 0 vs 3
# Create expansion1_time variable (Pre = 2006-2016, Post = 2017-2019)
seer_final <- seer_final %>%
  mutate(expansion3_time = case_when(
    dx.yr_grouped %in% c("2006-2010", "2011-2013","2014-2016") ~ "Pre",
    dx.yr_grouped %in% c( "2017-2019") ~ "Post",
    TRUE ~ NA_character_  # Assign NA for years outside the study period
  ))

# Create expansion1_group variable (1 = expansion3 states, 0 = non-expansion states)
seer_final <- seer_final %>%
  mutate(expansion3_group = case_when(
    Medicaid_Expansion_Status == 3 ~ 1,  # Expansion 3 (late)
    Medicaid_Expansion_Status == 0 ~ 0,  # Non-expansion (control)
    TRUE ~ NA_real_  # Assign NA for other groups
  ))



# Subset to only include subjects from Expansion Group 0 and 1
seer_final_did <- seer_final %>%
  filter(Medicaid_Expansion_Status %in% c(0, 3)) %>%
  # Assign treatment variable: 1 for late expansion (3), 0 for non-expansion (0)
  mutate(did_unadj_treat_3 = ifelse(Medicaid_Expansion_Status == 3, 1, 0),
         # Define Pre (2006-2016) and Post (2017+)
         did_unadj_time_3 = ifelse(dx.yr_grouped %in% c("2006-2010", "2011-2013", "2014-2016"), "Pre", "Post"))

# Convert did_unadj_time and did_unadj_treat to factors with explicit ordering:
seer_final_did <- seer_final_did %>%
  mutate(
    did_unadj_time_3 = factor(did_unadj_time_3, levels = c("Pre", "Post")),
    did_unadj_treat_3 = factor(did_unadj_treat_3, levels = c(0, 1))
  )

seer_final_did <- seer_final_did %>%
  filter(!grepl("Unknown Race", trimws(Race.and.origin.recode..NHW..NHB..NHAIAN..NHAPI..Hispanic.)))


capped_cox_final_did_unadj_3 <- coxph(Surv(capped_time, capped_event) ~ did_unadj_time_3 * did_unadj_treat_3 + cluster(Medicaid.Expansion.Status),
                                data = seer_final_did)

summary(capped_cox_final_did_unadj_3)

capped_cox_final_did_unadj_3 <- tidy(capped_cox_final_did_unadj_3, conf.int = TRUE) %>%
  mutate(HR = exp(estimate),  # Convert coefficients to hazard ratios
         CI_lower = exp(conf.low),
         CI_upper = exp(conf.high)) %>%
  dplyr::select(term, HR, CI_lower, CI_upper, p.value)

# # # view(capped_cox_final_did_unadj_3)

write.csv(capped_cox_final_did_unadj_3, "coxdata/unadj_0v3.csv")






##### 5. Propensity Score Matching ######


# Ensure survival time is capped at 24 months
cleaned_data <- cleaned_data %>%
  mutate(capped_time = pmin(Survival.months, 24))
 

cleaned_data$binary_group_1 <- ifelse(cleaned_data$Medicaid_Expansion_Status == 0, 0, 
                                    ifelse(cleaned_data$Medicaid_Expansion_Status == 1, 1, NA))
cleaned_data$binary_group_2 <- ifelse(cleaned_data$Medicaid_Expansion_Status == 0, 0, 
                                    ifelse(cleaned_data$Medicaid_Expansion_Status == 2, 1, NA))
cleaned_data$binary_group_3 <- ifelse(cleaned_data$Medicaid_Expansion_Status == 0, 0, 
                                    ifelse(cleaned_data$Medicaid_Expansion_Status == 3, 1, NA))



cleaned_data_cleaned_1 <- cleaned_data[!is.na(cleaned_data$binary_group_1), ]
cleaned_data_cleaned_2 <- cleaned_data[!is.na(cleaned_data$binary_group_2), ]
cleaned_data_cleaned_3 <- cleaned_data[!is.na(cleaned_data$binary_group_3), ]


# Create propensity score model using logistic regression



# Group 1 vs 0
# Step 0: Filter data to complete cases for modeling
cleaned_data_cleaned_1_complete <- cleaned_data_cleaned_1[complete.cases(
  cleaned_data_cleaned_1[, c("age_category", "Sex", 
                           "Race.and.origin.recode..NHW..NHB..NHAIAN..NHAPI..Hispanic.",
                           "rad_group", "stage_group", "surgery_status", 
                           "income_quintile", "Rural.Urban.Continuum.Code", 
                           "Marital_Grouped")]), ]


# Step 1: Fit Propensity Score Model
ps_model_1 <- glm(binary_group_1 ~ age_category + Sex + Combined.Stage +
                    Race.and.origin.recode..NHW..NHB..NHAIAN..NHAPI..Hispanic. + 
                    Rural.Urban.Continuum.Code, 
                  data = cleaned_data_cleaned_1_complete, 
                  family = "binomial")

# Step 2: Attach predicted propensity scores to this filtered dataset
cleaned_data_cleaned_1_complete$pscore <- predict(ps_model_1, type = "response")

# Step 3: Plot Before Matching
ggplot(cleaned_data_cleaned_1_complete, aes(x = pscore, fill = as.factor(binary_group_1))) +
  geom_density(alpha = 0.4) +
  labs(title = "Propensity Score Distributions Before Matching for Group 1", 
       x = "Propensity Score", y = "Density") +
  theme_minimal() + theme(text = element_text(family = "Times"))
 ggsave("stats_results/propensity_scores_before_1.png", width = 8, height = 6, dpi = 300)


cleaned_data_cleaned_1_complete$logit_pscore <- log(cleaned_data_cleaned_1_complete$pscore / (1 - cleaned_data_cleaned_1_complete$pscore))

match_1 <- matchit(
  binary_group_1 ~ 1,
  data = cleaned_data_cleaned_1_complete,
  method = "nearest",
  distance = cleaned_data_cleaned_1_complete$logit_pscore,
  caliper = 0.2,  # on logit scale
  ratio = 2,
  replace = FALSE
)
# # Step 4: Perform Matching
# match_1 <- matchit(binary_group_1 ~ 1, 
#                    data = cleaned_data_cleaned_1_complete, 
#                    method = "nearest", 
#                    distance = cleaned_data_cleaned_1_complete$pscore, 
#                    ratio = 2, 
#                    caliper = 0.2, 
#                    replace = FALSE)

# Step 5: Extract Matched Data
matched_data_1 <- match.data(match_1)

# Step 6: Plot After Matching
ggplot(matched_data_1, aes(x = pscore, fill = as.factor(binary_group_1))) +
  geom_density(alpha = 0.4) +
  labs(title = "Propensity Score Distributions After Matching for Group 1", 
       x = "Propensity Score", y = "Density") +
  theme_minimal() + theme(text = element_text(family = "Times"))
 ggsave("stats_results/propensity_scores_after_1.png", width = 8, height = 6, dpi = 300)

# Step 7: Covariate Balance Plot
plot_1 <- love.plot(match_1, 
                    binary = "raw", 
                    threshold = 0.1, 
                    title = "Covariate Balance: Absolute Standardized Mean Differences for Group 1") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    theme(text = element_text(family = "Times"))
  ) +
  labs(x = "Covariates", y = "Absolute Standardized Mean Difference")

print(plot_1)
# ggsave("stats_results/loveplot_1.png", width = 16, height = 4, dpi = 300)


# Step 8: Show Table of Balance
bal.tab(match_1)

# Step 5: Perform Matching with all covariates in the model
match_1 <- matchit(
  binary_group_1 ~ age_category +
    Sex + 
    Combined.Stage +
    Race.and.origin.recode..NHW..NHB..NHAIAN..NHAPI..Hispanic. +
     Rural.Urban.Continuum.Code,  # Include all covariates used in the model
  data = cleaned_data_cleaned_1_complete,
  method = "nearest",
  distance = cleaned_data_cleaned_1_complete$logit_pscore,
  caliper = 0.2,  # On logit scale
  ratio = 2,
  replace = FALSE
)

# Step 6: Extract Matched Data
matched_data_1 <- match.data(match_1)

# Step 7: Plot Propensity Score Distributions After Matching
ggplot(matched_data_1, aes(x = pscore, fill = as.factor(binary_group_1))) +
  geom_density(alpha = 0.4) +
  labs(title = "Propensity Score Distributions After Matching for Group 1", 
       x = "Propensity Score", y = "Density") +
  theme_minimal()
# # ggsave("stats_results/propensity_scores_after_1.png", width = 8, height = 6, dpi = 300)


# Check balance again
bal.tab(match_1)
# Step 8: Covariate Balance Plot
plot_1 <- love.plot(match_1, 
                    binary = "raw", 
                    threshold = 0.1,  # Set threshold for standardized mean differences
                    title = "Covariate Balance: Absolute Standardized Mean Differences for Group 1") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(x = "Covariates", y = "Absolute Standardized Mean Difference")


#Print the love plot
print(plot_1)
# ggsave("stats_results/full_loveplot_1.png", width = 16, height = 4, dpi = 300)

bal.tab(match_1)


# Group 2 vs 0
# Step 0: Filter data to complete cases for modeling
cleaned_data_cleaned_2_complete <- cleaned_data_cleaned_2[complete.cases(
  cleaned_data_cleaned_2[, c("age_category", "Sex", 
                           "Race.and.origin.recode..NHW..NHB..NHAIAN..NHAPI..Hispanic.",
                           "rad_group", "surgery_status", 
                           "income_quintile", "Rural.Urban.Continuum.Code", 
                           "Marital_Grouped")]), ]

# Step 1: Fit Propensity Score Model
ps_model_2 <- glm(binary_group_2 ~ age_category + Sex + Combined.Stage +
                    Race.and.origin.recode..NHW..NHB..NHAIAN..NHAPI..Hispanic. + 
                   Rural.Urban.Continuum.Code
                  , 
                  data = cleaned_data_cleaned_2_complete, 
                  family = "binomial")

# Step 2: Attach predicted propensity scores to this filtered dataset
cleaned_data_cleaned_2_complete$pscore <- predict(ps_model_2, type = "response")

# Step 3: Plot Before Matching
ggplot(cleaned_data_cleaned_2_complete, aes(x = pscore, fill = as.factor(binary_group_2))) +
  geom_density(alpha = 0.4) +
  labs(title = "Propensity Score Distributions Before Matching for Group 2", 
       x = "Propensity Score", y = "Density") +
  theme_minimal() + theme(text = element_text(family = "Times"))
 ggsave("stats_results/propensity_scores_before_2.png", width = 8, height = 6, dpi = 300)

# Step 4: Perform Matching
cleaned_data_cleaned_2_complete$logit_pscore <- log(cleaned_data_cleaned_2_complete$pscore / (1 - cleaned_data_cleaned_2_complete$pscore))

match_2 <- matchit(
  binary_group_2 ~ 1,
  data = cleaned_data_cleaned_2_complete,
  method = "nearest",
  distance = cleaned_data_cleaned_2_complete$logit_pscore,
  caliper = 0.2,  # on logit scale
  ratio = 2,
  replace = FALSE
)

# match_2 <- matchit(binary_group_2 ~ 1, 
#                    data = cleaned_data_cleaned_2_complete, 
#                    method = "nearest", 
#                    distance = cleaned_data_cleaned_2_complete$pscore, 
#                    ratio = 2, 
#                    caliper = 0.2, 
#                    replace = FALSE)

# Step 5: Extract Matched Data
matched_data_2 <- match.data(match_2)

# Step 6: Plot After Matching
ggplot(matched_data_2, aes(x = pscore, fill = as.factor(binary_group_2))) +
  geom_density(alpha = 0.4) +
  labs(title = "Propensity Score Distributions After Matching for Group 2", 
       x = "Propensity Score", y = "Density") +
  theme_minimal() + theme(text = element_text(family = "Times"))
 ggsave("stats_results/propensity_scores_after_2.png", width = 8, height = 6, dpi = 300)

# Step 7: Covariate Balance Plot
plot_2 <- love.plot(match_2, 
                    binary = "raw", 
                    threshold = 0.1, 
                    title = "Covariate Balance: Absolute Standardized Mean Differences for Group 1") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(x = "Covariates", y = "Absolute Standardized Mean Difference")

print(plot_2)
# ggsave("stats_results/loveplot_2.png", width = 16, height = 4, dpi = 300)


# Step 8: Show Table of Balance
bal.tab(match_2)


# Step 5: Perform Matching with all covariates in the model
match_2 <- matchit(
  binary_group_2 ~ age_category + Sex + Combined.Stage +
    Race.and.origin.recode..NHW..NHB..NHAIAN..NHAPI..Hispanic. +  Rural.Urban.Continuum.Code,  # Include all covariates used in the model
  data = cleaned_data_cleaned_2_complete,
  method = "nearest",
  distance = cleaned_data_cleaned_2_complete$logit_pscore,
  caliper = 0.2,  # On logit scale
  ratio = 2,
  replace = FALSE
)

# Step 6: Extract Matched Data
matched_data_2 <- match.data(match_2)

# Step 7: Plot Propensity Score Distributions After Matching
ggplot(matched_data_2, aes(x = pscore, fill = as.factor(binary_group_2))) +
  geom_density(alpha = 0.4) +
  labs(title = "Propensity Score Distributions After Matching for Group 2", 
       x = "Propensity Score", y = "Density") +
  theme_minimal() + theme(text = element_text(family = "Times"))
 ggsave("stats_results/propensity_scores_after_2.png", width = 8, height = 6, dpi = 300)

# Step 8: Covariate Balance Plot
plot_2 <- love.plot(match_2, 
                    binary = "raw", 
                    threshold = 0.1,  # Set threshold for standardized mean differences
                    title = "Covariate Balance: Absolute Standardized Mean Differences for Group 2") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(x = "Covariates", y = "Absolute Standardized Mean Difference")

#Print the love plot
print(plot_2)
# ggsave("stats_results/full_loveplot_2.png", width = 16, height = 4, dpi = 300)

# Step 8: Show Table of Balance
bal.tab(match_2)



# Group 3 vs 0
# Step 0: Filter data to complete cases for modeling
cleaned_data_cleaned_3_complete <- cleaned_data_cleaned_3[complete.cases(
  cleaned_data_cleaned_3[, c("age_category", "Sex", 
                           "Race.and.origin.recode..NHW..NHB..NHAIAN..NHAPI..Hispanic.",
                           "rad_group", "stage_group", "surgery_status", 
                           "income_quintile", "Rural.Urban.Continuum.Code", 
                           "Marital_Grouped")]), ]

# Step 1: Fit Propensity Score Model
ps_model_3 <- glm(binary_group_3 ~ age_category + Sex + Combined.Stage +
                    Race.and.origin.recode..NHW..NHB..NHAIAN..NHAPI..Hispanic. + 
                    Rural.Urban.Continuum.Code
                    , 
                  data = cleaned_data_cleaned_3_complete, 
                  family = "binomial")

# Step 2: Attach predicted propensity scores to this filtered dataset
cleaned_data_cleaned_3_complete$pscore <- predict(ps_model_3, type = "response")

# Step 3: Plot Before Matching
ggplot(cleaned_data_cleaned_3_complete, aes(x = pscore, fill = as.factor(binary_group_3))) +
  geom_density(alpha = 0.4) +
  labs(title = "Propensity Score Distributions Before Matching for Group 3", 
       x = "Propensity Score", y = "Density") +
  theme_minimal() + theme(text = element_text(family = "Times"))
 ggsave("stats_results/propensity_scores_before_3.png", width = 8, height = 6, dpi = 300)

# Step 4: Perform Matching
cleaned_data_cleaned_3_complete$logit_pscore <- log(cleaned_data_cleaned_3_complete$pscore / (1 - cleaned_data_cleaned_3_complete$pscore))

match_3 <- matchit(
  binary_group_3 ~ 1,
  data = cleaned_data_cleaned_3_complete,
  method = "nearest",
  distance = cleaned_data_cleaned_3_complete$logit_pscore,
  caliper = 0.2,  # on logit scale
  ratio = 2,
  replace = FALSE
)
# 
# 
# match_3 <- matchit(binary_group_3 ~ 1, 
#                    data = cleaned_data_cleaned_3_complete, 
#                    method = "nearest", 
#                    distance = cleaned_data_cleaned_3_complete$pscore, 
#                    ratio = 2, 
#                    caliper = 0.2, 
#                    replace = FALSE)

# Step 5: Extract Matched Data
matched_data_3 <- match.data(match_3)

# Step 6: Plot After Matching
ggplot(matched_data_3, aes(x = pscore, fill = as.factor(binary_group_3))) +
  geom_density(alpha = 0.4) +
  labs(title = "Propensity Score Distributions After Matching for Group 3", 
       x = "Propensity Score", y = "Density") +
  theme_minimal() + theme(text = element_text(family = "Times"))
ggsave("stats_results/propensity_scores_after_3.png", width = 8, height = 6, dpi = 300)

# Step 7: Covariate Balance Plot
plot_3 <- love.plot(match_3, 
                    binary = "raw", 
                    threshold = 0.1, 
                    title = "Covariate Balance: Absolute Standardized Mean Differences for Group 1") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(x = "Covariates", y = "Absolute Standardized Mean Difference")

print(plot_3)


print(plot_3)
# ggsave("stats_results/loveplot_3.png", width = 16, height = 4, dpi = 300)


# Step 8: Show Table of Balance
bal.tab(match_3)


# Step 5: Perform Matching with all covariates in the model
match_3 <- matchit(
  binary_group_3 ~ age_category + Sex + Combined.Stage +
    Race.and.origin.recode..NHW..NHB..NHAIAN..NHAPI..Hispanic. + Rural.Urban.Continuum.Code,  # Include all covariates used in the model
  data = cleaned_data_cleaned_3_complete,
  method = "nearest",
  distance = cleaned_data_cleaned_3_complete$logit_pscore,
  caliper = 0.2,  # On logit scale
  ratio = 2,
  replace = FALSE
)

# Step 6: Extract Matched Data
matched_data_3 <- match.data(match_3)

# Step 7: Plot Propensity Score Distributions After Matching
ggplot(matched_data_3, aes(x = pscore, fill = as.factor(binary_group_3))) +
  geom_density(alpha = 0.4) +
  labs(title = "Propensity Score Distributions After Matching for Group 3", 
       x = "Propensity Score", y = "Density") +
  theme_minimal()+ theme(text = element_text(family = "Times"))
 ggsave("stats_results/propensity_scores_after_3.png", width = 8, height = 6, dpi = 300)

# Step 8: Covariate Balance Plot
plot_3 <- love.plot(match_3, 
                    binary = "raw", 
                    threshold = 0.1,  # Set threshold for standardized mean differences
                    title = "Covariate Balance: Absolute Standardized Mean Differences for Group ") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(x = "Covariates", y = "Absolute Standardized Mean Difference")

#Print the love plot
print(plot_3)
# ggsave("stats_results/full_loveplot_3.png", width = 16, height = 4, dpi = 300)

# Step 8: Show Table of Balance
bal.tab(match_3)


bal.tab(match_1)  # For Group 1 vs 0
bal.tab(match_2)  # For Group 2 vs 0
bal.tab(match_3)  # For Group 3 vs 0

###### 6. Cox Proportional Hazards Model for PSM groups ##### 

# overall (not inclusive of expansion status) does not include Status + dx.yr_grouped
cleaned_data$capped_time <- pmin(cleaned_data$Survival.months, 24)
cleaned_data$capped_event <- ifelse(cleaned_data$Survival.months > 24, 0, cleaned_data$Event)
cleaned_data <- cleaned_data %>%
  filter(!grepl("Unknown Race", trimws(Race.and.origin.recode..NHW..NHB..NHAIAN..NHAPI..Hispanic.)))

#Combine previous 3 Matched Data tables into 1 group
# Combine the matched datasets and add a group identifier
matched_data_1$group <- "Group 1"
matched_data_2$group <- "Group 2"
matched_data_3$group <- "Group 3"

# Combine all matched datasets into one
combined_matched_data <- bind_rows(matched_data_1, matched_data_2, matched_data_3)

# Check the structure of the combined data
str(combined_matched_data)

combined_matched_data$capped_time <- pmin(combined_matched_data$Survival.months, 24)
combined_matched_data$capped_event <- ifelse(combined_matched_data$Survival.months > 24, 0, combined_matched_data$Event)
combined_matched_data <- combined_matched_data[combined_matched_data$capped_time > 0, ]

cox_psm_final <- coxph(
  Surv(capped_time, capped_event) ~ age_category + Sex + Marital_Grouped +
    Rural.Urban.Continuum.Code + Race.and.origin.recode..NHW..NHB..NHAIAN..NHAPI..Hispanic. +
    cluster(Medicaid.Expansion.Status), 
  data = combined_matched_data,
  robust = TRUE  # Adds robust standard errors
)

capped_fin_cox_psm_final_sum <- tidy(cox_psm_final, conf.int = TRUE) %>%
  mutate(HR = exp(estimate),  # Convert coefficients to hazard ratios
         CI_lower = exp(conf.low),
         CI_upper = exp(conf.high)) %>%
  dplyr::select(term, HR, CI_lower, CI_upper, p.value)

# # view(capped_fin_cox_psm_final_sum)

write.csv(capped_fin_cox_psm_final_sum, "coxdata/psm_covariates.csv")

# Group by Group Difference in Differences analysis
# 0 vs 1
# Create expansion1_time variable (Pre = 2006-2010, Post = 2011-2019)
matched_data_1 <- matched_data_1 %>%
  mutate(expansion1_time = case_when(
    dx.yr_grouped %in% c("2006-2010") ~ "Pre",
    dx.yr_grouped %in% c("2011-2013", "2014-2016", "2017-2019") ~ "Post",
    TRUE ~ NA_character_  # Assign NA for years outside the study period
  ))

# Create expansion1_group variable (1 = expansion1 states, 0 = non-expansion states)
matched_data_1 <- matched_data_1 %>%
  mutate(expansion1_group = case_when(
    Medicaid_Expansion_Status == 1 ~ 1,  # Expansion 1 (early)
    Medicaid_Expansion_Status == 0 ~ 0,  # Non-expansion (control)
    TRUE ~ NA_real_  # Assign NA for other groups
  ))

# Subset to only include subjects from Expansion Group 0 and 1
matched_data_1_did <- matched_data_1 %>%
  filter(Medicaid_Expansion_Status %in% c(0, 1)) %>%
  # Create a binary treatment variable:
  # did_treat = 1 if Expansion Group 1, 0 if Expansion Group 0
  mutate(did_treat = ifelse(Medicaid_Expansion_Status == 1, 1, 0),
         # Create a time variable: "Pre" for 2006-2010, "Post" for later years.
         did_time = ifelse(dx.yr_grouped == "2006-2010", "Pre", "Post"))

# Convert did_time and did_treat to factors with explicit ordering:
matched_data_1_did <- matched_data_1_did %>%
  mutate(
    did_time = factor(did_time, levels = c("Pre", "Post")),
    did_treat = factor(did_treat, levels = c(0, 1))
  )

matched_data_1_did$capped_time <- pmin(matched_data_1_did$Survival.months, 24)
matched_data_1_did$capped_event <- ifelse(matched_data_1_did$Survival.months > 24, 0, matched_data_1_did$Event)



capped_fin_cox_psm_final_did_1 <- coxph(Surv(capped_time, capped_event) ~ did_time * did_treat + cluster(Medicaid.Expansion.Status),
                                          data = matched_data_1_did, robust = TRUE)


summary(capped_fin_cox_psm_final_did_1)

capped_fin_cox_psm_final_did_1 <- tidy(capped_fin_cox_psm_final_did_1, conf.int = TRUE) %>%
  mutate(HR = exp(estimate),  # Convert coefficients to hazard ratios
         CI_lower = exp(conf.low),
         CI_upper = exp(conf.high)) %>%
  dplyr::select(term, HR, CI_lower, CI_upper, p.value)

# # # view(capped_fin_cox_psm_final_did_1)

write.csv(capped_fin_cox_psm_final_did_1, "coxdata/psm_0v1.csv")



# EARLY - consider early months (post-implementation periods)
# Step 1: Subset data to focus on 2014-2019
matched_data_1_did_early_1 <- matched_data_1 %>%
  filter(dx.yr_grouped %in% c("2006-2010", "2011-2013"))

# Step 2: Update the time variable for 2014-2016 and 2017-2019 periods
matched_data_1_did_early_1 <- matched_data_1_did_early_1 %>%
  mutate(did_time = ifelse(dx.yr_grouped %in% c("2011-2013"), "Post", "Pre"))

# Create treatment variable (same as before)
matched_data_1_did_early_1 <- matched_data_1_did_early_1 %>%
  mutate(did_treat = ifelse(Medicaid_Expansion_Status == 1, 1, 0))

# Step 3: Ensure 'did_time' is a factor with 'Pre' as reference
matched_data_1_did_early_1 <- matched_data_1_did_early_1 %>%
  mutate(
    did_time = factor(did_time, levels = c("Pre", "Post")),
    did_treat = factor(did_treat, levels = c(0, 1))
    
  )

# Check reference levels
contrasts(matched_data_1_did_early_1$did_time)
contrasts(matched_data_1_did_early_1$did_treat)

# Step 4: Run Cox regression model
capped_fin_cox_psm_final_early_1 <- coxph(Surv(capped_time, capped_event) ~ did_time * did_treat+ cluster(Medicaid.Expansion.Status),
                                          data = matched_data_1_did_early_1, robust = TRUE)

# Step 5: Summary of the model
summary(capped_fin_cox_psm_final_early_1)

capped_fin_cox_psm_final_did_early_1 <- tidy(capped_fin_cox_psm_final_early_1, conf.int = TRUE) %>%
  mutate(HR = exp(estimate),  # Convert coefficients to hazard ratios
         CI_lower = exp(conf.low),
         CI_upper = exp(conf.high)) %>%
  dplyr::select(term, HR, CI_lower, CI_upper, p.value)

# # # view(capped_fin_cox_psm_final_did_early_1)

write.csv(capped_fin_cox_psm_final_did_early_1, "coxdata/psm_0v1_early.csv")


# LATER
# Step 1: Subset data to focus on 2014-2019
matched_data_1_did_later_1 <- matched_data_1 %>%
  filter(dx.yr_grouped %in% c("2006-2010", "2014-2016", "2017-2019"))

# Step 2: Update the time variable for 2014-2016 and 2017-2019 periods
matched_data_1_did_later_1 <- matched_data_1_did_later_1 %>%
  mutate(did_time = ifelse(dx.yr_grouped %in% c("2014-2016", "2017-2019"), "Post", "Pre"))

# Create treatment variable (same as before)
matched_data_1_did_later_1 <- matched_data_1_did_later_1 %>%
  mutate(did_treat = ifelse(Medicaid_Expansion_Status == 1, 1, 0))

# Step 3: Ensure 'did_time' is a factor with 'Pre' as reference
matched_data_1_did_later_1 <- matched_data_1_did_later_1 %>%
  mutate(
    did_time = factor(did_time, levels = c("Pre", "Post")),
    did_treat = factor(did_treat, levels = c(0, 1))
    
  )

# Check reference levels
contrasts(matched_data_1_did_later_1$did_time)
contrasts(matched_data_1_did_later_1$did_treat)

# Step 4: Run Cox regression model
capped_fin_cox_psm_final_later_1 <- coxph(Surv(capped_time, capped_event) ~ did_time * did_treat + cluster(Medicaid.Expansion.Status),
                                            data = matched_data_1_did_later_1, robust = TRUE)

# Step 5: Summary of the model
summary(capped_fin_cox_psm_final_later_1)

capped_fin_cox_psm_final_did_later_1 <- tidy(capped_fin_cox_psm_final_later_1, conf.int = TRUE) %>%
  mutate(HR = exp(estimate),  # Convert coefficients to hazard ratios
         CI_lower = exp(conf.low),
         CI_upper = exp(conf.high)) %>%
  dplyr::select(term, HR, CI_lower, CI_upper, p.value)

# # # view(capped_fin_cox_psm_final_did_later_1)

write.csv(capped_fin_cox_psm_final_did_later_1, "coxdata/psm_0v1_later.csv")

### LAST ####

# Step 1: Subset data to focus on 2014-2019
matched_data_1_did_later_1 <- matched_data_1 %>%
  filter(dx.yr_grouped %in% c("2006-2010", "2017-2019"))

# Step 2: Update the time variable for 2014-2016 and 2017-2019 periods
matched_data_1_did_later_1 <- matched_data_1_did_later_1 %>%
  mutate(did_time = ifelse(dx.yr_grouped %in% c("2017-2019"), "Post", "Pre"))

# Create treatment variable (same as before)
matched_data_1_did_later_1 <- matched_data_1_did_later_1 %>%
  mutate(did_treat = ifelse(Medicaid_Expansion_Status == 1, 1, 0))

# Step 3: Ensure 'did_time' is a factor with 'Pre' as reference
matched_data_1_did_later_1 <- matched_data_1_did_later_1 %>%
  mutate(
    did_time = factor(did_time, levels = c("Pre", "Post")),
    did_treat = factor(did_treat, levels = c(0, 1))
    
  )

# Check reference levels
contrasts(matched_data_1_did_later_1$did_time)
contrasts(matched_data_1_did_later_1$did_treat)

# Step 4: Run Cox regression model
capped_fin_cox_psm_final_later_1 <- coxph(Surv(capped_time, capped_event) ~ did_time * did_treat + cluster(Medicaid.Expansion.Status),
                                          data = matched_data_1_did_later_1, robust = TRUE)

# Step 5: Summary of the model
summary(capped_fin_cox_psm_final_later_1)

capped_fin_cox_psm_final_did_later_1 <- tidy(capped_fin_cox_psm_final_later_1, conf.int = TRUE) %>%
  mutate(HR = exp(estimate),  # Convert coefficients to hazard ratios
         CI_lower = exp(conf.low),
         CI_upper = exp(conf.high)) %>%
  dplyr::select(term, HR, CI_lower, CI_upper, p.value)

# # # view(capped_fin_cox_psm_final_did_later_1)

write.csv(capped_fin_cox_psm_final_did_later_1, "coxdata/psm_0v1_LAST.csv")


#### Mid Group ####

# Step 1: Subset data to focus on 2014-2019
matched_data_1_did_later_1 <- matched_data_1 %>%
  filter(dx.yr_grouped %in% c("2006-2010", "2014-2016"))

# Step 2: Update the time variable for 2014-2016 and 2017-2019 periods
matched_data_1_did_later_1 <- matched_data_1_did_later_1 %>%
  mutate(did_time = ifelse(dx.yr_grouped %in% c("2014-2016"), "Post", "Pre"))

# Create treatment variable (same as before)
matched_data_1_did_later_1 <- matched_data_1_did_later_1 %>%
  mutate(did_treat = ifelse(Medicaid_Expansion_Status == 1, 1, 0))

# Step 3: Ensure 'did_time' is a factor with 'Pre' as reference
matched_data_1_did_later_1 <- matched_data_1_did_later_1 %>%
  mutate(
    did_time = factor(did_time, levels = c("Pre", "Post")),
    did_treat = factor(did_treat, levels = c(0, 1))
    
  )

# Check reference levels
contrasts(matched_data_1_did_later_1$did_time)
contrasts(matched_data_1_did_later_1$did_treat)

# Step 4: Run Cox regression model
capped_fin_cox_psm_final_later_1 <- coxph(Surv(capped_time, capped_event) ~ did_time * did_treat + cluster(Medicaid.Expansion.Status),
                                          data = matched_data_1_did_later_1, robust = TRUE)

# Step 5: Summary of the model
summary(capped_fin_cox_psm_final_later_1)

capped_fin_cox_psm_final_did_later_1 <- tidy(capped_fin_cox_psm_final_later_1, conf.int = TRUE) %>%
  mutate(HR = exp(estimate),  # Convert coefficients to hazard ratios
         CI_lower = exp(conf.low),
         CI_upper = exp(conf.high)) %>%
  dplyr::select(term, HR, CI_lower, CI_upper, p.value)

# # # view(capped_fin_cox_psm_final_did_later_1)

write.csv(capped_fin_cox_psm_final_did_later_1, "coxdata/psm_0v1_MID.csv")



##### LATE + MID #
# Step 1: Subset data to focus on 2011-2016
matched_data_1_did_later_1 <- matched_data_1 %>%
  filter(dx.yr_grouped %in% c("2006-2010", "2011-2013", "2014-2016"))

# Step 2: Update the time variable for 2014-2016 and 2017-2019 periods
matched_data_1_did_later_1 <- matched_data_1_did_later_1 %>%
  mutate(did_time = ifelse(dx.yr_grouped %in% c("2011-2013", "2014-2016"), "Post", "Pre"))

# Create treatment variable (same as before)
matched_data_1_did_later_1 <- matched_data_1_did_later_1 %>%
  mutate(did_treat = ifelse(Medicaid_Expansion_Status == 1, 1, 0))

# Step 3: Ensure 'did_time' is a factor with 'Pre' as reference
matched_data_1_did_later_1 <- matched_data_1_did_later_1 %>%
  mutate(
    did_time = factor(did_time, levels = c("Pre", "Post")),
    did_treat = factor(did_treat, levels = c(0, 1))
    
  )

# Check reference levels
contrasts(matched_data_1_did_later_1$did_time)
contrasts(matched_data_1_did_later_1$did_treat)

# Step 4: Run Cox regression model
capped_fin_cox_psm_final_later_1 <- coxph(Surv(capped_time, capped_event) ~ did_time * did_treat + cluster(Medicaid.Expansion.Status),
                                          data = matched_data_1_did_later_1, robust = TRUE)

# Step 5: Summary of the model
summary(capped_fin_cox_psm_final_later_1)

capped_fin_cox_psm_final_did_later_1 <- tidy(capped_fin_cox_psm_final_later_1, conf.int = TRUE) %>%
  mutate(HR = exp(estimate),  # Convert coefficients to hazard ratios
         CI_lower = exp(conf.low),
         CI_upper = exp(conf.high)) %>%
  dplyr::select(term, HR, CI_lower, CI_upper, p.value)

# # # view(capped_fin_cox_psm_final_did_later_1)

write.csv(capped_fin_cox_psm_final_did_later_1, "coxdata/psm_0v1_2011-2016.csv")



# 0 vs 2
# Create expansion1_time variable (Pre = 2006-2013, Post = 2014-2019)
matched_data_2 <- matched_data_2 %>%
  mutate(expansion2_time = case_when(
    dx.yr_grouped %in% c("2006-2010", "2011-2013") ~ "Pre",
    dx.yr_grouped %in% c("2014-2016", "2017-2019") ~ "Post",
    TRUE ~ NA_character_  # Assign NA for years outside the study period
  ))

# Create expansion1_group variable (1 = expansion2 states, 0 = non-expansion states)
matched_data_2 <- matched_data_2 %>%
  mutate(expansion2_group = case_when(
    Medicaid_Expansion_Status == 2 ~ 1,  # Expansion 2 (early)
    Medicaid_Expansion_Status == 0 ~ 0,  # Non-expansion (control)
    TRUE ~ NA_real_  # Assign NA for other groups
  ))



# Subset to only include subjects from Expansion Group 0 and 2
matched_data_2_did <- matched_data_2 %>%
  filter(Medicaid_Expansion_Status %in% c(0, 2)) %>%
  # Assign treatment variable: 1 for expansion (2), 0 for non-expansion (0)
  mutate(did_treat_2 = ifelse(Medicaid_Expansion_Status == 2, 1, 0),
         # Define Pre (2006-2013) and Post (2014+)
         did_time_2 = ifelse(dx.yr_grouped %in% c("2006-2010", "2011-2013"), "Pre", "Post"))

# Convert did_time and did_treat to factors with explicit ordering:
matched_data_2_did <- matched_data_2_did %>%
  mutate(
    did_time_2 = factor(did_time_2, levels = c("Pre", "Post")),
    did_treat_2 = factor(did_treat_2, levels = c(0, 1))
  )

matched_data_2_did <- matched_data_2_did %>%
  filter(!grepl("Unknown Race", trimws(Race.and.origin.recode..NHW..NHB..NHAIAN..NHAPI..Hispanic.)))


capped_fin_cox_psm_final_did_2 <- coxph(Surv(capped_time, capped_event) ~ did_time_2 * did_treat_2 + cluster(Medicaid.Expansion.Status),
                                          data = matched_data_2_did, robust = TRUE)

summary(capped_fin_cox_psm_final_did_2)

capped_fin_cox_psm_final_did_2 <- tidy(capped_fin_cox_psm_final_did_2, conf.int = TRUE) %>%
  mutate(HR = exp(estimate),  # Convert coefficients to hazard ratios
         CI_lower = exp(conf.low),
         CI_upper = exp(conf.high)) %>%
  dplyr::select(term, HR, CI_lower, CI_upper, p.value)

# # # view(capped_fin_cox_psm_final_did_2)

write.csv(capped_fin_cox_psm_final_did_2, "coxdata/psm_0v2.csv")

# EARLY - consider early months (post-implementation periods)
# Step1: Subset data to focus on 2017-2019
matched_data_2_did_early_2 <- matched_data_2 %>%
  filter(dx.yr_grouped %in% c("2006-2010", "2011-2013", "2014-2016"))

# Step 2: Update the time variable for 2014-2016 and 2017-2019 periods
matched_data_2_did_early_2 <- matched_data_2_did_early_2 %>%
  mutate(did_time_2 = ifelse(dx.yr_grouped %in% c("2014-2016"), "Post", "Pre"))

# Create treatment variable (same as before)
matched_data_2_did_early_2 <- matched_data_2_did_early_2 %>%
  mutate(did_treat_2 = ifelse(Medicaid_Expansion_Status == 2, 1, 0))

# Step 3: Ensure 'did_time' is a factor with 'Pre' as reference
matched_data_2_did_early_2 <- matched_data_2_did_early_2 %>%
  mutate(
    did_time_2 = factor(did_time_2, levels = c("Pre", "Post")),
    did_treat_2 = factor(did_treat_2, levels = c(0, 1))
  )

# Check reference levels
contrasts(matched_data_2_did_early_2$did_time_2)
contrasts(matched_data_2_did_early_2$did_treat_2)



# Step 4: Run Cox regression model
capped_fin_cox_psm_final_early_2 <- coxph(Surv(capped_time, capped_event) ~ did_time_2 * did_treat_2 + cluster(Medicaid.Expansion.Status),
                                          data = matched_data_2_did_early_2, robust = TRUE)

# Step 5: Summary of the model
summary(capped_fin_cox_psm_final_early_2)

capped_fin_cox_psm_final_did_early_2 <- tidy(capped_fin_cox_psm_final_early_2, conf.int = TRUE) %>%
  mutate(HR = exp(estimate),  # Convert coefficients to hazard ratios
         CI_lower = exp(conf.low),
         CI_upper = exp(conf.high)) %>%
  dplyr::select(term, HR, CI_lower, CI_upper, p.value)

# # # view(capped_fin_cox_psm_final_did_early_2)

write.csv(capped_fin_cox_psm_final_did_early_2, "coxdata/psm_0v2_early.csv")




# LATER
# Step1: Subset data to focus on 2017-2019
matched_data_2_did_later_2 <- matched_data_2 %>%
  filter(dx.yr_grouped %in% c("2006-2010", "2011-2013", "2017-2019"))

# Step 2: Update the time variable for 2014-2016 and 2017-2019 periods
matched_data_2_did_later_2 <- matched_data_2_did_later_2 %>%
  mutate(did_time_2 = ifelse(dx.yr_grouped %in% c("2017-2019"), "Post", "Pre"))

# Create treatment variable (same as before)
matched_data_2_did_later_2 <- matched_data_2_did_later_2 %>%
  mutate(did_treat_2 = ifelse(Medicaid_Expansion_Status == 2, 1, 0))

# Step 3: Ensure 'did_time' is a factor with 'Pre' as reference
matched_data_2_did_later_2 <- matched_data_2_did_later_2 %>%
  mutate(
    did_time_2 = factor(did_time_2, levels = c("Pre", "Post")),
    did_treat_2 = factor(did_treat_2, levels = c(0, 1))
  )

# Check reference levels
contrasts(matched_data_2_did_later_2$did_time_2)
contrasts(matched_data_2_did_later_2$did_treat_2)




# Step 4: Run Cox regression model
capped_fin_cox_psm_final_later_2 <- coxph(Surv(capped_time, capped_event) ~ did_time_2 * did_treat_2 + cluster(Medicaid.Expansion.Status),
                                            data = matched_data_2_did_later_2, robust = TRUE)

# Step 5: Summary of the model
summary(capped_fin_cox_psm_final_later_2)

capped_fin_cox_psm_final_did_later_2 <- tidy(capped_fin_cox_psm_final_later_2, conf.int = TRUE) %>%
  mutate(HR = exp(estimate),  # Convert coefficients to hazard ratios
         CI_lower = exp(conf.low),
         CI_upper = exp(conf.high)) %>%
  dplyr::select(term, HR, CI_lower, CI_upper, p.value)

# # # view(capped_fin_cox_psm_final_did_later_2)

write.csv(capped_fin_cox_psm_final_did_later_2, "coxdata/psm_0v2_later.csv")



# 0 vs 3
# Subset to only include subjects from Expansion Group 0 and 3
matched_data_3_did <- matched_data_3 %>%
  filter(Medicaid_Expansion_Status %in% c(0, 3)) %>%
  # Assign treatment variable: 1 for expansion (2), 0 for non-expansion (0)
  mutate(did_treat_3 = ifelse(Medicaid_Expansion_Status == 3, 1, 0),
         # Define Pre (2006-2016) and Post (2017+)
         did_time_3 = ifelse(dx.yr_grouped %in% c("2006-2010", "2011-2013", "2014-2016"), "Pre", "Post"))

# Convert did_time and did_treat to factors with explicit ordering:
matched_data_3_did <- matched_data_3_did %>%
  mutate(
    did_time_3 = factor(did_time_3, levels = c("Pre", "Post")),
    did_treat_3 = factor(did_treat_3, levels = c(0, 1))
  )

matched_data_3_did <- matched_data_3_did %>%
  filter(!grepl("Unknown Race", trimws(Race.and.origin.recode..NHW..NHB..NHAIAN..NHAPI..Hispanic.)))


capped_fin_cox_psm_final_did_3 <- coxph(Surv(capped_time, capped_event) ~ did_time_3 * did_treat_3 + cluster(Medicaid.Expansion.Status),
                                        data = matched_data_3_did, robust = TRUE)

summary(capped_fin_cox_psm_final_did_3)

capped_fin_cox_psm_final_did_3 <- tidy(capped_fin_cox_psm_final_did_3, conf.int = TRUE) %>%
  mutate(HR = exp(estimate),  # Convert coefficients to hazard ratios
         CI_lower = exp(conf.low),
         CI_upper = exp(conf.high)) %>%
  dplyr::select(term, HR, CI_lower, CI_upper, p.value)

# # # view(capped_fin_cox_psm_final_did_3)

write.csv(capped_fin_cox_psm_final_did_3, "coxdata/psm_0v3.csv")



### 7. Placebo Falsification Modeling  ### 

# Calculate average outcome by diagnosis year and expansion status
avg_cleaned_data <- cleaned_data %>%
  group_by(dx.yr, Medicaid_Expansion_Status) %>%
  summarize(avg_survival = mean(capped_time, na.rm = TRUE))



# Shaded time periods with labels
shaded_eras <- data.frame(
  xmin = c(2006, 2011, 2014, 2017),
  xmax = c(2011, 2014, 2017, 2020),
  fill = c("Pre-Policy (2006–2010)", "Early Expansion (2011–2013)", 
           "Main Expansion (2014–2016)", "Late Expansion (2017–2019)")
)


avg_cleaned_data <- cleaned_data %>%
  group_by(dx.yr, Medicaid_Expansion_Status) %>%
  summarise(
    avg_survival = mean(capped_time, na.rm = TRUE),
    n = n(),
    sd_survival = sd(capped_time, na.rm = TRUE),
    se_survival = sd_survival / sqrt(n),
    avg_survival_lower = avg_survival - 1.96 * se_survival,
    avg_survival_upper = avg_survival + 1.96 * se_survival
  ) %>%
  ungroup()

# Base plot
exsurv_plot <- ggplot(avg_cleaned_data, aes(x = dx.yr, y = avg_survival, 
                                          color = factor(Medicaid_Expansion_Status))) +
  # Background shading
  geom_rect(data = shaded_eras, 
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill),
            inherit.aes = FALSE, alpha = 0.08, show.legend = FALSE) +
  # Confidence interval ribbons
  # # CI error bars
  # geom_errorbar(aes(ymin = avg_survival_lower, ymax = avg_survival_upper), 
  #               width = 0.2, linewidth = 0.7, alpha = 0.7) +
  # Survival trend lines and points
  geom_line(aes(group = Medicaid_Expansion_Status), linewidth = 1.2) +
  geom_point(size = 3, shape = 21, fill = "white", stroke = 1.2) +
  
  # Color and labels
  scale_color_manual(
    values = jama_colors,
    labels = c("Non-Expansion", "Early Expansion", "2014 Expansion", "Late Expansion")
  ) +
  
  # Axis labels and title
  labs(
    title = "Average 24-Month Survival of Whole Cohort by Medicaid Expansion Status",
    subtitle = "Trends in NSCLC survival by year of diagnosis, SEER 2006–2019",
    x = "Diagnosis Year",
    y = "Average Survival (months)",
    color = "Expansion Status"
  ) +
  
  # Scales and limits
  coord_cartesian(xlim = c(2006, 2020), ylim = c(17, 23)) +
  
  # Theme tweaks
  theme_minimal(base_size = 14) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "gray85"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 13, margin = margin(b = 10)),
    axis.title = element_text(size = 14),
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 12),
    text = element_text(family = "Times")
  )

# Save the plot
ggsave("figures/exsurv_plot.png", plot = exsurv_plot, width = 10, height = 6.5, dpi = 300)

print(exsurv_plot)



# Falsification test for parallel trends assumption
pre_expansion_data <- subset(avg_cleaned_data, dx.yr < 2014)

# Regression model
falsification_model <- lm(avg_survival ~ dx.yr, data = pre_expansion_data)


# Check the interaction term for significance
summary(falsification_model) # shows falsification is not violated! :D 

sink("coxdata/falsification.txt")
print(summary(falsification_model))
sink()


### Placebo Falsification for SURVIVAL ###

# Subset pre-expansion years only
placebo_data <- cleaned_data %>%
  filter(dx.yr_grouped %in% c("2006-2010"),   # Ensure dx.yr_grouped is treated as a string
         Medicaid_Expansion_Status %in% c("Non-Expansion", "Early Expansion")) %>%
  mutate(
    fake_time = ifelse(dx.yr == 2010, 1, 0),  # Pretend 2010 is "post"
    treat = ifelse(Medicaid_Expansion_Status == "Early Expansion", 1, 0)
  )

# Run placebo DiD model
placebo_model <- lm(capped_survival ~ treat * fake_time, data = placebo_data)
summary(placebo_model)


modelsummary(placebo_model, output = "figures/placebo_survival_model_output.docx")

# Output
tidy(placebo_model, conf.int = TRUE, exponentiate = TRUE)


# Assuming you have already run tidy() on the placebo model
placebo_tidy <- tidy(placebo_model, conf.int = TRUE, exponentiate = TRUE)

# Save the tidy result as a CSV
write.csv(placebo_tidy, "low_income/placebofalse_survival.csv", row.names = FALSE)



library(dplyr)

# 1. Filter to pre-expansion era and the two groups of interest
placebo_grouped_data <- cleaned_data %>%
  filter(dx.yr_grouped == "2006-2010",
         Medicaid_Expansion_Status %in% c("Non-Expansion", "Early Expansion")) %>%
  group_by(Medicaid_Expansion_Status, dx.yr) %>%
  summarize(
    mean_survival = mean(capped_survival, na.rm = TRUE),
    n = n(),
    .groups = 'drop'
  ) %>%
  mutate(
    fake_time = ifelse(dx.yr == 2010, 1, 0),  # Pretend 2010 is post
    treat = ifelse(Medicaid_Expansion_Status == "Early Expansion", 1, 0)
  )

# 2. Fit the DiD model by group
placebo_model_group <- lm(mean_survival ~ treat * fake_time, data = placebo_grouped_data)
summary(placebo_model_group)



# 1. Subset to pre-expansion years and relevant groups
placebo_data <- cleaned_data %>%
  filter(dx.yr_grouped == "2006-2010",
         Medicaid_Expansion_Status %in% c("Non-Expansion", "Early Expansion")) %>%
  mutate(
    fake_time = ifelse(dx.yr == 2010, 1, 0),  # Treat 2010 as "post"
    treat = ifelse(Medicaid_Expansion_Status == "Early Expansion", 1, 0)
  )

# 2. Aggregate to group-level proportion of early-stage diagnosis
grouped_placebo <- placebo_data %>%
  group_by(dx.yr, treat, fake_time) %>%
  summarise(
    early_stage_rate = mean(early_stage, na.rm = TRUE),
    n = n()
  ) %>%
  ungroup()

# 3. Run linear model on group-level early-stage diagnosis rates
placebo_model_grouped <- lm(early_stage_rate ~ treat * fake_time, data = grouped_placebo)

# 4. Output summary
summary(placebo_model_grouped)



### 8. Kaplan-Meier Curves of Propensity Score-Matched Populations ### 

# Get the union of all column names from the three datasets
all_columns <- union(union(colnames(matched_data_1), colnames(matched_data_2)), colnames(matched_data_3))

# Add missing columns to each dataset
add_missing_columns <- function(df, all_cols) {
  missing_cols <- setdiff(all_cols, colnames(df))
  for (col in missing_cols) {
    df[[col]] <- NA
  }
  return(df)
}

# Add missing columns to each dataset
matched_data_1 <- add_missing_columns(matched_data_1, all_columns)
matched_data_2 <- add_missing_columns(matched_data_2, all_columns)
matched_data_3 <- add_missing_columns(matched_data_3, all_columns)

# Reorder columns in each dataset to match the same order
matched_data_1 <- matched_data_1[, all_columns]
matched_data_2 <- matched_data_2[, all_columns]
matched_data_3 <- matched_data_3[, all_columns]

# Now combine the datasets
combined_psm_data <- dplyr::bind_rows(matched_data_1, matched_data_2, matched_data_3)


# Check the result
head(combined_psm_data)

# Turn each variable into factors
combined_psm_data$Race.and.origin.recode..NHW..NHB..NHAIAN..NHAPI..Hispanic. <- as.factor(combined_psm_data$Race.and.origin.recode..NHW..NHB..NHAIAN..NHAPI..Hispanic.)
combined_psm_data$age_category <- as.factor(combined_psm_data$age_category)
combined_psm_data$Sex <- as.factor(combined_psm_data$Sex)
combined_psm_data$income_quintile <- as.factor(combined_psm_data$income_quintile)
combined_psm_data$Marital_Grouped <- as.factor(combined_psm_data$Marital_Grouped)
combined_psm_data$surgery_status <- as.factor(combined_psm_data$surgery_status)
combined_psm_data$rad_group <- as.factor(combined_psm_data$rad_group)
combined_psm_data$stage_group <- as.factor(combined_psm_data$stage_group)
combined_psm_data$Rural.Urban.Continuum.Code <- as.factor(combined_psm_data$Rural.Urban.Continuum.Code)
combined_psm_data$Medicaid_Expansion_Status <- as.factor(combined_psm_data$Medicaid_Expansion_Status)
combined_psm_data$Medicaid.Expansion.Status <- as.factor(combined_psm_data$Medicaid.Expansion.Status)



# # Split PSM dataset by era
 psm_split <- split(combined_psm_data, combined_psm_data$dx.yr_grouped)


# Generate KM plots only (no risk table) and store the plots
psm_plot_ggs <- map(names(psm_split), function(era) {
  data_subset <- psm_split[[era]]
  surv_obj <- survfit(Surv(capped_time, capped_event) ~ Medicaid_Expansion_Status, data = data_subset)
  
  plot_obj <- ggsurvplot(
    surv_obj,
    data = data_subset,
    title = era,
    conf.int = TRUE,
    pval = TRUE,
    risk.table = FALSE,  # No risk table for cowplot layout
    palette = jama_colors,
    xlim = c(0, 24),
    ylim = c(0.5, 1),
    break.time.by = 6,
    #censor.shape = 21,
    #censor.size = 2.5,
    #censor.fill = "white",
    xlab = "Months Since Diagnosis",
    ylab = "Overall Survival Probability",
    legend = "none",
    ggtheme = theme_classic(base_size = 10) +
      theme(
        legend.position = "bottom",  # Extractable
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 14),
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"),
        theme(text = element_text(family = "Times"))
      )
  )
  plot_obj$plot  # only return the ggplot
})

# Temporarily recreate a clean test plot with wide margins
test_plot <- ggsurvplot(
  survfit(Surv(capped_time, capped_event) ~ Medicaid_Expansion_Status, data = psm_split[[1]]),
  data = psm_split[[1]],
  palette = jama_colors,
  conf.int = FALSE,
  risk.table = FALSE,
  legend.title = "Medicaid Expansion Group",
  ggtheme = theme_minimal(base_size = 14) + 
    theme(
      legend.position = "right",
      text = element_text(family = "Times")
    )
)

# Extract full, uncropped legend
jama_legend <- get_legend(
  test_plot$plot + 
    theme(
      legend.position = "right",
      legend.box = "vertical",
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 12, face = "bold"),
      text = element_text(family = "Times")
    )
)


# Save the legend as a standalone PDF
ggsave("figures/km_legend_jama.pdf", jama_legend, width = 4.5, height = 2)


final_plot

ggsave("figures/km_plots_by_era_PSM_JAMA.pdf", final_plot, width = 12, height = 10)



# Log-rank test for comparing survival between groups (Medicaid expansion status)
log_rank_psm <- survdiff(Surv(capped_time, capped_event) ~ Medicaid_Expansion_Status + dx.yr_grouped, data = combined_psm_data)

log_rank_psm <- data.frame(
  dx.yr_grouped = unique(combined_psm_data$dx.yr_grouped),
  p_value = sapply(unique(combined_psm_data$dx.yr_grouped), function(period) {
    period_data <- combined_psm_data %>% filter(dx.yr_grouped == period)
    log_rank_test <- survdiff(Surv(capped_time, capped_event) ~ Medicaid_Expansion_Status, data = period_data)
    log_rank_test$pvalue
  })
)

# Save the results to a CSV file
write.csv(log_rank_psm, "figures/log_rank_test_PSM.csv", row.names = FALSE)


# Define the expansion groups
expansion_groups <- unique(combined_psm_data$Medicaid_Expansion_Status)

# Create an empty data frame to store pairwise log-rank test results
pairwise_log_rank_psm <- data.frame(dx.yr_grouped = character(),
                                    group_comparison = character(),
                                    p_value = numeric(),
                                    stringsAsFactors = FALSE)

for (period in unique(combined_psm_data$dx.yr_grouped)) {
  
  period_data <- combined_psm_data %>% filter(dx.yr_grouped == period)
  
  for (i in 1:(length(expansion_groups) - 1)) {
    for (j in (i + 1):length(expansion_groups)) {
      group1 <- expansion_groups[i]
      group2 <- expansion_groups[j]
      
      pairwise_data <- period_data %>%
        filter(Medicaid_Expansion_Status %in% c(group1, group2))
      
      if (length(unique(pairwise_data$Medicaid_Expansion_Status)) == 2) {
        
        # Optional: check group sizes to debug
        cat("\nPeriod:", period, " | Comparing:", group1, "vs", group2, "\n")
        cat("Group counts:\n")
        print(table(pairwise_data$Medicaid_Expansion_Status))
        
        log_rank_test <- tryCatch({
          survdiff(Surv(capped_time, capped_event) ~ Medicaid_Expansion_Status, data = pairwise_data)
        }, error = function(e) {
          message("Error in survdiff:", e)
          return(NULL)
        })
        
        if (!is.null(log_rank_test)) {
          p_value <- 1 - pchisq(log_rank_test$chisq, df = length(log_rank_test$n) - 1)
          
          pairwise_log_rank_psm <- rbind(pairwise_log_rank_psm,
                                         data.frame(dx.yr_grouped = period,
                                                    group_comparison = paste("Group", group1, "vs Group", group2),
                                                    p_value = p_value))
        }
      }
    }
  }
}

# Add multiple comparison corrections to the results
pairwise_log_rank_psm$holm <- p.adjust(pairwise_log_rank_psm$p_value, method = "holm")
pairwise_log_rank_psm$BH <- p.adjust(pairwise_log_rank_psm$p_value, method = "BH")  # Benjamini-Hochberg
pairwise_log_rank_psm$bonferroni <- p.adjust(pairwise_log_rank_psm$p_value, method = "bonferroni")

# Save the updated results with corrections
write.csv(pairwise_log_rank_psm, "figures/pairwise_log_rank_psm_with_corrections.csv", row.names = FALSE)


# Step 1: Cap survival at 24 months
combined_matched_data <- combined_matched_data %>%
  mutate(capped_survival = pmin(Survival.months, 24))

# Step 2: Group and summarize
survival_summary <- combined_matched_data %>%
  group_by(dx.yr_grouped, Medicaid_Expansion_Status) %>%
  summarise(
    mean_survival_months = mean(capped_survival, na.rm = TRUE),
    sd_survival_months = sd(capped_survival, na.rm = TRUE),
    n = n()
  ) %>%
  ungroup()

# Step 3: Print the table
print(survival_summary)

# Save the table as a CSV file
write.csv(survival_summary, "figures/PSM_months_survived_capped.csv", row.names = FALSE)


# Step 1: Cap survival time at 24 months (if not already capped elsewhere)
combined_matched_data <- combined_matched_data %>%
  mutate(capped_time = pmin(Survival.months, 24))

# Step 2: Create 24-month survival indicator
combined_matched_data <- combined_matched_data %>%
  mutate(alive_at_24mo = if_else(Survival.months >= 24, 1, 0))

# Step 3: Group and summarize percent alive
survival_24mo_summary <- combined_matched_data %>%
  group_by(dx.yr_grouped, Medicaid_Expansion_Status) %>%
  summarise(
    percent_alive_24mo = mean(alive_at_24mo, na.rm = TRUE) * 100,
    n = n()
  ) %>%
  ungroup()

# Step 4: Print and save
print(survival_24mo_summary)

write.csv(survival_24mo_summary, "figures/PSM_percent_alive_24mo.csv", row.names = FALSE)

# Calculate average outcome by diagnosis year and expansion status
avg_psm_data <- combined_psm_data %>%
  group_by(dx.yr, Medicaid_Expansion_Status) %>%
  summarize(avg_survival = mean(capped_time, na.rm = TRUE), .groups = 'drop')

avg_psm_data <- combined_psm_data %>%
  group_by(dx.yr, Medicaid_Expansion_Status) %>%
  summarise(
    avg_survival = mean(capped_time, na.rm = TRUE),
    n = n(),
    sd_survival = sd(capped_time, na.rm = TRUE),
    se_survival = sd_survival / sqrt(n),
    avg_survival_lower = avg_survival - 1.96 * se_survival,
    avg_survival_upper = avg_survival + 1.96 * se_survival
  ) %>%
  ungroup()

jama_colors <- c(
  "Non-Expansion"    = "black",
  "Early Expansion"  = "#003f5c",
  "2014 Expansion"   = "#bc5090",
  "Late Expansion"   = "#ffa600"
)

avg_psm_data$Medicaid_Expansion_Status <- factor(
  avg_psm_data$Medicaid_Expansion_Status,
  levels = c(0, 1, 2, 3),
  labels = c("Non-Expansion", "Early Expansion", "2014 Expansion", "Late Expansion")
)
# Create plot
exsurv_plot_psm <- ggplot(avg_psm_data, aes(x = dx.yr, y = avg_survival, color = factor(Medicaid_Expansion_Status))) +
  # Shaded eras
  annotate("rect", xmin = 2006, xmax = 2010.99, ymin = -Inf, ymax = Inf, fill = "#f4cccc", alpha = 0.2) +  # Light pink
  annotate("rect", xmin = 2011, xmax = 2013.99, ymin = -Inf, ymax = Inf, fill = "#cfe2f3", alpha = 0.2) +  # Light blue
  annotate("rect", xmin = 2014, xmax = 2016.99, ymin = -Inf, ymax = Inf, fill = "#d9ead3", alpha = 0.2) +  # Light green
  annotate("rect", xmin = 2017, xmax = 2019.99, ymin = -Inf, ymax = Inf, fill = "#fff2cc", alpha = 0.2) +  # Light yellow
  
  # Survival trend lines and points
  geom_line(aes(group = Medicaid_Expansion_Status), linewidth = 1.2) +
  geom_point(size = 3, shape = 21, fill = "white", stroke = 1.2) +
  
  
  # # CI error bars
  # geom_errorbar(aes(ymin = avg_survival_lower, ymax = avg_survival_upper), 
  #               width = 0.2, linewidth = 0.7, alpha = 0.7) +
  # Color and labels
  scale_color_manual(
    values = jama_colors,
    labels = c("Non-Expansion", "Early Expansion", "2014 Expansion", "Late Expansion")
  ) +
  labs(
    title = "Trends in Average Survival by Medicaid Expansion Group (Matched Cohorts)",
    x = "Diagnosis Year",
    y = "Average Survival (Months)",
    color = "Expansion Status"
  ) +
  
  # Limits
  scale_x_continuous(breaks = seq(2006, 2019, 1)) +
  coord_cartesian(xlim = c(2006, 2020), ylim = c(17, 23)) +
  
  # JAMA-style theme
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    legend.position = "top",
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    text = element_text(family = "Times")
  )

# Save
ggsave("figures/exsurv_plot_PSM_JAMAstyle.png", plot = exsurv_plot_psm, width = 10, height = 6, dpi = 300)


print(exsurv_plot_psm)



##### 9. Stage at Diagnosis of Matched Populations with Placebo Falsification Modeling ######


# Create a binary variable for early-stage diagnosis
cleaned_data <- cleaned_data %>%
  mutate(early_stage = ifelse(stage_group == "I", 1, 0))

# Summarize proportions
early_stage_summary <- cleaned_data %>%
  group_by(dx.yr_grouped, Medicaid_Expansion_Status) %>%
  summarise(
    n = n(),
    early_stage_n = sum(early_stage),
    early_stage_prop = mean(early_stage)
  ) %>%
  ungroup()

# View table
print(early_stage_summary)


# Calculate proportion and 95% CI (Wald method)
early_stage_summary <- cleaned_data %>%
  mutate(early_stage = ifelse(stage_group == "I", 1, 0)) %>%
  group_by(dx.yr_grouped, Medicaid_Expansion_Status) %>%
  summarise(
    n = n(),
    early_stage_n = sum(early_stage),
    early_stage_prop = mean(early_stage),
    se = sqrt((early_stage_prop * (1 - early_stage_prop)) / n),
    .groups = "drop"
  ) %>%
  mutate(
    ci_lower = pmax(0, early_stage_prop - 1.96 * se),
    ci_upper = pmin(1, early_stage_prop + 1.96 * se)
  )

# Optional: Set a colorblind-friendly palette
jama_colors <- c(
  "Non-Expansion"    = "black",
  "Early Expansion"  = "#003f5c",
  "2014 Expansion"   = "#bc5090",
  "Late Expansion"   = "#ffa600"
)

# Plot
plot_all = ggplot(early_stage_summary, 
       aes(x = dx.yr_grouped, 
           y = early_stage_prop, 
           fill = as.factor(Medicaid_Expansion_Status))) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  # Error bars (CI lines)
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), 
                position = position_dodge(width = 0.75), 
                width = 0.25, color = "black") +
  scale_fill_manual(
    values = jama_colors,
    name = "Medicaid Expansion Group",
    labels = c("0", "1", "2", "3")
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.05)),
    limits = c(0, 1)
  ) +
  labs(
    x = "Expansion Status by Diagnosis Era",
    y = "Proportion Diagnosed at Stage I"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 12),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(family = "Times")
  )

# Save the plot
ggsave("figures/early_stage_barplot_allgroups.png", 
       plot = plot_all, 
       width = 9, height = 5, dpi = 300)

# Optionally show the plot
print(plot_all)



# Create treatment and time binary variables (you may already have these)
cleaned_data <- cleaned_data %>%
  mutate(
    did_treat = ifelse(Medicaid_Expansion_Status == 1, 1, 0),
    did_time = ifelse(dx.yr_grouped %in% c("2014-2016", "2017-2019"), 1, 0)  # Adjust based on your study
  )

# Logistic regression for early stage diagnosis
early_stage_did_model <- glm(early_stage ~ did_treat * did_time, data = cleaned_data, family = binomial())

coeftest(early_stage_did_model, vcov = sandwich)
summary(early_stage_did_model)

# Tidy output
tidy(early_stage_did_model, conf.int = TRUE, exponentiate = TRUE)
write.csv(tidy(early_stage_did_model, conf.int = TRUE, exponentiate = TRUE), 
          "low_income/all_model_earlystage.csv", 
          row.names = FALSE)


ols_model <- lm(early_stage ~ did_treat * did_time, data = cleaned_data)
summary(ols_model)
tidy(ols_model, conf.int = TRUE, exponentiate = TRUE)
write.csv(tidy(ols_model, conf.int = TRUE, exponentiate = TRUE), 
          "low_income/OLSmodel_earlystage.csv", 
          row.names = FALSE)

# Subset to include only Medicaid Expansion group 1 and control (non-expansion) states
expansion_vs_control_data <- cleaned_data %>%
  filter(Medicaid_Expansion_Status %in% c(0, 1))  # 0 for control, 1 for expansion

# Create treatment and time binary variables for Expansion Group 1 vs Control
expansion_vs_control_data <- expansion_vs_control_data %>%
  mutate(
    did_treat = ifelse(Medicaid_Expansion_Status == 1, 1, 0),   # Expansion group (1) vs Control (0)
    did_time = ifelse(dx.yr_grouped %in% c("2014-2016", "2017-2019"), 1, 0)  # Post-expansion (2014-2016 & 2017-2019)
  )
# View summary of the model
summary(early_stage_did_model_expansion_vs_control)



# Logistic regression for early stage diagnosis
ADJUSTED_stage_did_model <- glm(early_stage ~ did_treat * did_time + age_category + Sex +
                               Race.and.origin.recode..NHW..NHB..NHAIAN..NHAPI..Hispanic. + 
                               Rural.Urban.Continuum.Code, data = cleaned_data, family = binomial())

coeftest(ADJUSTED_stage_did_model, vcov = sandwich)
summary(ADJUSTED_stage_did_model)

# Tidy output
tidy(ADJUSTED_stage_did_model, conf.int = TRUE, exponentiate = TRUE)
write.csv(tidy(ADJUSTED_stage_did_model, conf.int = TRUE, exponentiate = TRUE), 
          "low_income/ADJUSTED_model_earlystage.csv", 
          row.names = FALSE)



# First, create the tidy output
early_stage_did_model_expansion_vs_control_tidy <- tidy(early_stage_did_model_expansion_vs_control, conf.int = TRUE, exponentiate = TRUE)

# Then write it to CSV
write.csv(early_stage_did_model_expansion_vs_control_tidy, "low_income/earlystage_did.csv", row.names = FALSE)




# Subset pre-expansion years only
placebo_data <- cleaned_data %>%
  filter(dx.yr_grouped %in% c("2006-2010"),   # Ensure dx.yr_grouped is treated as a string
         Medicaid_Expansion_Status %in% c("Non-Expansion", "Early Expansion")) %>%
  mutate(
    fake_time = ifelse(dx.yr == 2010, 1, 0),  # Pretend 2010 is "post"
    treat = ifelse(Medicaid_Expansion_Status == "Early Expansion", 1, 0)
  ) %>%
  drop_na(early_stage)

# Run placebo DiD model
placebo_model <- glm(early_stage ~ treat * fake_time, data = placebo_data, family = binomial())

modelsummary(placebo_model, output = "low_income/placebo_stage_model_output.docx")

# Output
tidy(placebo_model, conf.int = TRUE, exponentiate = TRUE)


# Assuming you have already run tidy() on the placebo model
placebo_tidy <- tidy(placebo_model, conf.int = TRUE, exponentiate = TRUE)

# Save the tidy result as a CSV
write.csv(placebo_tidy, "low_income/placebofalse_stage.csv", row.names = FALSE)


#### Plot of Stage 1 survival, ends up kinda ugly so gonna just ignore this ### Bar plot is nicer 
# ggplot(cleaned_data,
#        aes(x = dx.yr, y = early_stage, color = factor(Medicaid_Expansion_Status))) +
#   stat_summary(fun = mean, geom = "line") +
#   labs(title = "Pre-treatment trends in Stage I diagnosis", y = "Proportion Stage I")

stage1dx <- ggplot(cleaned_data,
                   aes(x = dx.yr, y = (early_stage*100), color = factor(Medicaid_Expansion_Status))) +
                  stat_summary(fun = mean, geom = "line") +
  
  # Shaded eras
  annotate("rect", xmin = 2006, xmax = 2010.99, ymin = -Inf, ymax = Inf, fill = "#f4cccc", alpha = 0.2) +  # Light pink
  annotate("rect", xmin = 2011, xmax = 2013.99, ymin = -Inf, ymax = Inf, fill = "#cfe2f3", alpha = 0.2) +  # Light blue
  annotate("rect", xmin = 2014, xmax = 2016.99, ymin = -Inf, ymax = Inf, fill = "#d9ead3", alpha = 0.2) +  # Light green
  annotate("rect", xmin = 2017, xmax = 2019.99, ymin = -Inf, ymax = Inf, fill = "#fff2cc", alpha = 0.2) +  # Light yellow
  
  # Survival trend lines and points
  #geom_line(aes(), linewidth = 1.2) +
  geom_point(linewidth = 1.2, size = 3, shape = 21, fill = "white", stroke = 1.2) +
  
  # Color and labels
  scale_color_manual(
    values = jama_colors,
    labels = c("Non-Expansion", "Early Expansion", "2014 Expansion", "Late Expansion")
  ) +
  labs(
    title = "Trends in Proportion of Diagnoses at Stage I for Full Cohort",
    x = "Diagnosis Year",
    y = "% Population Diagnosed with Stage 1 Disease",
    color = "Expansion Status"
  ) +
  
  # Limits
  scale_x_continuous(breaks = seq(2006, 2019, 1)) +
  coord_cartesian(xlim = c(2006, 2020), ylim = c(40, 70)) +
  
  # JAMA-style theme
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    legend.position = "top",
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    text = element_text(family = "Times")
  )
print(stage1dx
      )

ggsave("figures/stage1_yrlyDiD_JAMAstyle.png", plot = stage1dx, width = 10, height = 6, dpi = 300)



stage1_summary_table <- cleaned_data %>%
  group_by(dx.yr, Medicaid_Expansion_Status) %>%
  summarise(
    n = n(),
    n_stage1 = sum(early_stage == 1),
    prop_stage1 = mean(early_stage),
    .groups = "drop"
  ) %>%
  mutate(
    percent_stage1 = round(prop_stage1 * 100, 1)
  )
stage1_summary_table <- stage1_summary_table %>%
  mutate(
    Expansion_Group = case_when(
      Medicaid_Expansion_Status == 0 ~ "Non-Expansion",
      Medicaid_Expansion_Status == 1 ~ "Early Expansion",
      Medicaid_Expansion_Status == 2 ~ "2014 Expansion",
      Medicaid_Expansion_Status == 3 ~ "Late Expansion",
      TRUE ~ as.character(Medicaid_Expansion_Status)
    )
  )

stage1_table_wide <- stage1_summary_table %>%
  dplyr::select(dx.yr, Expansion_Group, percent_stage1) %>%
  pivot_wider(names_from = Expansion_Group, values_from = percent_stage1)

library(knitr)

kable(stage1_table_wide, caption = "Proportion Diagnosed with Stage I NSCLC by Year and Expansion Group (%)")


group_counts <- cleaned_data %>%
  group_by(dx.yr, Medicaid_Expansion_Status) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(
    Expansion_Group = case_when(
      Medicaid_Expansion_Status == 0 ~ "Non-Expansion",
      Medicaid_Expansion_Status == 1 ~ "Early Expansion",
      Medicaid_Expansion_Status == 2 ~ "2014 Expansion",
      Medicaid_Expansion_Status == 3 ~ "Late Expansion"
    )
  ) %>%
  dplyr::select(dx.yr, Expansion_Group, n)

group_counts_wide <- group_counts %>%
  pivot_wider(names_from = Expansion_Group, values_from = n)

write.csv(group_counts_wide, "figures/group_counts_by_year.csv", row.names = FALSE)
#### EARLY VS CONTROL (1 v 0) #####

# Logistic regression for early-stage diagnosis (Expansion vs Control)
early_stage_did_model_expansion_vs_control <- glm(early_stage ~ did_treat * did_time, 
                                                  data = expansion_vs_control_data, 
                                                  family = binomial())

# Summarize proportions for Expansion vs Control
early_stage_summary_expansion_vs_control <- expansion_vs_control_data %>%
  group_by(dx.yr_grouped, Medicaid_Expansion_Status) %>%
  summarise(
    n = n(),
    early_stage_n = sum(early_stage),
    early_stage_prop = mean(early_stage)
  ) %>%
  ungroup()

# Plot early-stage proportions for Expansion vs Control
ggplot(early_stage_summary_expansion_vs_control, aes(x = dx.yr_grouped, y = early_stage_prop, fill = as.factor(Medicaid_Expansion_Status))) +
  geom_col(position = "dodge") +
  labs(x = "Diagnosis Year Group", y = "Proportion Stage I", fill = "Expansion Group") +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal() + theme(text = element_text(family = "Times"))






# Subset to include only Medicaid Expansion group 2 and control (non-expansion) states
expansion_group2_vs_control_data <- cleaned_data %>%
  filter(Medicaid_Expansion_Status %in% c(0, 2))  # 0 for control, 2 for expansion group 2


# Create treatment and time binary variables for Expansion Group 2 vs Control
expansion_group2_vs_control_data <- expansion_group2_vs_control_data %>%
  mutate(
    did_treat = ifelse(Medicaid_Expansion_Status == 2, 1, 0),   # Expansion group 2 (treated) vs Control (non-expansion)
    did_time = ifelse(dx.yr_grouped %in% c("2014-2016", "2017-2019"), 1, 0)  # Post-expansion (2014-2016 & 2017-2019)
  )


# Logistic regression for early-stage diagnosis (Expansion Group 2 vs Control)
early_stage_did_model_expansion_group2_vs_control <- glm(early_stage ~ did_treat * did_time, 
                                                         data = expansion_group2_vs_control_data, 
                                                         family = binomial())

# View summary of the model
summary(early_stage_did_model_expansion_group2_vs_control)

# First, create the tidy output
early_stage_did_model_expansion_group2_vs_control_tidy <- tidy(early_stage_did_model_expansion_group2_vs_control, conf.int = TRUE, exponentiate = TRUE)

# Then write it to CSV
write.csv(early_stage_did_model_expansion_group2_vs_control_tidy, "low_income/earlystage_did_2.csv", row.names = FALSE)

# Summarize proportions for Expansion Group 2 vs Control
early_stage_summary_expansion_group2_vs_control <- expansion_group2_vs_control_data %>%
  group_by(dx.yr_grouped, Medicaid_Expansion_Status) %>%
  summarise(
    n = n(),
    early_stage_n = sum(early_stage),
    early_stage_prop = mean(early_stage)
  ) %>%
  ungroup()

# Plot early-stage proportions for Expansion Group 2 vs Control
ggplot(early_stage_summary_expansion_group2_vs_control, aes(x = dx.yr_grouped, y = early_stage_prop, fill = as.factor(Medicaid_Expansion_Status))) +
  geom_col(position = "dodge") +
  labs(x = "Diagnosis Year Group", y = "Proportion Stage I", fill = "Expansion Group 2") +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal() + theme(text = element_text(family = "Times"))




avg_psm_by_stage <- combined_psm_data %>%
  group_by(dx.yr, Medicaid_Expansion_Status, stage_group) %>%  # Assuming you have `stage_group` as I vs II/IIIA
  summarize(avg_survival = mean(capped_time, na.rm = TRUE), .groups = 'drop')

gg_stage_plot <- ggplot(avg_psm_by_stage, aes(x = dx.yr, y = avg_survival, color = factor(Medicaid_Expansion_Status))) +
  geom_line(aes(group = Medicaid_Expansion_Status), linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~stage_group, nrow = 1) +  # Facet by stage
  scale_color_manual(
    values = c("black", "#bc5090"),
    labels = c("0", "2")
  ) +
  labs(
    title = "Average Survival by Stage and Medicaid Expansion (PSM)",
    x = "Diagnosis Year", y = "Avg Survival (Months)",
    color = "Expansion Status"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    text = element_text(family = "Times")
  )

# ggsave("figures/exsurv_plot_PSM_byStage.png", plot = gg_stage_plot, width = 12, height = 6, dpi = 300)


# placebo_model_pre <- glm(early_stage ~ did_treat * did_time, 
#                          data = cleaned_data %>% filter(dx.yr_grouped %in% c("2011-2013")),
#                          family = binomial())
# 
# coeftest(placebo_model_pre, vcov = sandwich)
# summary(placebo_model_pre)
# 
# # Tidy output
# tidy(placebo_model_pre, conf.int = TRUE, exponentiate = TRUE)
# write.csv(tidy(placebo_model_pre, conf.int = TRUE, exponentiate = TRUE), 
#           "low_income/placeboModel_pre.csv", 
#           row.names = FALSE)







#### 11. NSCLC Differential Benefit Analysis ####
### Changes in Healthcare Disparities for total Pre/Post populations ### 

# set Lung Income quintiles to be equal 

seer_final <- seer_final %>%
  mutate(income_quintile_equal = ntile(midpoint, 5)) %>%
  mutate(income_quintile_equal = factor(income_quintile_equal,
                                        levels = 1:5,
                                        labels = c("Q1", "Q2", "Q3", "Q4", "Q5")))

# seer_final$income_quintile <- relevel(seer_final$income_quintile, ref = "Q3")

seer_final <- seer_final %>%
  filter(!grepl("<NA>", trimws(income_quintile_equal)))

# # View result
table(seer_final$income_quintile_equal, useNA = "ifany")

# Step 1: Define expansion timing for each group # already done
seer_final <- seer_final %>%
  mutate(expansion1_time = if_else(
    Medicaid_Expansion_Status == 1,
    case_when(
      dx.yr_grouped %in% c("2006-2010") ~ "Pre",
      dx.yr_grouped %in% c("2011-2013", "2014-2016", "2017-2019") ~ "Post",
      TRUE ~ NA_character_
    ),
    NA_character_
  ))

seer_final <- seer_final %>%
  mutate(expansion2_time = if_else(
    Medicaid_Expansion_Status == 2,
    case_when(
      dx.yr_grouped %in% c("2006-2010", "2011-2013") ~ "Pre",
      dx.yr_grouped %in% c("2014-2016", "2017-2019") ~ "Post",
      TRUE ~ NA_character_
    ),
    NA_character_
  ))

seer_final <- seer_final %>%
  mutate(expansion3_time = if_else(
    Medicaid_Expansion_Status == 3,
    case_when(
      dx.yr_grouped %in% c("2006-2010", "2011-2013", "2014-2016") ~ "Pre",
      dx.yr_grouped %in% c("2017-2019") ~ "Post",
      TRUE ~ NA_character_
    ),
    NA_character_
  ))

# Step 2: Define Pre/Post periods across all expansion groups
seer_disp <- seer_final %>%
  mutate(
    Time_Period = case_when(
      (Medicaid_Expansion_Status == "1" & expansion1_time == "Pre") ~ "Pre",
      (Medicaid_Expansion_Status == "1" & expansion1_time == "Post") ~ "Post",
      (Medicaid_Expansion_Status == "2" & expansion2_time == "Pre") ~ "Pre",
      (Medicaid_Expansion_Status == "2" & expansion2_time == "Post") ~ "Post",
      (Medicaid_Expansion_Status == "3" & expansion3_time == "Pre") ~ "Pre",
      (Medicaid_Expansion_Status == "3" & expansion3_time == "Post") ~ "Post",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Time_Period))

seer_disp$Time_Period <- factor(seer_disp$Time_Period)

seer_disp$Time_Period <- relevel(seer_disp$Time_Period, ref = "Pre")


# Set references
seer_disp <- seer_disp %>%
  mutate(
    income_quintile_equal = relevel(factor(income_quintile_equal), ref = "Q5"),
    Race.and.origin.recode..NHW..NHB..NHAIAN..NHAPI..Hispanic.= relevel(factor(
      Race.and.origin.recode..NHW..NHB..NHAIAN..NHAPI..Hispanic.),ref = "Non-Hispanic White"),
    stage_group      = relevel(factor(stage_group),      ref = "II+III"),
    age_category     = relevel(factor(age_category),     ref = "Older Adult"),
    Marital_Grouped  = relevel(factor(Marital_Grouped),  ref = "Single/Unknown"),
    Sex = relevel(factor(Sex), ref = "Female"),
    Rural.Urban.Continuum.Code = relevel(factor(Rural.Urban.Continuum.Code), ref = "Counties in metropolitan areas ge 1 million pop")
  )

# Step 3: Create a Cox model with interaction terms to test for disparity change
cox_model <- coxph(
  Surv(capped_time, capped_event) ~ 
    (age_category + Marital_Grouped + Sex
     + Race.and.origin.recode..NHW..NHB..NHAIAN..NHAPI..Hispanic. + Rural.Urban.Continuum.Code
     + income_quintile_equal 
     + surgery_status
     + stage_group
     + rad_group
     ) * Time_Period,
  data = seer_disp,
  robust = TRUE
)

print(cox_model)

# Step 4: Extract interaction terms (i.e., change in disparities pre/post)
interaction_results <- tidy(cox_model, conf.int = TRUE) %>%
  filter(grepl(":Time_PeriodPost", term)) %>%
  mutate(
    HR = exp(estimate),
    CI_lower = exp(conf.low),
    CI_upper = exp(conf.high),
    sig = ifelse(p.value < 0.05, "Yes", "No")
  ) %>%
  dplyr::select(term, HR, CI_lower, CI_upper, p.value, sig)

# View interaction effects
print(interaction_results)

View(interaction_results)

write.csv(interaction_results, "coxdata/pre-postDisparities.csv")


#### 12. Number and Percent Censored at Each Time X ####
# 48 months
# Define 4 years in months
four_years_months <- 12 * 4

# Add eligibility column
combined_matched_data <- combined_matched_data %>%
  mutate(eligible_4yr = Survival.months >= four_years_months)

# Get counts and percentages
eligible_summary_48 <- combined_matched_data %>%
  group_by(Medicaid_Expansion_Status) %>%
  summarize(
    N_total = n(),
    N_eligible_4yr = sum(eligible_4yr),
    Pct_eligible_4yr = round(100 * mean(eligible_4yr), 1)
  )
print(eligible_summary_48
      )
write.csv(eligible_summary_48, "coxdata/Percent_Eligible_48.csv", row.names = FALSE)

# 60 Months
# Define 5 years in months
five_years_months <- 12 * 5

# Add eligibility column
combined_matched_data <- combined_matched_data %>%
  mutate(eligible_5yr = Survival.months >= five_years_months)

# Get counts and percentages
eligible_summary_60 <- combined_matched_data %>%
  group_by(Medicaid_Expansion_Status) %>%
  summarize(
    N_total = n(),
    N_eligible_5yr = sum(eligible_5yr),
    Pct_eligible_5yr = round(100 * mean(eligible_5yr), 1)
  )
print(eligible_summary_60
)
write.csv(eligible_summary_60, "coxdata/Percent_Eligible_60.csv", row.names = FALSE)

# 24 months
# Define 2 years in months
two_years_months <- 12 * 2

# Add eligibility column
combined_matched_data <- combined_matched_data %>%
  mutate(eligible_2yr = Survival.months >= two_years_months)

# Get counts and percentages
eligible_summary_24 <- combined_matched_data %>%
  group_by(Medicaid_Expansion_Status) %>%
  summarize(
    N_total = n(),
    N_eligible_2yr = sum(eligible_2yr),
    Pct_eligible_2yr = round(100 * mean(eligible_2yr), 1)
  )
print(eligible_summary_24
)
write.csv(eligible_summary_24, "coxdata/Percent_Eligible_24.csv", row.names = FALSE)


#### 10. STATE MAP ####
install.packages("usmap")
install.packages("ggplot2")


# Your Medicaid expansion classification
state_expansion <- data.frame(
  state = c(state.name),  # built-in state names
  status = NA_character_
)

# Define statuses manually based on your DiD design
early_states <- c("California", "Connecticut", "New Jersey", "Washington")
mid_states <- c("Kentucky", "New Mexico", "Hawaii", "Iowa")
late_states <- c("Louisiana", "Alaska")
non_expanded <- c("Utah", "Georgia")
not_considered <- setdiff(state.name, c(early_states, mid_states, late_states, non_expanded))

# Assign categories
state_expansion$status[state_expansion$state %in% early_states] <- "Early Expansion (2011)"
state_expansion$status[state_expansion$state %in% mid_states] <- "On-Time Expansion (2014)"
state_expansion$status[state_expansion$state %in% late_states] <- "Late Expansion (2016)"
state_expansion$status[state_expansion$state %in% non_expanded] <- "Non-Expansion"
state_expansion$status[state_expansion$state %in% not_considered] <- "N/A"

# Plot
state_map <- usmap::plot_usmap(data = state_expansion, values = "status") +
  scale_fill_manual(
    values = c(
      "Non-Expansion" = "gray2" ,          
      "Early Expansion (2011)" = "#003f5c",    
      "On-Time Expansion (2014)" = "#bc5090",      
      "Late Expansion (2016)" =  "#ffa600",
     "N/A" = "lightgrey"
    ),
    name = "Medicaid Expansion"
  ) +
  labs(title = "Medicaid Expansion Status by State") +
  theme(legend.position = "bottom", text = element_text(family = "Times")) 

ggsave("figures/StatePlot.pdf", 
       plot = state_map, 
       width = 15, height = 12)



#### 13. 4 YEAR SURVIVAL ####
# Refit Cox Proportions Hazard with Propensity Score Matching #

# overall (not inclusive of expansion status) does not include Status + dx.yr_grouped
cleaned_data$capped_time_48 <- pmin(cleaned_data$Survival.months, 48)
cleaned_data$capped_event_48 <- ifelse(cleaned_data$Survival.months > 48, 0, cleaned_data$Event)

#Combine previous 3 Matched Data tables into 1 group
# Combine the matched datasets and add a group identifier
matched_data_1$group <- "Group 1"
matched_data_2$group <- "Group 2"
matched_data_3$group <- "Group 3"

# Combine all matched datasets into one
combined_matched_data <- bind_rows(matched_data_1, matched_data_2, matched_data_3)


combined_matched_data$capped_time_48 <- pmin(combined_matched_data$Survival.months, 48)
combined_matched_data$capped_event_48 <- ifelse(combined_matched_data$Survival.months > 48, 0, combined_matched_data$Event)
combined_matched_data <- combined_matched_data[combined_matched_data$capped_time_48 > 0, ]

cox_psm_final <- coxph(
  Surv(capped_time_48, capped_event_48) ~ age_category + Sex + Marital_Grouped +
    Rural.Urban.Continuum.Code + Race.and.origin.recode..NHW..NHB..NHAIAN..NHAPI..Hispanic. + Medicaid_Expansion_Status,
  data = combined_matched_data,
  robust = TRUE  # Adds robust standard errors
)

capped_fin_cox_psm_final_sum <- tidy(cox_psm_final, conf.int = TRUE) %>%
  mutate(HR = exp(estimate),  # Convert coefficients to hazard ratios
         CI_lower = exp(conf.low),
         CI_upper = exp(conf.high)) %>%
  dplyr::select(term, HR, CI_lower, CI_upper, p.value)

print(capped_fin_cox_psm_final_sum)

write.csv(capped_fin_cox_psm_final_sum, "coxdata/48/psm_covariates.csv")

# Group by Group Difference in Differences analysis
# 0 vs 1
# Create expansion1_time variable (Pre = 2006-2010, Post = 2011-2016)
matched_data_1 <- matched_data_1 %>%
  mutate(expansion1_time = case_when(
    dx.yr_grouped %in% c("2006-2010") ~ "Pre",
    dx.yr_grouped %in% c("2011-2013", "2014-2016", "2017-2019") ~ "Post",
    TRUE ~ NA_character_  # Assign NA for years outside the study period
  ))

# Create expansion1_group variable (1 = expansion1 states, 0 = non-expansion states)
matched_data_1 <- matched_data_1 %>%
  mutate(expansion1_group = case_when(
    Medicaid_Expansion_Status == 1 ~ 1,  # Expansion 1 (early)
    Medicaid_Expansion_Status == 0 ~ 0,  # Non-expansion (control)
    TRUE ~ NA_real_  # Assign NA for other groups
  ))

# Subset to only include subjects from Expansion Group 0 and 1
matched_data_1_did <- matched_data_1 %>%
  filter(Medicaid_Expansion_Status %in% c(0, 1)) %>%
  # Create a binary treatment variable:
  mutate(did_treat = ifelse(Medicaid_Expansion_Status == 1, 1, 0),
         # Create a time variable: "Pre" for 2006-2010, "Post" for later years.
         did_time = ifelse(dx.yr_grouped == "2006-2010", "Pre", "Post"))

table(matched_data_1_did$dx.yr_grouped)


# Convert did_time and did_treat to factors with explicit ordering:
matched_data_1_did <- matched_data_1_did %>%
  mutate(
    did_time = factor(did_time, levels = c("Pre", "Post")),
    did_treat = factor(did_treat, levels = c(0, 1))
  )

matched_data_1_did$capped_time_48 <- pmin(matched_data_1_did$Survival.months, 48)
matched_data_1_did$capped_event_48 <- ifelse(matched_data_1_did$Survival.months > 48, 0, matched_data_1_did$Event)



capped_fin_cox_psm_final_did_1 <- coxph(Surv(capped_time_48, capped_event_48) ~ did_time * did_treat + cluster(Medicaid.Expansion.Status),
                                        data = matched_data_1_did, robust = TRUE)


summary(capped_fin_cox_psm_final_did_1)

capped_fin_cox_psm_final_did_1 <- tidy(capped_fin_cox_psm_final_did_1, conf.int = TRUE) %>%
  mutate(HR = exp(estimate),  # Convert coefficients to hazard ratios
         CI_lower = exp(conf.low),
         CI_upper = exp(conf.high)) %>%
  dplyr::select(term, HR, CI_lower, CI_upper, p.value)

print(capped_fin_cox_psm_final_did_1)

write.csv(capped_fin_cox_psm_final_did_1, "coxdata/48/psm_0v1.csv")



# 0 vs 2
# Create expansion1_time variable (Pre = 2006-2013, Post = 2014-2019)
matched_data_2 <- matched_data_2 %>%
  mutate(expansion2_time = case_when(
    dx.yr_grouped %in% c("2006-2010", "2011-2013") ~ "Pre",
    dx.yr_grouped %in% c("2014-2016", "2017-2019") ~ "Post",
    TRUE ~ NA_character_  # Assign NA for years outside the study period
  ))

# Create expansion1_group variable (1 = expansion2 states, 0 = non-expansion states)
matched_data_2 <- matched_data_2 %>%
  mutate(expansion2_group = case_when(
    Medicaid_Expansion_Status == 2 ~ 1,  # Expansion 2 (early)
    Medicaid_Expansion_Status == 0 ~ 0,  # Non-expansion (control)
    TRUE ~ NA_real_  # Assign NA for other groups
  ))



# Subset to only include subjects from Expansion Group 0 and 2
matched_data_2_did <- matched_data_2 %>%
  filter(Medicaid_Expansion_Status %in% c(0, 2)) %>%
  # Assign treatment variable: 1 for expansion (2), 0 for non-expansion (0)
  mutate(did_treat_2 = ifelse(Medicaid_Expansion_Status == 2, 1, 0),
         # Define Pre (2006-2013) and Post (2014+)
         did_time_2 = ifelse(dx.yr_grouped %in% c("2006-2010", "2011-2013"), "Pre", "Post"))

# Convert did_time and did_treat to factors with explicit ordering:
matched_data_2_did <- matched_data_2_did %>%
  mutate(
    did_time_2 = factor(did_time_2, levels = c("Pre", "Post")),
    did_treat_2 = factor(did_treat_2, levels = c(0, 1))
  )


matched_data_2_did <- matched_data_2_did %>%
  filter(!grepl("Unknown Race", trimws(Race.and.origin.recode..NHW..NHB..NHAIAN..NHAPI..Hispanic.)))

matched_data_2_did$capped_time_48 <- pmin(matched_data_2_did$Survival.months, 48)
matched_data_2_did$capped_event_48 <- ifelse(matched_data_2_did$Survival.months > 48, 0, matched_data_2_did$Event)



capped_fin_cox_psm_final_did_2 <- coxph(Surv(capped_time_48, capped_event_48) ~ did_time_2 * did_treat_2 + cluster(Medicaid.Expansion.Status),
                                        data = matched_data_2_did, robust = TRUE)

summary(capped_fin_cox_psm_final_did_2)

capped_fin_cox_psm_final_did_2 <- tidy(capped_fin_cox_psm_final_did_2, conf.int = TRUE) %>%
  mutate(HR = exp(estimate),  # Convert coefficients to hazard ratios
         CI_lower = exp(conf.low),
         CI_upper = exp(conf.high)) %>%
  dplyr::select(term, HR, CI_lower, CI_upper, p.value)

print(capped_fin_cox_psm_final_did_2)

write.csv(capped_fin_cox_psm_final_did_2, "coxdata/48/psm_0v2.csv")


# 0 vs 3
# Create expansion1_time variable (Pre = 2006-2013, Post = 2014-2019)
matched_data_3 <- matched_data_3 %>%
  mutate(expansion2_time = case_when(
    dx.yr_grouped %in% c("2006-2010", "2011-2013") ~ "Pre",
    dx.yr_grouped %in% c("2014-2016", "2017-2019") ~ "Post",
    TRUE ~ NA_character_  # Assign NA for years outside the study period
  ))

# Create expansion1_group variable (1 = expansion2 states, 0 = non-expansion states)
matched_data_3 <- matched_data_3 %>%
  mutate(expansion2_group = case_when(
    Medicaid_Expansion_Status == 3 ~ 1,  # Expansion 2 (early)
    Medicaid_Expansion_Status == 0 ~ 0,  # Non-expansion (control)
    TRUE ~ NA_real_  # Assign NA for other groups
  ))



# Subset to only include subjects from Expansion Group 0 and 2
matched_data_3_did <- matched_data_3 %>%
  filter(Medicaid_Expansion_Status %in% c(0, 3)) %>%
  # Assign treatment variable: 1 for expansion (2), 0 for non-expansion (0)
  mutate(did_treat_3 = ifelse(Medicaid_Expansion_Status == 3, 1, 0),
         # Define Pre (2006-2013) and Post (2014+)
         did_time_3 = ifelse(dx.yr_grouped %in% c("2006-2010", "2011-2013"), "Pre", "Post"))

# Convert did_time and did_treat to factors with explicit ordering:
matched_data_3_did <- matched_data_3_did %>%
  mutate(
    did_time_3 = factor(did_time_3, levels = c("Pre", "Post")),
    did_treat_3 = factor(did_treat_3, levels = c(0, 1))
  )

matched_data_3_did <- matched_data_3_did %>%
  filter(!grepl("Unknown Race", trimws(Race.and.origin.recode..NHW..NHB..NHAIAN..NHAPI..Hispanic.)))

matched_data_3_did$capped_time_48 <- pmin(matched_data_3_did$Survival.months, 48)
matched_data_3_did$capped_event_48 <- ifelse(matched_data_3_did$Survival.months > 48, 0, matched_data_3_did$Event)



capped_fin_cox_psm_final_did_3 <- coxph(Surv(capped_time_48, capped_event_48) ~ did_time_3 * did_treat_3 + cluster(Medicaid.Expansion.Status),
                                        data = matched_data_3_did, robust = TRUE)

summary(capped_fin_cox_psm_final_did_3)

capped_fin_cox_psm_final_did_3 <- tidy(capped_fin_cox_psm_final_did_3, conf.int = TRUE) %>%
  mutate(HR = exp(estimate),  # Convert coefficients to hazard ratios
         CI_lower = exp(conf.low),
         CI_upper = exp(conf.high)) %>%
  dplyr::select(term, HR, CI_lower, CI_upper, p.value)

print(capped_fin_cox_psm_final_did_3)

write.csv(capped_fin_cox_psm_final_did_3, "coxdata/48/psm_0v3.csv")



### EARLY EXPANSION IMPLEMENTATION/POST-IMPLEMENTATION 48 MONTH SURVIVAL

# EARLY - consider early months (post-implementation periods)
# Step 1: Subset data to focus on 2014-2019
matched_data_1_did_early_1 <- matched_data_1 %>%
  filter(dx.yr_grouped %in% c("2006-2010", "2011-2013"))

# Step 2: Update the time variable for 2014-2016 and 2017-2019 periods
matched_data_1_did_early_1 <- matched_data_1_did_early_1 %>%
  mutate(did_time = ifelse(dx.yr_grouped %in% c("2011-2013"), "Post", "Pre"))

# Create treatment variable (same as before)
matched_data_1_did_early_1 <- matched_data_1_did_early_1 %>%
  mutate(did_treat = ifelse(Medicaid_Expansion_Status == 1, 1, 0))

# Step 3: Ensure 'did_time' is a factor with 'Pre' as reference
matched_data_1_did_early_1 <- matched_data_1_did_early_1 %>%
  mutate(
    did_time = factor(did_time, levels = c("Pre", "Post")),
    did_treat = factor(did_treat, levels = c(0, 1))
    
  )

# Check reference levels
contrasts(matched_data_1_did_early_1$did_time)
contrasts(matched_data_1_did_early_1$did_treat)

matched_data_1_did_early_1$capped_time_48 <- pmin(matched_data_1_did_early_1$Survival.months, 48)
matched_data_1_did_early_1$capped_event_48 <- ifelse(matched_data_1_did_early_1$Survival.months > 48, 0, matched_data_1_did_early_1$Event)

# Step 4: Run Cox regression model
capped_fin_cox_psm_final_early_1 <- coxph(Surv(capped_time_48, capped_event_48) ~ did_time * did_treat+ cluster(Medicaid.Expansion.Status),
                                          data = matched_data_1_did_early_1, robust = TRUE)

# Step 5: Summary of the model
summary(capped_fin_cox_psm_final_early_1)

capped_fin_cox_psm_final_did_early_1 <- tidy(capped_fin_cox_psm_final_early_1, conf.int = TRUE) %>%
  mutate(HR = exp(estimate),  # Convert coefficients to hazard ratios
         CI_lower = exp(conf.low),
         CI_upper = exp(conf.high)) %>%
  dplyr::select(term, HR, CI_lower, CI_upper, p.value)

print(capped_fin_cox_psm_final_did_early_1)

write.csv(capped_fin_cox_psm_final_did_early_1, "coxdata/48/psm_0v1_early.csv")


# LATER
# Step 1: Subset data to focus on 2014-2019
matched_data_1_did_later_1 <- matched_data_1 %>%
  filter(dx.yr_grouped %in% c("2006-2010", "2014-2016", "2017-2019"))

# Step 2: Update the time variable for 2014-2016 and 2017-2019 periods
matched_data_1_did_later_1 <- matched_data_1_did_later_1 %>%
  mutate(did_time = ifelse(dx.yr_grouped %in% c("2014-2016", "2017-2019"), "Post", "Pre"))

# Create treatment variable (same as before)
matched_data_1_did_later_1 <- matched_data_1_did_later_1 %>%
  mutate(did_treat = ifelse(Medicaid_Expansion_Status == 1, 1, 0))

# Step 3: Ensure 'did_time' is a factor with 'Pre' as reference
matched_data_1_did_later_1 <- matched_data_1_did_later_1 %>%
  mutate(
    did_time = factor(did_time, levels = c("Pre", "Post")),
    did_treat = factor(did_treat, levels = c(0, 1))
    
  )

# Check reference levels
contrasts(matched_data_1_did_later_1$did_time)
contrasts(matched_data_1_did_later_1$did_treat)

matched_data_1_did_later_1$capped_time_48 <- pmin(matched_data_1_did_later_1$Survival.months, 48)
matched_data_1_did_later_1$capped_event_48 <- ifelse(matched_data_1_did_later_1$Survival.months > 48, 0, matched_data_1_did_later_1$Event)

# Step 4: Run Cox regression model
capped_fin_cox_psm_final_later_1 <- coxph(Surv(capped_time_48, capped_event_48) ~ did_time * did_treat + cluster(Medicaid.Expansion.Status),
                                          data = matched_data_1_did_later_1, robust = TRUE)

# Step 5: Summary of the model
summary(capped_fin_cox_psm_final_later_1)

capped_fin_cox_psm_final_did_later_1 <- tidy(capped_fin_cox_psm_final_later_1, conf.int = TRUE) %>%
  mutate(HR = exp(estimate),  # Convert coefficients to hazard ratios
         CI_lower = exp(conf.low),
         CI_upper = exp(conf.high)) %>%
  dplyr::select(term, HR, CI_lower, CI_upper, p.value)

print(capped_fin_cox_psm_final_did_later_1)
print(capped_fin_cox_psm_final_did_early_1)
write.csv(capped_fin_cox_psm_final_did_later_1, "coxdata/48/psm_0v1_later.csv")



### 48 month survival for implementation/postimplemntation periods for on-time expansion group##

# EARLY - consider early months (post-implementation periods)
# Step1: Subset data to focus on 2017-2019
matched_data_2_did_early_2 <- matched_data_2 %>%
  filter(dx.yr_grouped %in% c("2006-2010", "2011-2013", "2014-2016"))

# Step 2: Update the time variable for 2014-2016 and 2017-2019 periods
matched_data_2_did_early_2 <- matched_data_2_did_early_2 %>%
  mutate(did_time_2 = ifelse(dx.yr_grouped %in% c("2014-2016"), "Post", "Pre"))

# Create treatment variable (same as before)
matched_data_2_did_early_2 <- matched_data_2_did_early_2 %>%
  mutate(did_treat_2 = ifelse(Medicaid_Expansion_Status == 2, 1, 0))

# Step 3: Ensure 'did_time' is a factor with 'Pre' as reference
matched_data_2_did_early_2 <- matched_data_2_did_early_2 %>%
  mutate(
    did_time_2 = factor(did_time_2, levels = c("Pre", "Post")),
    did_treat_2 = factor(did_treat_2, levels = c(0, 1))
  )

# Check reference levels
contrasts(matched_data_2_did_early_2$did_time_2)
contrasts(matched_data_2_did_early_2$did_treat_2)

matched_data_2_did_early_2$capped_time_48 <- pmin(matched_data_2_did_early_2$Survival.months, 48)
matched_data_2_did_early_2$capped_event_48 <- ifelse(matched_data_2_did_early_2$Survival.months > 48, 0, matched_data_2_did_early_2$Event)


# Step 4: Run Cox regression model
capped_fin_cox_psm_final_early_2 <- coxph(Surv(capped_time_48, capped_event_48) ~ did_time_2 * did_treat_2 + cluster(Medicaid.Expansion.Status),
                                          data = matched_data_2_did_early_2, robust = TRUE)

# Step 5: Summary of the model
summary(capped_fin_cox_psm_final_early_2)

capped_fin_cox_psm_final_did_early_2 <- tidy(capped_fin_cox_psm_final_early_2, conf.int = TRUE) %>%
  mutate(HR = exp(estimate),  # Convert coefficients to hazard ratios
         CI_lower = exp(conf.low),
         CI_upper = exp(conf.high)) %>%
  dplyr::select(term, HR, CI_lower, CI_upper, p.value)

print(capped_fin_cox_psm_final_did_early_2)

write.csv(capped_fin_cox_psm_final_did_early_2, "coxdata/48/psm_0v2_early.csv")




# LATER
# Step1: Subset data to focus on 2017-2019
matched_data_2_did_later_2 <- matched_data_2 %>%
  filter(dx.yr_grouped %in% c("2006-2010", "2011-2013", "2017-2019"))

# Step 2: Update the time variable for 2014-2016 and 2017-2019 periods
matched_data_2_did_later_2 <- matched_data_2_did_later_2 %>%
  mutate(did_time_2 = ifelse(dx.yr_grouped %in% c("2017-2019"), "Post", "Pre"))

# Create treatment variable (same as before)
matched_data_2_did_later_2 <- matched_data_2_did_later_2 %>%
  mutate(did_treat_2 = ifelse(Medicaid_Expansion_Status == 2, 1, 0))

# Step 3: Ensure 'did_time' is a factor with 'Pre' as reference
matched_data_2_did_later_2 <- matched_data_2_did_later_2 %>%
  mutate(
    did_time_2 = factor(did_time_2, levels = c("Pre", "Post")),
    did_treat_2 = factor(did_treat_2, levels = c(0, 1))
  )

# Check reference levels
contrasts(matched_data_2_did_later_2$did_time_2)
contrasts(matched_data_2_did_later_2$did_treat_2)

matched_data_2_did_later_2$capped_time_48 <- pmin(matched_data_2_did_later_2$Survival.months, 48)
matched_data_2_did_later_2$capped_event_48 <- ifelse(matched_data_2_did_later_2$Survival.months > 48, 0, matched_data_2_did_later_2$Event)



# Step 4: Run Cox regression model
capped_fin_cox_psm_final_later_2 <- coxph(Surv(capped_time_48, capped_event_48) ~ did_time_2 * did_treat_2 + cluster(Medicaid.Expansion.Status),
                                          data = matched_data_2_did_later_2, robust = TRUE)

# Step 5: Summary of the model
summary(capped_fin_cox_psm_final_later_2)

capped_fin_cox_psm_final_did_later_2 <- tidy(capped_fin_cox_psm_final_later_2, conf.int = TRUE) %>%
  mutate(HR = exp(estimate),  # Convert coefficients to hazard ratios
         CI_lower = exp(conf.low),
         CI_upper = exp(conf.high)) %>%
  dplyr::select(term, HR, CI_lower, CI_upper, p.value)

print(capped_fin_cox_psm_final_did_later_2)

write.csv(capped_fin_cox_psm_final_did_later_2, "coxdata/48/psm_0v2_later.csv")



#### 14. 5 YEAR SURVIVAL ####
# Refit Cox Proportions Hazard with Propensity Score Matching #

# overall (not inclusive of expansion status) does not include Status + dx.yr_grouped
cleaned_data$capped_time_60 <- pmin(cleaned_data$Survival.months, 60)
cleaned_data$capped_event_60 <- ifelse(cleaned_data$Survival.months > 60, 0, cleaned_data$Event)

#Combine previous 3 Matched Data tables into 1 group
# Combine the matched datasets and add a group identifier
matched_data_1$group <- "Group 1"
matched_data_2$group <- "Group 2"
matched_data_3$group <- "Group 3"

# Combine all matched datasets into one
combined_matched_data <- bind_rows(matched_data_1, matched_data_2, matched_data_3)


combined_matched_data$capped_time_60 <- pmin(combined_matched_data$Survival.months, 60)
combined_matched_data$capped_event_60 <- ifelse(combined_matched_data$Survival.months > 60, 0, combined_matched_data$Event)
combined_matched_data <- combined_matched_data[combined_matched_data$capped_time_60 > 0, ]

cox_psm_final <- coxph(
  Surv(capped_time_60, capped_event_60) ~ age_category + Sex + Marital_Grouped +
    Rural.Urban.Continuum.Code + Race.and.origin.recode..NHW..NHB..NHAIAN..NHAPI..Hispanic. + Medicaid_Expansion_Status,
  # cluster(Medicaid.Expansion.Status), 
  data = combined_matched_data,
  robust = TRUE  # Adds robust standard errors
)

capped_fin_cox_psm_final_sum <- tidy(cox_psm_final, conf.int = TRUE) %>%
  mutate(HR = exp(estimate),  # Convert coefficients to hazard ratios
         CI_lower = exp(conf.low),
         CI_upper = exp(conf.high)) %>%
  dplyr::select(term, HR, CI_lower, CI_upper, p.value)

print(capped_fin_cox_psm_final_sum)

write.csv(capped_fin_cox_psm_final_sum, "coxdata/60/psm_covariates.csv")

# Group by Group Difference in Differences analysis
# 0 vs 1
# Create expansion1_time variable (Pre = 2006-2010, Post = 2011-2016)
matched_data_1 <- matched_data_1 %>%
  mutate(expansion1_time = case_when(
    dx.yr_grouped %in% c("2006-2010") ~ "Pre",
    dx.yr_grouped %in% c("2011-2013", "2014-2016", "2017-2019") ~ "Post",
    TRUE ~ NA_character_  # Assign NA for years outside the study period
  ))

# Create expansion1_group variable (1 = expansion1 states, 0 = non-expansion states)
matched_data_1 <- matched_data_1 %>%
  mutate(expansion1_group = case_when(
    Medicaid_Expansion_Status == 1 ~ 1,  # Expansion 1 (early)
    Medicaid_Expansion_Status == 0 ~ 0,  # Non-expansion (control)
    TRUE ~ NA_real_  # Assign NA for other groups
  ))

# Subset to only include subjects from Expansion Group 0 and 1
matched_data_1_did <- matched_data_1 %>%
  filter(Medicaid_Expansion_Status %in% c(0, 1)) %>%
  # Create a binary treatment variable:
  # did_treat = 1 if Expansion Group 1, 0 if Expansion Group 0
  mutate(did_treat = ifelse(Medicaid_Expansion_Status == 1, 1, 0),
         # Create a time variable: "Pre" for 2006-2010, "Post" for later years.
         did_time = ifelse(dx.yr_grouped == "2006-2010", "Pre", "Post"))

table(matched_data_1_did$dx.yr_grouped)


# Convert did_time and did_treat to factors with explicit ordering:
matched_data_1_did <- matched_data_1_did %>%
  mutate(
    did_time = factor(did_time, levels = c("Pre", "Post")),
    did_treat = factor(did_treat, levels = c(0, 1))
  )

matched_data_1_did$capped_time_60 <- pmin(matched_data_1_did$Survival.months, 60)
matched_data_1_did$capped_event_60 <- ifelse(matched_data_1_did$Survival.months > 60, 0, matched_data_1_did$Event)



capped_fin_cox_psm_final_did_1 <- coxph(Surv(capped_time_60, capped_event_60) ~ did_time * did_treat + cluster(Medicaid.Expansion.Status),
                                        data = matched_data_1_did, robust = TRUE)


summary(capped_fin_cox_psm_final_did_1)

capped_fin_cox_psm_final_did_1 <- tidy(capped_fin_cox_psm_final_did_1, conf.int = TRUE) %>%
  mutate(HR = exp(estimate),  # Convert coefficients to hazard ratios
         CI_lower = exp(conf.low),
         CI_upper = exp(conf.high)) %>%
  dplyr::select(term, HR, CI_lower, CI_upper, p.value)

print(capped_fin_cox_psm_final_did_1)

write.csv(capped_fin_cox_psm_final_did_1, "coxdata/60/psm_0v1.csv")



# 0 vs 2
# Create expansion1_time variable (Pre = 2006-2013, Post = 2014-2019)
matched_data_2 <- matched_data_2 %>%
  mutate(expansion2_time = case_when(
    dx.yr_grouped %in% c("2006-2010", "2011-2013") ~ "Pre",
    dx.yr_grouped %in% c("2014-2016", "2017-2019") ~ "Post",
    TRUE ~ NA_character_  # Assign NA for years outside the study period
  ))

# Create expansion1_group variable (1 = expansion2 states, 0 = non-expansion states)
matched_data_2 <- matched_data_2 %>%
  mutate(expansion2_group = case_when(
    Medicaid_Expansion_Status == 2 ~ 1,  # Expansion 2 (early)
    Medicaid_Expansion_Status == 0 ~ 0,  # Non-expansion (control)
    TRUE ~ NA_real_  # Assign NA for other groups
  ))



# Subset to only include subjects from Expansion Group 0 and 2
matched_data_2_did <- matched_data_2 %>%
  filter(Medicaid_Expansion_Status %in% c(0, 2)) %>%
  # Assign treatment variable: 1 for expansion (2), 0 for non-expansion (0)
  mutate(did_treat_2 = ifelse(Medicaid_Expansion_Status == 2, 1, 0),
         # Define Pre (2006-2013) and Post (2014+)
         did_time_2 = ifelse(dx.yr_grouped %in% c("2006-2010", "2011-2013"), "Pre", "Post"))

# Convert did_time and did_treat to factors with explicit ordering:
matched_data_2_did <- matched_data_2_did %>%
  mutate(
    did_time_2 = factor(did_time_2, levels = c("Pre", "Post")),
    did_treat_2 = factor(did_treat_2, levels = c(0, 1))
  )

# matched_data_2_did <- matched_data_2_did %>%
#   filter(!grepl("2017-2019", trimws(dx.yr_grouped))
#   )


matched_data_2_did <- matched_data_2_did %>%
  filter(!grepl("Unknown Race", trimws(Race.and.origin.recode..NHW..NHB..NHAIAN..NHAPI..Hispanic.)))

matched_data_2_did$capped_time_60 <- pmin(matched_data_2_did$Survival.months, 60)
matched_data_2_did$capped_event_60 <- ifelse(matched_data_2_did$Survival.months > 60, 0, matched_data_2_did$Event)



capped_fin_cox_psm_final_did_2 <- coxph(Surv(capped_time_60, capped_event_60) ~ did_time_2 * did_treat_2 + cluster(Medicaid.Expansion.Status),
                                        data = matched_data_2_did, robust = TRUE)

summary(capped_fin_cox_psm_final_did_2)

capped_fin_cox_psm_final_did_2 <- tidy(capped_fin_cox_psm_final_did_2, conf.int = TRUE) %>%
  mutate(HR = exp(estimate),  # Convert coefficients to hazard ratios
         CI_lower = exp(conf.low),
         CI_upper = exp(conf.high)) %>%
  dplyr::select(term, HR, CI_lower, CI_upper, p.value)

print(capped_fin_cox_psm_final_did_2)

write.csv(capped_fin_cox_psm_final_did_2, "coxdata/60/psm_0v2.csv")


# 0 vs 3
# Create expansion1_time variable (Pre = 2006-2013, Post = 2014-2019)
matched_data_3 <- matched_data_3 %>%
  mutate(expansion2_time = case_when(
    dx.yr_grouped %in% c("2006-2010", "2011-2013") ~ "Pre",
    dx.yr_grouped %in% c("2014-2016", "2017-2019") ~ "Post",
    TRUE ~ NA_character_  # Assign NA for years outside the study period
  ))

# Create expansion1_group variable (1 = expansion2 states, 0 = non-expansion states)
matched_data_3 <- matched_data_3 %>%
  mutate(expansion2_group = case_when(
    Medicaid_Expansion_Status == 3 ~ 1,  # Expansion 2 (early)
    Medicaid_Expansion_Status == 0 ~ 0,  # Non-expansion (control)
    TRUE ~ NA_real_  # Assign NA for other groups
  ))



# Subset to only include subjects from Expansion Group 0 and 2
matched_data_3_did <- matched_data_3 %>%
  filter(Medicaid_Expansion_Status %in% c(0, 3)) %>%
  # Assign treatment variable: 1 for expansion (2), 0 for non-expansion (0)
  mutate(did_treat_3 = ifelse(Medicaid_Expansion_Status == 3, 1, 0),
         # Define Pre (2006-2013) and Post (2014+)
         did_time_3 = ifelse(dx.yr_grouped %in% c("2006-2010", "2011-2013"), "Pre", "Post"))

# Convert did_time and did_treat to factors with explicit ordering:
matched_data_3_did <- matched_data_3_did %>%
  mutate(
    did_time_3 = factor(did_time_3, levels = c("Pre", "Post")),
    did_treat_3 = factor(did_treat_3, levels = c(0, 1))
  )

matched_data_3_did <- matched_data_3_did %>%
  filter(!grepl("Unknown Race", trimws(Race.and.origin.recode..NHW..NHB..NHAIAN..NHAPI..Hispanic.)))

matched_data_3_did$capped_time_60 <- pmin(matched_data_3_did$Survival.months, 60)
matched_data_3_did$capped_event_60 <- ifelse(matched_data_3_did$Survival.months > 60, 0, matched_data_3_did$Event)



capped_fin_cox_psm_final_did_3 <- coxph(Surv(capped_time_60, capped_event_60) ~ did_time_3 * did_treat_3 + cluster(Medicaid.Expansion.Status),
                                        data = matched_data_3_did, robust = TRUE)

summary(capped_fin_cox_psm_final_did_3)

capped_fin_cox_psm_final_did_3 <- tidy(capped_fin_cox_psm_final_did_3, conf.int = TRUE) %>%
  mutate(HR = exp(estimate),  # Convert coefficients to hazard ratios
         CI_lower = exp(conf.low),
         CI_upper = exp(conf.high)) %>%
  dplyr::select(term, HR, CI_lower, CI_upper, p.value)

print(capped_fin_cox_psm_final_did_3)

write.csv(capped_fin_cox_psm_final_did_3, "coxdata/60/psm_0v3.csv")



### EARLY EXPANSION IMPLEMENTATION/POST-IMPLEMENTATION 60 MONTH SURVIVAL

# EARLY - consider early months (post-implementation periods)
# Step 1: Subset data to focus on 2014-2019
matched_data_1_did_early_1 <- matched_data_1 %>%
  filter(dx.yr_grouped %in% c("2006-2010", "2011-2013"))

# Step 2: Update the time variable for 2014-2016 and 2017-2019 periods
matched_data_1_did_early_1 <- matched_data_1_did_early_1 %>%
  mutate(did_time = ifelse(dx.yr_grouped %in% c("2011-2013"), "Post", "Pre"))

# Create treatment variable (same as before)
matched_data_1_did_early_1 <- matched_data_1_did_early_1 %>%
  mutate(did_treat = ifelse(Medicaid_Expansion_Status == 1, 1, 0))

# Step 3: Ensure 'did_time' is a factor with 'Pre' as reference
matched_data_1_did_early_1 <- matched_data_1_did_early_1 %>%
  mutate(
    did_time = factor(did_time, levels = c("Pre", "Post")),
    did_treat = factor(did_treat, levels = c(0, 1))
    
  )

# Check reference levels
contrasts(matched_data_1_did_early_1$did_time)
contrasts(matched_data_1_did_early_1$did_treat)

matched_data_1_did_early_1$capped_time_60 <- pmin(matched_data_1_did_early_1$Survival.months, 60)
matched_data_1_did_early_1$capped_event_60 <- ifelse(matched_data_1_did_early_1$Survival.months > 60, 0, matched_data_1_did_early_1$Event)

# Step 4: Run Cox regression model
capped_fin_cox_psm_final_early_1 <- coxph(Surv(capped_time_60, capped_event_60) ~ did_time * did_treat+ cluster(Medicaid.Expansion.Status),
                                          data = matched_data_1_did_early_1, robust = TRUE)

# Step 5: Summary of the model
summary(capped_fin_cox_psm_final_early_1)

capped_fin_cox_psm_final_did_early_1 <- tidy(capped_fin_cox_psm_final_early_1, conf.int = TRUE) %>%
  mutate(HR = exp(estimate),  # Convert coefficients to hazard ratios
         CI_lower = exp(conf.low),
         CI_upper = exp(conf.high)) %>%
  dplyr::select(term, HR, CI_lower, CI_upper, p.value)

print(capped_fin_cox_psm_final_did_early_1)

write.csv(capped_fin_cox_psm_final_did_early_1, "coxdata/60/psm_0v1_early.csv")


# LATER
# Step 1: Subset data to focus on 2014-2019
matched_data_1_did_later_1 <- matched_data_1 %>%
  filter(dx.yr_grouped %in% c("2006-2010", "2014-2016", "2017-2019"))

# Step 2: Update the time variable for 2014-2016 and 2017-2019 periods
matched_data_1_did_later_1 <- matched_data_1_did_later_1 %>%
  mutate(did_time = ifelse(dx.yr_grouped %in% c("2014-2016", "2017-2019"), "Post", "Pre"))

# Create treatment variable (same as before)
matched_data_1_did_later_1 <- matched_data_1_did_later_1 %>%
  mutate(did_treat = ifelse(Medicaid_Expansion_Status == 1, 1, 0))

# Step 3: Ensure 'did_time' is a factor with 'Pre' as reference
matched_data_1_did_later_1 <- matched_data_1_did_later_1 %>%
  mutate(
    did_time = factor(did_time, levels = c("Pre", "Post")),
    did_treat = factor(did_treat, levels = c(0, 1))
    
  )

# Check reference levels
contrasts(matched_data_1_did_later_1$did_time)
contrasts(matched_data_1_did_later_1$did_treat)

matched_data_1_did_later_1$capped_time_60 <- pmin(matched_data_1_did_later_1$Survival.months, 60)
matched_data_1_did_later_1$capped_event_60 <- ifelse(matched_data_1_did_later_1$Survival.months > 60, 0, matched_data_1_did_later_1$Event)

# Step 4: Run Cox regression model
capped_fin_cox_psm_final_later_1 <- coxph(Surv(capped_time_60, capped_event_60) ~ did_time * did_treat + cluster(Medicaid.Expansion.Status),
                                          data = matched_data_1_did_later_1, robust = TRUE)

# Step 5: Summary of the model
summary(capped_fin_cox_psm_final_later_1)

capped_fin_cox_psm_final_did_later_1 <- tidy(capped_fin_cox_psm_final_later_1, conf.int = TRUE) %>%
  mutate(HR = exp(estimate),  # Convert coefficients to hazard ratios
         CI_lower = exp(conf.low),
         CI_upper = exp(conf.high)) %>%
  dplyr::select(term, HR, CI_lower, CI_upper, p.value)

print(capped_fin_cox_psm_final_did_later_1)
print(capped_fin_cox_psm_final_did_early_1)
write.csv(capped_fin_cox_psm_final_did_later_1, "coxdata/60/psm_0v1_later.csv")



### 60 month survival for implementation/postimplemntation periods for on-time expansion group ##

# EARLY - consider early months (post-implementation periods)
# Step1: Subset data to focus on 2017-2019
matched_data_2_did_early_2 <- matched_data_2 %>%
  filter(dx.yr_grouped %in% c("2006-2010", "2011-2013", "2014-2016"))

# Step 2: Update the time variable for 2014-2016 and 2017-2019 periods
matched_data_2_did_early_2 <- matched_data_2_did_early_2 %>%
  mutate(did_time_2 = ifelse(dx.yr_grouped %in% c("2014-2016"), "Post", "Pre"))

# Create treatment variable (same as before)
matched_data_2_did_early_2 <- matched_data_2_did_early_2 %>%
  mutate(did_treat_2 = ifelse(Medicaid_Expansion_Status == 2, 1, 0))

# Step 3: Ensure 'did_time' is a factor with 'Pre' as reference
matched_data_2_did_early_2 <- matched_data_2_did_early_2 %>%
  mutate(
    did_time_2 = factor(did_time_2, levels = c("Pre", "Post")),
    did_treat_2 = factor(did_treat_2, levels = c(0, 1))
  )

# Check reference levels
contrasts(matched_data_2_did_early_2$did_time_2)
contrasts(matched_data_2_did_early_2$did_treat_2)

matched_data_2_did_early_2$capped_time_60 <- pmin(matched_data_2_did_early_2$Survival.months, 60)
matched_data_2_did_early_2$capped_event_60 <- ifelse(matched_data_2_did_early_2$Survival.months > 60, 0, matched_data_2_did_early_2$Event)


# Step 4: Run Cox regression model
capped_fin_cox_psm_final_early_2 <- coxph(Surv(capped_time_60, capped_event_60) ~ did_time_2 * did_treat_2 + cluster(Medicaid.Expansion.Status),
                                          data = matched_data_2_did_early_2, robust = TRUE)

# Step 5: Summary of the model
summary(capped_fin_cox_psm_final_early_2)

capped_fin_cox_psm_final_did_early_2 <- tidy(capped_fin_cox_psm_final_early_2, conf.int = TRUE) %>%
  mutate(HR = exp(estimate),  # Convert coefficients to hazard ratios
         CI_lower = exp(conf.low),
         CI_upper = exp(conf.high)) %>%
  dplyr::select(term, HR, CI_lower, CI_upper, p.value)

print(capped_fin_cox_psm_final_did_early_2)

write.csv(capped_fin_cox_psm_final_did_early_2, "coxdata/60/psm_0v2_early.csv")




# LATER
# Step1: Subset data to focus on 2017-2019
matched_data_2_did_later_2 <- matched_data_2 %>%
  filter(dx.yr_grouped %in% c("2006-2010", "2011-2013", "2017-2019"))

# Step 2: Update the time variable for 2014-2016 and 2017-2019 periods
matched_data_2_did_later_2 <- matched_data_2_did_later_2 %>%
  mutate(did_time_2 = ifelse(dx.yr_grouped %in% c("2017-2019"), "Post", "Pre"))

# Create treatment variable (same as before)
matched_data_2_did_later_2 <- matched_data_2_did_later_2 %>%
  mutate(did_treat_2 = ifelse(Medicaid_Expansion_Status == 2, 1, 0))

# Step 3: Ensure 'did_time' is a factor with 'Pre' as reference
matched_data_2_did_later_2 <- matched_data_2_did_later_2 %>%
  mutate(
    did_time_2 = factor(did_time_2, levels = c("Pre", "Post")),
    did_treat_2 = factor(did_treat_2, levels = c(0, 1))
  )

# Check reference levels
contrasts(matched_data_2_did_later_2$did_time_2)
contrasts(matched_data_2_did_later_2$did_treat_2)

matched_data_2_did_later_2$capped_time_60 <- pmin(matched_data_2_did_later_2$Survival.months, 60)
matched_data_2_did_later_2$capped_event_60 <- ifelse(matched_data_2_did_later_2$Survival.months > 60, 0, matched_data_2_did_later_2$Event)



# Step 4: Run Cox regression model
capped_fin_cox_psm_final_later_2 <- coxph(Surv(capped_time_60, capped_event_60) ~ did_time_2 * did_treat_2 + cluster(Medicaid.Expansion.Status),
                                          data = matched_data_2_did_later_2, robust = TRUE)

# Step 5: Summary of the model
summary(capped_fin_cox_psm_final_later_2)

capped_fin_cox_psm_final_did_later_2 <- tidy(capped_fin_cox_psm_final_later_2, conf.int = TRUE) %>%
  mutate(HR = exp(estimate),  # Convert coefficients to hazard ratios
         CI_lower = exp(conf.low),
         CI_upper = exp(conf.high)) %>%
  dplyr::select(term, HR, CI_lower, CI_upper, p.value)

print(capped_fin_cox_psm_final_did_later_2)

write.csv(capped_fin_cox_psm_final_did_later_2, "coxdata/60/psm_0v2_later.csv")

