# Load libraries
library(tidyverse)
library(ggplot2)
library(tidyr)

# Import CKD dataset
ckd_data <- read.csv("kidney_disease.csv")

# Drop 'id' column
ckd_data <- ckd_data %>% select(-id)

# Clean classification column
ckd_data$classification <- ckd_data$classification %>% 
  str_trim() %>% 
  tolower()

# Verify unique values
unique(ckd_data$classification)

# Define numerical columns
numerical_columns <- c('age', 'bp', 'sg', 'al', 'su', 'bgr', 'bu', 'sc', 'sod', 'pot', 'hemo', 'pcv', 'wc', 'rc')

# Handle missing values in numerical columns and replace with mean value
ckd_data[numerical_columns] <- lapply(ckd_data[numerical_columns], function(x) {
  x[x == "" | x == "\t?"] <- NA
  x <- as.numeric(as.character(x))
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  return(x)
})

# Define categorical columns
categorical_columns <- c('rbc', 'pc', 'pcc', 'ba', 'htn', 'dm', 'cad', 'appet', 'pe', 'ane', 'classification')

# Clean and convert categorical columns
clean_categorical <- function(x) {
  x <- tolower(x) %>%                
    str_squish() %>%                 
    str_replace_all("[^a-z0-9]", "") 
  return(x)
}

# Calculate the mode values in categorical columns
get_mode <- function(v) {
  uniq_v <- unique(v)
  uniq_v[which.max(tabulate(match(v, uniq_v)))]
}

# Handle missing values in categorical columns and convert to factors
for (col in categorical_columns) {
  # Clean it
  ckd_data[[col]] <- clean_categorical(ckd_data[[col]])
  
  # Replace NA values with mode values
  if (all(is.na(ckd_data[[col]]) | ckd_data[[col]] == "")) {
    ckd_data[[col]] <- "unknown"
  } else {
    mode_value <- as.character(get_mode(ckd_data[[col]][!is.na(ckd_data[[col]]) & ckd_data[[col]] != ""]))
    ckd_data[[col]][is.na(ckd_data[[col]]) | ckd_data[[col]] == ""] <- mode_value
  }
  
  # Convert to factor
  ckd_data[[col]] <- factor(ckd_data[[col]])
}

# Print the summary statistics for categorical columns
summary_stats_cat <- summary(ckd_data %>% select(all_of(categorical_columns)))
print(summary_stats_cat)

# Print the summary statistics for numerical columns
summary_stats_cat <- summary(ckd_data %>% select(all_of(numerical_columns)))
print(summary_stats_cat)

# Create a correlation matrix and convert it to a data frame
cor_matrix <- cor(ckd_data %>% select_if(is.numeric), use = "complete.obs")
cor_matrix_df <- as.data.frame(as.table(cor_matrix))

# Rename the columns for easy reading
colnames(cor_matrix_df) <- c("Variable1", "Variable2", "Correlation")

# Filter out any self-correlations and repeat pairs
cor_filtered <- cor_matrix_df %>%
  filter(Variable1 != Variable2) %>% # Remove rows where Variable1 is the same as Variable2.Get rid of self correlations.
  # Create a new column called pair that combines Variable1 and Variable2 into a single string, This makes sure the order doesn’t matter
  # (e.g., “A_B” is the same as “B_A”). The sort function ensures that the smaller variable name comes first. This makes each pair unique.
  mutate(pair = pmap_chr(list(Variable1, Variable2), ~ paste(sort(c(...)), collapse = "_"))) %>%
  # This line keeps the first occurrence of each pair to make sure there are no repeat correlations in our list.
  distinct(pair, .keep_all = TRUE) %>%
  arrange(desc(abs(Correlation)))

# Print the top 10 correlations in order of strongest to weakest to console
print(head(cor_filtered, n = 10))

# Create heatmap of numerical variables with decimal values to represent percentages
ggplot(cor_matrix_df, aes(Variable1, Variable2, fill = Correlation)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "lightblue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-1, 1), space = "Lab",
                       name = "Correlation") +
  geom_text(aes(label = sprintf("%.2f", Correlation)), color = "black") +
  theme_minimal() +
  labs(title = "Heatmap of Numerical Correlations", x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5))

# Scatter plot of PCV and Hemoglobin - strongest correlation
ggplot(ckd_data, aes(x = pcv, y = hemo, color = classification)) +
  geom_point() +
  labs(title = "Scatter Plot of PCV vs Hemoglobin by Classification", x = "Packed Cell Volume (PCV)", y = "Hemoglobin (Hemo)") +
  theme_minimal() +
  scale_color_manual(values = c("blue", "orange"), name = "Classification")

# Scatter plot of RC vs PCV
ggplot(ckd_data, aes(x = rc, y = pcv, color = classification)) +
  geom_point() +
  labs(title = "Scatter Plot of RC vs PCV by Classification", x = "Red Blood Cell Count (RC)", y = "Packed Cell Volume (PCV)") +
  theme_minimal() +
  scale_color_manual(values = c("red", "yellow"), name = "Classification")

# Scatter plot of RC vs Hemo
ggplot(ckd_data, aes(x = rc, y = hemo, color = classification)) +
  geom_point() +
  labs(title = "Scatter Plot of RC vs Hemoglobin (Hemo) by Classification", x = "Red Blood Cell Count (RC)", y = "Hemoglobin (Hemo)") +
  theme_minimal() +
  scale_color_manual(values = c("purple", "green"), name = "Classification")

# Histogram of age distribution in the dataset
ggplot(ckd_data, aes(x = age)) +
  geom_histogram(binwidth = 5, fill = "blue", color = "black") +
  labs(title = "Age Distribution", x = "Age", y = "Frequency") +
  theme_minimal()

# Grouped bar chart for hypertension by classification
ggplot(ckd_data, aes(x = htn, fill = classification)) +
  geom_bar(position = "dodge") +
  geom_text(stat = 'count', aes(label = after_stat(count)), 
            position = position_dodge(0.9), 
            vjust = -0.5) +
  labs(title = "Hypertension (htn)", x = "Hypertension", y = "Count") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

# Grouped bar chart for diabetes by classification
ggplot(ckd_data, aes(x = dm, fill = classification)) +
  geom_bar(position = "dodge") +
  geom_text(stat = 'count', aes(label = after_stat(count)), 
            position = position_dodge(0.9), 
            vjust = -0.5) +
  labs(title = "Diabetes (dm)", x = "Diabetes Mellitus", y = "Count") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

# Grouped bar chart for coronary artery disease by classification
ggplot(ckd_data, aes(x = cad, fill = classification)) +
  geom_bar(position = "dodge") +
  geom_text(stat = 'count', aes(label = after_stat(count)), 
            position = position_dodge(0.9), 
            vjust = -0.5) +
  labs(title = "Coronary Artery Disease (cad)", x = "Coronary Artery Disease", y = "Count") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

# Grouped bar chart for red blood cell count by classification
ggplot(ckd_data, aes(x = rbc, fill = classification)) +
  geom_bar(position = "dodge") +
  geom_text(stat = 'count', aes(label = after_stat(count)), 
            position = position_dodge(0.9), 
            vjust = -0.5) +
  labs(title = "Red Blood Cell Count (rbc)", x = "RBC", y = "Count") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

# Grouped bar chart for pus cell by classification
ggplot(ckd_data, aes(x = pc, fill = classification)) +
  geom_bar(position = "dodge") +
  geom_text(stat = 'count', aes(label = after_stat(count)), 
            position = position_dodge(0.9), 
            vjust = -0.5) +
  labs(title = "Pus Cell Count (pc)", x = "Pus Cell", y = "Count") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

# Grouped bar chart for pus cell clumps by classification
ggplot(ckd_data, aes(x = pcc, fill = classification)) +
  geom_bar(position = "dodge") +
  geom_text(stat = 'count', aes(label = after_stat(count)), 
            position = position_dodge(0.9), 
            vjust = -0.5) +
  labs(title = "Pus Cell Clumps (pcc)", x = "Pus Cell Clumps", y = "Count") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

# Grouped bar chart for blood analysis by classification
ggplot(ckd_data, aes(x = ba, fill = classification)) +
  geom_bar(position = "dodge") +
  geom_text(stat = 'count', aes(label = after_stat(count)), 
            position = position_dodge(0.9), 
            vjust = -0.5) +
  labs(title = "Blood Analysis (ba)", x = "Blood Analysis", y = "Count") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

# Grouped bar chart for appetite by classification
ggplot(ckd_data, aes(x = appet, fill = classification)) +
  geom_bar(position = "dodge") +
  geom_text(stat = 'count', aes(label = after_stat(count)), 
            position = position_dodge(0.9), 
            vjust = -0.5) +
  labs(title = "Appetite Status (appet)", x = "Appetite", y = "Count") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

# Grouped bar chart for pedal edema by classification
ggplot(ckd_data, aes(x = pe, fill = classification)) +
  geom_bar(position = "dodge") +
  geom_text(stat = 'count', aes(label = after_stat(count)), 
            position = position_dodge(0.9), 
            vjust = -0.5) +
  labs(title = "Pedal Edema Status (pe)", x = "Pedal Edema", y = "Count") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

# Grouped bar chart for anemia by classification
ggplot(ckd_data, aes(x = ane, fill = classification)) +
  geom_bar(position = "dodge") +
  geom_text(stat = 'count', aes(label = after_stat(count)), 
            position = position_dodge(0.9), 
            vjust = -0.5) +
  labs(title = "Anemia Status (ane)", x = "Anemia", y = "Count") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5)
  )
