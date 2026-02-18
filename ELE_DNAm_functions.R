# Functions for regression analyses

#-------------------------------------------------------------------------------

## DEFINE FUNCTIONS 

#-------------------------------------------------------------------------------


# Create combinations of all possible outcome and exposure relationships and save in new df
func_var_combinations <- function(exposure_vars, outcome_vars) {
  expand.grid(outcome = outcome_vars, exposure = exposure_vars)
}
#creates a 2 column df (n=30): column 1 corresponds to outcome; column 2 corresponds to exposure


# Run linear regression models
run_linear_models <- function(func_var_combinations, data, covars = NULL, family = gaussian(link = "identity")) {
  apply(func_var_combinations, 1, function(vars) {
    if (is.null(covars)) {
      formula <- as.formula(paste(vars["outcome"], vars["exposure"], sep = " ~ "))
    } else {
      covars_formula <- paste(covars, collapse = " + ")
      formula <- as.formula(paste(vars["outcome"], paste(vars["exposure"], covars_formula, sep = " + "), sep = " ~ "))
    }
    model <- glm(formula, data = data, family = family)
    list(formula = formula, model = model)
  })
}
#the 2nd line apply function iterates over each row of the 2 column df created above
#the 3rd line specifies that the model we are running is not using any covariates
#the 4th line creates a formula to iterate the func_var_combinations models: outcome ~ exposure
#the 5th line specifies what to do if we are running a model with covariates
# the 6th and 7th line specifies another formula to iterate over  the func_var_combinations models: outcome ~ exposure + the covariates as well
#the 8th line creates a df called model that fits a glm/the generalized linear model using the formula created above 
#last line returns a list storing the iterated glm models 



# Get tidy model summaries
# Add label to each model summary with the respective outcome-exposure it is testing
# Get model summary statistics, and include confidence intervals + FDR adjusted p-values
models_summarize <- function(func_var_combinations, models_list) {
  labeled_models_list <- models_list
  names(labeled_models_list) <- apply(func_var_combinations, 1, function(vars) {
    paste(vars["outcome"], vars["exposure"], sep = "_vs_")
  })
  
  tidy_summaries_df <- bind_rows(
    lapply(names(labeled_models_list), function(model_name) {
      extract_summary(labeled_models_list[[model_name]]$model, model_name)
    })
  )
  tidy_summaries_df
}
#apply creates a unique name for each model summary we get from df unadjusted_m1 using the two columns of the func_var_combinations df putting "vs" to separate each name
#lapply iterates over each item in the list unadjusted_m1 and stores the below statistics in a new list




# Organize model summaries and add FDR adjusted p-values
extract_summary <- function(model, model_name) {
  tidy_model <- tidy(model)
  conf_intervals <- confint(model)
  tidy_model <- cbind(tidy_model, conf_intervals)
  tidy_model$FDR_p.value <- p.adjust(tidy_model$p.value, method = "fdr")
  tidy_model$model_name <- model_name
  return(tidy_model)
}



# Perform diagnostics on model results
perform_diagnostics <- function(model, model_name) {
  # Extract the model diagnostics 
  model_data <- augment(model)
  #augment creates a df with the model predictions, residuals, and other diagnostic criteria from each model 
  
  
  # Check Linearity Assumption (Residuals vs Fitted)
  residuals_fitted_plot <- ggplot(model_data, aes(.fitted, .resid)) +
    geom_point() +
    geom_smooth(method = "loess", se = FALSE, color = "blue") +
    labs(title = paste("Residuals vs Fitted:", model_name), x = "Fitted values", y = "Residuals") +
    theme(plot.title = element_text(size = 8))  
  
  # Check Homoscedasticity Assumption (Scale-Location plot)
  scale_location_plot <- ggplot(model_data, aes(.fitted, abs(.std.resid))) +
    geom_point() +
    geom_smooth(method = "loess", se = FALSE, color = "blue") +
    labs(title = paste("Scale-Location Plot:", model_name), x = "Fitted values", y = "Standardized Residuals") +
    theme(plot.title = element_text(size = 8))
  
  # Check Normality of residuals Assumption (Q-Q plot & Shapiro-Wilk test)
  shapiro_test <- shapiro.test(model_data$.resid) # Performing the test on regular residuals
  shapiro_pvalue <- shapiro_test$p.value
  qq_plot <- ggplot(model_data, aes(sample = .std.resid)) + # Using standardized residuals for the Q-Q plot
    stat_qq() +
    stat_qq_line() +
    labs(title = paste("Normal Q-Q Plot:", model_name, sprintf("(Shapiro-Wilk p = %.2e)", shapiro_pvalue)),
         x = "Theoretical Quantiles", y = "Standardized Residuals") +
    theme(plot.title = element_text(size = 6))
  
  # Histogram of residuals (Check of Normality)
  hist_plot <- ggplot(model_data, aes(x = .resid)) +
    geom_histogram(binwidth = 0.5, fill = "blue", alpha = 0.7) +
    labs(title = paste("Histogram of Residuals:", model_name), x = "Residuals", y = "Frequency") +
    theme(plot.title = element_text(size = 8))
  
  # Combine the diagnostoc plots from above using patchwork
  combined_plot <- (residuals_fitted_plot + scale_location_plot) / (qq_plot + hist_plot)
  return(combined_plot)
}


# Format and filter forest plot data
# Filter out the intercept rows and include only the main exposure terms
filtered_forest_plot_data <- function(models_summaries_df, terms_to_exclude) {
  filtered_df <- models_summaries_df %>%
    filter(!term %in% terms_to_exclude & !grepl("Sample_PlateThyagarajan_Sample_", term)) %>%
    select(term, estimate, std.error, statistic, p.value, `2.5 %`, `97.5 %`, FDR_p.value, model_name)
  
  # Add outcome and exposure columns by splitting model names
  filtered_df <- filtered_df %>%
    mutate(outcome = sub("_vs_.*", "", model_name),
           exposure = sub(".*_vs_", "", model_name))
  
  filtered_df
}


# Create forest plots
forest_plotunadj <- function(filtered_df) {
  ggplot(filtered_df, aes(x = estimate, y = exposure, color = outcome)) +
    geom_point(position = position_dodge(width = 0.9), size = 3) +
    scale_color_jama(palette = c("default"), alpha = 1) +
    geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.4, color = "gray40") +
    geom_errorbarh(aes(xmin = `2.5 %`, xmax = `97.5 %`), height = 0.3, position = position_dodge(width = 0.9)) +
    geom_text(aes(label = sprintf("%0.2f [%0.2f, %0.2f]%s", 
                                  round(estimate, 2), 
                                  round(`2.5 %`, 2), 
                                  round(`97.5 %`, 2),
                                  ifelse(`FDR_p.value` < 0.05, "*", ""))), 
              vjust = -0.5, hjust = -0.04, size = 6,fontface = "bold", position = position_dodge(width = 0.9), show.legend = FALSE) +
    labs(x = "Unadjusted difference in DNAm age acceleration (years)", y = "Exposure Variables", color = "Outcome") +
    theme_classic() +
    guides(color = guide_legend(reverse = TRUE)) +
    theme(axis.text = element_text(size = 18, color = "black"), axis.title = element_text(size = 22), legend.position = "right", legend.text = element_text(size = 18), legend.title = element_text(size = 18))
}

forest_plotunadjGM <- function(filtered_df) {
  ggplot(filtered_df, aes(x = estimate, y = exposure, color = outcome)) +
    geom_point(position = position_dodge(width = 0.9), size = 3) +
    scale_color_jama(palette = c("default"), alpha = 1) +
    geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.4, color = "gray40") +
    geom_errorbarh(aes(xmin = `2.5 %`, xmax = `97.5 %`), height = 0.3, position = position_dodge(width = 0.9)) +
    geom_text(aes(label = sprintf("%0.2f [%0.2f, %0.2f]%s", 
                                  round(estimate, 2), 
                                  round(`2.5 %`, 2), 
                                  round(`97.5 %`, 2),
                                  ifelse(`FDR_p.value` < 0.05, "*", ""))), 
              vjust = -0.5, hjust = -0.04, size = 6,fontface = "bold", position = position_dodge(width = 0.9), show.legend = FALSE) +
    labs(x = "Unadjusted difference in methylation (percent)", y = "Exposure Measures", color = "Outcome") +
    theme_classic() +
    guides(color = guide_legend(reverse = TRUE)) +
    theme(axis.text = element_text(size = 18, color = "black"), axis.title = element_text(size = 22), legend.position = "right", legend.text = element_text(size = 18), legend.title = element_text(size = 18))
}


forest_plot <- function(filtered_df) {
  outcome_levels <- c("GrimAge", "MPOA", "Levine", "Horvath", "Hannum")
  blues <- c("GrimAge" = "#191970", "MPOA" = "#0047AB", "Levine" = "#1F51FF", "Horvath" = "#00688B", "Hannum" = "#1874CD")
  ggplot(filtered_df, aes(x = estimate, y = exposure, color = outcome)) +
    geom_point(position = position_dodge(width = 0.9), size = 3) +
    scale_color_manual(values = blues) +
    geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.4, color = "gray40") +
    geom_errorbarh(aes(xmin = `2.5 %`, xmax = `97.5 %`), height = 0.3, position = position_dodge(width = 0.9)) +
    geom_text(aes(label = sprintf("%0.2f [%0.2f, %0.2f]%s", 
                                  round(estimate, 2), 
                                  round(`2.5 %`, 2), 
                                  round(`97.5 %`, 2),
                                  ifelse(`FDR_p.value` < 0.05, "*", ""))), 
              vjust = -0.7, hjust = -0.04, size = 6,fontface = "bold", position = position_dodge(width = 0.9), show.legend = FALSE) +
    labs(x = "Adjusted difference in DNAm age acceleration (years)", y = "Exposure measures", color = "Outcome") +
    theme_classic() +
    guides(color = guide_legend(reverse = TRUE)) +
    theme(axis.text = element_text(size = 18, color = "black"), axis.title = element_text(size = 22), legend.position = "right", legend.text = element_text(size = 18), legend.title = element_text(size = 18))
}

forest_plotGM <- function(filtered_df) {
  outcome_levels <- c("Global methylation", "Island methylation", "Shelf methylation", "Shore methylation", "Open Sea methylation")
  reds <- c("Global methylation" = "#870411", "Island methylation" = "#b01020", "Shelf methylation" = "#e81e32", "Shore methylation" = "#f05663", "Open Sea methylation" = "#ff6975")
  ggplot(filtered_df, aes(x = estimate, y = exposure, color = outcome)) +
    geom_point(position = position_dodge(width = 0.9), size = 3) +
    scale_color_manual(values = reds) +
    geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.4, color = "gray40") +
    geom_errorbarh(aes(xmin = `2.5 %`, xmax = `97.5 %`), height = 0.3, position = position_dodge(width = 0.9)) +
    geom_text(aes(label = sprintf("%0.2f [%0.2f, %0.2f]%s", 
                                  round(estimate, 2), 
                                  round(`2.5 %`, 2), 
                                  round(`97.5 %`, 2),
                                  ifelse(`FDR_p.value` < 0.05, "*", ""))), 
              vjust = -0.7, hjust = -0.04, size = 6,fontface = "bold", position = position_dodge(width = 0.9), show.legend = FALSE) +
    labs(x = "Adjusted difference in methylation (percent)", y = "Exposure measures", color = "Outcome") +
    theme_classic() +
    guides(color = guide_legend(reverse = TRUE)) +
    theme(axis.text = element_text(size = 18, color = "black"), axis.title = element_text(size = 22), legend.position = "right", legend.text = element_text(size = 18), legend.title = element_text(size = 18))
}
