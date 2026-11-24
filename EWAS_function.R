# Functions for EWAS plots

#-------------------------------------------------------------------------------

## DEFINE FUNCTIONS

#-------------------------------------------------------------------------------

# Lambda function - this is a summary measure of genomic inflation  
calculate_lambda <- function(p_values) {
  observed_chisq <- qchisq(1 - p_values, df = 1)
  lambda <- median(observed_chisq) / qchisq(0.5, df = 1)
  return(lambda)
}

#Function to create qq plots (original)
qq_plot <- function(p_values, lambda_value, title) {
  observed <- -log10(sort(p_values))
  expected <- -log10(ppoints(length(p_values)))
  ggplot(data.frame(expected, observed), aes(x = expected, y = observed)) +
    geom_point(size = 0.5, alpha = 0.6) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
    labs(title = title, subtitle = paste("Lambda: ", round(lambda_value, 2)),
         x = "Expected -log10(p)",
         y = "Observed -log10(p)") +
    theme_minimal() +
    theme(plot.title = element_text(size = 22, face = "bold"), plot.subtitle = element_text(size = 20, color = "black"), axis.text = element_text(size = 20, color = "black"), axis.title = element_text(size = 20))
}

#Function to create qq plots (years in school)
qq_plotschlyrs <- function(p_values, lambda_value, title) {
  observed <- -log10(sort(p_values))
  expected <- -log10(ppoints(length(p_values)))
  ggplot(data.frame(expected, observed), aes(x = expected, y = observed)) +
    geom_point(size = 0.5, alpha = 0.6) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
    labs(title = title, subtitle = paste("Lambda: ", round(lambda_value, 2)),
         x = "Expected -log10(p)",
         y = "Observed -log10(p)") +
    scale_y_continuous(limits = c(0, 60)) +
    theme_minimal() +
    theme(plot.title = element_text(size = 20, face = "bold"), plot.subtitle = element_text(size = 20, color = "black"), axis.text = element_text(size = 20, color = "black"), axis.title = element_text(size = 20))
}

#Function to create qq plots (other exposures)
qq_plot2 <- function(p_values, lambda_value, title) {
  observed <- -log10(sort(p_values))
  expected <- -log10(ppoints(length(p_values)))
  ggplot(data.frame(expected, observed), aes(x = expected, y = observed)) +
    geom_point(size = 0.5, alpha = 0.6) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
    labs(title = title, subtitle = paste("Lambda: ", round(lambda_value, 2)),
         x = "Expected -log10(p)",
         y = "Observed -log10(p)") +
    scale_y_continuous(limits = c(0, 8)) +
    theme_minimal() +
    theme(plot.title = element_text(size = 22, face = "bold"), plot.subtitle = element_text(size = 20, color = "black"), axis.text = element_text(size = 20, color = "black"), axis.title = element_text(size = 20))
}



#Create tbl_sumamry function  
top_tbl <- function(data, title) {
  select_data <- data %>%
    select(CpG, logFCpercent, P.Value, adj.P.Val, AveExpr, CHR, BP, UCSC_RefGene_Name, Relation_to_Island)
  tbl <- gt(select_data) %>%
    tab_header(title = title) %>%
    cols_label(CpG = "CpG",
               logFCpercent = "Estimate",
               P.Value = "P-Value",
               adj.P.Val = "FDR P",
               AveExpr = "Average DNA Methylation",
               CHR = "Chromosome",
               BP = "Position",
               UCSC_RefGene_Name = "UCSC Reference Gene",
               Relation_to_Island = "CpG Island Region"
    ) %>%
    fmt_scientific(columns = c(P.Value, adj.P.Val), decimals = 2) %>%
    tab_options(table.font.size = "small") %>%
    opt_align_table_header(align = "left")
  return(tbl)}


# wrapper function to input results into methylRAA
do_gsea <- function(df){
  input <- df$P.Value
  names(input) <- df$CpG
  res <- methylRRA(cpg.pval = input, array.type = "EPIC", method = "GSEA")
  return(res[,c(1,2,3,4,5,6,7)])
}

  