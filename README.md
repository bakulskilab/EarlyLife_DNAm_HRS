# Early life experiences are associated with later life DNA methylation signatures in the Health and Retirement Study

The repository contains the datasets and R files used in the analysis of early life exposures and later life DNA methylation using data from the Health and Retirement Study (HRS).

Original HRS datasets are available at: https://hrsdata.isr.umich.edu/data-products

DNA methylation data for the single-site analysis requires an application to the National Institute on Aging Genetics of Alzheimer's Disease Data Storage Site (NIAGADS) for use (Accession Number: [NG00153.v1](https://dss.niagads.org/datasets/ng00153/)).




## Citation
Kanney N, Cockell S, Wang H, Fu M, Dou J, Hicken MT, Payne-Sturges D, Needham B, Ware EB, Bakulski KM. Early life experiences are associated with later life DNA methylation signatures in the Health and Retirement Study. Journal. Year/Volume. doi/PMID




## Repository File Descriptions

### Descriptive Statistics

* __01_DataManage_Exclusion.qmd__ – imports and merges all files, data cleaning, management, and exclusion to get filtered, complete analytic sample

* __02_Bivariate.qmd__ – data visualizations on exposures, outcomes, and covariate variables

  

### DNA Methylation Age Acceleration and Global Methylation Files

* __ELE_DNAm_functions.R__ - functions to run code files 03, 04, Model_Diagnostics, and Forest_Plots

* __03_Primary.qmd__ – linear regression models for unadjusted, health behavior SES, and mutually adjusted models for DNAm age acceleration analysis

* __04_SensitivityModels.qmd__ – linear regression models for core and health only sensitivity models for DNAm age acceleration analysis

* __Model_Diagnostics.qmd__ – diagnostic plots from DNAm age acceleration model results

* __Forest_Plots.qmd__ – forest plots for DNAm age acceleration model results

  

### Epigenome-Wide Association Study Files

* __EWAS_function.R__ - functions for running code files 06 and EWAS_Tables_Plots

* __05_Primary_Sensitvity_EWAS.qmd__ – epigenome-wide association study (EWAS) for core, health only, and health behavior SES models 

* __EWAS_Tables_Plots__ - plots and tables for epigenome-wide association study (EWAS) for health behavior SES and sensitivity from __05_Primary_Sensitvity_EWAS.qmd__ DNAm analysis

* __EWAS_Compare.qmd__ - table output to compare health behavior SES EWAS model results across exposures and two prior EWAS studies on maternal educational attainment (PMID: 38052982) and cigarette smoking (PMID: 27651444)

### Gene Set Enrichment Analysis

* __06_GeneSet.qmd__ – methylGSA gene ontology pathway enrichment analysis


