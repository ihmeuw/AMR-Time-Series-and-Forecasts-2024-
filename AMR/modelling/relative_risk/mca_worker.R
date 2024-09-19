## Purpose: Runs MCA - Not a necessary part of the pipeline but useful to run from time to time as we add more data 

username <- Sys.getenv("USER")
amr_repo <- sprintf("/ihme/homes/%s/amr/", username)
library(tidyverse)
library(readr)
library(dplyr)
library("FactoMineR")
library("factoextra")
library(ggplot2)

print("Loading arguments")
args <- commandArgs(trailingOnly = T)
print(args)
print(args[1])
description <- as.character(args[1])

coeff_folder <- paste0("/mnt/team/amr/priv/intermediate_files/05_rr/coefficients/", description)
mca_folder <- paste0("/mnt/team/amr/priv/intermediate_files/05_rr/mca/", description)

## Helper Functions --------------------------------------------------------------------------------

dimensions_of_variation <- function(df){
  # This function finds, of the possible descriptive columns, which have variation (IE: different column values)
  var_dims <- c()
  for (col in names(df)){
    if ((length(unique(df[[col]])) > 1) & (col != "adjusted") &  (col != "pathogen") & (col != "abx_class")) {
      var_dims <- c(var_dims, col)
    }
  }

  return(var_dims)
}

get_col_type <- function(df_dims){ 
  # MCA requires that pertinent dimensions be classified as categorical, quantitative, etc. 
  col_types <- c()
  for (col in df_dims){
    if (col %in% c("year_id", "age_group", "loc_id", "infs", "high_income")){
      col_type <- "n"
      col_types <- c(col_types, col_type)
    }else{
      col_type <- "c"
      col_types <- c(col_types, col_type)
    }
  }

  return(col_types)
}

## Read in Coefficients -------------------------------------------------------------------------------- 

# read in all sources together 
files <- list.files(path = coeff_folder, pattern = "\\.csv$", full.names = TRUE)
data_list <- lapply(files, read_csv)
df <- bind_rows(data_list)

# subset to possible relevant columns - 
mca_prep_df <- df[c("pathogen", "abx_class", "year_id", "age_group", "loc_id", "infs", "high_income", "adjusted", "coeff_adjusted")]
var_dims <- dimensions_of_variation(mca_prep_df)
col_types <- get_col_type(var_dims)

if (length(var_dims) < 2){
  stop("There are not enough dimensions to perform MCA - Script Exiting - Try again after running Stage 1 w/ more splits")
}else{
  print(paste0("Performing MCA on the following dimensions: ", var_dims))
  mca_df <- mca_prep_df[c("pathogen", "abx_class", "adjusted", var_dims)]
  
  factor_dims <- Filter(function(x) x != "coeff_adjusted", var_dims)
  mca_df[factor_dims] <- lapply(mca_df[factor_dims], factor)
}


## MCA --------------------------------------------------------------------------------
# Perform MCA on each pathogen-drug combination by the given dimension
# Perform one on:
### unadjusted RR coefficient
### the adjusted by Level 2 UCoD coefficient ('adjusted')
### the adjusted by origin coefficient ("hosp_adjusted")
# additional info and source code: http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/116-mfa-multiple-factor-analysis-in-r-essentials/

#  "adjusted", "hosp_adjusted" - can run for other adjustments if that's important
for (adjustment in c("unadjusted")){
  print(paste0("Working on ", adjustment))
  mca <- mca_df %>%
    filter(adjusted == adjustment) %>%
    filter(!is.na(coeff_adjusted)) %>%
    filter(!is.infinite(coeff_adjusted))
  
  mca <- mca %>%
    mutate(combo = paste(pathogen, abx_class, sep = "-"))

  for (drug_bug_combo in unique(mca$combo)){
    mca_combo <- mca %>%
      filter(combo == drug_bug_combo)
    
    mca_combo_small <- mca_combo[var_dims]
    sample_size <- nrow(mca_combo_small)

    mca_result <- try({
      mca_combo_small.mfa <- MFA(mca_combo_small, 
               group = rep(1, length(var_dims)), # We do not have multiple columns describing the same category
               type = col_types,
               name.group = var_dims,
               graph = FALSE)

    #eig.val <- get_eigenvalue(mca_combo_small.mfa)
    plot_scree <- fviz_screeplot(mca_combo_small.mfa) + 
      ggtitle(paste0("Scree Plot - ", drug_bug_combo, " - ", adjustment, "RR - SS: ", sample_size))
    
    group <- get_mfa_var(mca_combo_small.mfa, "group")
    plot_mfa_var <- fviz_mfa_var(mca_combo_small.mfa, choice = "group") +
      ggtitle(paste0("MFA Vars - ", drug_bug_combo, " - ", adjustment, "RR - SS: ", sample_size))

    scree_name <- paste0("scree_", drug_bug_combo, "_", adjustment, ".pdf")
    var_name <- paste0("mfa_var_", drug_bug_combo, "_", adjustment, ".pdf")
    
    scree_full_path <- paste0(mca_folder, "/", scree_name)
    var_full_path <- paste0(mca_folder, "/", var_name)
    
    ggsave(scree_full_path, plot_scree, device = "pdf", width = 11, height = 8.5)
    ggsave(var_full_path, plot_mfa_var, device = "pdf", width = 11, height = 8.5)
    }, silent=FALSE
  )

  if(inherits(mca_result, "try-error")) {
    print(paste0("Dimensions were too low for ", drug_bug_combo, " with a sample size of ", sample_size))
  } else {
    print(paste0("Execution successful. Outputs have been saved for ", drug_bug_combo))
  }
  }
}
