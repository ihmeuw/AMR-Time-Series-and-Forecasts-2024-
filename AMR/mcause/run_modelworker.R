Sys.setenv(MKL_VERBOSE=0)

.libPaths(c(.libPaths(), "FILEPATH"))
library(data.table)
library(boot)
library(arm)
library(tidyr)
library(ggplot2)
library(dplyr)
library(lme4)
library(caret)
library(lmerTest)
library(DescTools)

print("Loading arguments")
args <- commandArgs(trailingOnly = T)
int_cause <- as.character(args[1])
data_dir <- "FILEPATH"
diag_dir <- "FILEPATH"
release_id <- as.integer(args[4])
in_sample_dir <- "FILEPATH"
out_of_sample_dir <- "FILEPATH"
oosv <- as.integer(args[5])

source("FILEPATH")
source("FILEPATH")


get_formula <- function(covariates) {
  base_covariates <- c("sex_id", "age_group_id", "(1|level_1/level_2)")
  if (grepl("by_age", data_dir)) {base_covariates <- base_covariates[base_covariates != "age_group_id"]}
  exp_vars <- c(covariates, base_covariates)
  resp_vars <- c("successes", "failures")
  formula <- reformulate(exp_vars, parse(text = sprintf("cbind(%s)", toString(resp_vars)))[[1]])
  return(formula)
}

print("Loading in Input Data")
data <- fread("FILEPATH")
covariates <- get_covariates(int_cause)
data <- convert_factor_variables(data)
loc_df <- get_location_metadata(location_set_id=35, release_id=release_id)[, c("location_id", "ihme_loc_id")]
data <- merge(loc_df, data, by=c("location_id"))
data$iso3 <- substr(data$ihme_loc_id, 1, 3)
set.seed(52)

print("Get Formula")
save_pre_model_diagnostics(data, in_sample_dir, covariates)
formula <- get_formula(covariates)

if (oosv == 1){
  print("Out of Sample Validation")
  data <- cbind(ID = rownames(data), data)
  rownames(data) <- 1:nrow(data)

  folds <- floor(runif(nrow(data), min = 0, max = 5)) + 1
  folds <- data %>% split(folds)

  preds <- c()
  for (i in 1:5){
    test <- folds[[i]]
    train <- setDT(data[!data$ID %in% test$ID,])

    stopifnot(nrow(test) + nrow(train) == nrow(data))
    
    test = subset(test, select = -c(ID))
    train = subset(train, select = -c(ID))
    
    model <- run_model(train, formula, int_cause, data_dir, save=FALSE, fold=i)
    save_betas(model, train, out_of_sample_dir, int_cause, i)
    
    preds[[i]] <- make_predictions(test, model)
  }

  oos_preds <- bind_rows(preds, .id = "column_label")
  stopifnot(nrow(oos_preds) == nrow(data))
  oos_preds <- data.table(oos_preds)
  save_diagnostics(oos_preds, out_of_sample_dir, int_cause)
}

print("Saving In Sample")
model <- run_model(data, formula, int_cause, data_dir, save=TRUE)
save_betas(model, data, in_sample_dir, int_cause, fold="all_data")
in_sample_preds <- make_predictions(data, model)
save_diagnostics(in_sample_preds, in_sample_dir, int_cause)
