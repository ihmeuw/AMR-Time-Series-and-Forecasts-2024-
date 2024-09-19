rm(list = ls())
source('/FILEPATH')
pacman::p_load(data.table, openxlsx, dplyr, ggplot2, gridExtra, argparse)
options(stringsAsFactors = FALSE)
central_root <- 'FILEPATH'
setwd(central_root)

get_confusion_matrix <- function(dt) {
  # Get a confusion matrix from the predictions
  # See doc: https://en.wikipedia.org/wiki/Evaluation_of_binary_classifiers
  
  # Number of resistant as predicted from model
  dt[, pred_resistant := gpr_mean * sample_size]
  # True positives
  dt[, TP := pmin(resistant, pred_resistant)]
  # True negatives
  dt[, TN := pmin(susceptible, sample_size - pred_resistant)]
  # False positives
  dt[, FP := pmax(0, pred_resistant - resistant)]
  # False negatives
  dt[, FN := pmax(0, sample_size - pred_resistant - susceptible)]
  
  confusion_matrix <- dt[, lapply(.SD, sum), .SDcols = c("TP", "FP", "FN", "TN")]
  return(confusion_matrix)
}

get_bin_class_metrics <- function(cm) {
  cm <- copy(cm)
  cm[, sensitivity := TP / (TP + FN)]  # AKA TPR
  cm[, specificity := TN / (TN + FP)]  # AKA 1 - FPR
  cm[, accuracy := (TP + TN) / (TP + TN + FP + FN)]
  return(cm)
}

finalconfig <- read.csv("FILEPATH")

finalconfig <- finalconfig[!finalconfig$description %in% c(
  'escherichia_coli-sulfa', 'enterococcus_spp-fluoroquinolone'),]

combos <- unique(finalconfig$description)

oos_summary <- c()
oos_metrics <- c()
for (combo in combos){
  bug <- str_split(combo, '-')[[1]][1]
  drug <- str_split(combo, '-')[[1]][2]
  input <- read.csv(paste0('FILEPATH', bug, '/', drug, '/input_w_ho.csv'))
  
  # drop large India CRAB datapoint
  if (combo == 'acinetobacter_baumanii-carbapenem'){
    input <- input[input$sample_size < 100000, ]
  }
  
  allpreds <- c()
  
  for (holdout in unique(finalconfig$holdout)){
    holdouts <- input[input$master_fold_id == holdout,]
    gpr <- model_load(finalconfig$run_id[finalconfig$holdout == holdout & finalconfig$description == combo], 'gpr')
    pred_set <- merge(holdouts, gpr, by = c('location_id', 'year_id', 'sex_id'))
    stg1 <- model_load(finalconfig$run_id[finalconfig$holdout == holdout & finalconfig$description == combo], 'stage1')
    stg1$stage1 <- inv.logit(stg1$stage1)
    pred_set <- merge(pred_set, stg1[, c('location_id', 'year_id', 'sex_id', 'stage1')], by = c('location_id', 'year_id', 'sex_id'))
    pred_set$pathogen <- bug
    pred_set$abx_class <- drug
    allpreds <- bind_rows(allpreds, pred_set)
  }
  
  allpreds$gpr_var <- ((allpreds$gpr_mean-allpreds$gpr_upper)/1.96)^2
  allpreds$adj_var <- allpreds$gpr_var+allpreds$variance
  allpreds$gpr_lower_adj <- pmax(0, allpreds$gpr_mean - 1.96*sqrt(allpreds$adj_var))
  allpreds$gpr_upper_adj <- pmin(1, allpreds$gpr_mean + 1.96*sqrt(allpreds$adj_var))
  
  bugdrugsum <- allpreds %>% group_by(pathogen, abx_class) %>%
    summarize(
      bias = sum(gpr_mean-val)/length(val),
      bias_stg1 = sum(stage1-val)/length(val),
      rmse = sqrt(sum((gpr_mean-val)^2)/length(val)),
      rmse_stg1 = sqrt(sum((stage1-val)^2)/length(val)),
      coverage = length(val[val <= gpr_upper_adj & val >= gpr_lower_adj])/length(val),
      old_coverage = length(val[val <= gpr_upper & val >= gpr_lower])/length(val)
    )
  
  bugdrugsum$data_points <- nrow(input)
  bugdrugsum$n_locations <- length(unique(input$location_id))
  
  oos_summary <- bind_rows(oos_summary, bugdrugsum)
  
  allpreds <- as.data.table(allpreds)
  allpreds$resistant <- (allpreds$val*allpreds$sample_size)
  allpreds$susceptible <- allpreds$sample_size-allpreds$resistant
  cm <- get_confusion_matrix(allpreds)
  metrics <- get_bin_class_metrics(cm)
  metrics$pathogen <- bug
  metrics$abx_class <- drug
  oos_metrics <- bind_rows(oos_metrics, metrics)
}

oos_summary <- merge(oos_metrics[, c('pathogen', 'abx_class', 'sensitivity', 'specificity', 'accuracy')],
                     oos_summary, by = c('pathogen', 'abx_class'))

write.csv(oos_summary, 'FILEPATH', row.names = F)

ggplot(oos_summary) +
  geom_point(aes(x = coverage, y = bias))

ggplot(oos_summary) +
  geom_point(aes(x = n_locations, y = coverage))


testconf <- finalconfig[finalconfig$success == 1,]

stg1 <- model_load(195166, 'stage1')
stg1$stage1 <- inv.logit(stg1$stage1)


allpreds$pathogen <- 'escherichia_coli'
allpreds$abx_class <- 'beta_lactamase_inhibitor'



oosdir <- 'FILEPATH'
bugs <- list.files(oosdir)[!list.files(oosdir) %in% 'run_logs']

stg1nas <- c()
for (bug in bugs){
  drugs = list.files(paste0(oosdir, bug, '/'))
  for (drug in drugs){
    na_total = 0
    stg1s = list.files(paste0(oosdir, bug, '/', drug, '/'), recursive = TRUE, pattern = 'custom_stage1_df.csv', full.names = T)
    for (stg1 in stg1s){
       result = read.csv(stg1)
       ho_na = anyNA(result$cv_custom_stg1)
       na_total = na_total + ho_na
    }
    na_result <- data.frame(bug, drug, na_total)
    names(na_result) <- c('pathogen', 'abx_class', 'hos_w_na_stg1')
    stg1nas <- bind_rows(stg1nas, na_result)
  }
}




metrics <- get_bin_class_metrics(cm)

library(DescTools, lib.loc = "FILEPATH")

get_roc_auc <- function(dt) {
  # Calculate ROC curve and AUC score AKA C-Statistic
  # https://developers.google.com/machine-learning/crash-course/classification/thresholding
  # https://developers.google.com/machine-learning/crash-course/classification/roc-and-auc
  
  # For a bunch of thresholds from 0 - 1
  thresholds <- 0:1000 / 1000
  results <- as.list(rep(NA, length(thresholds)))
  counter <- 1
  for (i in thresholds) {
    # Initialize 4 cols (TP, TN, FP, FN) as all 0
    dt[, (c("TP", "TN", "FP", "FN")) := 0]
    # if pred > threshold, resistant = TP, susceptible = FP
    # If pred < threshold, resistant = FN, susceptible = TN
    dt[gpr_mean >= i, TP := resistant]
    dt[gpr_mean >= i, FP := susceptible]
    dt[gpr_mean < i, FN := resistant]
    dt[gpr_mean < i, TN := susceptible]
    result <- dt[, lapply(.SD, sum), .SDcols = c("TP", "FP", "FN", "TN")]
    result[, threshold := i]
    results[[counter]] <- result
    counter <- counter + 1
  }
  results <- rbindlist(results)
  # After done, calculate true positive rate/false positive rate for each threshold value
  results[, tpr := TP / (TP + FN)]
  results[, fpr := FP / (FP + TN)]
  # Plot the roc curve
  roc_curve <- ggplot(results, aes(fpr, tpr)) + geom_line() + geom_point() +
    ggtitle("Receiver operating characteristic (ROC) curve")
  # Calculate AUC using trapezoids
  auc_score <- AUC(results$fpr, results$tpr, method = "trapezoid")
  return(list(roc_curve = roc_curve, auc_score = auc_score))
}

get_roc_auc(allpreds)

ggplot(allpreds) +
  geom_point(aes(x = gpr_mean, y = val, size = sample_size)) +
  facet_wrap(~super_region_id)

adjdat <- model_load(finalconfig$run_id[finalconfig$holdout == holdout & finalconfig$description == combo], 'adj_data')



allpreds$old_covered <- ifelse(allpreds$val >= allpreds$gpr_lower & allpreds$val <= allpreds$gpr_upper, 1, 0)
allpreds$new_covered <- ifelse(allpreds$val >= allpreds$gpr_lower_adj & allpreds$val <= allpreds$gpr_upper_adj, 1, 0)


ggplot(oos_summary) +
  geom_histogram(aes(x = coverage))


fold5 <- list.files(path=paste0("FILEPATH", 195145, '/draws_temp_0'), full.names = TRUE)
combined_fold5<- rbindlist(lapply(fold5, fread), use.names=TRUE)
combined_fold5_long<- data.table(tidyr::gather(combined_fold5, draw, val, draw_0: draw_999))

vars <- combined_fold5_long %>% group_by(location_id, year_id, sex_id, age_group_id) %>%
  summarize(
    var = var(val)
  )

gpr <- model_load(195145, 'gpr')
gpr$gpr_var <- (gpr$gpr_mean-gpr$gpr_upper)^2/1.96
gpr$other_var <- ((gpr$gpr_mean-gpr$gpr_upper)/1.96)^2

chk <- merge(gpr, vars, by = c('location_id', 'year_id', 'sex_id', 'age_group_id'))

adj1 <- model_load(195145, 'adj_data')
adj2 <- model_load(195146, 'adj_data')

chk <- merge(adj1[!is.na(adj1$original_data), c('original_data', 'original_variance', 'sample_size',
                                                'year_id', 'location_id', 'nsv')],
             adj2[!is.na(adj2$original_data), c('original_data', 'original_variance', 'sample_size',
                                                'year_id', 'location_id', 'nsv')],
             by = c('original_data', 'original_variance', 'sample_size',
                    'year_id', 'location_id'), suffixes = c('_1', '_2'))

ggplot(chk) +
  geom_point(aes(x = nsv_1, y = nsv_2))
