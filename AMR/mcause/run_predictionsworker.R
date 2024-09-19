Sys.setenv(MKL_VERBOSE=0)

.libPaths(c(.libPaths(), "FILEPATH"))
library(data.table)
library(boot)
library(tidyr)
library(arm)
library(merTools)
library(utils)
library(MASS)
library(dplyr)
library(arrow)

source("FILEPATH")

args <- commandArgs(trailingOnly = T)
description <- as.character(args[1])
int_cause <- as.character(args[2])
dir <- "FILEPATH"
year_id <- as.integer(args[4])

cause_id <- as.integer(args[5])


draw_num <- 100



print(paste0("year id ", year_id))
print(paste0("cause id ", cause_id))
print(paste0("model description ", description))

get_prediction_interval <- function(model, template, by_age=FALSE) {
    if ((nrow(template) == 0) & (by_age)) {
        return(data.frame())
    } else {
        template <- convert_factor_variables(template)
        print(paste0(Sys.time(), " Getting predictions for ", draw_num, " simluations"))
        df_pred <- predictInterval(merMod=model, newdata=template, level=0.95, n.sims=draw_num, stat="mean",
                                   type="linear.prediction", returnSims=TRUE, include.resid.var=FALSE,
                                   which="fixed")
        draws <- attr(df_pred, "sim.results")

        lvl1_cause <- as.character(unique(template$level_1))
        lvl2_cause <- as.character(unique(template$level_2))
        stopifnot(length(lvl1_cause) == 1)
        stopifnot(length(lvl2_cause) == 1)
        cause_effects <- ranef(model)
        lvl1_effect <- cause_effects$level_1[lvl1_cause, ]
        lvl2_effect <- cause_effects[['level_2:level_1']][paste(lvl2_cause, lvl1_cause, sep = ":"), ]
        if (!any(is.na(c(lvl1_effect, lvl2_effect)))) {
            draws <- draws + lvl1_effect + lvl2_effect
            point_predict_vars <- NULL
        } else if (!is.na(lvl1_effect)) {
            print("No conditional level 2 effect available, using only level 1")
            draws <- draws + lvl1_effect
            point_predict_vars <- NA
        } else if (!is.na(lvl2_effect)) {
            stop("You have a null level 1 effect but a not null level 2, how is this possible?")
        } else {
            print("No conditional random effects available, using only fixed effects")
            point_predict_vars <- NA
        }

        draws <- inv.logit(draws)
        template <- cbind(template, draws %>% as.data.table %>% setnames(paste0("draw_", 0:(draw_num-1))))

        fixed_effects_pred <- predict(model, newdata=template, re.form=point_predict_vars, type="response")

        template$point_estimate <- fixed_effects_pred

        return(template)
    }
}

get_prediction_ci <- function(model) {
    draws <- MASS::mvrnorm(n=draw_num, mu=coef(model), Sigma=vcov(model))
    targets <- read_parquet("FILEPATH")
    template <- data.frame()
    for (cause_id in targets$keep_causes) {
        df <- read_parquet("FILEPATH")
        template <- rbind(template, df)
    }
    template <- convert_factor_variables(template)
    df_pred <- model.matrix(formula(model)[-2], data = template)
    y_hat <- (df_pred %>% as.matrix) %*% (draws %>% as.matrix %>% t)
    y_hat <- inv.logit(y_hat)
    cbind(template, y_hat %>% as.data.table %>% setnames(paste0("draw_", 0:(draw_num-1))))
}

set_sepsis_fraction_to_1 <- function(predictions, dir, year_id, cause_id) {
    if (is.element(cause_id, c(368, 383))) {
        write_parquet(predictions, "FILEPATH")
        predictions[, paste0("draw_", 0:(draw_num-1))] = 1
    }
    return(predictions)
}

if (grepl("simple_glm", description)) {
    model <- readRDS("FILEPATH")
    predictions <- get_prediction_ci(model)
    print(paste0(Sys.time(), " Saving output"))
    for (cause in unique(predictions$cause_id)) {
        print(paste0(Sys.time(), " Saving output for cause_id ", cause))
        write_parquet(predictions[predictions$cause_id == cause], "FILEPATH")
    }
} else {
    if (grepl("by_age", description)) {
        sub_dirs <- "FILEPATHS"
        age_group_ids <- "FILEPATHS"
        predictions <- data.frame()
        for (age in age_group_ids) {
            print(paste0("Working on age_group_id: ", age))
            age_dir <- "FILEPATH"
            model <- readRDS("FILEPATH")
            template <- read_parquet("FILEPATH")
            template <- as.data.table(template)
            if (dim(template)[1] != 0) {
                age_predictions <- get_prediction_interval(model, template, by_age=TRUE)
                if (dim(predictions)[1] == 0){
                    predictions <- age_predictions
                }else{
                    predictions <- rbind(predictions, age_predictions, fill=TRUE)
                }
            }
        }
    } else {
        model <- readRDS("FILEPATH")
        template <- read_parquet("FILEPATH")
        predictions <- get_prediction_interval(model, template)
    }

    if ("detailed_age_group_id" %in% colnames(predictions)) { 
        predictions$detailed_age_group_id <- as.factor(predictions$detailed_age_group_id)
        predictions[is.na(detailed_age_group_id), detailed_age_group_id := age_group_id]
        
        predictions <- predictions %>%
            dplyr::rename(
                agg_age_group_id = age_group_id,
                age_group_id = detailed_age_group_id
                )
        predictions <- subset(predictions, select = -c(agg_age_group_id) )
    }

    print(paste0(Sys.time(), " Saving output"))
    if ((int_cause == "sepsis") & grepl("mortality", dir)) {predictions <- set_sepsis_fraction_to_1(predictions, dir, year_id, cause_id)}
    write_parquet(predictions, "FILEPATH")
}
