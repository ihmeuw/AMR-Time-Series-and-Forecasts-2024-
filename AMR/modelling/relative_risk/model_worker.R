
library(reticulate)
reticulate::use_python("FILEPATH")
mr <- import("mrtool")
library(dplyr)
library(readr)
library(data.table)
library(tidyverse)
library(assertthat)
library(metafor)
library(lme4)

source("FILEPATH")
source("FILEPATH")
source("FILEPATH")

source("FILEPATH")

amr_repo <- "FILEPATH"

print("Loading arguments")
args <- commandArgs(trailingOnly = T)
input_description <- as.character(args[1])
output_description <- as.character(args[2])
run_isv <- as.integer(args[3])
run_oosv <- as.integer(args[4])
syn_type <- as.character(args[5])

inpath <- "FILEPATH"
outpath <- "FILEPATH"
data_path <- "FILEPATH"

abx_meta <- read.csv("FILEPATH", stringsAsFactors = F) 

gamma_SD <- c(0.01, 0.000001)
intercept_sd <- 0.01
intercept_sd_strong <- 0.000001
syndrome_sd <- 1
n_samples <- 100L


get_pval <- function(beta, beta_sd, one_sided = FALSE) {
  zscore <- abs(beta/beta_sd)
  if (one_sided) 1 - pnorm(zscore)
  if (!one_sided) (1 - pnorm(zscore))*2
}
rmse <- function(m, o){
  sqrt(mean(abs(m - o)^2, na.rm = T))
}
mae <- function(m, o) {
  mean(abs(m - o), na.rm = T)
}

solvable_se <- function(data, min_se=1e-3){
  for (row in 1:nrow(data)){
    if (data$se[row] < min_se){
      bad_source <- data$source[row]
      bad_combo <- data$combo[row]
      print(paste0("Going to unadjusted coeff for: ", bad_source, " & ", bad_combo))
      data$coeff[row] <- data$coeff_unadjusted[row]
      data$se[row] <- data$se_unadjusted[row]
    }
  }
  return(data)
}

generate_crude_rr <- function(df, filter_cols){
  rr_df <- df %>% 
    dplyr::group_by_at(vars(filter_cols)) %>%
    summarize(deaths = sum(deaths), cases = sum(cases))
  
  rr_df_r <- rr_df %>%
    filter(resistant == 1) %>%
    rename(deaths_r = deaths) %>%
    rename(cases_r = cases)
  
  rr_df_s <- rr_df %>%
    filter(resistant == 0) %>%
    rename(deaths_s = deaths) %>%
    rename(cases_s = cases)
  
  combine_cols <- setdiff(filter_cols, "resistant")
  
  rr_df <- rr_df_r %>% left_join(rr_df_s, by=combine_cols)
  rr_df <- rr_df %>% 
    mutate(crude_rr = (deaths_r/cases_r)/(deaths_s/cases_s)) %>%
    mutate(log_crude_rr = log(crude_rr)) %>%
    mutate(lnrr_var = (1/deaths_r) - (1/cases_r) + (1/deaths_s) - (1/cases_s)) %>%
    mutate(lnrr_se = sqrt(lnrr_var)) %>%
    filter(!is.na(log_crude_rr)) %>%
    filter(!is.infinite(log_crude_rr)) %>%
    filter(lnrr_se != 0)
  
  return(rr_df)
}

generate_prior <- function(df, stage){
  fix_model <- c()
  mix_model <- c()
  for (j in unique(df$abx_class)){
    print(paste0("Working on ", j))
    
    dataslice <- df %>%
      filter(abx_class == j)
    
    if (j == "carbapenem" & stage == "stage_1"){
      dataslice <- dataslice %>%
        filter(source != "SOURCES")
    }
    
    m2 <- glmer(cbind(deaths, cases-deaths) ~ resistant + (1|source), data = dataslice, family=binomial(link="log"))
    
    OR_m <- exp(fixef(m2))[[2]]
    suscep <- dataslice %>% 
      filter(resistant == 0) %>% 
      dplyr::group_by(abx_class) %>%
      summarize(deaths = sum(deaths), cases=sum(cases))
    p0 <- suscep$deaths/suscep$cases
    RR_m <- OR_m / ((1 - p0) + (p0 * OR_m))
    
    
    mix_coeff <- c(j, p0, OR_m, RR_m, log(RR_m))
    mix_model <- rbind(mix_model, mix_coeff)
  }
  
  colnames(mix_model) <- c('abx_class', 'cfr_susceptible', 'OR', 'RR', 'log_RR')
  
  write.csv(mix_model, "FILEPATH", row.names=FALSE)
}


antibiotic_group <- c('third_gen_ceph','carbapenem','fluoroquinolone','penicillin','aminopenicillin','beta_lactamase_inhibitor','anti_pseudomonal_penicillin','methicillin','vancomycin','fourth_gen_ceph','sulfa','aminoglycoside','macrolide')

ir2 <- read.csv("FILEPATH", stringsAsFactors = F) %>% filter(include==1) %>% dplyr::select(abx_class,pathogen) %>% distinct()
ir2$pathogen[ir2$pathogen == "acinetobacter_baumanii"] <- 'acinetobacter_baumannii'
ir <- read.csv("FILEPATH", stringsAsFactors = F) %>% dplyr::select(abx_class, atc) %>% distinct()
ir$atc[ir$abx_class == 'vancomycin'] <- 'J01D'
ir<-ir[!ir$abx_class %in% c('polymixin','chloramphenicol','nitrofurantoin','second_gen_ceph','rifamycin'),]
irall<-c('abx_class' = 'all_resistant', 'atc' = 'ddd_per_1000')
ir <- rbind(ir,irall)

bdc <- read.csv("FILEPATH")
bdc <- bdc %>%
  filter(kevin_catrin_exclude == 0) %>%
  mutate(pathogen = case_when(
    pathogen == "group_b_strep" ~ 'streptococcus_group_b',
    pathogen == "enterococcus_spp" ~ 'enterococcus_others',
    pathogen == "group_a_strep" ~ 'streptococcus_group_a',
    pathogen == "non_typhoidal_salmonellae" ~ 'salmonella_ints_non_typhi/paratyphi',
    TRUE ~ pathogen
  )) %>%
  select(pathogen, abx_class) %>%
  distinct()

bdc <- bdc %>%
  mutate(combo = paste(pathogen, abx_class, sep = "-"))


files <- list.files(path = inpath, pattern = "\\.csv$", full.names = TRUE)
data_list <- lapply(files, read_csv)
data_list <- lapply(data_list, function(df) {
  df %>%
    mutate(sample_id = as.character(sample_id))
})
df <- bind_rows(data_list) %>%
  mutate(combo = paste(pathogen, abx_class, sep = "-")) %>%
  filter(combo %in% unique(bdc$combo))

if (syn_type == "bsi"){
  df <- df %>%
    mutate(infectious_syndrome = if_else(infectious_syndrome %in% c("L2_lower_respiratory_infection", "L2_urinary_tract_infection"), 
                                         "other_syndromes", 
                                         infectious_syndrome))
}

model_dat_1 <- generate_crude_rr(df, c("abx_class", "source", "pathogen", "resistant"))
model_dat_2 <- generate_crude_rr(df, c("abx_class", "source", "pathogen", "infectious_syndrome", "resistant"))

if (syn_type != "none"){
  cf <- model_dat_2
}else{
  cf <- model_dat_1
}
cf$coeff <- cf$log_crude_rr
cf$se <- cf$lnrr_se

cf <- cf %>%
  mutate(combo = paste(pathogen, abx_class, sep = "-"))

cf <- cf %>%
    filter(!is.na(coeff)) %>%
    filter(!is.infinite(coeff)) %>%
    filter(se != 0)

cf <- cf %>% 
  filter(combo %in% unique(bdc$combo))

og_data <- cf
og_data$key <- apply(og_data, 1, paste, collapse = "-")

cf <- cf %>%
  filter(!(source %in% c("SOURCES") & abx_class == "carbapenem" & pathogen == "acinetobacter_baumannii" & crude_rr >= 2)) %>%
  filter(!(source == "SOURCES" & abx_class == "carbapenem" & pathogen == "citrobacter_spp" & crude_rr >= 2)) %>%
  filter(!(source == "SOURCES" & abx_class == "third_gen_ceph" & pathogen == "citrobacter_spp" & crude_rr >= 2)) %>%
  filter(!(source %in% c("SOURCES") & abx_class == "carbapenem" & pathogen == "escherichia_coli" & crude_rr >= 2)) %>%
  filter(!(source %in% c("SOURCES") & abx_class == "fourth_gen_ceph" & pathogen == "acinetobacter_baumannii" & crude_rr >= 2)) %>%
  filter(!(source %in% c("SOURCES") & abx_class == "fourth_gen_ceph" & pathogen == "citrobacter_spp" & crude_rr >= 2)) %>%
  filter(!(source %in% c("SOURCES") & abx_class == "fourth_gen_ceph" & pathogen == "serratia_spp"& crude_rr >= 2)) %>%
  filter(!(source %in% c("SOURCES") & abx_class == "third_gen_ceph" & pathogen == "serratia_spp"& crude_rr >= 2)) %>%
  filter(!(source == "SOURCES" & abx_class %in% c("aminoglycoside", "fourth_gen_ceph") & pathogen == "pseudomonas_aeruginosa")) %>%
  filter(!(source %in% c("SOURCES") & abx_class == "carbapenem" & pathogen == "pseudomonas_aeruginosa" & crude_rr >= 2)) %>%
  filter(!(source %in% c("SOURCES") & abx_class == "third_gen_ceph" & pathogen == "pseudomonas_aeruginosa"& crude_rr >= 2)) %>%
  filter(!(source == "SOURCES" & abx_class == "third_gen_ceph" & pathogen == "pseudomonas_aeruginosa"& crude_rr >= 2)) %>%
  filter(!(source == "SOURCES" & abx_class == "vancomycin" & pathogen %in% c("staphylococcus_aureus", "enterococcus_faecium") & crude_rr >= 2)) %>%
  filter(!(source %in% c("SOURCES") & abx_class == "fluoroquinolone" & pathogen == "citrobacter_spp" & crude_rr >= 2)) %>%
  filter(!(source == "SOURCES" & abx_class == "fluoroquinolone" & pathogen %in% c("enterococcus_faecalis", "escherichia_coli") & crude_rr >= 2)) %>%
  filter(!(source %in% c("SOURCES") & abx_class == "fluoroquinolone" & pathogen == "morganella_spp" & crude_rr <= 0.5)) %>%
  filter(!(source == "SOURCES" & abx_class == "sulfa" & pathogen == "streptococcus_pneumoniae" & crude_rr >= 2)) %>%
  filter(!(source %in% c("SOURCES") & abx_class == "fluoroquinolone" & pathogen == "enterobacter_spp" & crude_rr >= 2)) %>%
  filter(!(source %in% c("SOURCES") & abx_class == "fluoroquinolone" & pathogen == "proteus_spp" & crude_rr >= 2)) %>%
  filter(!(source %in% c("SOURCES") & abx_class == "aminoglycoside" & pathogen == "acinetobacter_baumannii" & crude_rr >= 2)) %>%
  filter(!(source %in% c("SOURCES") & abx_class == "aminoglycoside" & pathogen == "pseudomonas_aeruginosa" & crude_rr >= 2)) %>%
  filter(!(source %in% c("SOURCES") & abx_class == "fluoroquinolone" & pathogen == "enterococcus_faecium" & crude_rr >= 2)) %>%
  filter(!(source %in% c("SOURCES") & abx_class == "fluoroquinolone" & pathogen == "enterococcus_faecalis" & crude_rr >= 2)) %>%
  filter(!(source %in% c("SOURCES") & abx_class == "fluoroquinolone" & pathogen == "enterococcus_others" & crude_rr >= 2))
  
cf$key <- apply(cf, 1, paste, collapse = "-")
  
outlier <- og_data %>%
  rowwise() %>%
  mutate(outlier = if_else(key %in% cf$key, 0, 1)) %>%
  ungroup() %>%
  select(-key)

write.csv(outlier, "FILEPATH", row.names = FALSE)

if (syn_type == "none"){
  merge_cols = c("pathogen", "abx_class", "source")
  
  stage_1_outlier <- og_data %>%
    filter(!(abx_class %in% c("macrolide", "fluoroquinolone") & crude_rr >= 2)) %>%
    filter(!(abx_class %in% c("macrolide", "fluoroquinolone") & crude_rr > 1.2 & pathogen == "staphylococcus_aureus")) %>%
    filter(!(abx_class == "methicillin" & crude_rr <1)) %>%
    mutate(key = paste(pathogen, abx_class, source, sep = "-"))
  
  stage_1 <- df %>% 
    dplyr::group_by(abx_class, source, pathogen, resistant) %>%
    summarize(deaths = sum(deaths), cases=sum(cases)) %>%
    mutate(combo = paste(pathogen, abx_class, sep = "-")) %>%
    mutate(key = paste(pathogen, abx_class, source, sep = "-")) %>%
    filter(key %in% unique(stage_1_outlier$key))
}else{
  merge_cols = c("pathogen", "abx_class", "source", "infectious_syndrome")
  
  stage_1_outlier <- og_data %>%
    mutate(key = paste(pathogen, abx_class, source, infectious_syndrome, sep = "-"))
  
  stage_1 <- df %>% 
    dplyr::group_by(abx_class, source, pathogen, infectious_syndrome, resistant) %>%
    summarize(deaths = sum(deaths), cases=sum(cases)) %>%
    mutate(combo = paste(pathogen, abx_class, sep = "-")) %>%
    mutate(key = paste(pathogen, abx_class, source, infectious_syndrome, sep = "-")) %>%
    filter(key %in% unique(stage_1_outlier$key))
}

stage_1_outlier <- og_data %>%
  filter(!(abx_class %in% c("macrolide", "fluoroquinolone") & crude_rr >= 2)) %>%
  filter(!(abx_class %in% c("macrolide", "fluoroquinolone") & crude_rr > 1.2 & pathogen == "staphylococcus_aureus")) %>%
  filter(!(abx_class == "methicillin" & crude_rr <1)) %>%
  mutate(key = paste(pathogen, abx_class, source, sep = "-"))

stage_1 <- df %>% 
  dplyr::group_by(abx_class, source, pathogen, resistant) %>%
  summarize(deaths = sum(deaths), cases=sum(cases)) %>%
  mutate(combo = paste(pathogen, abx_class, sep = "-")) %>%
  mutate(key = paste(pathogen, abx_class, source, sep = "-")) %>%
  filter(key %in% unique(stage_1_outlier$key))

generate_prior(stage_1, "stage_1")



make_predict_template <- function(bdc, j, syn_type, ref_syn=NULL, syndromes=NULL, template_data=NULL){
    print("Making template")
  if (is.null(template_data)){
    if (syn_type != "none"){
      unique_pathogens <- bdc %>%
        filter(abx_class == j) %>%
        distinct(pathogen)
      
      abx_template <- expand.grid(abx_class = j,
                                  pathogen = unique_pathogens$pathogen,
                                  infectious_syndrome = syndromes)
      for (i in syndromes) {
        if (i != ref_syn){
          abx_template$new <- ifelse(abx_template$infectious_syndrome == i,1,0)
          colnames(abx_template)[colnames(abx_template) == 'new'] <- paste0(i)
        }
      }
      
      abx_template <- abx_template[order(abx_template$pathogen), ]
    }else{
      print("Making template")
      unique_pathogens <- bdc %>%
        filter(abx_class == j) %>%
        distinct(pathogen)
      
      abx_template <- expand.grid(abx_class = j,
                                  pathogen = unique_pathogens$pathogen)
    }
  }else{
    if (syn_type != "none"){
      abx_template <- template_data[,c("pathogen", "abx_class", "infectious_syndrome", "coeff", "se")] %>%
        rename("obs" = "coeff") %>%
        rename("obs_se" = "se")
      for (i in syndromes){
        if (i != ref_syn){
          abx_template$new <- ifelse(abx_template$infectious_syndrome == i,1,0)
          colnames(abx_template)[colnames(abx_template) == 'new'] <- paste0(i)
        }
      }
    }else{
      abx_template <- template_data[,c("pathogen", "abx_class", "coeff", "se")] %>%
        rename("obs" = "coeff") %>%
        rename("obs_se" = "se")
    }
  }
    
  return(abx_template)
}

run_mr_brt <- function(abx_prior, syn_priors_list, intercept_sd, syndrome_sd, gamma_sd, abx_template, rr_dat, syndromes, ref_syn){
  if (syn_type != "none"){
    syn_model <- paste0("mr$LinearCovModel(\"",syn_priors_list[,1], "\", prior_beta_gaussian = array(c(0,", syndrome_sd, ")))")
    model_definition <- paste(syn_model, collapse = ",\n")
  }
  mr_pred_dat <- mr$MRData()
  if(syn_type == "all"){
    if("obs" %in% colnames(abx_template)){
      mr_pred_dat$load_df(
        data = abx_template,
        col_obs = "obs",
        col_obs_se = "obs_se",
        col_covs = c(syn_priors_list),
        col_study_id = "pathogen" )
    }else{
      mr_pred_dat$load_df(
        data = abx_template,
        col_covs = c(syn_priors_list),
        col_study_id = "pathogen" )
    }
  }else{
    if ("obs" %in% colnames(abx_template)){
      mr_pred_dat$load_df(
        data = abx_template,
        col_obs = "obs",
        col_obs_se = "obs_se",
        col_study_id = "pathogen" )
    }else{
      mr_pred_dat$load_df(
        data = abx_template,
        col_study_id = "pathogen" )
    }
  }
  
  mr_pred <- mr_pred_dat$to_df()

  if (syn_type != "none"){
    mrbrtschema <- paste0(
      "mr$MRBRT(
  data = rr_dat,
  cov_models = list(mr$LinearCovModel('intercept', use_re = TRUE, prior_beta_gaussian = array(c(", abx_prior[1],",", intercept_sd,")), prior_gamma_gaussian = array(c(0,", gamma_sd, "))),",
      model_definition,
      "))"
    ) 
  }else{
    mrbrtschema <- paste0(
      "mr$MRBRT(
  data = rr_dat,
  cov_models = list(mr$LinearCovModel('intercept', use_re = TRUE, prior_beta_gaussian = array(c(", abx_prior[1],",",intercept_sd,")), prior_gamma_gaussian = array(c(0,", gamma_sd,")))",
      "))"
    )
  }
  
  print("Predicting!")
  rr_mod <- eval(parse(text = paste(mrbrtschema)))
  rr_mod$fit_model(inner_print_level = 5L, inner_max_iter = 1000L)
  
  mr_pred$predict <- rr_mod$predict(data = mr_pred_dat, predict_for_study = TRUE)
  
  if (syn_type != "none"){
    fe_syns <- syndromes[syndromes != ref_syn]
    mr_pred <- mr_pred %>%
      mutate(ref_syn = if_else(rowSums(mr_pred[,fe_syns]) == 0, 1, 0))
    names(mr_pred)[names(mr_pred) == "ref_syn"] <- ref_syn
  }
  if (syn_type == "all"){
    if ("L2_blood_stream_infection" %in% colnames(mr_pred)){
      mr_pred <- mr_pred %>% 
        mutate(infectious_syndrome = case_when(
          L2_blood_stream_infection == 1 ~ "L2_blood_stream_infection"))
    }
    if ("L2_urinary_tract_infection" %in% colnames(mr_pred)){
      mr_pred <- mr_pred %>% 
        mutate(infectious_syndrome = case_when(
          L2_urinary_tract_infection == 1 ~ "L2_urinary_tract_infection", 
          TRUE ~ infectious_syndrome))
    }
    if ("L2_lower_respiratory_infection" %in% colnames(mr_pred)){
      mr_pred <- mr_pred %>% 
        mutate(infectious_syndrome = case_when(
          L2_lower_respiratory_infection == 1  ~ "L2_lower_respiratory_infection", 
          TRUE ~ infectious_syndrome))
    }
    if ("other_syndromes" %in% colnames(mr_pred)){
      mr_pred <- mr_pred %>% 
        mutate(infectious_syndrome = case_when(
          other_syndromes == 1  ~ "other_syndromes", 
          TRUE ~ infectious_syndrome))
    }
    mr_pred <- mr_pred %>% 
      rename(pathogen = study_id)
    }else{
      mr_pred <- mr_pred %>% 
        rename(pathogen = study_id)
    }
  
  mr_pred$abx_class <- j
  
  samples <- rr_mod$sample_soln(sample_size = n_samples)
  
  draws <- rr_mod$create_draws(
    data = mr_pred_dat,
    beta_samples = samples[[1]], 
    gamma_samples = samples[[2]],
    random_study = TRUE)
  
  if (syn_type != "none"){
    pred_w_draws <- cbind(mr_pred[, colnames(mr_pred) %in% c("abx_class", "pathogen", "infectious_syndrome", "predict")],
                          draws) %>%
      rename(point_estimate = predict)
  }else{
    pred_w_draws <- cbind(mr_pred[, colnames(mr_pred) %in% c("abx_class", "pathogen", "predict")],
                          draws) %>%
      rename(point_estimate = predict)
  }
  
  return(list(pred = mr_pred, mod = rr_mod, draws = pred_w_draws, samp = samples))
}

get_summary <- function(df, SD, samples=NA){
  if (syn_type != "none"){
    sum_tab <- df[, c("abx_class", "pathogen", "infectious_syndrome")]
  }else{
    sum_tab <- df[, c("abx_class", "pathogen")]
  }
  sum_tab <- sum_tab %>%
    mutate(pval = get_pval(df$predict,SD)) %>%
    mutate(pval05 = ifelse(pval<0.05,"*","")) %>%
    mutate(log_point_estimate = df$predict) %>%
    mutate(log_lower = log_point_estimate - (1.96*SD)) %>%
    mutate(log_upper = log_point_estimate + (1.96*SD)) %>%
    mutate(point_estimate = exp(log_point_estimate)) %>%
    mutate(lower = exp(log_lower)) %>%
    mutate(upper = exp(log_upper))
  
  if (length(samples) > 1){
    sum_tab <- sum_tab %>%
      mutate(sd = sd(samples[[1]]))
  }
  
  return(sum_tab)
}

filter_data <- function(cf, ir2, j){
    nonintrinsic <- unique(ir2$pathogen[ir2$abx_class==j])

    data <- cf[cf$abx_class == paste0(j),] %>%
      filter(pathogen %in% nonintrinsic)

  return(data)
}

filter_priors <- function(table_priors, j){
  if (j %in% c("beta_lactamase_inhibitor")){
    abx_prior <- c(table_priors[table_priors$abx_class == j,]$log_RR, intercept_sd*10, intercept_sd_strong) 
  }else if (j == "fluoroquinolone"){
    abx_prior <- c(table_priors[table_priors$abx_class == j,]$log_RR, intercept_sd/10, intercept_sd_strong)
    }else{
    abx_prior <- c(table_priors[table_priors$abx_class == j,]$log_RR, intercept_sd, intercept_sd_strong) 
  }

  return(abx_prior)
}

run_model <- function(j, data, ir2, table_priors, template_data=NULL){
  print(paste0("Working on ", j))
  abx_prior <- filter_priors(table_priors, j)
  
  syn_priors_list <- c()
  if (syn_type != "none"){
    data <- data %>%
      filter(!infectious_syndrome == "L2_unspecified")
    
    syndromes <- unique(data$infectious_syndrome)
    ref_syn <- "other_syndromes"
    
    for (i in syndromes) {
      if (i != ref_syn){
        data$new <- ifelse(data$infectious_syndrome == i,1,0)
        colnames(data)[colnames(data) == 'new'] <- paste0(i)
        syn_effs <- c(i)
        syn_priors_list <- rbind(syn_priors_list, syn_effs)
      }
    }
    
    abx_template <- make_predict_template(bdc, j, syn_type, ref_syn=ref_syn, syndromes=syndromes, template_data=template_data)
    
    data <- data[,c("pathogen", "abx_class", syn_priors_list, "coeff", "se")]
  }else{
    abx_template <- make_predict_template(bdc, j, syn_type, template_data=template_data)
    syndromes <- c()
    ref_syn <- "none"
  }

  if (syn_type != "none"){
    rr_dat <- mr$MRData()
    rr_dat$load_df(
      data = data,  col_obs = "coeff", col_obs_se = "se",
      col_covs = c(syn_priors_list),
      col_study_id = "pathogen" )
  }else{
    rr_dat <- mr$MRData()
    rr_dat$load_df(
      data = data,  col_obs = "coeff", col_obs_se = "se",
      col_study_id = "pathogen" )
  }
  
  mr_detailed = run_mr_brt(abx_prior, syn_priors_list, abx_prior[2], syndrome_sd, gamma_SD[1], abx_template, rr_dat, syndromes, ref_syn)
  mr_strong_prior = run_mr_brt(abx_prior, syn_priors_list, abx_prior[3], syndrome_sd, gamma_SD[2], abx_template, rr_dat, syndromes, ref_syn)
  
  print("Splicing in result from prior where there is protective!")
  
  mr_pred_detailed <- mr_detailed$pred %>%
    mutate(se = abx_prior[2])
  mr_pred_sp <- mr_strong_prior$pred %>%
    mutate(se = abx_prior[3])
  
  mr_pred_chim <- mr_pred_detailed %>%
    mutate(predict = if_else(predict < 0, mr_pred_sp$predict, predict))
  
  mr_pred_det_check <- mr_pred_detailed %>%
    mutate(strong_prior = if_else(predict < 0, 1, 0))
  
  mr_draws_detailed <- mr_detailed$draws
  mr_draws_sp <- mr_strong_prior$draws
  
  if (syn_type != "none"){
    mr_draws_chim_det <- mr_draws_detailed %>%
      left_join(mr_pred_det_check[, c("abx_class", "pathogen", "infectious_syndrome", "strong_prior")], by=c("abx_class", "pathogen", "infectious_syndrome")) %>%
      filter(strong_prior == 0)
    mr_draws_chim_sp <- mr_draws_detailed %>%
      left_join(mr_pred_det_check[, c("abx_class", "pathogen", "infectious_syndrome", "strong_prior")], by=c("abx_class", "pathogen", "infectious_syndrome")) %>%
      filter(strong_prior == 1)
    mr_draws_chim <- rbind(mr_draws_chim_det, mr_draws_chim_sp)
  }else{
    mr_draws_chim_det <- mr_draws_detailed %>%
      left_join(mr_pred_det_check[, c("abx_class", "pathogen", "strong_prior")], by=c("abx_class", "pathogen")) %>%
      filter(strong_prior == 0)
    mr_draws_chim_sp <- mr_draws_detailed %>%
      left_join(mr_pred_det_check[, c("abx_class", "pathogen", "strong_prior")], by=c("abx_class", "pathogen")) %>%
      filter(strong_prior == 1)
    mr_draws_chim <- rbind(mr_draws_chim_det, mr_draws_chim_sp)
  }
  return(list(mr_draws_detailed = mr_draws_detailed,
    mr_draws_sp = mr_draws_sp,
    mr_draws_chim = mr_draws_chim, 
    mr_pred_detailed = mr_pred_detailed,
    mr_pred_sp = mr_pred_sp, 
    mr_pred_chim = mr_pred_chim, 
    abx_prior = abx_prior, 
    mr_detailed = mr_detailed, 
    mr_strong_prior = mr_strong_prior))
}


table_priors <- read.csv("FILEPATH")
syn_priors <- c()
for (j in antibiotic_group) {
  data <- filter_data(cf, ir2, j)
  model_run <- run_model(j, data, ir2, table_priors)
  mr_draws_sp <- model_run$mr_draws_sp
  mr_draws_detailed <- model_run$mr_draws_detailed
  mr_draws_chim <- model_run$mr_draws_chim 
  mr_pred_detailed <- model_run$mr_pred_detailed
  mr_pred_sp <- model_run$mr_pred_sp
  mr_pred_chim <- model_run$mr_pred_chim 
  abx_prior <- model_run$abx_prior
  mr_detailed <- model_run$mr_detailed
  mr_strong_prior <- model_run$mr_strong_prior

  
  print("Saving!")
  write.csv(mr_draws_detailed,"FILEPATH", row.names = FALSE)
  write.csv(mr_draws_sp,"FILEPATH", row.names = FALSE)
  write.csv(mr_draws_chim,"FILEPATH", row.names = FALSE)
  
  print("Creating summary files")
  sum_tab_det <- get_summary(mr_pred_detailed, abx_prior[2], mr_detailed$samp)
  sum_tab_sp <- get_summary(mr_pred_sp, abx_prior[3], mr_strong_prior$samp)
  sum_tab <- get_summary(mr_pred_chim, 0.1)
  
  write.csv(sum_tab_det,"FILEPATH", row.names = FALSE)
  write.csv(sum_tab_sp,"FILEPATH", row.names = FALSE)
  write.csv(sum_tab,"FILEPATH", row.names = FALSE)
}


print("Saving draws")
draws <- c() 
for (j in antibiotic_group) {
  abx_draw <-read.csv("FILEPATH", stringsAsFactors = FALSE)
  draws <- bind_rows(draws, abx_draw)
}
print(outpath)
write.csv(draws,"FILEPATH", row.names = FALSE)

print("Saving summary tables")
summary_table <- c()
for (j in antibiotic_group) {
  abx_sum<-read.csv("FILEPATH", stringsAsFactors = FALSE)
  summary_table <- bind_rows(summary_table,abx_sum)
}

write.csv(summary_table,"FILEPATH", row.names = FALSE)

write.csv(table_priors,"FILEPATH", row.names = FALSE)


if (run_isv == 1){
  isv <- c()
  for (j in antibiotic_group){
    data <- filter_data(cf, ir2, j)
    print(paste0("Working on ISV for ", j))
    model_run <- run_model(j, data, ir2, table_priors, template_data=data)
    mr_pred_chim <- model_run$mr_pred_chim
    
    isv <- rbind(isv, mr_pred_chim)
  }
    
  stopifnot(nrow(isv) == nrow(cf))
  isv <- data.table(isv) %>%
    rename(obs_log = obs) %>%
    mutate(obs = exp(obs_log)) %>%
    mutate(pred = exp(predict)) %>%
    mutate(obs_var = ((exp(obs_log) - exp(obs_log+(1.96*obs_se)))/1.96)^2) %>%
    mutate(pred_var = ((exp(predict) - exp(predict+(1.96*se)))/1.96)^2)
  
  isv <- isv %>%
    filter(!(is.na(pred)))
  
  isv_table <- c()
  
  for (j in antibiotic_group){
    isv_j <- isv %>%
      filter(abx_class == j)
    
    pj <- ggplot(isv_j, aes(x = obs, y = pred)) + 
      geom_point() +
      geom_smooth(method = lm, se = FALSE, color = "blue") +
      labs(x = "Observed", y = "Predicted", title = "Predicted vs Observed") +
      theme_minimal()+
      ggtitle(paste0("Obs vs Preds, ", j))
    
    ggsave(paste0(j,"_obs_pred.png"), plot = pj, path = "FILEPATH", width = 10, height = 6, dpi = 300)
    
    rmse_is <- round(rmse(isv_j$obs, isv_j$pred),2)
    mae_is <- round(mae(isv_j$obs, isv_j$pred),2)
    
    rss <- sum((isv_j$obs - isv_j$pred) ^ 2)
    tss <- sum((isv_j$obs - mean(isv_j$obs)) ^ 2)
    rsq_is <- 1 - (rss/tss)
    
    isv_j$var_adj_is <- isv_j$obs_var + isv_j$pred_var
    isv_j$pred_lb_adj_is <- isv_j$pred - 1.96*sqrt(isv_j$var_adj_is)
    isv_j$pred_ub_adj_is <- isv_j$pred + 1.96*sqrt(isv_j$var_adj_is)
    
    
    isv_j$covered <- ifelse((isv_j$obs >= isv_j$pred_lb_adj_is) & (isv_j$obs <= isv_j$pred_ub_adj_is),1,0)
    coverage_is <- nrow(isv_j[isv_j$covered == 1,])/length(isv_j$covered)
    
    is <- c(paste0(j),rsq_is,rmse_is,mae_is, coverage_is)
    isv_table <- rbind(isv_table,is)
  }
  
  print("Saving ISV")
  isv_table <- as.data.frame(isv_table)
  colnames(isv_table) <- c('abx_class','R2','RMSE','MAE', 'coverage')
  isv_table <- left_join(isv_table, abx_meta) 
  
  write.csv(isv_table, "FILEPATH", row.names = F)
  }


if (run_oosv == 1){
  print("Out of Sample Validation")
  
  table_oos <- c()
  
  for (j in antibiotic_group){
    data <- filter_data(cf, ir2, j)
    data <- cbind(ID = rownames(data), data)
    rownames(data) <- 1:nrow(data)
    
    folds <- floor(runif(nrow(data), min = 0, max = 4)) + 1
    folds <- data %>% split(folds)

    preds <- c()
    for (i in 1:4){
      test <- folds[[i]]
      train <- setDT(data[!data$ID %in% test$ID,])

      stopifnot(nrow(test) + nrow(train) == nrow(data))
      
      test = subset(test, select = -c(ID))
      train = subset(train, select = -c(ID))
      

      print(paste0("Working on OOSV for ", j))
      model_run <- run_model(j, data, ir2, table_priors, template_data=test)
      mr_pred_chim <- model_run$mr_pred_chim
      

      preds <- rbind(preds, mr_pred_chim)
    }
    oos_preds <- rbind(preds)
    stopifnot(nrow(oos_preds) == nrow(data))
    oos_preds <- data.table(oos_preds) %>%
      rename(obs_log = obs) %>%
      mutate(obs = exp(obs_log)) %>%
      mutate(pred = exp(predict)) %>%
      mutate(obs_var = ((exp(obs_log) - exp(obs_log+(1.96*obs_se)))/1.96)^2) %>%
      mutate(pred_var = ((exp(predict) - exp(predict+(1.96*se)))/1.96)^2)
    
    oos_preds_j <- oos_preds %>%
      filter(abx_class == j)
    
    pj <- ggplot(oos_preds_j, aes(x = obs, y = pred)) + 
      geom_point() +
      geom_smooth(method = lm, se = FALSE, color = "blue") +
      labs(x = "Observed", y = "Predicted", title = "Predicted vs Observed") +
      theme_minimal()+
      ggtitle(paste0("OOSV Obs vs Preds, ", j))
    
    ggsave(paste0(j,"_obs_pred.png"), plot = pj, path = "FILEPATH", width = 10, height = 6, dpi = 300)
    
    rmse_oos <- rmse(oos_preds_j$obs, oos_preds_j$predict)
    
    rss <- sum((oos_preds_j$obs - oos_preds_j$pred)^ 2)
    tss <- sum((oos_preds_j$obs - mean(oos_preds_j$obs) ^ 2))
    rsq <- 1 - (rss/tss)
    
    mae_oos <- mae(oos_preds_j$obs, oos_preds_j$pred)
    
    oos_preds_j$var_adj_is <- oos_preds_j$obs_var + oos_preds_j$pred_var
    oos_preds_j$pred_lb_adj_is <- oos_preds_j$pred - 1.96*sqrt(oos_preds_j$var_adj_is)
    oos_preds_j$pred_ub_adj_is <- oos_preds_j$pred + 1.96*sqrt(oos_preds_j$var_adj_is)
    
    
    oos_preds_j$covered <- ifelse((oos_preds_j$obs >= oos_preds_j$pred_lb_adj_is) & (oos_preds_j$obs <= oos_preds_j$pred_ub_adj_is),1,0)
    coverage_oos <- nrow(oos_preds_j[oos_preds_j$covered == 1,])/length(oos_preds_j$covered)
    
    table <- c(paste0(j),rsq,round(rmse_oos,2),round(mae_oos,2), coverage_oos)
    table_oos <- rbind(table_oos, table)
  }
  
  print("Saving OOSV")
  table_oos <- as.data.frame(table_oos)
  colnames(table_oos) <- c('abx_class','R2','RMSE','MAE', 'coverage')
  table_oos <- left_join(table_oos, abx_meta) 
  
  write.csv(table_oos, "FILEPATH", row.names = F)

  }
