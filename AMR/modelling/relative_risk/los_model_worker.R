library(reticulate)
reticulate::use_python("FILEPATH")
library(dplyr)
library(readr)
library(data.table)
library(tidyverse)
library(assertthat)
library(metafor)
library(lme4)
library(lmtest)
library(sandwich)

source("FILEPATH")
source("FILEPATH")
source("FILEPATH")

mr <- import("mrtool")
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

gamma_SD <- c(0.0001, 0.000001)
intercept_sd <- 0.0001
intercept_sd_strong <- 0.00001
syndrome_sd <- c(0.01, 0.00001)
n_samples <- 1000L


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

generate_crude_rr <- function(df, filter_cols){
  rr_df <- df %>% 
    filter(is.na(mean_LOS_r) == TRUE) %>%
    dplyr::group_by_at(vars(filter_cols)) %>%
    summarize(days_infection = sum(days_infection), cases = sum(cases))
  
  if (!("infectious_syndrome" %in% filter_cols)){
    rr_df$infectious_syndrome <- "all"
  }
  if (!("source" %in% filter_cols)){
    rr_df$source <- "all"
  }
  
  
  rr_nolit <- c()
  for (j in unique(rr_df$abx_class)){
    for (p in unique(rr_df$pathogen)){
      for (s in unique(rr_df$source)){
        for (i in unique(rr_df$infectious_syndrome)){
          print(paste0("Working on ", j, " ", p, " ", s, " ", i))
          data <- rr_df %>%
            filter(abx_class == j) %>%
            filter(pathogen == p) %>%
            filter(source == s) %>%
            filter(infectious_syndrome == i)
          total <- data %>% dplyr::group_by(abx_class, resistant) %>% summarize(admit = sum(cases), los = mean(days_infection,na.rm = TRUE))
          if(!is.na(total[1,3]) & !is.na(total[2,3]) & !is.na(total[2,4]) & !is.na(total[1,4])) {
            rr_poisson_mod1 <- glm(data = data, formula = days_infection ~ resistant, family = 'poisson'(link='log'))
            
            r1 <- c(s, p, j, i, summary.glm(rr_poisson_mod1)$coefficients[2,1:2],total[1,3],total[1,4],total[2,3],total[2,4])
            rr_nolit <- rbind(rr_nolit, r1)
          }else{
          print("Moving onto next combo")
        }
        }
      }
  }
  }
  rr_nolit2 <- as.data.frame(rr_nolit)
  colnames(rr_nolit2) <- c('source','pathogen','abx_class', "infectious_syndrome", 'log_crude_rr','lnrr_se', 'cases_s','mean_LOS_s','cases_r','mean_LOS_r')
  rr_nolit2$log_crude_rr <- sapply(rr_nolit2$log_crude_rr, function(x) as.numeric(x))
  rr_nolit2$lnrr_se <- sapply(rr_nolit2$lnrr_se, function(x) as.numeric(x))
  rr_nolit2$cases_s <- sapply(rr_nolit2$cases_s, function(x) as.numeric(x))
  rr_nolit2$cases_r <- sapply(rr_nolit2$cases_r, function(x) as.numeric(x))
  rr_nolit2$mean_LOS_r <- sapply(rr_nolit2$mean_LOS_r, function(x) as.numeric(x))
  rr_nolit2$mean_LOS_s <- sapply(rr_nolit2$mean_LOS_s, function(x) as.numeric(x))
  
  
  rr_lit <- df %>%
    filter(is.na(mean_LOS_r) == FALSE) %>%
    mutate(ri = 0)
  
  rr_lit <- escalc(data=rr_lit, measure="ROMC", m1i=mean_LOS_r,sd1i=sd_r,m2i=mean_LOS_s,sd2i = sd_s,ni=sample_size,ri= ri, append = TRUE)
  rr_lit <- rr_lit %>% 
    mutate(lnrr_se = sqrt(vi)) %>%
    mutate(log_crude_rr = yi)
  rr_lit <- rr_lit[,colnames(rr_nolit2)]
  
  rr_fin <- rr_nolit2 %>%
    rbind(rr_lit) %>%
    mutate(combo = paste(pathogen, abx_class, sep = "-")) %>%
    filter(!is.na(log_crude_rr)) %>%
    filter(!is.infinite(log_crude_rr)) %>%
    filter(lnrr_se != 0)
  
  rr_fin <- as.data.frame(rr_fin)
  
  print("Fixing list issue")
  list_columns <- sapply(rr_fin , function(x) class(x) == "list")
  list_column_names <- names(rr_fin )[list_columns]
  print(list_column_names)
  
  for (col_name in list_column_names) {
    rr_fin[[col_name]] <- sapply(rr_fin [[col_name]], function(x) x[1])
  }
  
  return(rr_fin)
}


antibiotic_group <- c('third_gen_ceph','carbapenem','fluoroquinolone','penicillin','aminopenicillin','beta_lactamase_inhibitor','anti_pseudomonal_penicillin','methicillin','vancomycin','fourth_gen_ceph','sulfa','aminoglycoside','macrolide')

ir2 <- read.csv(paste0(data_path,'intrinsic_resistance.csv'), stringsAsFactors = F) %>% filter(include==1) %>% dplyr::select(abx_class,pathogen) %>% distinct()
ir2$pathogen[ir2$pathogen == "acinetobacter_baumanii"] <- 'acinetobacter_baumannii'
ir <- read.csv(paste0(data_path,'intrinsic_resistance.csv'), stringsAsFactors = F) %>% dplyr::select(abx_class, atc) %>% distinct()
ir$atc[ir$abx_class == 'vancomycin'] <- 'J01D'
ir<-ir[!ir$abx_class %in% c('polymixin','chloramphenicol','nitrofurantoin','second_gen_ceph','rifamycin'),]
irall<-c('abx_class' = 'all_resistant', 'atc' = 'ddd_per_1000')
ir <- rbind(ir,irall)

bdc <- read.csv("FILEPATH")
bdc <- bdc %>%
  filter(NAMES_exclude == 0) %>%
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
df <- bind_rows(data_list) %>%
  mutate(combo = paste(pathogen, abx_class, sep = "-")) %>%
  filter(combo %in% unique(bdc$combo))

stage_1 <- generate_crude_rr(df, c("abx_class", "source", "pathogen", "resistant"))
write.csv(stage_1, "FILEPATH", row.names=FALSE)

stage_1 <- stage_1 %>%
  filter(!(abx_class == "fluoroquinolone" & log_crude_rr > 24)) %>%
  filter(!(abx_class == "fluoroquinolone" & log_crude_rr < 0 & lnrr_se > 0.05)) %>%
  filter(!(abx_class %in% c("beta_lactamase_inhibitor",
                            "methicillin") & log_crude_rr < 0 & lnrr_se > 0.1)) %>%
  filter(!(abx_class %in% c("vancomycin") & log_crude_rr < 0 & lnrr_se > 0.05)) %>%
  filter(!(abx_class %in% c("anti_pseudomonal_penicillin", "fourth_gen_ceph", 
                            "sulfa", "macrolide") & log_crude_rr < 0 & lnrr_se > 0.01)) %>%
  filter(!(abx_class %in% c("aminoglycoside") & log_crude_rr < 0 & lnrr_se > 0.001)) %>%
  filter(!(abx_class %in% c("penicillin") & log_crude_rr < -0.15)) %>%
  filter(!(abx_class == "aminopenicillin" & pathogen == "escherichia_coli" & log_crude_rr > 2)) %>%
  filter(!(abx_class == "methicillin" & log_crude_rr < 0)) %>%
  filter(!(abx_class %in% c("aminoglycoside", "aminopenicillin",
                            "beta_lactamase_inhibitor", 
                            "macrolide", "methicillin") & log_crude_rr > 2)) %>%
  filter(!(abx_class %in% c("fluoroquinolone") & log_crude_rr > 3)) %>%
  filter(!(abx_class %in% c("aminopenicillin", "macrolide", "fourth_gen_ceph",  "anti_pseudomonal_penicillin") & log_crude_rr > 1.4)) %>%
  filter(!(abx_class %in% c("aminoglycoside") & log_crude_rr > 0.7)) %>%
  filter(!(abx_class %in% c("carbapenem", "third_gen_ceph", "vancomycin") & log_crude_rr < -1))
  
if (syn_type != "none"){
  stage_2 <- generate_crude_rr(df, c("abx_class", "pathogen", "infectious_syndrome", "resistant"))
}else{
  stage_2 <- generate_crude_rr(df, c("abx_class", "pathogen", "resistant"))
}
write.csv(stage_2, "FILEPATH", row.names=FALSE)

og_data <- stage_2
og_data$key <- apply(og_data, 1, paste, collapse = "-")

stage_2_data <- stage_2 %>%
  filter(!(abx_class %in% c("methicillin") & log_crude_rr < 0 & pathogen == "streptococcus_pneumoniae")) %>%
  filter(!(abx_class %in% c("third_gen_ceph", "carbapenem") & log_crude_rr < 0 & pathogen == "pseudomonas_aeruginosa")) %>%
  filter(!(pathogen == "acinetobacter_baumannii" & log_crude_rr > 1 & abx_class %in% c("fluoroquinolone", "third_gen_ceph","carbapenem", "fourth_gen_ceph", "anti_pseudomonal_penicillin"))) %>%
  filter(!(pathogen %in% c("morganella_spp", "proteus_spp", "pseudomonas_aeruginosa", "staphylococcus_aureus") & log_crude_rr < 0 & abx_class %in% c("fluoroquinolone"))) %>%
  filter(!(pathogen %in% c("enterococcus_faecalis", "enterococcus_faecium") & abx_class %in% c("vancomycin"))) %>%
  filter(!(pathogen %in% c("escherichia_coli", "klebsiella_pneumoniae") & log_crude_rr < 0 & abx_class %in% c("third_gen_ceph", "carbapenem", "fluoroquinolone", "beta_lactamase_inhibitor", "sulfa", "aminoglycoside"))) %>%
  filter(!(pathogen %in% c("proteus_spp", "pseudomonas_aeruginosa", "serratia_spp",
                           "enterobacter_spp", "citrobacter_spp", "morganella_spp", "staphylococcus_aureus") & log_crude_rr < -0.6 & !(abx_class %in% c("anti_pseudomonal_penicillin", "fourth_gen_ceph")))) %>%
  filter(!(pathogen %in% c("enterococcus_faecium") & log_crude_rr > 2)) %>%
  filter(!(pathogen %in% c("streptococcus_pneumoniae") & log_crude_rr < -0.5)) %>%
  filter(!(abx_class %in% c("anti_pseudomonal_penicillin", "fourth_gen_ceph") & log_crude_rr < -0.8))

stage_2_data <- stage_2 %>%
  mutate(log_crude_rr = ifelse(abx_class == "anti_pseudomonal_penicillin" & pathogen %in% c("enterobacter_spp", "pseudomonas_aeruginosa") & log_crude_rr < 0, 0, log_crude_rr)) %>%
  mutate(log_crude_rr = ifelse(abx_class == "fourth_gen_ceph" & pathogen %in% c("pseudomonas_aeruginosa") & log_crude_rr < 0, 0, log_crude_rr)) %>%
  mutate(log_crude_rr = ifelse(abx_class == "sulfa" & pathogen %in% c("proteus_spp") & log_crude_rr < 0, 0, log_crude_rr))

stage_2_data$key <- apply(stage_2_data, 1, paste, collapse = "-")
  
outlier <- og_data %>%
  rowwise() %>%
  mutate(outlier = if_else(key %in% stage_2_data$key, 0, 1)) %>%
  ungroup() %>%
  select(-key)

write.csv(outlier,"FILEPATH", row.names = FALSE)
  



stage1_priors <- c()
for (j in antibiotic_group) { 
  
  stage1_data <- stage_1 %>%
    filter(abx_class == j)
  
  stage1_dat_prior <- mr$MRData()
  stage1_dat_prior$load_df(
    data = stage1_data,  col_obs = "log_crude_rr", col_obs_se = "lnrr_se",
    col_study_id = "source" )
  
  stage1_mrbrtprior <- mr$MRBRT(
    data = stage1_dat_prior,
    cov_models = list(mr$LinearCovModel('intercept', use_re = TRUE)
                      )
    )
  
  stage1_mrbrtprior$fit_model(inner_print_level = 5L, inner_max_iter = 1000L)
  stage1_abx_prior <- stage1_mrbrtprior$beta_soln[1]
  
  stage1_priors_row <- c(j, stage1_abx_prior, exp(stage1_abx_prior), intercept_sd, intercept_sd_strong)
  stage1_priors <- rbind(stage1_priors, stage1_priors_row)
}
colnames(stage1_priors) <- c('abx_class', 'log_RR', 'RR', 'intercept_sd', 'intercept_sd_strong')




run_mr_brt <- function(abx_prior, syn_priors_list, intercept_sd, syndrome_sd, gamma_sd, abx_template, rr_dat, syndromes, ref_syn){
  if (syn_type != "none"){
    syn_model <- paste0("mr$LinearCovModel(\"",syn_priors_list[,1], "\", prior_beta_gaussian = array(c(0,", syndrome_sd, ")))")
    model_definition <- paste(syn_model, collapse = ",\n")
  }
  mr_pred_dat <- mr$MRData()
  if(syn_type == "all"){
    mr_pred_dat$load_df(
      data = abx_template,
      col_covs = c(syn_priors_list),
      col_study_id = "pathogen" )
  }else{
    mr_pred_dat$load_df(
      data = abx_template,
      col_study_id = "pathogen" )
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
    mr_pred <- mr_pred %>% 
      mutate(infectious_syndrome = case_when(
        L2_blood_stream_infection == 1 ~ "L2_blood_stream_infection",
        L2_urinary_tract_infection == 1 ~ "L2_urinary_tract_infection",
        other_syndromes == 1 ~ "other_syndromes",
        L2_lower_respiratory_infection == 1  ~ "L2_lower_respiratory_infection"
      )) %>%
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


table_priors <- as.data.frame(stage1_priors)
syn_priors <- c()
cf <- stage_2_data %>%
  mutate(coeff = log_crude_rr) %>%
  mutate(se = lnrr_se)

for (j in antibiotic_group) {
  print(paste0("Working on ", j))
  print("Making priors")
  nonintrinsic <- unique(ir2$pathogen[ir2$abx_class==j])
  
  data <- cf[cf$abx_class == paste0(j),] %>%
    filter(pathogen %in% nonintrinsic)
  
  abx_prior <- c(table_priors[table_priors$abx_class == j,]$log_RR, intercept_sd, intercept_sd_strong)
  
  if (syn_type != "none"){
    data <- data %>%
      filter(!infectious_syndrome == "L2_unspecified") %>%
      filter(!infectious_syndrome == "L2_unspecified_infection")
    
    syndromes <- unique(data$infectious_syndrome)
    ref_syn <- "other_syndromes"
  
    syn_priors_list <- c()
    for (i in syndromes) {
      if (i != ref_syn){
        data$new <- ifelse(data$infectious_syndrome == i,1,0)
        colnames(data)[colnames(data) == 'new'] <- paste0(i)
        syn_effs <- c(i)
        syn_priors_list <- rbind(syn_priors_list, syn_effs)
      }
    }
    
    print("Making template")
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
    
    
    data <- data[,c("pathogen", "abx_class", syn_priors_list, "coeff", "se")]
  }else{
    syn_priors_list <- c()
    print("Making template")
    unique_pathogens <- bdc %>%
      filter(abx_class == j) %>%
      distinct(pathogen)
    
    abx_template <- expand.grid(abx_class = j,
                                pathogen = unique_pathogens$pathogen)
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
  
  mr_detailed = run_mr_brt(abx_prior, syn_priors_list, abx_prior[2], syndrome_sd[1], gamma_SD[1], abx_template, rr_dat, syndromes, ref_syn)
  mr_strong_prior = run_mr_brt(abx_prior, syn_priors_list, abx_prior[3], syndrome_sd[2], gamma_SD[2], abx_template, rr_dat, syndromes, ref_syn)
  
  print("Splicing in result from prior where there is protective!")
  
  mr_pred_detailed <- mr_detailed$pred
  mr_pred_sp <- mr_strong_prior$pred
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
  
  print("Saving!")
  write.csv(mr_draws_detailed,"FILEPATH", row.names = FALSE)
  write.csv(mr_draws_sp,"FILEPATH", row.names = FALSE)
  write.csv(mr_draws_chim,"FILEPATH", row.names = FALSE)
  
  print("Creating summary files")
  sum_tab_det <- get_summary(mr_pred_detailed, as.numeric(abx_prior[2]), mr_detailed$samp)
  sum_tab_sp <- get_summary(mr_pred_sp, as.numeric(abx_prior[3]), mr_strong_prior$samp)
  sum_tab <- get_summary(mr_pred_chim, 0.1)
  
  write.csv(sum_tab_det,"FILEPATH", row.names = FALSE)
  write.csv(sum_tab_sp,"FILEPATH", row.names = FALSE)
  write.csv(sum_tab,"FILEPATH", row.names = FALSE)
}


draws <- c() 
for (j in antibiotic_group) {
  abx_draw <-read.csv("FILEPATH", stringsAsFactors = FALSE)
  draws <- bind_rows(draws, abx_draw)
}
write.csv(draws,"FILEPATH", row.names = FALSE)

summary <- c()
for (j in antibiotic_group) {
  abx_sum<-read.csv("FILEPATH", stringsAsFactors = FALSE)
  summary <- bind_rows(summary,abx_sum)
}

write.csv(summary,"FILEPATH", row.names = FALSE)

write.csv(table_priors,"FILEPATH", row.names = FALSE)

