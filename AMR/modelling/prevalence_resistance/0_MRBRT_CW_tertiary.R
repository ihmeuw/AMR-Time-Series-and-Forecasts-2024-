##################################################
## Project: Antimicrobial Resistance - CW Resistance Input Data
## Script purpose: Crosswalk tertiary resistance levels to
##                 non-tertiary prior to running Ensemble-STGPR
##################################################
##### R Initialization and Functions #####
Sys.setenv("RETICULATE_PYTHON" = 'FILEPATH')
reticulate::use_python("FILEPATH")
library(reticulate)
mrtool <- NULL
xwalk <- NULL
.onLoad <- function(libname, pkgname) {
  use_python(python = "FILEPATH")
  mrtool <<- import("mrtool")
  for (nm in names(mrtool)) assign(nm, mrtool[[nm]], parent.env(environment()))
  xwalk <<- import("crosswalk")
  for (nm2 in c("linear_to_log", "linear_to_logit", "log_to_linear", "logit_to_linear")) {
    assign(nm2, xwalk[["utils"]][[nm2]], parent.env(environment()))
  }
}
.onLoad()
exists("core")
###
username <- Sys.getenv("USER")
hdrive <- sprintf("FILEPATH", username)
library(dplyr)
library(ggplot2)
source("FILEPATH")
library(reticulate)
reticulate::use_python("FILEPATH")
add_loc_info <- function(df){
  locs = get_location_metadata(location_set_id = 35, 
                               release_id =  16)[, c("location_id", "location_name", "region_name", "super_region_name")]
  df <- merge(df, locs, by = 'location_id', all.x = TRUE)
  return(df)
}

# function to convert logits to probabilities
logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}

delta_transform <- function(mean, sd, transformation) {
  
  if (transformation == "linear_to_log") f <- xwalk$utils$linear_to_log
  if (transformation == "log_to_linear") f <- xwalk$utils$log_to_linear
  if (transformation == "linear_to_logit") f <- xwalk$utils$linear_to_logit
  if (transformation == "logit_to_linear") f <- xwalk$utils$logit_to_linear
  
  out <- do.call("cbind", f(mean = array(mean), sd = array(sd)))
  colnames(out) <- paste0(c("mean", "sd"), "_", strsplit(transformation, "_")[[1]][3])
  return(out)
}

calculate_diff <- function(df, alt_mean, alt_sd, ref_mean, ref_sd) {
  df <- as.data.frame(df)
  out <- data.frame(
    diff_mean =  df[, alt_mean] - df[, ref_mean],
    diff_sd = sqrt(df[, alt_sd]^2 + df[, ref_sd]^2)
  )
  return(out)
}

# set seed for reproducibility
set.seed(123)

modelpath <- "FILEPATH"
output_path <- 'FILEPATH'

# READ combos
combos<-  read.csv('FILEPATH', stringsAsFactors = FALSE) %>% filter(COLUMN==0) %>% dplyr::select(pathogen, abx_class) 
combos$pathogen[combos$pathogen == "acinetobacter_baumanii"] <- 'acinetobacter_baumannii'
combos$pathogen[grepl('gonorrheae', combos$pathogen)] <- 'neisseria_gonorrhoeae'
combos$abx_class[combos$abx_class=='rifampicin_new']<-'mdr'
combos$abx_class[combos$abx_class=='isoniazid_new']<-'xdr'
combos <- combos[!grepl('prop',combos$abx_class),]
combos$abx_class <- sub('_retreated',"",combos$abx_class)
combos$combo<-paste0(combos$pathogen,"-",combos$abx_class)
location_md2 <- subset(get_location_metadata(location_set_id=35, release_id = 16),select = c(location_id,level,location_name,location_name_short,super_region_id,super_region_name, region_name,region_id,_loc_id,parent_id), droplevels = TRUE) 

oxc <- c('escherichia_coli-third_gen_ceph','escherichia_coli-fluoroquinolone','klebsiella_pneumoniae-third_gen_ceph','klebsiella_pneumoniae-carbapenem',
         'staphylococcus_aureus-methicillin','streptococcus_pneumoniae-penicillin',"salmonella_paratyphi-fluoroquinolone","salmonella_typhi-fluoroquinolone",
         'shigella-fluoroquinolone',"salmonella_ints_non_typhi/paratyphi-fluoroquinolone","salmonella_paratyphi-mdr","salmonella_typhi-mdr",
         "neisseria_gonorrhoeae-third_gen_ceph","non_typhoidal_salmonellae-fluoroquinolone")
active4 <- read.csv(paste0('FILEPATH')) %>% filter(active_4=="1") %>% dplyr::select(source)

##### Read in Data #####
unadjusted_data <- read.csv('FILEPATH') %>%
  mutate(combo = paste0(pathogen,'-',abx_class)) %>% filter(combo %in% unique(combos$combo))
unadjusted_data <- unadjusted_data[unadjusted_data$source %in% c('SOURCE'),]
unadjusted_data[,c('level_0','index')] <- NULL

resdat <- fread('FILEPATH')
resdat <- resdat %>% mutate(combo = paste0(pathogen,'-',abx_class)) %>% filter(combo %in% unique(combos$combo))
resdat <- resdat[!((resdat$combo %in% oxc) & (resdat$source %in% c('SOURCE','SOURCE'))),]
resdat <- resdat %>% filter(cases != 0) %>% filter(cases < 10000000) %>% filter(cases >= resistant) 
resdat[,c('index')] <- NULL
resdat <- rbind(resdat,unadjusted_data)
resdat$nid <- as.character(resdat$nid)

unique(active4$source[!(active4$source %in% unique(resdat$source))])
unique(resdat$source[!(resdat$source %in% unique(active4$source))])

# merge in oxford data
file <- list.files(path = paste0("FILEPATH"),
                   pattern = '*.csv', full.names = T) 
file <- file[!grepl("FILEPATH",file)]

oxdata<-data.frame()
for(i in file){
  raw_df <- read.csv(i)
  if(!('tertiary' %in% colnames(raw_df))){
    raw_df$tertiary <- 2L
  }
  if(!('is_outlier' %in% colnames(raw_df))){
    raw_df$is_outlier <- 0L
  }
  if(!('variance' %in% colnames(raw_df))){
    raw_df$variance <- raw_df$val*(1-raw_df$val)*raw_df$sample_size
  }
  if(!('n_resistant' %in% colnames(raw_df))){
    raw_df$n_resistant <- raw_df$val*raw_df$sample_size
  }
  if(!('pathogen' %in% colnames(raw_df))){
    raw_df$pathogen <- ifelse(grepl("paratyphi",i),'salmonella_paratyphi',ifelse(grepl("nts",i),'non_typhoidal_salmonellae','salmonella_typhi'))
  }
  if(!('antimicrobial' %in% colnames(raw_df))){
    raw_df$antimicrobial <- ifelse(grepl("fqn",i),'fluoroquinolone',ifelse(grepl("FQN",i),'fluoroquinolone','mdr'))
  }  
  raw_df <- raw_df %>% dplyr::select("country", "year_id", "location_id","pathogen","antimicrobial","nid","n_resistant","val","sample_size",
                                   "tertiary", "variance","is_outlier")
  raw_df$filename <- paste0(substring(paste0(i),87,101))
  oxdata<-rbind(oxdata,raw_df)
}
table(oxdata$filename,oxdata$pathogen)
oxdata$source <- oxdata$nid
oxdata <- oxdata[!is.na(oxdata$val),]
oxdata$location_id[(oxdata$country == 'Kosovo')] <- 53
oxdata$country <- NULL

duplicated_sources <- c('SOURCE') 
oxdata_to_dedup <- oxdata
oxdata_to_dedup <- oxdata_to_dedup[!(oxdata_to_dedup$source %in% duplicated_sources),]

# PUT DATA TOGETHER
oxdata <- oxdata_to_dedup
oxdata$pathogen <- tolower(sub(" ","_",oxdata$pathogen))
oxdata$pathogen[oxdata$pathogen == "shigella"] <- 'shigella_spp'
oxdata$pathogen[oxdata$pathogen == "non_typhoidal_salmonellae"] <- "salmonella_ints_non_typhi/paratyphi"
oxdata$antimicrobial <- tolower(oxdata$antimicrobial)
oxdata$abx_class <- case_when(oxdata$antimicrobial %in% c('3gc','third-generation cephalosporins','3rd gen cephalosporin','3rd gen ceph',
                                      'ceftriaxone','cefixime','ceftazidime','cefoperazone','ceftriaxone/cefotaxime','cefotaxime') ~ 'third_gen_ceph',
                              oxdata$antimicrobial %in% c('fluoroquinolones','ciprofloxacin','ofloxacin','norfloxacin','levofloxacin', 
                                                          'moxifloxacin','nalidixic acid') ~ 'fluoroquinolone',
                              oxdata$antimicrobial %in% c('carbapenems','carbapenem') ~ 'carbapenem',
                              oxdata$antimicrobial %in% c('penicillin','methicillin','mdr') ~ oxdata$antimicrobial)
oxdata$antimicrobial <- NULL

oxdata$hospital_type <- ifelse(oxdata$tertiary == 1, 'tertiary', 'nontert-mixed')
oxdata$cases <- oxdata$sample_size
oxdata$resistant <- oxdata$n_resistant
oxdata$susceptible <- oxdata$sample_size - oxdata$n_resistant
oxdata$combo <- paste0(oxdata$pathogen, '-', oxdata$abx_class)
oxdata <- oxdata[,colnames(oxdata) %in% c("location_id","year_id","source","nid","pathogen","abx_class","resistant","susceptible","cases","combo","hospital_type","variance","is_outlier")]
write.csv(oxdata, paste0(output_path,'FILEPATH'), row.names = F)

resdat <- bind_rows(resdat, oxdata) %>% dplyr::select(location_id,year_id,source,nid,pathogen,abx_class,combo,resistant,susceptible,cases,variance,is_outlier,hospital_type)

# reasigning subnational units to national location
resdat <- left_join(resdat, location_md2[,c('location_id','level','parent_id')], by = c('location_id'))
resdat$location_id <- ifelse((resdat$level == 5), resdat$parent_id, resdat$location_id)
resdat$level<-NULL
resdat$parent_id<-NULL
resdat <- left_join(resdat, location_md2[,c('location_id','level','parent_id')], by = c('location_id'))
resdat$location_id <- ifelse((resdat$level == 4), resdat$parent_id, resdat$location_id)
resdat$level<-NULL
resdat$parent_id<-NULL

resdat <- resdat %>% group_by(location_id,year_id,source,nid,pathogen,abx_class,combo,hospital_type) %>% summarise(resistant = sum(resistant, na.rm = T), 
                                                                                                 susceptible = sum(susceptible, na.rm = T), 
                                                                                                 cases = sum(cases, na.rm = T))
resdat <- resdat %>% left_join(location_md2[,c("location_id",'level','super_region_id','super_region_name','region_id','region_name','_loc_id')])

resdat$val <- resdat$resistant/resdat$cases
resdat$val[resdat$val == 0] <- 0.001
resdat$val[resdat$val == 1] <- 0.999
# calculate SEs
resdat$se <- sqrt(((resdat$val)*(1-resdat$val))/resdat$cases)
resdat$variance <- ((resdat$val)*(1-resdat$val))/resdat$cases
resdat$variance[is.na(resdat$variance)] <- resdat$se[is.na(resdat$variance)]*resdat$se[is.na(resdat$variance)]

resdat$sample_size <-  resdat$cases
resdat <- resdat[,colnames(resdat)[!colnames(resdat) %in% c('val','se')]]


# drop overlaps 
tertdat <- resdat[resdat$source != 'SOURCE' | resdat$location_id != 95,]

# add location and hospital group
tertdat <- tertdat %>% left_join(location_md2[,c('location_id','region_name','super_region_name')]) 
tertdat$tert_bin <- ifelse(tertdat$hospital_type %in% 'tertiary', 'tertiary', 'nontert-mixed')

# limit to data points with sample size >= 5
tertdat <- tertdat[tertdat$cases >= 5,]

# calculate prevalence and SEs
tertdat$prev <- tertdat$resistant/tertdat$cases
tertdat$prev[tertdat$prev == 0] <- 0.001
tertdat$prev[tertdat$prev == 1] <- 0.999
# calculate SEs
tertdat$se <- sqrt((tertdat$prev*(1-tertdat$prev))/tertdat$cases)


# hospital type counts by super region
htcounts <- tertdat[!is.na(tertdat$super_region_name),] %>% dplyr::group_by(super_region_name) %>%
  dplyr::summarize(
    non_tertiary = sum(cases[hospital_type == 'non-tertiary'])/sum(cases),
    tertiary = sum(cases[hospital_type == 'tertiary'])/sum(cases),
    mixed_unknown = sum(cases[hospital_type == 'm/u'])/sum(cases)
  )

#aggregate to super region
tertdat <- tertdat %>% dplyr::group_by(pathogen, abx_class, source,nid, tert_bin, super_region_name, year_id) %>%
  dplyr::summarize(
    resistant = sum(resistant),
    susceptible = sum(susceptible),
    cases = sum(cases)
  )

# limit to data points with sample size >= 5
tertdat <- tertdat[tertdat$cases >= 5,]

# calculate prevalence and SEs
tertdat$prev <- tertdat$resistant/tertdat$cases
tertdat$prev[tertdat$prev == 0] <- 0.001
tertdat$prev[tertdat$prev == 1] <- 0.999
# calculate SEs
tertdat$se <- sqrt((tertdat$prev*(1-tertdat$prev))/tertdat$cases)

# create index for investigating matches later
tertdat$index <- seq.int(nrow(tertdat))

# create matching dataframe
keepcols <- c('location_id', 'region_name', 'super_region_name', 'pathogen', 'abx_class',
              'year_id', 'source', 'nid','hospital_type','is_outlier' ,'prev', 'se', 'index', 'cases')

df_ref <- tertdat[tertdat$tert_bin == 'nontert-mixed', colnames(tertdat) %in% keepcols]
df_alt <- tertdat[tertdat$tert_bin == 'tertiary', colnames(tertdat) %in% keepcols]

df_matched <- merge(df_ref, df_alt, by = c('super_region_name', 'pathogen', 'abx_class'), suffixes = c('_ref', '_alt'), allow.cartesian=TRUE)

rm(df_ref, df_alt, iord, oxdata, tertdat)

# calculate logit_diffs and SEs using delta method
dat_diff <- as.data.frame(cbind(
  delta_transform(
    mean = df_matched$prev_alt, 
    sd = df_matched$se_alt,
    transformation = "linear_to_logit" ),
  delta_transform(
    mean = df_matched$prev_ref, 
    sd = df_matched$se_ref,
    transformation = "linear_to_logit")
))
names(dat_diff) <- c("mean_alt", "mean_se_alt", "mean_ref", "mean_se_ref")

df_matched[, c("logit_diff", "logit_diff_se")] <- calculate_diff(
  df = dat_diff, 
  alt_mean = "mean_alt", alt_sd = "mean_se_alt",
  ref_mean = "mean_ref", ref_sd = "mean_se_ref" )

df_matched$altvar <- 'tertiary'
df_matched$refvar <- 'nontert-mixed'

# limit to matches at most 5 years apart
match5yr <- df_matched[abs(df_matched$year_id_alt - df_matched$year_id_ref) <= 5,]

# abbreviate super region
match5yr$super_region_name <- case_when(match5yr$super_region_name == "Central Europe, Eastern Europe, and Central Asia" ~ "CEurCAsia",
                                        match5yr$super_region_name == "High-income" ~ "HI",
                                        match5yr$super_region_name == "Latin America and Caribbean" ~ "LatAm",
                                        match5yr$super_region_name == "North Africa and Middle East" ~ "NAfrME",
                                        match5yr$super_region_name == "South Asia" ~ "SAsia",
                                        match5yr$super_region_name == "Southeast Asia, East Asia, and Oceania" ~ "SEAsiaOce",
                                        match5yr$super_region_name == "Sub-Saharan Africa" ~ "SSAfr")

# assign bug group and drug group for "supercombos"
match5yr$bug_group <- case_when(match5yr$pathogen %in% c('enterococcus_faecalis', 'enterococcus_faecium', 
                                                         'enterococcus_spp', 'staphylococcus_aureus',
                                                         'group_b_strep', 'group_a_strep', 'streptococcus_pneumoniae') ~ 'strep_group',
                                match5yr$pathogen %in% c('citrobacter_spp', 'enterobacter_spp', 'haemophilus_influenzae', 
                                                         'klebsiella_pneumoniae', 'proteus_spp', 'serratia_spp','morganella_spp') ~ 'kpneu_group',
                                match5yr$pathogen %in% c("escherichia_coli",'salmonella_typhi','salmonella_paratyphi',
                                                       'shigella_spp','non_typhoidal_salmonellae','neisseria_gonorrhoeae') ~ 'ecoli_group',
                                match5yr$pathogen %in% c("pseudomonas_aeruginosa", "acinetobacter_baumannii") ~ 'acinetobacter_group')


match5yr$drug_group <- case_when(match5yr$abx_class %in% c('anti_pseudomonal_penicillin', 'beta_lactamase_inhibitor',
                                                           'carbapenem', 'fourth_gen_ceph', 'third_gen_ceph', 'aminopenicillin',
                                                           'methicillin', 'penicillin') ~ 'beta_lactam',
                                 TRUE ~ match5yr$abx_class)

# assign regular combo for RE
match5yr$combo <- paste0(match5yr$pathogen, ' - ', match5yr$abx_class)

##############
# RUN MODELS
##############
# global model by supercombo
glob_mod_combos <- match5yr %>% group_by(bug_group, drug_group) %>%
  summarize(
    n = length(logit_diff),
    nspr = length(unique(super_region_name))
  )

glob_mod_combos$mean_eff <- NA

for (i in c(1:4,7:nrow(glob_mod_combos))){
  bug = glob_mod_combos$bug_group[i]
  drug = glob_mod_combos$drug_group[i]
  data <- match5yr[match5yr$bug_group == bug & match5yr$drug_group == drug,]
  
  mrdat <- mr$MRData()
  mrdat$load_df(
    data = data,  col_obs = "logit_diff", col_obs_se = "logit_diff_se",
    col_study_id = "combo")
  
  combomod <- mr$MRBRT(
    data = mrdat,
    cov_models = list(
      LinearCovModel("intercept", use_re = TRUE, prior_beta_uniform = array(c(0,Inf)))))
  
  combomod$fit_model(inner_print_level = 5L, inner_max_iter = 10000L)
  
  glob_mod_combos$mean_eff[i] <- combomod$summary()[[1]][[1]]
  
  py_save_object(object = combomod, filename = paste0(modelpath, "global/", bug, ":", drug, ".pkl"), pickle = "dill")
}

# super-region models by supercombo
spr_mod_combos <- match5yr %>% group_by(bug_group, drug_group, super_region_name) %>%
  summarize(
    n = length(logit_diff)
  )

spr_mod_combos$mean_eff <- NA

for (i in c(1:34,40:41,43:nrow(spr_mod_combos))){
  print(i)
  bug = spr_mod_combos$bug_group[i]
  drug = spr_mod_combos$drug_group[i]
  spr = spr_mod_combos$super_region_name[i]
  data <- match5yr[match5yr$bug_group == bug & match5yr$drug_group == drug & match5yr$super_region_name == spr,]
  
  mrdat <- mr$MRData()
  mrdat$load_df(
    data = data,  col_obs = "logit_diff", col_obs_se = "logit_diff_se",
    col_study_id = "combo")
  
  combomod <- mr$MRBRT(
    data = mrdat,
    cov_models = list(
      LinearCovModel("intercept", use_re = TRUE, prior_beta_uniform = array(c(0,Inf)))))
  
  combomod$fit_model(inner_print_level = 5L, inner_max_iter = 10000L)
  
  spr_mod_combos$mean_eff[i] <- combomod$summary()[[1]][[1]]
  
  if(!dir.exists(paste0(modelpath, spr, '/'))){
    dir.create(paste0(modelpath, spr, '/'), recursive = TRUE)
  }
  py_save_object(object = combomod, filename = paste0(modelpath, spr, "/", bug, ":", drug, ".pkl"), pickle = "dill")
}

##############
# GENERATE PREDICTIONS
##############
# estimate proportion of data tertiary
resdat <- resdat %>% left_join(location_md2[,c('location_id','location_name')])

tertprop <- resdat %>% group_by(super_region_name) %>%
  summarize(
    prop_n_tertiary = sum(cases[hospital_type == 'tertiary'])/sum(cases),
    prop_dp_tertiary = length(cases[hospital_type == 'tertiary' & cases >= 5])/length(cases[cases >= 5])
  )
tertprop <- tertprop[!is.na(tertprop$super_region_name),]
write.csv(tertprop, paste0(modelpath,'tertiary_prop_by_SR.csv'), row.names = F)

# abbreviate super region
resdat$super_region_name <- case_when(resdat$super_region_name == "Central Europe, Eastern Europe, and Central Asia" ~ "CEurCAsia",
                                      resdat$super_region_name == "High-income" ~ "HI",
                                      resdat$super_region_name == "Latin America and Caribbean" ~ "LatAm",
                                      resdat$super_region_name == "North Africa and Middle East" ~ "NAfrME",
                                      resdat$super_region_name == "South Asia" ~ "SAsia",
                                      resdat$super_region_name == "Southeast Asia, East Asia, and Oceania" ~ "SEAsiaOce",
                                      resdat$super_region_name == "Sub-Saharan Africa" ~ "SSAfr")

# calculate resistance, standard error, and apply cutoff
resdat$prev <- resdat$resistant/resdat$cases
resdat$COprev <- resdat$prev
resdat$COprev[resdat$prev == 0] <- 0.001
resdat$COprev[resdat$prev == 1] <- 0.999
resdat$se <- sqrt((resdat$COprev*(1-resdat$COprev))/resdat$cases)

# remove handful of rows that got pulled that have 0 cases
resdat <- resdat[!is.na(resdat$COprev),]

# calculate logit resistance
resdat <- as.data.table(resdat) %>% 
  cbind(delta_transform(resdat$COprev, resdat$se, "linear_to_logit"))

# assign hospital type, and supercombo bug and drug group
resdat$hospital_type <- ifelse(resdat$hospital_type == 'tertiary', 'tertiary', 'nontert-mixed')

resdat$bug_group <- case_when(resdat$pathogen %in% c('enterococcus_faecalis', 'enterococcus_faecium', 
                                                     'enterococcus_spp', 'staphylococcus_aureus',
                                                     'group_b_strep', 'group_a_strep', 'streptococcus_pneumoniae') ~ 'strep_group',
                              resdat$pathogen %in% c('citrobacter_spp', 'enterobacter_spp','haemophilus_influenzae', 
                                                     'klebsiella_pneumoniae', 'proteus_spp', 'serratia_spp') ~ 'kpneu_group',
                              resdat$pathogen == "acinetobacter_baumannii" ~ 'acinetobacter_group',
                              resdat$pathogen == "escherichia_coli" ~ 'ecoli_group')


resdat$drug_group <- case_when(resdat$abx_class %in% c('anti_pseudomonal_penicillin', 'beta_lactamase_inhibitor',
                                                       'carbapenem', 'fourth_gen_ceph', 'third_gen_ceph', 'aminopenicillin',
                                                       'methicillin', 'penicillin') ~ 'beta_lactam',
                               TRUE ~ resdat$abx_class)

# create an index to facilitate recombining later
resdat$index <- c(1:nrow(resdat))

#### Global Adjustment ####
globmods <- gsub('.pkl', '', list.files(paste0(modelpath, 'global/'), pattern = ':'))
if(pulliord){
  globmods <- 'acinetobacter_group:beta_lactam'
}
tertadjglob <- as.data.frame(c())

for (mod in globmods){
  bug = strsplit(mod, ':')[[1]][1]
  drug = strsplit(mod, ':')[[1]][2]
  
  result = py_load_object(filename = paste0(modelpath, 'global/', mod, ".pkl"), pickle = "dill")
  
  tertdata = resdat[resdat$hospital_type == 'tertiary' & resdat$bug_group == bug & resdat$drug_group == drug,]
  if(pulliord){
    tertdata <- resdat
  }
  tertpred <- mr$MRData()
  
  # have to feed a numeric covariate to get it to load the data, prev just a placeholder
  tertpred$load_df(
    data = tertdata,
    col_covs=list('prev')
  )
  
  # predict adjustment
  tertdata$adj <- result$predict(data = tertpred)
  # subtract adjustment from cutoff resistance
  tertdata$adj_logit <- tertdata$mean_logit - tertdata$adj
  # convert to prevalence
  tertdata$pred_prev <- logit2prob(tertdata$adj_logit)
  
  tertadjglob <- bind_rows(tertadjglob, tertdata)
}

adjdat <- bind_rows(resdat[resdat$hospital_type != 'tertiary',], tertadjglob)

adjdat$pred_prev[adjdat$hospital_type != 'tertiary'] <- adjdat$COprev[adjdat$hospital_type != 'tertiary']

#### Super Region Adjustment ####
spr_mod_combos_datarich <- spr_mod_combos[spr_mod_combos$n >= 250,]
if(pulliord){
  bug_group = 'acinetobacter_group'
  drug_group = 'beta_lactam'
  super_region_name = 'HI'
  spr_mod_combos_datarich <- data.frame(bug_group, drug_group, super_region_name)
}
tertadjspr <- as.data.frame(c())

for (i in 1:nrow(spr_mod_combos_datarich)){
  bug = spr_mod_combos_datarich$bug_group[i]
  drug = spr_mod_combos_datarich$drug_group[i]
  spr = spr_mod_combos_datarich$super_region_name[i]
  
  result = py_load_object(filename = paste0(modelpath, spr, '/', bug, ':', drug, ".pkl"), pickle = "dill")
  
  tertdata = resdat[resdat$hospital_type == 'tertiary' & resdat$bug_group == bug & resdat$drug_group == drug & resdat$super_region_name == spr,]
  if(length(tertdata$location_id) < 1){
    tertdata = resdat[resdat$bug_group == bug & resdat$drug_group == drug & resdat$super_region_name == spr,]
  }
  if(pulliord){
    tertdata <- resdat
  }
  
  # have to feed a numeric covariate to get it to load the data, prev just a placeholder
  tertpred <- mr$MRData()
  tertpred$load_df(
    data = tertdata,
    col_covs=list('prev')
  )
  
  tertdata$adj <- result$predict(data = tertpred)
  tertdata$adj_logit <- tertdata$mean_logit - tertdata$adj
  tertdata$pred_prev <- logit2prob(tertdata$adj_logit)
  
  tertadjspr <- bind_rows(tertadjspr, tertdata)
}

# combine with global adjustment
adjdat <- merge(adjdat, tertadjspr[,c('pathogen', 'abx_class', 'year_id', 'source', 'nid','hospital_type', 'super_region_name', 'cases', 'resistant', 'susceptible', 'index', 'pred_prev'),],
                  by = c('pathogen', 'abx_class', 'year_id', 'source','nid' ,'hospital_type', 'super_region_name', 'cases', 'resistant', 'susceptible', 'index'),
                  suffixes = c('_glob', '_spr'), all.x = TRUE)

# substitute global prevalence estimate wherever super region is unavailable (24 supercombo/superregions)
adjdat$pred_prev_spr[is.na(adjdat$pred_prev_spr)] <- adjdat$pred_prev_glob[is.na(adjdat$pred_prev_spr)]


#### Data ####

adjdat$resistant[adjdat$resistant != 0 & adjdat$hospital_type == 'tertiary'] <-
  adjdat$pred_prev_spr[adjdat$resistant != 0 & adjdat$hospital_type == 'tertiary'] *
  adjdat$cases[adjdat$resistant != 0 & adjdat$hospital_type == 'tertiary']

adjdat$susceptible <- adjdat$cases - adjdat$resistant

keepcols <- c('source', 'location_id', 'year_id', 'pathogen', 'abx_class', 'hospital_type', 'cases',
                  'resistant', 'susceptible','nid','region_id','super_region_id')

adjdat <- adjdat[, ..keepcols]

aggadj <- adjdat %>% dplyr::group_by(source, location_id, year_id, pathogen, abx_class,nid,region_id,super_region_id) %>%
  dplyr::summarize(
    cases = sum(cases),
    resistant = sum(resistant),
    susceptible = sum(susceptible),
    nsub_dp = length(hospital_type)
  )

aggadj <- aggadj[, !colnames(aggadj) %in% 'nsub_dp']
aggadj$prev <- aggadj$resistant/aggadj$cases
aggadj$prev[aggadj$prev == 0] <- 0.001
aggadj$prev[aggadj$prev == 1] <- 0.999
# calculate SEs
aggadj$se <- sqrt((aggadj$prev*(1-aggadj$prev))/aggadj$cases)
aggadj$val <- aggadj$resistant/aggadj$cases
aggadj$val[aggadj$val == 0] <- 0.001
aggadj$val[aggadj$val == 1] <- 0.999

if(pulliord){
  write.csv(aggadj, 'FILEPATH', row.names = F)
} else {
  write.csv(aggadj, 'FILEPATH', row.names = F)
}

##############
# PLOT RESULTS
##############
adjdat <- adjdat%>% left_join(location_md2[,c('location_id','super_region_name')])
# summarize the change in crude resistance by super region/combo
resdiff <- adjdat %>% group_by(pathogen, abx_class, super_region_name) %>%
  summarize(
    crude_res_base = sum(COprev*cases)/sum(cases),
    crude_res_glob_adj = sum(pred_prev_glob*cases)/sum(cases),
    crude_res_spr_adj = sum(pred_prev_spr*cases)/sum(cases)
  )

resdiff$base_glob_diff <- resdiff$crude_res_glob_adj - resdiff$crude_res_base
resdiff$base_spr_diff <- resdiff$crude_res_spr_adj - resdiff$crude_res_base
resdiff$glob_spr_diff <- resdiff$crude_res_glob_adj - resdiff$crude_res_spr_adj

resdiff$combo <- paste0(resdiff$pathogen, ':', resdiff$abx_class)

resdiff$base_glob_diff[resdiff$base_glob_diff > 0] <- 0
resdiff$base_spr_diff[resdiff$base_spr_diff > 0] <- 0

min <- round(min(c(resdiff$base_glob_diff, resdiff$base_spr_diff)), digits = 2)

resdiff$base_glob_diff[resdiff$base_glob_diff == 0] <- NA
resdiff$base_spr_diff[resdiff$base_spr_diff == 0] <- NA
resdiff$glob_spr_diff[resdiff$glob_spr_diff == 0] <- NA

ggplot(resdiff) +
  geom_tile(aes(x = super_region_name, y = combo, fill = base_glob_diff)) +
  scale_y_discrete(limits = rev) +
  scale_fill_viridis(option = 'turbo', na.value = '#ffffff')+
  #  scale_fill_gradient(low = '#4d004b', high = '#e0ecf4', na.value = '#ffffff', limits = c(min, 0)) +
  labs(x = 'Super Region', y = 'Combo', title = 'Change in Crude Resistance: Global Models', fill = 'Absolute change')
ggsave('FILEPATH', height = 9, width = 16)

ggplot(resdiff) +
  geom_tile(aes(x = super_region_name, y = combo, fill = base_spr_diff)) +
  scale_y_discrete(limits = rev) +
  scale_fill_gradient(low = '#4d004b', high = '#e0ecf4', na.value = '#ffffff', limits = c(-0.14, 0)) +
  labs(x = 'Super Region', y = 'Combo', title = 'Change in Crude Resistance: Super Region Stratified', fill = 'Absolute change')
ggsave('FILEPATH', height = 9, width = 16)

ggplot(resdiff) +
  geom_tile(aes(x = super_region_name, y = combo, fill = glob_spr_diff)) +
  scale_y_discrete(limits = rev) +
  scale_fill_gradient2(low = 'red', high = 'blue', mid = 'white', na.value = '#ffffff') +
  labs(x = 'Super Region', y = 'Combo', title = 'Crude Resistance: Global - Super Region', fill = 'Absolute difference')
ggsave('FILEPATH', height = 9, width = 16)


# illustrate adjusted data points by supercombo
ggplot(adjdat[adjdat$bug_group == 'acinetobacter_group' & adjdat$hospital_type == 'tertiary',]) +
  geom_abline(aes(slope = 1, intercept = 0), alpha = 0.7) +
  geom_point(aes(x = COprev, y = pred_prev_spr, color = super_region_name), size = 2) +
  theme_bw() +
  labs(x = 'Original Prevalence', y = 'Adjusted Prevalance', color = 'Super Region', title = 'Pseudomonadales Adjustment') +
  facet_wrap(~drug_group)

longadj <- melt(adjdat, id.vars = c('pathogen', 'abx_class', 'bug_group', 'drug_group', 'super_region_name', 'hospital_type'),
                measure.vars = c('COprev', 'pred_prev_glob', 'pred_prev_spr'), variable.name = 'estimate',
                value.name = 'prevalence')
longadj$estimate <- case_when(longadj$estimate == 'COprev' ~ 'Unadjusted',
                              longadj$estimate == 'pred_prev_glob' ~ 'Global Adjustment',
                              longadj$estimate == 'pred_prev_spr' ~ 'SupReg Hybrid Adjustment')

longadj$estimate <- factor(longadj$estimate, levels = c('Unadjusted', 'Global Adjustment', 'SupReg Hybrid Adjustment'))

for (bug in c('acinetobacter_group', 'ecoli_group', 'strep_group')) {
  plottitle = case_when(bug == 'acinetobacter_group' ~ 'Pseudomonadales Adjustment',
                    bug == 'ecoli_group' ~ 'Enterobacterales Adjustment',
                    bug == 'strep_group' ~ 'Gram-Positives Adjustment')
  
  ggplot(longadj[longadj$bug_group == bug & longadj$hospital_type == 'tertiary',]) +
    geom_boxplot(aes(x = drug_group, y = prevalence, fill = estimate), outlier.alpha = 0.2) +
    theme_bw() +
    labs(x = 'Antibiotic Super Group', y = 'Prevalence of Resistance', fill = 'Crude Resistance', title = plottitle) +
    facet_wrap(~super_region_name)
  ggsave(paste0('FILEPATH', bug, '240713.png'), height = 9, width = 16)
}