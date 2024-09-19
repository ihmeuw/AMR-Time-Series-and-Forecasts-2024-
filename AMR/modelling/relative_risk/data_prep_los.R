
username <- Sys.getenv("USER")
amr_repo <- "FILEPATH"
source("FILEPATH")
library(tidyverse)
library(lmtest)
library(sandwich)
library(assertthat)
library(lme4)
library(lubridate)
source("FILEPATH")
source("FILEPATH")
source("FILEPATH")
source("FILEPATH")
source("FILEPATH")
source("FILEPATH")

print("Loading arguments")
args <- commandArgs(trailingOnly = T)
source <- as.character(args[1])
description <- as.character(args[2])
burden <- as.integer(args[3])

if (burden == "fatal"){
  val = "deaths"
  out_dir = "coefficients"
  suffix = "coeff"
}else{
  val = "days_infection"
  out_dir = "los"
  suffix = "los"
}

outpath <- "FILEPATH"


infsyn_hierarchy <- read_csv("FILEPATH")
l2 <- infsyn_hierarchy %>% 
  filter(level == 2)
peri <- infsyn_hierarchy %>%
  filter(infectious_syndrome == "L1_peritoneal_and_intra_abdomen_infection")
card <- infsyn_hierarchy %>%
  filter(infectious_syndrome == "L1_myocarditis_pericarditis_carditis")
bone <- infsyn_hierarchy %>%
  filter(infectious_syndrome == "L1_bone_joint_infection")
l2 <- l2 %>%
  filter(!(infectious_syndrome %in% c("L2_myocarditis", "L2_pericarditis", "L2_carditis",
                                      "L2_bone_joint_bacterial_inf", "L2_bone_joint_non_bacterial_inf",
                                      "L2_peritonitis", "L2_gastrointestinal_infection_other",
                                      "L2_intra_abdomen_abscess")))
l2 <- rbind(l2, peri, card, bone)

available_parents <- l2$infectious_syndrome

organism <- c('acinetobacter_baumannii','citrobacter_spp','enterobacter_spp','enterococcus_faecalis','enterococcus_faecium','enterococcus_spp','escherichia_coli','group_a_strep','group_b_strep','haemophilus_influenzae','klebsiella_pneumoniae','moraxella_spp','morganella_spp','neisseria_meningitidis','salmonella_ints_non_typhi/paratyphi','proteus_spp','providencia_spp','pseudomonas_aeruginosa','pseudomonas_spp','neisseria_gonorrhoeae','shigella_spp','mycobacterium_tuberculosis','salmonella_paratyphi','salmonella_typhi','salmonella_typhi_paratyphi','serratia_spp','staphylococcus_aureus','streptococcus_pneumoniae','mycoplasma','listeria','legionella_spp')
antibiotic_group <- c('third_gen_ceph','carbapenem','fluoroquinolone','penicillin','aminopenicillin','beta_lactamase_inhibitor','anti_pseudomonal_penicillin','methicillin','vancomycin','fourth_gen_ceph','sulfa','aminoglycoside','macrolide')

level_amr_map <- get_cause_metadata(cause_set_id = 4, release_id = 10) %>% 
  dplyr::select(acause, cause_id, level, parent_id, acause_parent) %>% distinct()


clean_cause <- function(df, has_ucod, level_amr_map){
  if (has_ucod==1){
    df$cause_id[df$cause_id == 919] <- 743 
    
    df <- df %>% left_join(level_amr_map[level_amr_map$level %in% c(4,3),c('parent_id','cause_id')], by = 'cause_id') %>% dplyr::rename(levelup = parent_id)
    df <- df %>% left_join(level_amr_map[level_amr_map$level %in% c(3),c('parent_id','cause_id')], by = c("levelup" = "cause_id")) %>% rename(level2 = parent_id)
    df$level2[is.na(df$level2) & !is.na(df$levelup)] <- df$levelup[is.na(df$level2) & !is.na(df$levelup)]
    df$level2[is.na(df$level2) & !is.na(df$cause_id)] <- df$cause_id[is.na(df$level2) & !is.na(df$cause_id)]
    df <- df %>% left_join(level_amr_map[,c('level','cause_id','acause')], by = c("level2" = "cause_id")) %>% rename(level2_amr = acause)
    df$level2_amr[df$cause_id==743] <- '_gc'
    df[,c('level','levelup','level3','level2')] <- NULL
  }else{
    df$level2_amr  <- '_gc'
  }
  
  return(df)
}

clean_resistance <- function(df){
  df$resistant[df$resistance == 'resistant'] <- 1L
  df$resistant[df$resistance %in% c('susceptible', 'sensitive')] <- 0L
  
  if("hosp" %in% names(df)) {
    df$origin[(df$hosp == 'hospital')] <- 1L
    df$origin[!(df$hosp == 'hospital')] <- 0L
  }else{
    df$origin <- 0
  }
  
  return(df)
}

clean_syndrome <- function(df){
  df <- assign_l2_syndrome(df, available_parents, infsyn_hierarchy)
    df <- df %>%
      mutate(infectious_syndrome = replace_na(infectious_syndrome, "L2_unspecified")) %>%
      mutate(infectious_syndrome = case_when(
        infectious_syndrome %in% c("L2_OAMOIS", "L2_unspecified") ~ "L2_unspecified",
        TRUE ~ infectious_syndrome
        )) %>%
      mutate(infs = case_when(
        infectious_syndrome != "L2_blood_stream_infection" & 
          infectious_syndrome != "L2_lower_respiratory_infection" & 
          infectious_syndrome != "L2_unspecified" &
          infectious_syndrome != "L2_urinary_tract_infection" ~ "other_syndromes",
        TRUE ~ infectious_syndrome
      ))

  return(df)
}

clean_time <- function(df){
  df$days[df$days < 0] <- 0L
  df$days_prior[df$days_prior < 0] <- 0L
}


print(paste0("Working on ", source))

df <- get_amr_data(sources = source, redistribute_pathogens=TRUE)

if (source %in% c("SOURCES")){
  df$sample_id <- runif(n = nrow(df), min = 1, max = nrow(df))
}

df <- clean_resistance(df)
df <- clean_syndrome(df)
if (burden == "nonfatal"){
  df <- clean_time(df)
}
  df <- df %>% dplyr::filter(!is.na(!!sym(val)) & !is.na(abx_class) & !is.na(pathogen) & !is.na(resistant)) 

  df <- df %>%  dplyr::group_by(sample_id,abx_class,pathogen,source) %>% 
    dplyr::summarize(resistant = max(resistant),deaths = max(!!sym(val)), cases=max(cases))

write.csv("FILEPATH", row.names = F)