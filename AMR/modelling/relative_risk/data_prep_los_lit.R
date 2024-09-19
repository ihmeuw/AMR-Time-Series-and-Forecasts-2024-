
username <- Sys.getenv("USER")
amr_repo <- "FILEPATH"
library(tidyverse)
library(metafor)

print("Loading arguments")
args <- commandArgs(trailingOnly = T)
description <- as.character(args[1])

outpath <- "FILEPATH"
suffixdt <- "FILEPATH"


nidmeta <- read.csv("FILEPATH") %>% dplyr::select(source,file_path,active_5)


read_lit_data <- function(){
  lospath <- "FILEPATH"
  litrev <- as.data.frame(readxl::read_xlsx("FILEPATH")) 
  litrev$source <- litrev$author


  return(litrev)
}

correct_mean_iqr <- function(df){


  df$sd_r[(is.na(df$sd_r) & !is.na(df$half_IQR_r))] <- df$half_IQR_r[(is.na(df$sd_r) & !is.na(df$half_IQR_r))]
  df$sd_s[(is.na(df$sd_s) & !is.na(df$half_IQR_s))] <- df$half_IQR_s[(is.na(df$sd_s) & !is.na(df$half_IQR_s))]
  df$mean_LOS_r[(is.na(df$mean_LOS_r) & !is.na(df$median_LOS_r))] <- df$median_LOS_r[(is.na(df$mean_LOS_r) & !is.na(df$median_LOS_r))]
  df$mean_LOS_s[(is.na(df$mean_LOS_s) & !is.na(df$median_LOS_s))] <- df$median_LOS_s[(is.na(df$mean_LOS_s) & !is.na(df$median_LOS_s))]

  return(df)
}

clean_lit_syndrome <- function(df){
  df <- df %>%
    dplyr::mutate(infectious_syndrome = case_when(
      syndrome %in% c("Pneumonia", "LRI", "Ventilator-associated pneumonia", "CAPneumonia", "pneumonia") ~ "L2_lower_respiratory_infection",
      syndrome %in% c("blood", "BSI" ,"Neutropenia") ~ "L2_blood_stream_infection",
      syndrome %in% c("UTI") ~ "L2_urinary_tract_infection",
      syndrome %in% c("Prosthetic joint infection", "Septic arthritis", "tb") ~ "other_syndromes",
      TRUE ~ "L2_unspecified_infection"
    ))
  return(df)
}


df <- read_lit_data()
df <- correct_mean_iqr(df)

df <- df[!is.na(df$sd_r),]
df <- clean_lit_syndrome(df)

df$pathogen[df$pathogen == "acinetobacter_baumanii"] <- 'acinetobacter_baumannii'
df$pathogen[df$pathogen == "neisseria_gonorrheae"] <- 'neisseria_gonorrhoeae'

df <- df[,c('pathogen','abx_class','infectious_syndrome','source','cases_s','cases_r','mean_LOS_r',
                    'sd_r','mean_LOS_s','sd_s','sample_size')]

write.csv(df,"FILEPATh", row.names = F)

