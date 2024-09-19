library(dplyr)
library(data.table)
library(tidyverse)

main_path <- "FILEPATH"

all_syn <- "FILEPATH"
syn <- "FILEPATH"
moh_map <- read_csv("FILEPATH")

out_path <- "FILEPATH"

single_syn_keep <- moh_map %>% 
  filter(abx_class != "aminopenicillin") %>%
  mutate(combo = paste(pathogen, abx_class, sep = "-"))

syndrome_keep <- moh_map %>% 
  filter(abx_class %in% c("aminopenicillin")) %>%
  mutate(combo = paste(pathogen, abx_class, sep = "-"))

all_syn_sum <- read_csv("FILEPATH")
all_syn_sum <- all_syn_sum %>%
  mutate(combo = paste(pathogen, abx_class, sep = "-")) %>% 
  filter(combo %in% unique(single_syn_keep$combo))

syn_sum <- read_csv("FILEPATH")
syn_sum <- syn_sum %>%
  mutate(combo = paste(pathogen, abx_class, sep = "-")) %>% 
  filter(combo %in% unique(syndrome_keep$combo))

all_syn_draw <- read_csv("FILEPATH")
all_syn_draw <- all_syn_draw %>%
  mutate(combo = paste(pathogen, abx_class, sep = "-")) %>% 
  filter(combo %in% unique(single_syn_keep$combo))

syn_draw <- read_csv("FILEPATH")
syn_draw <- syn_draw %>%
  mutate(combo = paste(pathogen, abx_class, sep = "-")) %>% 
  filter(combo %in% unique(syndrome_keep$combo))


syndromes <- c("L2_blood_stream_infection", "L2_urinary_tract_infection", "L2_lower_respiratory_infection", "other_syndromes")

expanded_allsyn_sum <- all_syn_sum %>% 
  crossing(infectious_syndrome = syndromes)

expanded_allsyn_draw <- all_syn_draw %>% 
  crossing(infectious_syndrome = syndromes)

summary <- rbind(expanded_allsyn_sum, syn_sum)
draws <- rbind(expanded_allsyn_draw, syn_draw)

write.csv(summary, "FILEPATH", row.names=FALSE)
write.csv(draws, "FILEPATH", row.names=FALSE)


secondary <- read_csv("FILEPATH")
secondary <- secondary %>% 
  mutate(combo = paste(pathogen, abx_class, sep = "-")) %>%
  filter(!(combo %in% c("escherichia_coli-sulfa", "acinetobacter_baumannii-anti_pseudomonal_penicillin", 
                        "acinetobacter_baumannii-fourth_gen_ceph", "acinetobacter_baumannii-aminoglycoside")))

sulfa_1 <- read_csv("FILEPATH")
ap_1 <- read_csv("FILEPATH")
fgc_1 <- read_csv("FILEPATH")
ag_1 <- read_csv("FILEPATH")

spot_fix <- function(df, path){
  rep <- df %>%
    filter(pathogen == path)
  if (path == "enterobacter_spp"){
    rep$pathogen <- "escherichia_coli"
  }
  return(rep)
}

sulfa_2 <- spot_fix(sulfa_1, "enterobacter_spp")
ap_2 <- spot_fix(ap_1, "acinetobacter_baumannii")
fgc_2 <- spot_fix(fgc_1, "acinetobacter_baumannii")
ag_2 <- spot_fix(ag_1, "acinetobacter_baumannii")

secondary <- secondary[ , !(names(secondary) %in% c("combo", "strong_prior"))]
names(secondary) = gsub(pattern = "X*", replacement = "", x = names(secondary))

draws_fin <- rbind(secondary, sulfa_2, ap_2, fgc_2, ag_2)
write.csv(draws_fin, "FILEPATH", row.names=FALSE)
