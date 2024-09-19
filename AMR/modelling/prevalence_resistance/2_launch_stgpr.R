### Launching ST-GPR
rm(list = ls())
central_root <- 'FILEPATH'
setwd(central_root)
source('FILEPATH')
library(dplyr)

# subset input data to combos of interest
username <- Sys.getenv("USER")
hdrive <- sprintf("FILEPATH", username)
combos <- read.csv(paste0(hdrive,'amr/maps/bug_drug_combos.csv'))


path_config <- 'FILEPATH'
read_config <- 'FILEPATH'
config2  <- read.csv(paste0(read_config,'FILEPATH'), stringsAsFactors = F)
config <- rbind(config)
config$success <- 1L
config$gpr_amp_factor <- 1L
config <- config[(config$run_id %in% unique(combos$resistance_run_id)),]
length(unique(config$description))
length(config$description)
length(unique(combos$resistance_run_id))

config$path_to_data <- paste0('FILEPATH',
                              substring(config$description,1,config$notes-1),"/",
                              substring(config$description,config$notes+1,),"/input.csv")
config$path_to_custom_stage_1 <-paste0('FILEPATH',
                                       substring(config$description,1,config$notes-1),"/",
                                       substring(config$description,config$notes+1,),"/custom_stage1_df.csv")
config$run_id <- NULL
config$success <- NULL

finalconfig<-c()
for(i in 2:nrow(config)) {
  temp_config <- config[i,]
  write.csv(temp_config,paste0(path_config,'temp_config.csv'),row.names=F) 
  run_id <- register_stgpr_model(path_to_config = paste0(path_config,'temp_config.csv'))
  stgpr_sendoff(run_id, 
              log_path=paste0('FILEPATH'))
  temp_config$run_id <-run_id 
  finalconfig<-rbind(finalconfig,temp_config)
}

finalconfig2 <- as.data.frame(finalconfig)
write.csv(finalconfig2, paste0(path_config,'FILEPATH'), row.names = FALSE)

#To check each run_id
for(i in unique(finalconfig2$run_id)) {
  print(i)
  check_run(i)
}

#To check each run_id
for(i in 1:nrow(finalconfig2)) {
  finalconfig2$success[i] <- check_run(finalconfig2$run_id[i])
}
write.csv(finalconfig2, paste0(path_config,'FILEPATH'), row.names = FALSE)

reruns <- finalconfig2[finalconfig2$success != 1, !colnames(finalconfig2) %in% "success"]

rerunned <- c()
for(i in 1:nrow(reruns)) {
  temp_config <- reruns[i,]
  y <- substring(temp_config$description,1,stringr::str_locate(temp_config$description,"-")-1)[1]
  z <- substring(temp_config$description,stringr::str_locate(temp_config$description,"-")+1,)[1]
  data <-  read.csv(temp_config$path_to_data, stringsAsFactors = F)
  #temp_config$year_start <-min(data$year_id[data$is_outlier==0])
  rm(data,y,z)
  write.csv(temp_config,paste0(path_config,'temp_config.csv'),row.names=F) 
  run_id <- register_stgpr_model(path_to_config = paste0(path_config,'temp_config.csv'))
  stgpr_sendoff(run_id, 
                log_path=paste0('FILEPATH'))
  temp_config$run_id <-run_id 
  rerunned<-rbind(rerunned,temp_config)
}

#To check each run_id
for(i in 1:nrow(rerunned)) {
  rerunned$success[i] <- check_run(rerunned$run_id[i])
}

finalconfig2 <- finalconfig2[finalconfig2$success == 1,]
finalconfig2 <- rbind(finalconfig2, rerunned)

write.csv(finalconfig2,paste0(path_config, 'FILEPATH'), row.names = FALSE)

failures <- finalconfig[finalconfig$success == 0,]


config <- read.csv(paste0(path_config,'FILEPATH'))
config <- config[config$description != 'enterococcus_spp-fluoroquinolone',]

config <- bind_rows(config, finalconfig2)
write.csv(config, paste0(path_config,'FILEPATH'), row.names = FALSE)


# manually fix some 0's and 1's in Enterococcus_spp fluoroquinolone
for (i in 1:nrow(reruns)){
  dat <- read.csv(reruns$path_to_custom_stage_1[i])
  dat$cv_custom_stage_1[dat$cv_custom_stage_1 == 0] <- 0.00001
  dat$cv_custom_stage_1[dat$cv_custom_stage_1 == 1] <- 0.99999
  write.csv(dat, reruns$path_to_custom_stage_1[i], row.names = F)
}
