
get_adjusted_rr <- function(dataor = dataor, sources = sources){
  dataor = dataor
  org_adjusted_dataor <- organism[organism %in% unique(dataor$pathogen)]
  antibiotic_group_dataor <-antibiotic_group[antibiotic_group %in% unique(dataor$abx_class)]
  yearseries <- unique(dataor$year)
  agegroups <- unique(dataor$group_id)
  locseries <- unique(dataor$loc_id)
  incomeseries <- unique(dataor$high_income)
  infseries <- unique(dataor$infs)
  
  
  adjusted_coeff <-c()
  for(l in locseries) {
    for (i in org_adjusted_dataor) {
      for (j in antibiotic_group_dataor) {
        for (y in yearseries) {
          for (a in agegroups) {
            for (n in infseries){
              for (h in incomeseries){
                text_for_data <- paste0("dataor[dataor$pathogen == '",i,
                                        "' & dataor$abx_class == '",j,
                                        "' & dataor$loc_id==",l,
                                        " & dataor$year=='",y,
                                        "' & dataor$group_id=='",a,
                                        "' & dataor$infs=='",n, 
                                        "' & dataor$high_income=='",h,
                                        "' & !is.na(dataor$resistant) & !is.na(dataor$deaths),]")
                dataslice <- eval(parse(text = paste(text_for_data)))
                if(length(dataslice$sample_id) > 1 & (length(unique(dataslice$resistant))> 1)) {
                  total <- dataslice %>% dplyr::group_by(resistant) %>% dplyr::summarize(admit = n(), died = sum(deaths,na.rm = TRUE))
                  rr <- (total[2,3]/total[2,2])/(total[1,3]/total[1,2])
                  lnrr_var <- (1/total[2,3]) - (1/total[2,2]) + (1/total[1,3]) - (1/total[1,2])
                  lnrr_se <- sqrt(lnrr_var)
                  if(!is.na(total[1,3]) & !is.na(total[2,3]) & (total[1,3] > 0) & (total[2,3] > 0) & (total[1,3] < total[1,2]) & (total[2,3] < total[2,2])) {
                    r2 <- c(source, i, j, y, a, l, n, h, log(rr$died),lnrr_se$died,total$resistant[1],total$admit[1],total$died[1],total$resistant[2],total$admit[2],total$died[2])
                    if (length(unique(dataslice$level2_amr)) > 1) {
                      model2 <- glm(data = dataslice, formula = deaths ~ resistant + origin + level2_amr, family = 'poisson'(link='log'))
                      model2 <- coeftest(model2, vcov = vcovHC(model2,type="HC1"))
                      if((model2[2,1] > 0) & (model2[2,1] < 10) & (length(model2) > 8)){
                        r2 <- as.vector(c(r2,'adjusted',model2[2,1],model2[2,2], model2[3,1],model2[3,2]))
                      }
                      else {
                        if (length(unique(dataslice$origin)) > 1) {
                          model2 <- glm(data = dataslice, formula = deaths ~ resistant + origin , family = 'poisson'(link='log'))
                          if((model2$coefficients[2] > 0) & (model2$coefficients[2] < 10) & (length(unique(is.na(model2$coefficients))) == 1)){
                            model2 <- coeftest(model2, vcov = vcovHC(model2,type="HC1"))
                            r2 <- as.vector(c(r2,'hosp_adjusted',model2[2,1],model2[2,2], model2[3,1],model2[3,2]))
                          }
                        }
                      }
                    }
                    else {
                      if (length(unique(dataslice$origin)) > 1) {
                        model2 <- glm(data = dataslice, formula = deaths ~ resistant + origin , family = 'poisson'(link='log'))
                        if((model2$coefficients[2] > 0) & (model2$coefficients[2] < 10) & (length(unique(is.na(model2$coefficients))) == 1)){
                          model2 <- coeftest(model2, vcov = vcovHC(model2,type="HC1"))
                          r2 <- c(r2,'hosp_adjusted',model2[2,1],model2[2,2], model2[3,1],model2[3,2])
                        }
                        else {
                          r2 <- c(source, i, j, y, a, l, n, h, log(rr$died),lnrr_se$died,total$resistant[1],total$admit[1],total$died[1],total$resistant[2],total$admit[2],total$died[2],'unadjusted',log(rr$died),lnrr_se$died,log(rr$died),lnrr_se$died)
                        }
                        rm(list = objects()[grepl("model",objects())])
                      }
                    }
                  } else {
                    r2 <- c(source, i, j, y, a, l, n, h, log(rr$died),lnrr_se$died,total$resistant[1],total$admit[1],total$died[1],total$resistant[2],total$admit[2],total$died[2])
                  }
                  if(length(r2)<=17){r2 <- c(r2,'unadjusted',log(rr$died),lnrr_se$died,log(rr$died),lnrr_se$died)}
                  adjusted_coeff <- rbind(adjusted_coeff,r2)
                  rm(total,rr,lnrr_se,lnrr_var,r2)
                } 
                rm(dataslice,text_for_data)
              }
            }
          }
        }
      }
    }
  }
  
  adjusted_coeff <- as.data.frame(adjusted_coeff)
  colnames(adjusted_coeff) <- c('refvar','pathogen','abx_class','year_id','age_group','loc_id', 'infs', 'high_income', 'coeff_unadjusted','se_unadjusted',
                                'res1','admit_s','died_s','res2','admit_r','died_r','adjusted','coeff_adjusted','se_adjusted','coeff_hosp_adjusted','se_hosp_adjusted')
  adjusted_coeff$coeff_adjusted <- as.numeric(adjusted_coeff$coeff_adjusted)
  adjusted_coeff$se_adjusted <- as.numeric(adjusted_coeff$se_adjusted)
  adjusted_coeff$coeff_unadjusted <- as.numeric(adjusted_coeff$coeff_unadjusted)
  adjusted_coeff$se_unadjusted <- as.numeric(adjusted_coeff$se_unadjusted)
  adjusted_coeff$adjusted <- as.character(adjusted_coeff$adjusted)
  adjusted_coeff$source <- paste0(sources)
  return(adjusted_coeff)
}


under_over_5 <- function(dataor = dataor) {
  dataor = dataor
  dataor$group_id <- 'All ages'
  dataor$group_id[dataor$age_group_id %in% c(42,2,3,4,5,28,34,388,389,390,238,407,177,325)] <- 'under_5' 
  dataor$group_id[dataor$age_group_id %in% c(6,7,8,9,10,11,12,13,14,15,16,17,18,25,329,362,239,209,214,215,222,224,229,284,329,335,422,425)] <- '5 to 69'
  dataor$group_id[dataor$age_group_id %in% c(19,20,26,30,31,32,40,235)] <- '70 plus'
  dataor$group_id[dataor$age_group_id %in% c(322,323,324,385,39)] <- 'under_5'
  dataor$group_id[dataor$age_group_id %in% c(29,37,41,157,176,221,257,338)] <- '5 to 69'
  dataor$group_id[dataor$age_group_id %in% c(40,154)] <- '70 plus'
  return(dataor)
}

backtracking_pathogens <- function(dataor = dataor) {
  for(i in c('citrobacter','serratia','providencia','proteus','enterobacter','moraxella','morganella','shigella')) {
    dataor$pathogen[grepl(paste0(i), dataor$pathogen)] <- paste0(i,'_spp') 
  }
  for(i in c('campylobacter','listeria')) {
    dataor$pathogen[grepl(paste0(i), dataor$pathogen)] <- paste0(i) 
  }
  for(i in c('group_b','group_a')) {
    dataor$pathogen[grepl(paste0(i), dataor$pathogen)] <- paste0(i,'_strep') 
  }
  dataor$pathogen[grepl('enterococcus', dataor$pathogen) & 
                    !grepl('faec',dataor$pathogen)] <- 'enterococcus_spp'
  dataor$pathogen[grepl('ints', dataor$pathogen)|grepl('salmonellae', dataor$pathogen)|grepl('non_ty', dataor$pathogen)|dataor$pathogen=="salmonella"] <- 'non_typhoidal_salmonellae'
  dataor$pathogen[grepl('pseudomonas', dataor$pathogen) & 
                    !grepl('aeruginosa',dataor$pathogen)] <- 'pseudomonas_spp'
  dataor$pathogen[grepl('klebsiella', dataor$pathogen) & 
                    !grepl('pneumoniae',dataor$pathogen)] <- 'klebsiella_spp'
  dataor$pathogen[(dataor$pathogen == 'klebsiella')] <- 'klebsiella_spp'
  dataor$pathogen[grepl('streptococcus', dataor$pathogen) & 
                    !grepl('pneumoniae',dataor$pathogen)] <- 'streptococcus_spp'
  return(dataor)
}

assign_l2_syndrome <- function(df, available_parents, infsyn) {
  infsyn$path_to_top_parent <- strsplit(as.character(infsyn$path_to_top_parent), ",")
  infsyn$parent_candidates <- lapply(infsyn$path_to_top_parent, function(x) {
    intersect(x, available_parents)
  })
  
  infsyn <- infsyn[sapply(infsyn$parent_candidates, length) > 0, ]
  
  infsyn$parent <- sapply(infsyn$parent_candidates, function(x) { tail(x, n=1) })
  
  merged_df <- merge(df, infsyn[c("infectious_syndrome", "parent")], by.x="infectious_syndrome", by.y="infectious_syndrome", all.x=TRUE)
  merged_df$infectious_syndrome <- merged_df$parent

  final_df <- merged_df %>% select(-parent)

  return(final_df)
}


iscategory <- function(dataor = dataor) {
  for(i in c('blood','respiratory','peritonitis','mening','bone','diarrh','skin','urinary')) {
    dataor$infs[grepl(paste0(i), dataor$infectious_syndrome)] <- paste0(i) 
  }
  for(i in c('oral','eye','uri','otitis')) {
    dataor$infs[grepl(paste0(i), dataor$infectious_syndrome)] <- 'skin'
  }
  for(i in c('ditis','cardio')) {
    dataor$infs[grepl(paste0(i), dataor$infectious_syndrome)] <- 'card'
  }
  for(i in c('typh','ints','gi','gastro','enteric')) {
    dataor$infs[grepl(paste0(i), dataor$infectious_syndrome)] <- 'diarrh' 
  }
  dataor$infs[grepl('bsi',dataor$infectious_syndrome)] <- 'blood'
  dataor$infs[grepl('cns',dataor$infectious_syndrome)|grepl('encephal',dataor$infectious_syndrome)|grepl('cosis',dataor$infectious_syndrome)] <- 'mening'
  dataor$infs[grepl('intra_abdomen',dataor$infectious_syndrome)|grepl('hepatitis',dataor$infectious_syndrome)] <- 'peritonitis'
  dataor$infs[grepl('lri',dataor$infectious_syndrome)] <- 'respiratory'
  dataor$infs[grepl('genital', dataor$infectious_syndrome)|grepl('syphilis',dataor$infectious_syndrome)|grepl('sti',dataor$infectious_syndrome)|grepl('uti',dataor$infectious_syndrome)] <- 'urinary'
  dataor$infs[grepl('tb', dataor$infectious_syndrome)] <- 'tb'
  dataor$infs[grepl('unsp', dataor$infectious_syndrome)|dataor$infectious_syndrome=='L1_non_infectious_cause'|dataor$infectious_syndrome=='L1_other_infection'] <- 'unsp'
  return(dataor)
}

decade_bins <- function(dataor = dataor) {
  dataor$year <- ''
  for (t in   seq.int(from = 1990, to = 2022, by = 10)) {
    top <- t+9
    dataor$new <- ifelse(((dataor$year_id >= t) & (dataor$year_id <= top)),1,0)
    dataor$year[dataor$new == 1] <- t
    dataor$new<-NULL
  }
  dataor$year <- as.numeric(dataor$year)
  return(dataor)
}

hi_nhi <- function(dataor = dataor){
  dataor <- dataor %>% left_join(location_md2[,c('location_id','super_region_id')])
  dataor$super_region_id[is.na(dataor$super_region_id)] <- 1L
  dataor$high_income <- ifelse(dataor$super_region_id == 64,"High Income", ifelse(dataor$super_region_id == 1, "Global", "non-HighIncome"))
  dataor$high_income[is.na(dataor$super_region_id)] <- 'Global'
  return(dataor)
}

locid <- function(dataor = dataor) {
  litlocs <- read.csv(paste0(nidmeta$file_path[nidmeta$source == 'SOURCE'], suffixdt), stringsAsFactors = F) %>% dplyr::select(nid,location_id) %>% dplyr::rename(litlocs = location_id) %>% distinct()
  dataor <- dataor %>% left_join(litlocs, by = c('refvar' = 'nid')) 
  dataor$locid <- case_when(dataor$source == 'SOURCE' ~ 102, dataor$source %in% c('SOURCE','SOURCE') ~ 95,
                            dataor$source == 'SOURCE' ~ 18, dataor$source == 'SOURCE' ~ 12, dataor$source == 'SOURCE' ~ 196,
                            dataor$source == 'SOURCE' ~ 20, dataor$source == 'SOURCE' ~ 163,
                            dataor$source == 'SOURCE' ~ 10, dataor$source == 'SOURCE' ~ 45,
                            dataor$source == 'SOURCE' ~ 123, dataor$source == 'SOURCE' ~ 130,
                            dataor$source == 'SOURCE' ~ 22, dataor$source == 'SOURCE' ~ 8,
                            dataor$source == 'SOURCE' ~ 95, dataor$source == 'SOURCE' ~ 207,
                            dataor$source == 'SOURCE' ~ 161, grepl("SOURCE",dataor$source) ~ 81,
                            dataor$source == 'SOURCE' & dataor$super_region_id=='non-HighIncome' ~ 163, 
                            dataor$source == 'SOURCE' & dataor$super_region_id=='HighIncome' ~ 81,
                            dataor$source == 'SOURCE' ~ 214, dataor$source == 'SOURCE' ~ 180,
                            dataor$source == 'SOURCE' & dataor$super_region_id=='non-HighIncome' ~ 163, 
                            dataor$source == 'SOURCE' & dataor$super_region_id=='HighIncome' ~ 81,
                            dataor$source == 'SOURCE' ~ 125
                            
  )
  dataor$locid[is.na(dataor$locid)] <- dataor$litlocs[is.na(dataor$locid)]
  return(dataor)
}  

make_plot <- function(dat, size_multiplier = 1) {
  if (0) {
    dat <- data_year_age
    size_multiplier = 0.2
  }
  with(dat, plot(
    x = spliny, 
    y = coeff_unadjusted,
    cex = 1/se_unadjusted * size_multiplier
  ))
  grid()
}


plot_spline_results <- function(mod_tmp, show_knots = TRUE, data_year_age = data_year_age) {
  df_pred_tmp <- data.frame(spliny = seq(0, 100, by = 0.1))
  dat_pred_tmp <- mr$MRData()
  dat_pred_tmp$load_df(
    data = df_pred_tmp,
    col_covs=as.list('spliny')
  )
  df_pred_tmp$pred_tmp <- mod_tmp$predict(data = dat_pred_tmp)
  with(df_pred_tmp, lines(spliny, pred_tmp))
  if (show_knots) {
    for (k in as.numeric(mod_tmp$cov_models[[2]]$spline_knots)) abline(v = k, col = "gray")
  }
}