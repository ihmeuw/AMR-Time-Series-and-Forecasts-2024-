library(reticulate)
reticulate::use_python("FILEPATH")
mr <- import("mrtool")
username <- Sys.getenv("USER")
hdrive <- sprintf("FILEPATH", username)
setwd(hdrive)
amr_repo <- sprintf("FILEPATH", username)

library(dplyr)
library(ggforce)
library(ggplot2)
library(data.table)
library(gbm)
library(xgboost)
library(mgcv)
library(tidyverse)
library(glmnet)
library(matrixStats)
library(quadprog)
library(gtools)
library(nnet)
library(randomForest)
library(stats)
library(lme4)
library(Cubist)
library(caret, lib.loc="FILEPATH")
source("FILEPATH")

#~~~~~~~~~~~~~#
# i. Setup ####
#~~~~~~~~~~~~~#

oos <- FALSE
n_holdouts <- 10
oos_method <- 'country'
output_path <- 'FILEPATH'

# loading metadata
location_md2 <- subset(get_location_metadata(location_set_id=35, release_id = 16),select = c(location_id,
              level,location_name,location_name_short,super_region_id,super_region_name, region_name,region_id,ihme_loc_id,parent_id), droplevels = TRUE) 

datapath <- 'FILEPATH'
abxcovs <- read.csv(paste0(datapath,"FILEPATH")) %>% 
  dplyr::mutate(as.vector(scale(ddd_per_1000_fitted, center = TRUE, scale = TRUE))) %>%
  dplyr::mutate(as.vector(scale(prop_j01a_final, center = TRUE, scale = TRUE))) %>%
  dplyr::mutate(as.vector(scale(prop_j01b_final, center = TRUE, scale = TRUE))) %>%
  dplyr::mutate(as.vector(scale(prop_j01c_final, center = TRUE, scale = TRUE))) %>%
  dplyr::mutate(as.vector(scale(prop_j01d_final, center = TRUE, scale = TRUE))) %>%
  dplyr::mutate(as.vector(scale(prop_j01e_final, center = TRUE, scale = TRUE))) %>%
  dplyr::mutate(as.vector(scale(prop_j01f_final, center = TRUE, scale = TRUE))) %>%
  dplyr::mutate(as.vector(scale(prop_j01g_final, center = TRUE, scale = TRUE))) %>%
  dplyr::mutate(as.vector(scale(prop_j01m_final, center = TRUE, scale = TRUE)))
abxcovs <- cbind(abxcovs$location_id,abxcovs$year_id,abxcovs[,13:21])
colnames(abxcovs) <-c('location_id','year_id','ddd_per_1000','j01a','j01b','j01c','j01d','j01e','j01f','j01g','j01m')


covs_dict <- read.csv('FILEPATH', stringsAsFactors = F) %>% filter(!cov %in% c(2436,2437,319,320))
covs_dict <- covs_dict %>% filter(cov != 57) %>% distinct()
covs_list <- unique(covs_dict$cov)
covs_list <- c(covs_list, 118,	142,	145,	149,	160,	848,	1068,	1075,	1099,	1165,	1208,	1999,	2279,	2284,	2314,	2339,	2342,47, 2313)
covs_list <- unique(covs_list)

covbind <- tibble()
for(i in covs_list) { 
  j <- as.numeric(i) 
  if(is.numeric(j) & !is.na(j)){
    covtest <- get_covariate_estimates(eval(parse(text = paste0("covariate_id =", j))), release_id = 16, year_id = 1990:2022) %>% 
      dplyr::group_by(location_id,year_id,covariate_name_short) %>% dplyr::summarise(mean_value = mean(mean_value)) 
    z <- scale(covtest$mean_value, center = TRUE, scale = TRUE)
    covtest$val <- z[,1] 
    covtest <- covtest[,c('location_id','year_id','covariate_name_short','val')]
    covs_dict$covariate_name_short[covs_dict$cov == j] <- paste0(unique(covtest$covariate_name_short))[1]
  covbind <- rbind(covbind, covtest) 
}
}
covs_dict$covariate_name_short[is.na(covs_dict$covariate_name_short)] <- covs_dict$cov[is.na(covs_dict$covariate_name_short)]
covbind <- pivot_wider(covbind, id_cols = c('location_id','year_id'), names_from = 'covariate_name_short', values_from = 'val')
covbind <- covbind %>% left_join(abxcovs, by = c('location_id', 'year_id')) 
  covbind <- covbind %>% select(colnames(covbind)[!colnames(covbind) %in% c('LDI_pc','level')])
covbind <- data.frame(covbind)  

# Load dictionary of abx_class to ATCC class and plausible combinations
data_path <-'FILEPATH'
ir <- read.csv(paste0(data_path,'intrinsic_resistance.csv'), stringsAsFactors = F) %>% dplyr::select(abx_class, atc) %>% distinct()
ir$atc[ir$abx_class == 'vancomycin'] <- 'J01D'
ir<-ir[!ir$abx_class %in% c('polymixin','chloramphenicol','nitrofurantoin','second_gen_ceph','rifamycin'),]
irall<-c('abx_class' = 'all_resistant', 'atc' = 'ddd_per_1000')
ir <- rbind(ir,irall)
ir <- ir %>% add_row(abx_class = "fluoroquinolone", pathogen = "neisseria_gonorrheae",
                     atc = 'j01m', combo = 'neisseria_gonorrheae-fluoroquinolone')


# Load initial data
unadjusted_data <- fread(paste0('FILEPATH')) 
unadjusted_data[,c('index','hospital_type','nid')] <- NULL
unadjusted_data <- left_join(unadjusted_data, location_md2[,c('location_id','level','parent_id')], by = c('location_id'))
unadjusted_data$location_id <- ifelse((unadjusted_data$level == 5), unadjusted_data$parent_id, unadjusted_data$location_id)
unadjusted_data$level<-NULL
unadjusted_data$parent_id<-NULL
unadjusted_data <- left_join(unadjusted_data, location_md2[,c('location_id','level','parent_id')], by = c('location_id'))
unadjusted_data$location_id <- ifelse((unadjusted_data$level == 4), unadjusted_data$parent_id, unadjusted_data$location_id)
unadjusted_data$level<-NULL
unadjusted_data$parent_id<-NULL
unadjusted_data <- unadjusted_data %>% group_by(location_id,year_id,source,pathogen,abx_class) %>%
  summarise(cases = sum(cases),resistant_unadj=sum(resistant),susceptible_unadj=sum(susceptible))
tempids <- colnames(unadjusted_data)[!grepl('unadj',colnames(unadjusted_data))]
tempids <- colnames(unadjusted_data)[!grepl('unadj',colnames(unadjusted_data))]

# append un-adjusted combos into initial data
initial_data <- read.csv('FILEPATH') 
initial_data <- merge(initial_data, unadjusted_data, by = c(tempids), all = TRUE)
initial_data$resistant[is.na(initial_data$resistant)] <- initial_data$resistant_unadj[is.na(initial_data$resistant)]
initial_data$susceptible[is.na(initial_data$susceptible)] <- initial_data$susceptible_unadj[is.na(initial_data$susceptible)]
raw_data <- initial_data[,c(tempids,'nid','resistant','susceptible')] %>% distinct()
raw_data$sample_size<-raw_data$cases
raw_data$val<-raw_data$resistant / raw_data$cases
raw_data$prev<-raw_data$resistant / raw_data$cases
raw_data$se <- sqrt(((raw_data$val)*(1-raw_data$val))/raw_data$cases)
raw_data$variance <- ((raw_data$val)*(1-raw_data$val))/raw_data$cases
raw_data$variance[is.na(raw_data$variance)] <- raw_data$se[is.na(raw_data$variance)]*raw_data$se[is.na(raw_data$variance)]
raw_data$pathogen[grepl("salmonella_ints_non_typhi/paratyphi", raw_data$pathogen)] <- 'non_typhoidal_salmonellae'
raw_data$combo <- paste0(raw_data$pathogen,"-",raw_data$abx_class)
raw_data <- raw_data %>% left_join(location_md2[,c('location_id', 'ihme_loc_id','region_id','super_region_id')])
rm(unadjusted_data,initial_data)


# subset to the combos of interest
combos<-  read.csv(paste0(amr_repo,'FILEPATH'), stringsAsFactors = FALSE) %>% filter(COLUMN==0) %>% dplyr::select(pathogen, abx_class,modelable_entity_id) 
combos$pathogen[combos$pathogen == "acinetobacter_baumanii"] <- 'acinetobacter_baumannii'
combos$pathogen[grepl('gonorrheae', combos$pathogen)] <- 'neisseria_gonorrhoeae'
combos$abx_class[combos$abx_class=='rifampicin_new']<-'mdr'
combos$abx_class[combos$abx_class=='isoniazid_new']<-'xdr'
combos <- combos[!grepl('prop',combos$abx_class),]
combos$abx_class <- sub('_retreated',"",combos$abx_class)
combos$combo<-paste0(combos$pathogen,"-",combos$abx_class)

combos <- combos[!grepl("tuberculosis",combos$combo),]
problem_xgboost <- c('acinetobacter_baumannii-fourth_gen_ceph')
problem_gam <- c()
problem_rf <- c()
 problem_rfg <- c()
 problem_rfl <- c()


if (oos){
  args <- commandArgs(trailingOnly = T)
  list <- as.character(args[1])
}
outliering = TRUE
set.seed(123)
finder <- list()
covfound <- tibble()
# Initiate the loop for those plausible combinations 
for(x in unique(combos$combo)){

#components of combination
  y <- substring(x,1,stringr::str_locate(x,"-")-1)[1]
  z <- substring(x,stringr::str_locate(x,"-")+1,)[1]

  outputdir <-  paste0(output_path, y,'/', z,'/') 
  dir.create(outputdir, showWarnings = F, recursive = T)
  outputplots <-  paste0(output_path, '/plots/') 
  dir.create(outputplots, showWarnings = F, recursive = T)

  #any restrictions on the data. If you want to restrict to certain countries/regions add this here
  mydata <- raw_data[((raw_data$cases>=5) & (raw_data$pathogen == paste0(y)) & (raw_data$abx_class == paste0(z))),] 
    
  colids <- c('location_id','source','nid','year_id','resistant','susceptible','sample_size','super_region_id','region_id','val','se','combo')
  colcovs <- unique(covs_dict[covs_dict$combo == paste0(x), 'covariate_name_short'])
  abx_for_combo <- tolower(unique(ir$atc[ir$abx_class==paste0(z)]))
  
  if(length(colcovs) == 0){
    colcovs <- c(colnames(covbind)[!colnames(covbind) %in% c('location_id','year_id','level')], abx_for_combo) 
      }
  colcovs <- unique(colcovs)
  ncolours <- length(colcovs) + 1
  
 mydata <- mydata[,c(colids)] %>% left_join(covbind[,c('location_id','year_id',colcovs)], by = c('location_id','year_id')) 
 mydata$age_group_id <- 22
 mydata$sex_id<- 3
 mydata$measure_id<-1
 mydata$pathogen <- paste0(y)
 mydata$abx_class <- paste0(z)
 mydata$measure <- 'continuous'
 # calculate variance
 mydata$variance <- ((mydata$val)*(1-mydata$val))/mydata$sample_size
 mydata$variance[mydata$variance == 0] <- 0.0001
 
  if(outliering == TRUE) {
   sqoutliers <- mydata[(mydata$abx_class == 'aminoglycoside') & (mydata$val > 0.98),] 
  if(length(sqoutliers$nid)>0){
   sqoutliers$is_outlier <- 1
   sqoutliers <- sqoutliers[,c("location_id","ihme_loc_id","year_id","source","nid","super_region_id","region_id",
                               "age_group_id","sex_id","measure_id","measure","is_outlier","susceptible",
                               'resistant',"sample_size","val","se","variance","pathogen","abx_class")]
    check1 <-length(mydata$nid)
   mydata <- mydata[!((mydata$abx_class == 'aminoglycoside') & (mydata$val > 0.98)),]
   check2 <-length(mydata$nid)+length(sqoutliers$nid)
   stopifnot(check1==check2)
   rm(check1,check2)
  }
  }
 
 response = cbind(successes = mydata$resistant,
                  failures = mydata$susceptible)

 response = cbind(successes = mydata$resistant,
                  failures = mydata$susceptible)
 vars <- as.matrix(mydata[, c(colcovs,'location_id')])
 colnames(vars) <- c(colcovs,'location_id')
 mydata$wcv <- 1 
 family <- 'binomial'
 cv_lasso = cv.glmnet(x = vars , y= response, family = family, alpha = 1, weights = mydata$wcv, nfolds = 5, foldid = mydata$fold_id)
 
 #print out the model coefficients.
 cv_lasso$lambda.1se
   important_covs <- coef(cv_lasso, s = median(cv_lasso$lambda[cv_lasso$lambda>cv_lasso$lambda.min]))@i
   important_covs <- important_covs+1 
   all_covs <-  data.frame(coef(cv_lasso, s = "lambda.1se")@Dimnames[1])
  all_covs$row_id <- 1:nrow(all_covs)
  all_covs <- all_covs[all_covs$row_id %in% important_covs,]
 covs_to_include <- unique(all_covs[,1])
 covs_to_include <- covs_to_include[covs_to_include!="(Intercept)"]
 covs_to_include <- c("location_id", covs_to_include)

  mydatatoplot <-mydata[,c('val',covs_to_include)]

 png(filename=paste0(outputplots,x,"r.png"))
 tryCatch(pairs(mydatatoplot,                     
       labels = colnames(mydatatoplot), 
       pch = 19L,                
       bg = rainbow(ncolours), 
       main = paste0(x),    
       row1attop = TRUE,    
       gap = 1,             
       cex.labels = NULL,   
       font.labels = 1,     
       upper.panel = panel.cor, 
       lower.panel = panel.smooth),
       error = function(e) print('nograph'))
 dev.off()

 write.csv(covs_to_include, paste0(outputdir, 'covs_to_include.csv'), row.names = FALSE)
 covs_to_include <- na.omit(covs_to_include)
 covs_to_include <- covs_to_include[!is.na(covs_to_include)] 
 text <- paste(covs_to_include, collapse = " + ")
 textmodel <-paste0("glmer(response ~ 1 +",text,"+ (1 |super_region_id / region_id / location_id ) , data = mydata, family = 'binomial')")
  model1 <- eval(parse(text = (textmodel)))

    covs_to_include <- na.omit(covs_to_include)
    covs_to_include <- covs_to_include[!is.na(covs_to_include)] 
    covs <- covbind[,c('location_id','year_id', covs_to_include)] %>% right_join(location_md2[(location_md2$level > 2),c('location_id',
          'super_region_id','region_id')]) %>% left_join(location_md2[,c('location_id','ihme_loc_id')])
  covs <- data.table(covs)
  covs$pred <- predict(model1, newdata = covs, type = 'response', allow.new.levels = TRUE)
  other<- mydata[,c('location_id', 'year_id', 'val','se', 'sample_size')]
  covs <- covs %>% left_join(other, by = c('location_id','year_id')) 
  covs <- data.table(covs)
  mydata$upper_bound <-  NULL
  mydata$lower_bound <-  NULL
  mydata$is_outlier <- 0

  
  value_deviation <-2
  MADs <-  covs[,.(upper_bound = pred + value_deviation*mad(pred[!is.na(val)],val[!is.na(val)]),
                   lower_bound = pred - value_deviation*mad(pred[!is.na(val)],val[!is.na(val)])),
                by = c('location_id','year_id')]
  
  MADs <- MADs[,c('lower_bound','upper_bound')]
  covs <- cbind(covs, MADs)
  covs$upper_bound[covs$upper_bound>1] <- 1
  covs$lower_bound[covs$lower_bound<0] <- 0
  
  #define outliers in the dataset
  MADs <- covs[,.(location_id, year_id, lower_bound, upper_bound)]
  MADs <-  unique(MADs)
  mydata <- merge(mydata, MADs, by = c('location_id', 'year_id'))
  mydata$is_outlier[mydata$val<mydata$lower_bound |mydata$val>mydata$upper_bound] <- 1
  mydata$is_outlier[mydata$val>mydata$lower_bound & mydata$val<mydata$upper_bound] <- 0
  #remove value = 0 and 1
  mydata$val[mydata$val==0]<-0.02
  mydata$val[mydata$val==1]<-0.98
  
  # calculate variance
  mydata$variance <- ((mydata$val)*(1-mydata$val))/mydata$sample_size
  mydata$variance[mydata$variance == 0] <- 0.0001

# this is correction for smart outliers
  mydatainput <- mydata[,c("location_id","ihme_loc_id","year_id","source","nid","super_region_id","region_id","age_group_id","sex_id","measure_id","measure","is_outlier","susceptible",
                      'resistant',"sample_size","val","se","variance","pathogen","abx_class")]
  if(length(sqoutliers$nid)>0){
    mydatainput <- rbind(mydatainput, sqoutliers)
  }
  write.csv(mydatainput, paste0(outputdir,'input.csv'), row.names = F)
  rm(mydatainput)

# country fixed effects
    covs_to_include[covs_to_include == 'location_id'] <- 'ihme_loc_id'
    # simple fixed effect in linear prediction
    textl <- paste(covs_to_include2, collapse = " + ")
    textlmodel <-paste0("glm(response ~ 1 +",textl,", data = mydata, family = 'binomial')")
    response = cbind(successes = mydata$resistant,
                     failures = mydata$susceptible)
    linearmodel <- eval(parse(text = (textlmodel)))
    
	#getting fixed effects by country
    cf <- linearmodel$coefficients[grepl('ihme_loc_id',names(linearmodel$coefficients))]
    names(cf) = gsub('ihme_loc_id','', names(cf))
    cf<-data.frame(cf)
    cf$ihme_loc_id = row.names(cf)
    
    ldf <- covs[,c('ihme_loc_id','location_id','super_region_id')] %>% distinct()
    ldf$nodata <- ifelse(ldf$ihme_loc_id %in% unique(cf$ihme_loc_id),0,1)
    ldf <- ldf %>% left_join(cf, by = 'ihme_loc_id') 
    ldf <- ldf %>% group_by(super_region_id) %>% mutate(mcf = median(cf,na.rm  = T))
    ldf <- ldf %>% group_by(super_region_id) %>% mutate(lmcf1 = abs(cf-mcf))
    ldf<- ldf%>%group_by(super_region_id)%>%mutate(lmcf2=min(lmcf1, na.rm=T))
    sldf <- ldf[ldf$lmcf1==ldf$lmcf2, c('super_region_id','ihme_loc_id')] %>% distinct()
    gldf <- unique(ldf$ihme_loc_id[abs(ldf$cf-median(ldf$cf,na.rm=T))==min(abs(ldf$cf-median(ldf$cf,na.rm = T)), na.rm=T)])
    ldf$lmcf <- ifelse(ldf$nodata==1,sldf$ihme_loc_id[ldf$super_region_id==sldf$super_region_id],ldf$ihme_loc_id)
    ldf$lmcf[is.na(ldf$lmcf)]<-gldf[!is.na(gldf)]
    ldf$ihme_loc_id<-ldf$lmcf
    ldf$locid <- ifelse(ldf$nodata==1,location_md2$location_id[location_md2$ihme_loc_id==ldf$ihme_loc_id],ldf$location_id)

covs$ihme_loc_id<-NULL 
covs<-covs%>%left_join(ldf[,c('location_id','ihme_loc_id','locid')]) 
covs_to_include2 <- covs_to_include
covs_to_include2[covs_to_include2 == 'ihme_loc_id'] <- 'locid'
mydata$locid <- mydata$location_id
covs$linear_cv_pred <- predict(linearmodel, newdata = covs, type = 'response', allow.new.levels = TRUE)


if (oos) {
  if(oos_method == 'random'){
    mydata <- mydata[sample(nrow(mydata)),]
    mydata[,master_fold_id := cut(seq(1,nrow(mydata)),breaks=n_holdouts,labels=FALSE)]
  }
  
  if(oos_method == 'country'){
    country <- unique(mydata[, location_id])
    country <- country[sample(length(country))]
    master_fold_id <- cut(seq(1,length(country)),breaks=n_holdouts,labels=FALSE)
    folds <- as.data.table(cbind(country, master_fold_id))
    folds <- rename(folds, location_id = country)
    
    mydata <- merge(mydata, folds, by = c('location_id'))
    mydata$master_fold_id <- as.numeric(mydata$master_fold_id)
    rm(country, master_fold_id)
  }
}

#specify child models to include
if(x %in% c(problem_xgboost)) {
  child_models <- c('gam', 'ridge', 'lasso', 'enet', 'nnet', 'rf', 'cubist')  
} else {
  child_models <- c('xgboost','gam', 'ridge', 'lasso', 'enet', 'nnet', 'rf', 'cubist')  
}
  if(x %in% c(problem_gam)) {
    child_models <- c('xgboost','ridge', 'lasso', 'enet', 'nnet', 'rf', 'cubist')  
  }

  
#specify the stacker 
if(x %in% c(glm_stacker) ){
  stacker <- 'GLM'
} else {
  stacker <- 'RWM'
}

#specify the family you are modelling 
family <- 'binomial'

transformation <- NULL

centre_scale <- FALSE

include_year <-  TRUE

p <- 'val'        
n <- 'resistant'            
d <- 'sample_size'  
w <- NULL           

min_year <- 1990 
max_year <- 2022

colnames(mydata)[colnames(mydata)==d] <- 'd' 
if(!is.null(p)) {colnames(mydata)[colnames(mydata)==p] <- 'p'} 
if(!is.null(n)) {colnames(mydata)[colnames(mydata)==n] <- 'n'} 

if(is.null(n) &!is.null(p)&!is.null(d)){mydata$n <- mydata$p*mydata$d}

if(is.null(p) &!is.null(n)&!is.null(d)){mydata$p <- mydata$n/mydata$d}
mydata$p[mydata$p==0]<-0.02
mydata$p[mydata$p==1]<-0.98

#if you havent specified a weights column set to 1
if(is.null(w)){
  mydata$w <- 1
} else {
  colnames(mydata)[colnames(mydata)==w] <- 'w' 
}

#perform transformations as specified
if(is.null(transformation)){
} else if(transformation == 'log'){
  if(family == 'binomial'){
    mydata$n <- log(mydata$n)
    mydata$d <- log(mydata$d)
    mydata$p <- log(mydata$p)
  } else if(family == 'gaussian')
    mydata$n <- log(mydata$n)
} else if(transformation == 'logit'){  
  if(family == 'binomial'){
    mydata$p <- log(mydata$p/(1-mydata$p))
  } else if(family == 'gaussian'){
    message('should not be using logit transformation with gaussian data')
  }
}

#merge covs onto data
mydata <- merge(mydata, covs, by = c('location_id', 'year_id'))
mydata <- data.table(mydata)

## remove NAs
if(family == 'binomial'){
  mydata    <- na.omit(mydata, c('n', 'd', 'p', names(covs)))
}

if(family == 'gaussian'){
  mydata    <- na.omit(mydata, c('n', names(covs)))
}

if (oos){
  loopval <- n_holdouts
  # set fixed data that isn't disrupted by holdout procedure
  master_data <- mydata
  # similarly need a static output dir
  staticoutput <- paste0('FILEPATH',y,'/', z,'/')
  # and static covs
  staticcovs <- covs
  # and child models
  staticchildren <- child_models
} else {
  loopval <- 1
}

#specify within ensemble holdout method, currently can use random or country
holdout_method <- 'random'

## Run holdouts
for(h in 1:loopval){
  if (oos){
    message(paste0('Running holdout ', h))
    #limit data to a holdout
    mydata <- master_data[master_data$master_fold_id !=h,]
    mydata$master_fold_id <- NULL
    covs <- staticcovs
    child_models <- staticchildren
  }
 
  ## shuffle the data into five random folds
  if(holdout_method == 'random'){
    mydata <- mydata[sample(nrow(mydata)),]
    mydata[,fold_id := cut(seq(1,nrow(mydata)),breaks=5,labels=FALSE)]
  }
  
  if(holdout_method == 'country'){
    country <- unique(mydata[, country])
    country <- country[sample(length(country))]
    fold_id <- cut(seq(1,length(country)),breaks=5,labels=FALSE)
    folds <- as.data.table(cbind(country, fold_id))
    
    mydata <- merge(mydata, folds, by = c('country'))
    mydata$fold_id <- as.numeric(mydata$fold_id)
    rm(country, fold_id)
  }

  # Limit to required years
  mydata <- mydata[mydata$year_id >= min_year & mydata$year_id <= max_year,]
  covs <- covs[covs$year_id >= min_year & covs$year_id <= max_year,]
  
  ## add a row id column
  mydata[, a_rowid := seq(1:nrow(mydata))]
  
  # I exclude year
  if(include_year == TRUE){
    covs_to_include <- c('year_id', covs_to_include)
  }
  
  # create subfolder for holdout
  if (oos){
    outputdir <- paste0(staticoutput, '/holdout_', h, '/')
    dir.create(outputdir, showWarnings = F)
  }
  
#~~~~~~~~~~~~~~~~~~~~~#
# Fit child models ####
#~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~#
# 1. XGBoost ####
#~~~~~~~~~~~~~~~#

if('xgboost' %in% child_models){
  dir.create(paste0(outputdir, '/xgboost'), showWarnings = F)
  
  # Create model formula
  if(family == 'binomial'){
    form <- as.formula(paste0('p ~ ', paste(covs_to_include, collapse = " + ")))
  } else if(family == 'gaussian'){
    form <- as.formula(paste0('n ~ ', paste(covs_to_include, collapse = " + ")))
  }
  
  
  #tune the XGBoost
  # message("Model tuning xgboost")
  # #set the options for parameters to look at whilst tuning - can adjust these as required
  xg_grid <- expand.grid(nrounds = c(50, 100, 200),
                         max_depth = c(4, 6, 8),
                         eta = (3:6) / 100, #learning rate 
                         colsample_bytree = .5,
                         min_child_weight = 1,
                         subsample = 1,
                         gamma = 0)
  
  
  # Set cross validation options, default to 5 times repeated 5-fold cross validation
  # Selection function is "oneSE" to pick simplest model within one standard error of minimum
  # then imput this into the training model
  train_control <- trainControl(selectionFunction = "oneSE",
                                method = "repeatedcv",
                                number = 5,
                                repeats = 5,
                                index = list(mydata$a_rowid[mydata$fold_id!=1],
                                             mydata$a_rowid[mydata$fold_id!=2],
                                             mydata$a_rowid[mydata$fold_id!=3],
                                             mydata$a_rowid[mydata$fold_id!=4],
                                             mydata$a_rowid[mydata$fold_id!=5]),
                                indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                               mydata$a_rowid[mydata$fold_id==2],
                                               mydata$a_rowid[mydata$fold_id==3],
                                               mydata$a_rowid[mydata$fold_id==4],
                                               mydata$a_rowid[mydata$fold_id==5]))
  
  
  # Fit model
  xg_fit <- caret::train(form,
                  data = mydata,
                  trControl = train_control,
                  verbose = T,
                  tuneGrid = xg_grid,
                  metric = "RMSE",
                  method = "xgbTree",
                  objective = if(family == 'binomial'){"reg:logistic"}else if(family == 'gaussian'){"reg:linear"}else{message('Family of model not compatiable')},
                  weights = mydata$w)
  
  # Save model fit object
  saveRDS(xg_fit, paste0(outputdir, "/xgboost/xg_fit.RDS"))
  
  # Save the best parameters to csv file
  write.csv(xg_fit$bestTune, paste0(outputdir, 'xgboost/xgboost_best_tune_.csv'))
  xg_best_tune <- xg_fit$bestTune
  
  #set up final parameters based on the model tuning
  xg_grid_final <- expand.grid(nrounds = xg_best_tune$nrounds,
                               max_depth = xg_best_tune$max_depth,
                               eta = xg_best_tune$eta,
                               colsample_bytree = .5,
                               min_child_weight = 1,
                               subsample = 1,
                               gamma = 0)
  
  train_control_final <- trainControl(method = "cv",
                                      number = 5,
                                      savePredictions = "final",
                                      index = list(mydata$a_rowid[mydata$fold_id!=1],
                                                   mydata$a_rowid[mydata$fold_id!=2],
                                                   mydata$a_rowid[mydata$fold_id!=3],
                                                   mydata$a_rowid[mydata$fold_id!=4],
                                                   mydata$a_rowid[mydata$fold_id!=5]),
                                      indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                                     mydata$a_rowid[mydata$fold_id==2],
                                                     mydata$a_rowid[mydata$fold_id==3],
                                                     mydata$a_rowid[mydata$fold_id==4],
                                                     mydata$a_rowid[mydata$fold_id==5]))
  
  #Fit final model
  message("Fitting xgboost on final tuned hyperparameters")
  xg_fit_final <- train(form,
                        data = mydata,
                        trControl = train_control_final,
                        verbose = F,
                        tuneGrid = xg_grid_final,
                        metric = "RMSE",
                        method = "xgbTree",
                        objective = if(family == 'binomial'){"reg:logistic"}else if(family == 'gaussian'){"reg:linear"}else{message('Family of model not compatiable')},
                        weights = mydata$w)
  
  # Plot the covariate importance of final model
  cov_plot <-
    ggplot(varImp(xg_fit_final, scale = FALSE)) +
    labs(x = "Covariate", y = "Relative Importance") +
    theme_bw()
  ggsave(filename = paste0(outputdir, '/xgboost/_covariate_importance.png'),
         plot = cov_plot)
  
  # Extract out of sample and in sample predictions
  mydata[, 'xgboost_cv_pred'   := arrange(xg_fit_final$pred, rowIndex)[,"pred"]]
  mydata[, 'xgboost_full_pred' := predict(xg_fit_final, mydata)]
  
  #save model fit
  xg_fit_final$model_name <- "xgboost"
  saveRDS(xg_fit_final, paste0(outputdir, '/xgboost/full_xgboost.RDS'))
  
  #predict out for all locations
  covs[, 'xgboost' := predict(xg_fit_final, covs)]
  
  rm(form, cov_plot, train_control, train_control_final, xg_best_tune, xg_fit, xg_fit_final, xg_grid, xg_grid_final)
}

#~~~~~~~~~~~#
# 2. GAM ####
#~~~~~~~~~~~#
if('gam' %in% child_models){
  dir.create(paste0(outputdir, '/gam/'), showWarnings = F)
  
  #If there are any binary covariates then remove them from the cov list and add as additional terms
  
  #set response variable
  if(family == 'binomial'){
    response <- cbind(sucesses = mydata$n, 
                      failures = mydata$d - mydata$n)
  } else if (family == 'gaussian'){
    response <- mydata$n
  }
  
  #build the GAM formula:
  gam_formula <- paste0('response ~ 1+ s(', paste(covs_to_include2, collapse = ", bs = 'ts', k = 3) + s("), ", bs = 'ts', k = 3)")
  gam_formula <- as.formula(gam_formula)
  
  # Fit full model
  #has some sort of parrallelisation inbuilt - set using this 
  full_gam = mgcv::gam(gam_formula, 
                       data = mydata, 
                       family = if(family =='binomial'){'quasibinomial'}else if(family == 'gaussian'){'gaussian'}, 
                       weights = mydata$w, 
                       control = list(nthreads = 2))
  full_gam$model_name = 'GAM'
  
  mydata[,'gam_full_pred' := predict(full_gam, mydata, type = 'response')]
  
  for(i in 1:5){
    if(family == 'binomial'){
      response <- cbind(successes = mydata$n[mydata$fold_id!=i], 
                        failures = mydata$d[mydata$fold_id!=i] - mydata$n[mydata$fold_id!=i])
    } else if (family == 'gaussian'){
      response <- mydata$n[mydata$fold_id!=i] 
    }      
    
    baby_gam = mgcv::gam(gam_formula, 
                         data = mydata[mydata$fold_id!=i], 
                         family = if(family =='binomial'){'quasibinomial'}else if(family == 'gaussian'){'gaussian'}, 
                         weights = mydata$weight[mydata$fold_id!=i], 
                         control = list(nthreads = 2))
    
    #fill in the data
    mydata[fold_id==i, 'gam_cv_pred' := predict(baby_gam, mydata[fold_id==i,],type = 'response')] 
    
  }
  
  #save full model fit
  saveRDS(full_gam, paste0(outputdir, '/gam/full_gam.RDS'))
  
  #predict out for all locations
  covs[,'gam' := predict(full_gam, covs, type = 'response')]
  
  #plot out GAM results to analyse
  pdf(paste0(outputdir, '/gam/plots2.pdf'))
  gam.check(full_gam)
  dev.off()
  
  rm(baby_gam, full_gam, gam_formula, response)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 3. Penalised regression (E-net/Ridge/Lasso) ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
if('enet' %in% child_models | 'ridge' %in% child_models | 'lasso' %in% child_models){
  
  dir.create(paste0(outputdir, '/glmnet'),showWarnings = F)
  
  if(family == 'binomial'){
    response <- cbind(failures   = mydata$d - mydata$n, 
                      successes = mydata$n)
  }else if(family == 'gaussian'){
    response <- mydata$n
  }
  
  vars <- as.matrix(mydata[, covs_to_include2, with = F])
  colnames(vars) <- covs_to_include2
  
  cv_lambda0 = cv.glmnet(x = vars , y= response, family = family, alpha = 0, weights = mydata$w, nfolds = 5, foldid = mydata$fold_id)
  cv_lambda0.25 = cv.glmnet(x = vars , y= response, family = family, alpha = 0.25, weights = mydata$w, nfolds = 5, foldid = mydata$fold_id)
  cv_lambda0.5 = cv.glmnet(x = vars , y= response, family = family, alpha = 0.5, weights = mydata$w, nfolds = 5, foldid = mydata$fold_id)
  cv_lambda0.75 = cv.glmnet(x = vars , y= response, family = family, alpha = 0.75, weights = mydata$w, nfolds = 5, foldid = mydata$fold_id)
  cv_lambda1 = cv.glmnet(x = vars , y= response, family = family, alpha = 1, weights = mydata$w, nfolds = 5, foldid = mydata$fold_id)
  
  #plot out the lambda and alpha options
  pdf(paste0(outputdir, '/glmnet/parameter_selection.pdf'))
  par(mfrow=c(3,2))
  plot(cv_lambda0)
  plot(cv_lambda0.25)
  plot(cv_lambda0.5)
  plot(cv_lambda0.75)
  plot(cv_lambda1)
  plot(log(cv_lambda0$lambda),cv_lambda0$cvm,pch=19,col="red",xlab="log(Lambda)",ylab=cv_lambda0$name)
  points(log(cv_lambda0.25$lambda),cv_lambda0.25$cvm,pch=19,col="pink")
  points(log(cv_lambda0.5$lambda),cv_lambda0.5$cvm,pch=19,col="blue")
  points(log(cv_lambda0.75$lambda),cv_lambda0.75$cvm,pch=19,col="yellow")
  points(log(cv_lambda1$lambda),cv_lambda1$cvm,pch=19,col="green")
  legend("bottomright",legend=c("alpha= 1","alpha= .75", "alpha= .5", "alpha= .25","alpha 0"),pch=19,col=c("green","yellow","blue","pink","red"))
  dev.off()
  
  #fit the full model using selected lambda and alpha
  if('ridge' %in% child_models){full_ridge = glmnet(x = vars , y= response, family = family, alpha = 0, weights = mydata$w)}
  if('enet' %in% child_models){full_enet = glmnet(x = vars , y= response, family = family, alpha = 0.5, weights = mydata$w)}
  if('lasso' %in% child_models){full_lasso = glmnet(x = vars , y= response, family = family, alpha = 1, weights = mydata$w)}
  
  #predict full model results (requires matrix)
  if('ridge' %in% child_models){mydata[,'ridge_full_pred' := predict(full_ridge,newx = vars, s = cv_lambda0$lambda.1se, type = 'response')]}
  if('lasso' %in% child_models){mydata[,'lasso_full_pred' := predict(full_lasso,newx = vars, s = cv_lambda1$lambda.1se, type = 'response')]}
  if('enet' %in% child_models){mydata[,'enet_full_pred' := predict(full_enet,newx = vars, s = cv_lambda0.5$lambda.1se, type = 'response')]}
  
  #fit the model on the holdouts
  for(i in 1:5){
    if(family == 'binomial'){
      response <- cbind(failures = mydata$d[mydata$fold_id!=i] - mydata$n[mydata$fold_id!=i], 
                        successes = mydata$n[mydata$fold_id!=i])
    } else if(family == 'gaussian'){
      response <- mydata$n[mydata$fold_id!=i] 
    }
    
    vars <- as.matrix(mydata[fold_id != i, covs_to_include2, with = F])
    colnames(vars) <- covs_to_include2
    
    if('ridge' %in% child_models){baby_ridge = glmnet(x = vars , y= response, family = family, lambda = cv_lambda0$lambda.1se, alpha = 0, weights = mydata$w[mydata$fold_id!=i])}
    if('lasso' %in% child_models){baby_lasso = glmnet(x = vars , y= response, family = family, lambda = cv_lambda1$lambda.1se, alpha = 1, weights = mydata$w[mydata$fold_id!=i])}
    if('enet' %in% child_models){baby_enet = glmnet(x = vars , y= response, family = family, lambda = cv_lambda0.5$lambda.1se, alpha = 0.5, weights = mydata$w[mydata$fold_id!=i])}
    
    new_vars <- as.matrix(mydata[fold_id == i, covs_to_include2, with = F])
    
    #fill in the data
    if('ridge' %in% child_models){mydata[fold_id==i,'ridge_cv_pred' := predict(baby_ridge,newx = new_vars, s = cv_lambda0$lambda.1se, type = 'response')]}
    if('lasso' %in% child_models){mydata[fold_id==i,'lasso_cv_pred' := predict(baby_lasso,newx = new_vars, s = cv_lambda1$lambda.1se, type = 'response')]}
    if('enet' %in% child_models){mydata[fold_id==i,'enet_cv_pred' := predict(baby_enet,newx = new_vars, s = cv_lambda0.5$lambda.1se, type = 'response')]}
  }
  
  #save the model and relevent coefficients
  if('ridge' %in% child_models){saveRDS(cv_lambda0, paste0(outputdir, '/glmnet/full_ridge.rds'))}
  if('enet' %in% child_models){saveRDS(cv_lambda0.5, paste0(outputdir, '/glmnet/full_enet.rds'))}
  if('lasso' %in% child_models){saveRDS(cv_lambda1, paste0(outputdir, '/glmnet/full_lasso.rds'))}
  
  #predict out for all locations - can change the lambdas if required
  new_covs <- as.matrix(covs[, covs_to_include2, with = F])
  names(new_covs) <- covs_to_include2
  if('ridge' %in% child_models){covs[,'ridge' := predict(full_ridge,newx = new_covs[,rownames(full_ridge$beta)], s = cv_lambda0$lambda.1se, type = 'response')]}
  if('enet' %in% child_models){covs[,'enet' := predict(full_enet,newx = new_covs[,rownames(full_enet$beta)], s = cv_lambda0.5$lambda.1se, type = 'response')]}
  if('lasso' %in% child_models){covs[,'lasso' := predict(full_lasso,newx = new_covs[,rownames(full_lasso$beta)], s = cv_lambda1$lambda.1se, type = 'response')]}
  
  rm(cv_lambda1, cv_lambda0.5, cv_lambda0, cv_lambda0.25, cv_lambda0.75, full_lasso, full_enet, full_ridge, baby_lasso, baby_ridge, baby_enet, new_vars, response, vars, i, new_covs, all_names)    
}

#~~~~~~~~~~~~~~~~~~~~~#
# 4. Random forest ####
#~~~~~~~~~~~~~~~~~~~~~#
if('rf' %in% child_models){
  dir.create(paste0(outputdir, '/rf'), showWarnings = F)
  
  # Create model formula
  if(family == 'binomial'){
    form <- as.formula(paste0('p ~ ', paste(covs_to_include, collapse = " + ")))
  } else if(family == 'gaussian'){
    form <- as.formula(paste0('n ~ ', paste(covs_to_include, collapse = " + ")))
  }
  
  train_control <- trainControl(selectionFunction = "best",
                                method = "repeatedcv",
                                number = 5,
                                repeats = 5,
                                index = list(mydata$a_rowid[mydata$fold_id!=1],
                                             mydata$a_rowid[mydata$fold_id!=2],
                                             mydata$a_rowid[mydata$fold_id!=3],
                                             mydata$a_rowid[mydata$fold_id!=4],
                                             mydata$a_rowid[mydata$fold_id!=5]),
                                indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                               mydata$a_rowid[mydata$fold_id==2],
                                               mydata$a_rowid[mydata$fold_id==3],
                                               mydata$a_rowid[mydata$fold_id==4],
                                               mydata$a_rowid[mydata$fold_id==5]),
                                search = 'grid')
  
  tunegrid <- expand.grid(.mtry=c(1:length(covs_to_include)))
  
  # Fit model
  rf_fit <- train(form,
                  data = mydata,
                  trControl = train_control,
                  verbose = T,
                  tuneGrid = tunegrid,
                  metric = "RMSE",
                  method = "rf",
                  weights = mydata$w)
  
  # Save model fit object 
  saveRDS(rf_fit, paste0(outputdir, "/rf/rf_fit.RDS"))
  png(paste0(outputdir, '/rf/rf_fit.png'))
  plot(rf_fit)  
  dev.off()
  
  #specify the parameters
  mtry_tune <- rf_fit$bestTune$mtry
  write.csv(mtry_tune, paste0(outputdir, '/rf/rf_params.csv'), row.names = F)
  tunegrid_final <- expand.grid(.mtry=mtry_tune)
  
  #specify the folds in the train control section
  train_control_final <- trainControl(method = "cv",
                                      number = 5,
                                      savePredictions = "final",
                                      index = list(mydata$a_rowid[mydata$fold_id!=1],
                                                   mydata$a_rowid[mydata$fold_id!=2],
                                                   mydata$a_rowid[mydata$fold_id!=3],
                                                   mydata$a_rowid[mydata$fold_id!=4],
                                                   mydata$a_rowid[mydata$fold_id!=5]),
                                      indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                                     mydata$a_rowid[mydata$fold_id==2],
                                                     mydata$a_rowid[mydata$fold_id==3],
                                                     mydata$a_rowid[mydata$fold_id==4],
                                                     mydata$a_rowid[mydata$fold_id==5]))
  
  
  #fit the final model
  rf_fit_final <- train(form,
                        data = mydata,
                        trControl = train_control_final,
                        verbose = T,
                        tuneGrid = tunegrid_final,
                        metric = "RMSE",
                        method = "rf",
                        importance=T,
                        weights = mydata$w)
  
  # Extract out of sample and in sample predictions
  mydata[, 'rf_cv_pred'   := arrange(rf_fit_final$pred, rowIndex)[,"pred"]]
  mydata[, 'rf_full_pred' := predict(rf_fit_final, mydata)]
  
  #save model fit
  rf_fit_final$model_name <- "rf"
  saveRDS(rf_fit_final, paste0(outputdir, '/rf/full_rf.RDS'))
  
  cov_plot <-
    ggplot(varImp(rf_fit_final, scale = FALSE)) +
    labs(x = "Covariate", y = "Relative Importance") +
    theme_bw()
  ggsave(filename = paste0(outputdir, '/rf/rf_covariate_importance.png'),
         plot = cov_plot)
  
  #predict out for all locations
  covs[, 'rf' := predict(rf_fit_final, covs)]
  covs$rf <- ifelse(covs$rf < 0, 0, covs$rf)
  
  rm(form, cov_plot, train_control, train_control_final, rf_fit, rf_fit_final, tunegrid, tunegrid_final)
}

#~~~~~~~~~~~~~~~~~~~~~~~#
# 5. Neural networks ####
#~~~~~~~~~~~~~~~~~~~~~~~#
if('nnet' %in% child_models){
  dir.create(paste0(outputdir, '/nnet'), showWarnings = F)
  # Create model formula
  if(family == 'binomial'){
    form <- as.formula(paste0('p ~ ', paste(covs_to_include, collapse = " + ")))
    linoutval = FALSE
  } else if(family == 'gaussian'){
    form <- as.formula(paste0('n ~ ', paste(covs_to_include, collapse = " + ")))
    linoutval = TRUE
  }
  
  train_control <- trainControl(selectionFunction = "best",
                                method = "repeatedcv",
                                number = 5,
                                repeats = 5,
                                index = list(mydata$a_rowid[mydata$fold_id!=1],
                                             mydata$a_rowid[mydata$fold_id!=2],
                                             mydata$a_rowid[mydata$fold_id!=3],
                                             mydata$a_rowid[mydata$fold_id!=4],
                                             mydata$a_rowid[mydata$fold_id!=5]),
                                indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                               mydata$a_rowid[mydata$fold_id==2],
                                               mydata$a_rowid[mydata$fold_id==3],
                                               mydata$a_rowid[mydata$fold_id==4],
                                               mydata$a_rowid[mydata$fold_id==5]),
                                search = 'grid')
  
  tunegrid <- expand.grid(.decay = c(1, 0.5, 0.1, 0.01, 0.001, 0.0001, 0.00001), .size = c(4, 5, 6, 7, 8, 9))
  
  # Fit model
  nn_fit <- train(form,
                  data = mydata,
                  trControl = train_control,
                  verbose = T,
                  tuneGrid = tunegrid,
                  metric = "RMSE",
                  method = "nnet",
                  linout = linoutval,
                  maxit = 1000,
                  weights = mydata$w)
  
  # Save model fit object 
  saveRDS(nn_fit, paste0(outputdir, "/nnet/nn_fit.RDS"))
  png(paste0(outputdir, '/nnet/nn_fit.png'))
  plot(nn_fit)  
  dev.off()
  
  # Save the best parameters to csv file
  write.csv(nn_fit$bestTune, paste0(outputdir, '/nnet/nnet_best_tune.csv'))
  
  #specify the parameters
  tunegrid_final <- expand.grid(.decay=nn_fit$bestTune$decay, .size=nn_fit$bestTune$size)
  
  #specify the folds in the train control section
  train_control_final <- trainControl(method = "cv",
                                      number = 5,
                                      savePredictions = "final",
                                      index = list(mydata$a_rowid[mydata$fold_id!=1],
                                                   mydata$a_rowid[mydata$fold_id!=2],
                                                   mydata$a_rowid[mydata$fold_id!=3],
                                                   mydata$a_rowid[mydata$fold_id!=4],
                                                   mydata$a_rowid[mydata$fold_id!=5]),
                                      indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                                     mydata$a_rowid[mydata$fold_id==2],
                                                     mydata$a_rowid[mydata$fold_id==3],
                                                     mydata$a_rowid[mydata$fold_id==4],
                                                     mydata$a_rowid[mydata$fold_id==5]))
  
  
  nn_fit_final <- train(form,
                        data = mydata,
                        trControl = train_control_final,
                        verbose = T,
                        tuneGrid = tunegrid_final,
                        metric = "RMSE",
                        method = "nnet",
                        linout = linoutval,
                        maxit = 1000,
                        weights = mydata$w)
  
  # Extract out of sample and in sample predictions
  mydata[, 'nnet_cv_pred'   := arrange(nn_fit_final$pred, rowIndex)[,"pred"]]
  mydata[, 'nnet_full_pred' := predict(nn_fit_final, mydata)]
  
  #save model fit
  nn_fit_final$model_name <- "nn"
  saveRDS(nn_fit_final, paste0(outputdir, '/nnet/full_nn.RDS'))
  
  cov_plot <-
    ggplot(varImp(nn_fit_final, scale = FALSE)) +
    labs(x = "Covariate", y = "Relative Importance") +
    theme_bw()
  ggsave(filename = paste0(outputdir, '/nnet/nn_covariate_importance.png'),
         plot = cov_plot)
  
  #predict out for all locations
  covs[, 'nnet' := predict(nn_fit_final, covs)]
  
  rm(nn_fit, nn_fit_final, train_control, train_control_final, tunegrid, tunegrid_final, cov_plot)
}

#~~~~~~~~~~~~~~~~~~~~#
# 6. Cubist model ####
#~~~~~~~~~~~~~~~~~~~~#
if('cubist' %in% child_models){
  dir.create(paste0(outputdir, '/cubist'), showWarnings = F)
  
  # Create model formula
  if(family == 'binomial'){
    form <- as.formula(paste0('p ~ ', paste(covs_to_include, collapse = " + ")))
  } else if(family == 'gaussian'){
    form <- as.formula(paste0('n ~ ', paste(covs_to_include, collapse = " + ")))
  }
  
  train_control <- trainControl(selectionFunction = "best",
                                method = "repeatedcv",
                                number = 5,
                                repeats = 5,
                                index = list(mydata$a_rowid[mydata$fold_id!=1],
                                             mydata$a_rowid[mydata$fold_id!=2],
                                             mydata$a_rowid[mydata$fold_id!=3],
                                             mydata$a_rowid[mydata$fold_id!=4],
                                             mydata$a_rowid[mydata$fold_id!=5]),
                                indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                               mydata$a_rowid[mydata$fold_id==2],
                                               mydata$a_rowid[mydata$fold_id==3],
                                               mydata$a_rowid[mydata$fold_id==4],
                                               mydata$a_rowid[mydata$fold_id==5]),
                                search = 'grid')
  
  tunegrid <- expand.grid(.committees = c(seq(1, 40, 5)), 
                          .neighbors  = c(0, 3, 6, 9))
  
  # Fit model
  cubist_fit <- train(form,
                      data = mydata,
                      trControl = train_control,
                      verbose = T,
                      tuneGrid = tunegrid,
                      metric = "RMSE",
                      method = "cubist",
                      control = Cubist::cubistControl(),
                      weights = mydata$w)
  
  saveRDS(cubist_fit, paste0(outputdir, "/cubist/cubist_fit.RDS"))
  png(paste0(outputdir, '/cubist/cubist_fit.png'))
  plot(cubist_fit)  
  dev.off()
  
  write.csv(cubist_fit$bestTune, paste0(outputdir, '/cubist/cubist_best_tune.csv'))
  
  tunegrid_final <- expand.grid(.committees=cubist_fit$bestTune$committees, .neighbors=cubist_fit$bestTune$neighbors)
  
  train_control_final <- trainControl(method = "cv",
                                      number = 5,
                                      savePredictions = "final",
                                      index = list(mydata$a_rowid[mydata$fold_id!=1],
                                                   mydata$a_rowid[mydata$fold_id!=2],
                                                   mydata$a_rowid[mydata$fold_id!=3],
                                                   mydata$a_rowid[mydata$fold_id!=4],
                                                   mydata$a_rowid[mydata$fold_id!=5]),
                                      indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                                     mydata$a_rowid[mydata$fold_id==2],
                                                     mydata$a_rowid[mydata$fold_id==3],
                                                     mydata$a_rowid[mydata$fold_id==4],
                                                     mydata$a_rowid[mydata$fold_id==5]))
  
  
  cubist_fit_final <- train(form,
                            data = mydata,
                            trControl = train_control_final,
                            verbose = T,
                            tuneGrid = tunegrid_final,
                            metric = "RMSE",
                            method = "cubist",
                            weights = mydata$w)
  
  mydata[, 'cubist_cv_pred'   := arrange(cubist_fit_final$pred, rowIndex)[,"pred"]]
  mydata[, 'cubist_full_pred' := predict(cubist_fit_final, mydata)]
  
  mydata$cubist_cv_pred <- ifelse(mydata$cubist_cv_pred > 1, 1, mydata$cubist_cv_pred)
  mydata$cubist_full_pred <- ifelse(mydata$cubist_full_pred > 1, 1, mydata$cubist_full_pred)
  
  #save model fit
  cubist_fit_final$model_name <- "cubist"
  saveRDS(cubist_fit_final, paste0(outputdir, '/cubist/full_cubist.RDS'))
  
  cov_plot <-
    ggplot(varImp(cubist_fit_final, scale = FALSE)) +
    labs(x = "Covariate", y = "Relative Importance") +
    theme_bw()
  ggsave(filename = paste0(outputdir, '/cubist/cubist_covariate_importance.png'),
         plot = cov_plot)
  
  #predict out for all locations
  covs[, 'cubist' := predict(cubist_fit_final, covs)]
  covs$cubist <- ifelse(covs$cubist > 1, 1, covs$cubist)
  
  rm(form, cov_plot, train_control, train_control_final, cubist_fit, cubist_fit_final, tunegrid, tunegrid_final)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Print out correlations between data and predictions ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

table_r2 <-c()
for(i in 1:length(child_models)){
  if(family == 'binomial'){
    r2 <- c(child_models[i], formatC(cor(mydata$p, mydata[,get(paste0(child_models[i], '_cv_pred'))])^2), format = 'f', digits = 6)
    table_r2 <- rbind(table_r2,r2)
    message(paste0(child_models[i], ' correlation: ', round(cor(mydata$p, mydata[,get(paste0(child_models[i], '_cv_pred'))])^2,2)))
  }
  if(family == 'gaussian'){
    r2 <- c(child_models[i], formatC(cor(mydata$n, mydata[,get(paste0(child_models[i], '_cv_pred'))])^2), format = 'f', digits = 6)
    table_r2 <- rbind(table_r2,r2)
    message(paste0(child_models[i], ' correlation: ', round(cor(mydata$n, mydata[,get(paste0(child_models[i], '_cv_pred'))])^2,2)))
  }
}
colnames(table_r2) <- c('model','r2','format','digits')
table_r2 <- table_r2[,c('model','r2')]
table_r2 <- data.table(table_r2)
table_r2$r2 <- as.numeric(table_r2$r2)
table_r2 <- table_r2[order(-table_r2$r2),]
write.csv(table_r2, paste0(outputdir, '/table_r2.csv'), row.names = F)
table_r2 <- na.omit(table_r2)
table_r2 <- table_r2[order(-table_r2$r2),]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Check the correlation of the stackers ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

child_models <- unique(table_r2$model)
final_child_models <- unique(table_r2$model)
table_cor <- c()
for(i in 1:length(child_models)){
  for(j in 1:length(child_models)){
    if(i==j){
    } else{
      test1 <- sd(mydata[,get(paste0(child_models[i], '_cv_pred'))])
      test2 <- sd(mydata[,get(paste0(child_models[j], '_cv_pred'))])
      if(test1 > 0 & test2 > 0) {
        cor1 <- c(child_models[i],child_models[j],formatC(cor(mydata[,get(paste0(child_models[i], '_cv_pred'))],mydata[,get(paste0(child_models[j], '_cv_pred'))])^2, format = 'f', digits = 6))
        delete <- ifelse(cor1[[3]] < 0.8,"",ifelse(table_r2$r2[table_r2$model == child_models[i]] < table_r2$r2[table_r2$model == child_models[j]],child_models[i],child_models[j]))
        final_child_models <- final_child_models[!final_child_models == delete]
        row1 <- c(cor1,delete)
        table_cor <- rbind(table_cor,row1)
        rm(cor1, delete, row1)
        if(cor(mydata[,get(paste0(child_models[i], '_cv_pred'))],mydata[,get(paste0(child_models[j], '_cv_pred'))])^2>0.80){message(paste0(child_models[i],  ' and ', child_models[j], ' correlated, remove one'))
        }
      }
    }
  }
}

write.csv(table_cor, paste0(outputdir, '/table_correlation.csv'), row.names = F)
write.csv(final_child_models, paste0(outputdir, '/table_child_models.csv'), row.names = F)

final_child_columns <- c(paste(paste0(final_child_models,'_cv_pred'),sep = ","), paste(paste0(final_child_models,'_full_pred'),sep=","))
child_models <- final_child_models

if(centre_scale == TRUE){
  mydata$year <-  NULL
  covs$year <- NULL
}
write.csv(mydata, paste0(outputdir, '/fitted_child_models.csv'), row.names = F)
write.csv(covs, paste0(outputdir, '/child_model_preds.csv'), row.names = F)

mydata <- read.csv(paste0(outputdir, '/fitted_child_models.csv'), stringsAsFactors =  F)
covs <- read.csv(paste0(outputdir, '/child_model_preds.csv'), stringsAsFactors = F)
mydata <- data.table(mydata)
covs <- data.table(covs)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Combined the child model estimates ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
headcols <- colnames(mydata)[!colnames(mydata) %in% colnames(mydata)[grep('pred',colnames(mydata))]]
headcols <- c(headcols,final_child_columns)
stackers <- subset(mydata, select = headcols, drop.levels=TRUE)
stackers <- data.frame(stackers)
X <- as.matrix(stackers[colnames(stackers)[(grep('cv_pred', colnames(stackers)))]])
Y = if(family == 'binomial'){stackers$p}else if(family == 'gaussian'){stackers$n}

C <- data.frame(covs)
C <- as.matrix(C[c(child_models)])

#Select which stacker you are using and stack the estimates

if(stacker == 'CWM'){
  # code is from http://zoonek.free.fr/blosxom/R/2012-06-01_Optimization.html

  s <- solve.QP( 
    t(X) %*% X, t(Y) %*% X, 
    cbind(  
      matrix(1, nr=length(child_models), nc=1),
      diag(length(child_models)),
      -diag(length(child_models))
    ),
    c(1, 
      rep(0.000001, length(child_models)),
      rep(-1, length(child_models))), 
    meq = 1
	)
  
  mydata$stacked_preds <- rowWeightedMeans(X, w = s$solution)   
  
  #Calculate the stacked predictions
  covs$cv_custom_stage_1 <- rowWeightedMeans(C, w = s$solution)
}

if(stacker == 'RWM'){
  r2 <- rep(NA, length(child_models))
  for(i in 1:length(child_models)){
    if(family == 'binomial'){
      r2[i] <- round(cor(mydata$p, mydata[,get(paste0(child_models[i], '_cv_pred'))])^2,2)
    }
    if(family == 'gaussian'){
      r2[i] <- round(cor(mydata$n, mydata[,get(paste0(child_models[i], '_cv_pred'))])^2,2)
    }
  }
  
  total <-  sum(r2)
  weights <- r2/total
  mydata$stacked_preds <- rowWeightedMeans(X, w = weights)   
  covs$cv_custom_stage_1 <- rowWeightedMeans(C, w = weights)
}

if(stacker == 'GBM'){
  if(family == 'gaussian'){
    form <- as.formula(paste0('n ~ ', paste(colnames(stackers)[(grep('cv_pred', colnames(stackers)))], collapse = " + ")))
  }
  if(family == 'binomial'){
    form <- as.formula(paste0('logit(p) ~ ', paste(colnames(stackers)[(grep('cv_pred', colnames(stackers)))], collapse = " + ")))
  }
  train_control_final <- trainControl(method = "cv",
                                      number = 5,
                                      savePredictions = "final",
                                      index = list(mydata$a_rowid[mydata$fold_id!=1],
                                                   mydata$a_rowid[mydata$fold_id!=2],
                                                   mydata$a_rowid[mydata$fold_id!=3],
                                                   mydata$a_rowid[mydata$fold_id!=4],
                                                   mydata$a_rowid[mydata$fold_id!=5]),
                                      indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                                     mydata$a_rowid[mydata$fold_id==2],
                                                     mydata$a_rowid[mydata$fold_id==3],
                                                     mydata$a_rowid[mydata$fold_id==4],
                                                     mydata$a_rowid[mydata$fold_id==5]))
  
  model_gbm<- 
    train(form, data = mydata, method='gbm', trControl=train_control_final, tuneLength=3)
  
  mydata[, 'stacked_preds'   := inv.logit(predict(model_gbm, mydata))]
  covs <- data.frame(covs)
  colnames(covs)[colnames(covs) %in% child_models] <- colnames(stackers)[(grep('cv_pred', colnames(stackers)))]
  covs$cv_custom_stage_1 <- inv.logit(predict(model_gbm, covs[colnames(covs)[(grep('cv_pred', colnames(covs)))]])) 
  covs <- data.table(covs)
  
  #save the covariate importance plot
  jpeg(paste0(outputdir, 'stacker_cov_imporatance.jpeg'),
       height = 10, width = 10, units = 'cm', res = 150)
  ggplot(varImp(model_gbm, scale = FALSE)) +
    labs(x = "Covariate", y = "Relative Importance") +
    theme_bw()
  
  dev.off()
}
if(stacker == 'GLM'){
  if(family == 'gaussian'){
    form <- as.formula(paste0('n ~ ', paste(colnames(stackers)[(grep('cv_pred', colnames(stackers)))], collapse = " + ")))
  }
  if(family == 'binomial'){
    form <- as.formula(paste0('p ~ ', paste(colnames(stackers)[(grep('cv_pred', colnames(stackers)))], collapse = " + ")))
  }
  
  train_control_final <- trainControl(method = "cv",
                                      number = 5,
                                      savePredictions = "final",
                                      index = list(mydata$a_rowid[mydata$fold_id!=1],
                                                   mydata$a_rowid[mydata$fold_id!=2],
                                                   mydata$a_rowid[mydata$fold_id!=3],
                                                   mydata$a_rowid[mydata$fold_id!=4],
                                                   mydata$a_rowid[mydata$fold_id!=5]),
                                      indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                                     mydata$a_rowid[mydata$fold_id==2],
                                                     mydata$a_rowid[mydata$fold_id==3],
                                                     mydata$a_rowid[mydata$fold_id==4],
                                                     mydata$a_rowid[mydata$fold_id==5]))
  
  model_glm<- 
    train(form, data = mydata, method='glm', trControl=train_control_final, tuneLength=3)
  
  mydata[, 'stacked_preds'   := predict(model_glm, mydata)]
  
  covs <- data.frame(covs)
  colnames(covs)[colnames(covs) %in% child_models] <- colnames(stackers)[(grep('cv_pred', colnames(stackers)))]
  covs$cv_custom_stage_1 <- predict(model_glm, covs[colnames(covs)[(grep('cv_pred', colnames(covs)))]]) 
  covs <- data.table(covs)
}

if(stacker == 'nnet'){
  if(family == 'gaussian'){
    form <- as.formula(paste0('n ~ ', paste(colnames(stackers)[(grep('cv_pred', colnames(stackers)))], collapse = " + ")))
  }
  if(family == 'binomial'){
    form <- as.formula(paste0('p ~ ', paste(colnames(stackers)[(grep('cv_pred', colnames(stackers)))], collapse = " + ")))
  }  
  train_control_final <- trainControl(method = "cv",
                                      number = 5,
                                      savePredictions = "final",
                                      index = list(mydata$a_rowid[mydata$fold_id!=1],
                                                   mydata$a_rowid[mydata$fold_id!=2],
                                                   mydata$a_rowid[mydata$fold_id!=3],
                                                   mydata$a_rowid[mydata$fold_id!=4],
                                                   mydata$a_rowid[mydata$fold_id!=5]),
                                      indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                                     mydata$a_rowid[mydata$fold_id==2],
                                                     mydata$a_rowid[mydata$fold_id==3],
                                                     mydata$a_rowid[mydata$fold_id==4],
                                                     mydata$a_rowid[mydata$fold_id==5]))
  
  model_nnet<- 
    train(form, data = mydata, method='nnet', trControl=train_control_final, tuneLength=3)
  
  # mydata[, 'stacked_preds'   := arrange(model_nnet$pred, rowIndex)[,"pred"]]
  mydata[, 'stacked_preds'   := predict(model_nnet, mydata)]
  
  covs <- data.frame(covs)
  colnames(covs)[colnames(covs) %in% child_models] <- colnames(stackers)[(grep('cv_pred', colnames(stackers)))]
  covs$cv_custom_stage_1 <- predict(model_gbm, covs[colnames(covs)[(grep('cv_pred', colnames(covs)))]]) 
  covs <- data.table(covs)
}

stg1 <- covs[, .(location_id, year_id,
                 age_group_id = rep(22, length(covs$location_id)),
                 sex_id = rep(3, length(covs$location_id)),
                 cv_custom_stage_1)]


#remove covariates from the dataset
mydata[, colnames(mydata)[grep('^cv_', colnames(mydata))] := NULL]

#check that the estimates are within expected range
max(stg1$cv_custom_stage_1, na.rm = T)
stg1$cv_custom_stage_1[stg1$cv_custom_stage_1>=1] <- 0.9999

#save prediction
write.csv(mydata, paste0(outputdir, '/fitted_stackers.csv'), row.names = F)
write.csv(stg1, paste0(outputdir, '/custom_stage1_df.csv'), row.names = F)

# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# # Plot out the is and oos preds ####
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 
for(stack in c(child_models, 'stacked_preds')){
  png(paste0(outputdir, '/', stack, '_oos.png'),
      height = 12, width = 12, res = 350, unit = 'in')
  print(
    ggplot(mydata)+
      geom_point(aes(x = if(family == 'binomial'){p}else if(family == 'gaussian'){n},
                     y=if(stack == 'stacked_preds'){stacked_preds} else{get(paste0(stack, '_cv_pred'))}))+
      geom_abline(slope = 1, intercept = 0, colour = 'red')+
      theme_bw() +
      theme(strip.background = element_rect(fill = "white")) +
      labs(
        x = "Data Estimate",
        y = "Mean Prediction",
        size = "Weight",
        title = (paste0("Validation Plot for ", stack)),
        subtitle = "OOS: TRUE")
  )
  dev.off()
}

#in sample plots
for(stack in child_models){
  png(paste0(outputdir, '/', stack, '_is.png'),
      height = 12, width = 12, res = 350, unit = 'in')
  print(
    ggplot(mydata)+
      geom_point(aes(x = if(family == 'binomial'){p}else if(family == 'gaussian'){n},
                     y=get(paste0(stack, '_full_pred'))))+
      geom_abline(slope = 1, intercept = 0, colour = 'red')+
      theme_bw() +
      theme(strip.background = element_rect(fill = "white")) +
      labs(
        x = "Data Estimate",
        y = "Mean Prediction",
        size = "Weight",
        title = (paste0("Validation Plot for ", stack)),
        subtitle = "OOS: FALSE")
  )
  dev.off()
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Calculate metrics for all models ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

child_model_metrics <- data.frame(child_models)
child_model_metrics$r_sq <- NA

mydata <- data.frame(mydata)
for(i in 1:length(child_models)){
  cm <- paste0(child_models, '_cv_pred')[i]
  pred <- mydata[c(cm)]
  pred <-  unlist(pred)
  child_model_metrics$rmse[child_model_metrics$child_models == child_models[i]] <- round(RMSE(pred, if(family == 'binomial'){mydata$p}else if(family == 'gaussian'){mydata$n}),4)
  child_model_metrics$r_sq[child_model_metrics$child_models == child_models[i]] <- round(cor(pred, if(family == 'binomial'){mydata$p}else if(family == 'gaussian'){mydata$n})^2,2)
}

child_model_metrics$child_models <- as.character(child_model_metrics$child_models)
child_model_metrics <- rbind(child_model_metrics, c('Stackers', NA, NA) )
child_model_metrics$rmse[child_model_metrics$child_models=='Stackers'] <- round(RMSE(mydata$stacked_preds, if(family == 'binomial'){mydata$p}else if(family == 'gaussian'){mydata$n}),4)
child_model_metrics$r_sq[child_model_metrics$child_models=='Stackers'] <- round(cor(mydata$stacked_preds, if(family == 'binomial'){mydata$p}else if(family == 'gaussian'){mydata$n})^2,2)

write.csv(child_model_metrics, paste0(outputdir, '/national_stacker_metrics.csv'), row.names = F)

#~~~~~~~~~~~~~~~~~~~~~#
# Plot the results ####
#~~~~~~~~~~~~~~~~~~~~~#
library(foreign)

covs <- data.table(covs)

mydata <- data.frame(mydata)
input <- mydata[c('location_id', 'year_id', 'p')] 
covs <- merge(covs, input, by = c('location_id', 'year_id'), all.x = T, all.y = T, allow.cartesian = T)#, #, 'age_category'

colnames(covs) <-  gsub('_cv_pred', '', colnames(covs))
covs <- melt(covs, id.vars = c('location_id', 'year_id', 'p'), 
             measure.vars = c(child_models, 'cv_custom_stage_1'))

covs <- merge(covs, location_md2[,c('location_id','level','ihme_loc_id','region_id')], all.x = T, all.y = F)
covs <- covs[covs$level == 3,]

covs$lintype <- ifelse(covs$variable == "cv_custom_stage_1", 'a', 'b')
covs$size <- ifelse(covs$variable == "cv_custom_stage_1", 1.5, 1)

pdf(paste0(outputplots,'st1_',y,'-',z,'.pdf'),
    height = 8.3, width = 11.7)

#plot out a page for each region
for(i in 1:length(unique(covs$region_id))){
  subset <- covs[covs$region_id == unique(covs$region_id)[i],]
  print(
    ggplot(subset)+
      geom_line(aes(x=year_id, y = value, colour = variable, linetype = lintype, size = size), alpha = 0.8)+
      scale_size_continuous(range = c(1, 1.5)) +
      geom_point(aes(x=year_id, y = (p)))+
      facet_wrap(~ihme_lc_id, nrow = ceiling(sqrt(length(unique(subset$location_id)))))+
      ylim(0, max(covs$value))
  )
}
dev.off()
}
}


############################
# Config template for stgpr
###########################
source('FILEPATH')
  output_path <- 'FILEPATH'
if (!oos){
  config_table <- c()
  for(x in list) {
    #components of combination
    y <- substring(x,1,stringr::str_locate(x,"-")-1)[1]
    z <- substring(x,stringr::str_locate(x,"-")+1,)[1]
    description <- x #paste0("raw-data test for ",x)
    outputdir <-  paste0(output_path,y,'/', z,'/') 
    path_to_data <-  paste0(outputdir,'input.csv') 
    path_to_custom_stage_1 <-  paste0(outputdir,'custom_stage1_df.csv') 
    modelable_entity_id <- ifelse(length(combos$modelable_entity_id[combos$combo == paste0(x)]) == 0, 27429,combos$modelable_entity_id[combos$combo == paste0(x)])
    location_set_id <- 35
    release_id <- 16
    st_omega <- '1,1'
    st_zeta <- '0.075,0.001'
    st_lambda <- '0.4,0.4'
    data_transform <- 'logit'
    prediction_age_group_ids <- 22
    prediction_sex_ids <- 3

    prediction_units <-'proportion resistance'
    gpr_scale <- '20,20'
    gpr_amp_method <- 'broken_stick'
    gpr_draws <- 100
    metric_id <- 3
    dat <- read.csv(path_to_data, stringsAsFactors = FALSE)
    year_start <- 1990 
    year_end <- 2021 
    by_loc_sum = group_by(dat, location_id) %>%
      summarize(
        num_dp = length(val[is_outlier==0])
      )
    density_cutoffs <- max(2, ceiling(quantile(by_loc_sum$num_dp, 0.25))[[1]])
    
    config <- data.frame(modelable_entity_id,release_id,metric_id,data_transform,prediction_units,st_lambda,st_omega,st_zeta,gpr_scale,path_to_data,path_to_custom_stage_1,location_set_id,year_start,year_end,prediction_age_group_ids,prediction_sex_ids,description,gpr_amp_method,gpr_draws,density_cutoffs)
    write.csv(config,paste0(output_path,'config_table.csv'),row.names = F)
    stgpr_version_id <- register_stgpr_model(paste0(output_path,'config_table.csv'))
    stgpr_sendoff(stgpr_version_id, project = "proj_amr", nparallel = 50, log_path = paste0(output_path,'logs/'))
    config$run_id <- stgpr_version_id
    config_table <- rbind(config_table,config)
  }  
  config_table <- as.data.frame(config_table)
  config_table$model_index_id <- seq.int(nrow(config_table))
  write.csv(config_table,paste0(output_path,'FILEPATH'),row.names = F)
}
get_model_status(stgpr_version_id)
