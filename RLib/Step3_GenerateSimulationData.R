#################################################################################
# Project: Exposure measurement error                                           #
# Code: Step 3 - generate datasets for simulation                               #
# Machine: QNAP                                                                 #
#################################################################################

############################# 0. Setup ##############################
rm(list=ls())
gc()

library(utils)
library(data.table)
library(dplyr)

dir_uncertainty <- '~/data/Uncertainty/'
dir_simulation <- '~/data/SimulationData/'
dir_monitor_correlation <- '~/data/SpatialCorrelation/'



##################### 1. generate independent errors ######################
# load uncertainty data for all years
years <- 2000:2016
ls_data_uncertainty <- list()
for (i in 1:length(years)) {
  data_uncertainty_tmp <- readRDS(paste0(dir_uncertainty,'PredUncertainty_PM25_USZIP_',years[i],'.rds'))
  data_uncertainty_tmp <- data_uncertainty_tmp %>% group_by(ZIP) %>% summarise(pred_res_sd=mean(pred_res_sd,na.rm=TRUE))
  data_uncertainty_tmp <- data_uncertainty_tmp[complete.cases(data_uncertainty_tmp),]
  data_uncertainty_tmp$year <- years[i]
  ls_data_uncertainty[[i]] <- data_uncertainty_tmp
  rm(data_uncertainty_tmp)
  gc()
}
data_uncertainty <- rbindlist(ls_data_uncertainty)
data_uncertainty$year <- as.character(data_uncertainty$year)

# merge uncertainty data into covariate data 
data_counts <- readRDS(paste0(dir_simulation,'counts_weight_pred_sd_20211221.rds'))
data_counts$year <- as.character(data_counts$year)
data_counts <- left_join(data_counts,data_uncertainty,by=c('zip'='ZIP','year'='year'))

# generate 20*50 datasets for simulation
for (i in 1:20) {
  counts_covar_simu_set_50 <- data.frame()
  for (j in 1:50) {
    counts_covar <- data_counts
    counts_covar$annual_ind_err_sd <- rnorm(nrow(counts_covar),0,counts_covar$pred_res_sd)
    # annual_ind_err_sd <- counts_covar$annual_ind_err_sd[counts_covar$annual_ind_err_sd>=-2&counts_covar$annual_ind_err_sd<=2]
    # hist(annual_ind_err_sd, col=rgb(0,0,1,0.2), xlim=c(-2,2),ylim=c(0,4.5),breaks=seq(from=-2, to=2, by=0.1), freq=FALSE)
    counts_covar$annual_ind_err_2sd <- rnorm(nrow(counts_covar),0,2*counts_covar$pred_res_sd)
    # annual_ind_err_2sd <- counts_covar$annual_ind_err_2sd[counts_covar$annual_ind_err_2sd<=2&counts_covar$annual_ind_err_2sd>=-2]
    # hist(annual_ind_err_2sd, col=rgb(0,0,1,0.2), xlim=c(-2,2),ylim=c(0,4.5),breaks=seq(from=-2, to=2, by=0.1), freq=FALSE)
    counts_covar$annual_ind_err_3sd <- rnorm(nrow(counts_covar),0,3*counts_covar$pred_res_sd)
    # annual_ind_err_3sd <- counts_covar$annual_ind_err_3sd[counts_covar$annual_ind_err_3sd<=2&counts_covar$annual_ind_err_3sd>=-2]
    # hist(annual_ind_err_3sd, col=rgb(0,0,1,0.2), xlim=c(-2,2),ylim=c(0,4.5),breaks=seq(from=-2, to=2, by=0.1), freq=FALSE)
    counts_covar$set <- j
    counts_covar_simu_set_50 <- rbind.data.frame(counts_covar_simu_set_50,counts_covar)
    
    print(paste0("set ",i,'-',j," is done"))
    rm(counts_covar)
    gc()
  }
  counts_covar_simu_set_50 <- counts_covar_simu_set_50[,c('zip','year','pm25','annual_ind_err_sd','annual_ind_err_2sd','annual_ind_err_3sd')]
  saveRDS(counts_covar_simu_set_50,file = paste0(dir_simulation,'counts_covar_simu_set_50_p',i,'.rds'))
  rm(counts_covar_simu_set_50)
  gc()
}



##################### 2. generate spatially correlated errors ######################
# load paired zip code datasets
join_zip_less10km <- readRDS(paste0(dir_monitor_correlation,"join_zip_less10km_all.rds"))
join_zip_10_20km <- readRDS(paste0(dir_monitor_correlation,"join_zip_10_20km_all.rds"))
join_zip_20_40km <- readRDS(paste0(dir_monitor_correlation,"join_zip_20_40km_all.rds"))

# create a function to add spatial correlation adjustment 
spatial_corr_adj_add <- function(data_counts,join_zip_less10km,join_zip_10_20km,join_zip_20_40km,b0,b1,b2,b3){
  # merge in the error for <=10km. should have all the merged annual errors for nearest zipcodes for each zipcode in the dataset 
  join_zip_less10km <- merge(join_zip_less10km,data_counts[c("zip","year","pm25","annual_ind_err_sd","annual_ind_err_2sd","annual_ind_err_3sd")],
                             all.x=TRUE,by.x=c('zipcode_near','year'),by.y=c('zip','year')) 
  # compute the mean annual error for each zipcode
  output_less10km_sd <- join_zip_less10km %>% group_by(year,zipcode) %>% summarize(annual_ind_err_sd_avg_less10km=mean(annual_ind_err_sd,na.rm=TRUE))
  output_less10km_2sd <- join_zip_less10km %>% group_by(year,zipcode) %>% summarize(annual_ind_err_2sd_avg_less10km=mean(annual_ind_err_2sd,na.rm=TRUE))
  output_less10km_3sd <- join_zip_less10km %>% group_by(year,zipcode) %>% summarize(annual_ind_err_3sd_avg_less10km=mean(annual_ind_err_3sd,na.rm=TRUE))
  output_less10km <- left_join(output_less10km_sd,output_less10km_2sd,by=c('zipcode','year'))
  output_less10km <- left_join(output_less10km,output_less10km_3sd,by=c('zipcode','year'))
  # link back to original count data
  data_counts <- merge(data_counts,output_less10km,all.x=TRUE,by.x=c('zip','year'),by.y=c('zipcode','year')) # keep all rows of count data valid
  # clean up memory
  rm(join_zip_less10km,output_less10km_sd,output_less10km_2sd,output_less10km_3sd,output_less10km)
  gc()
  
  # do the same for 10-20km
  join_zip_10_20km <- merge(join_zip_10_20km,data_counts[c("zip","year","pm25","annual_ind_err_sd","annual_ind_err_2sd","annual_ind_err_3sd")],
                            all.x=TRUE,by.x=c('zipcode_near','year'),by.y=c('zip','year')) 
  # compute the mean annual error for each zipcode
  output_10_20km_sd <- join_zip_10_20km %>% group_by(year,zipcode) %>% summarize(annual_ind_err_sd_avg_10_20km=mean(annual_ind_err_sd,na.rm=TRUE))
  output_10_20km_2sd <- join_zip_10_20km %>% group_by(year,zipcode) %>% summarize(annual_ind_err_2sd_avg_10_20km=mean(annual_ind_err_2sd,na.rm=TRUE))
  output_10_20km_3sd <- join_zip_10_20km %>% group_by(year,zipcode) %>% summarize(annual_ind_err_3sd_avg_10_20km=mean(annual_ind_err_3sd,na.rm=TRUE))
  output_10_20km <- left_join(output_10_20km_sd,output_10_20km_2sd,by=c('zipcode','year'))
  output_10_20km <- left_join(output_10_20km,output_10_20km_3sd,by=c('zipcode','year'))
  # link back to original count data
  data_counts <- merge(data_counts,output_10_20km,all.x=TRUE,by.x=c('zip','year'),by.y=c('zipcode','year'))
  # clean up memory
  rm(join_zip_10_20km,output_10_20km_sd,output_10_20km_2sd,output_10_20km_3sd,output_10_20km)
  gc()
  
  # do the same for 20-40km
  join_zip_20_40km <- merge(join_zip_20_40km,data_counts[c("zip","year","pm25","annual_ind_err_sd","annual_ind_err_2sd","annual_ind_err_3sd")],
                            all.x=TRUE,by.x=c('zipcode_near','year'),by.y=c('zip','year')) 
  # compute the mean annual error for each zipcode
  output_20_40km_sd <- join_zip_20_40km %>% group_by(year,zipcode) %>% summarize(annual_ind_err_sd_avg_20_40km=mean(annual_ind_err_sd,na.rm=TRUE))
  output_20_40km_2sd <- join_zip_20_40km %>% group_by(year,zipcode) %>% summarize(annual_ind_err_2sd_avg_20_40km=mean(annual_ind_err_2sd,na.rm=TRUE))
  output_20_40km_3sd <- join_zip_20_40km %>% group_by(year,zipcode) %>% summarize(annual_ind_err_3sd_avg_20_40km=mean(annual_ind_err_3sd,na.rm=TRUE))
  output_20_40km <- left_join(output_20_40km_sd,output_20_40km_2sd,by=c('zipcode','year'))
  output_20_40km <- left_join(output_20_40km,output_20_40km_3sd,by=c('zipcode','year'))
  # link back to original count data
  data_counts <- merge(data_counts,output_20_40km,all.x=TRUE,by.x=c('zip','year'),by.y=c('zipcode','year'))
  # clean up memory
  rm(join_zip_20_40km,output_20_40km_sd,output_20_40km_2sd,output_20_40km_3sd,output_20_40km)
  gc()
  
  # replace NAs with 0s
  data_counts$annual_ind_err_sd_avg_less10km <- ifelse(is.na(data_counts$annual_ind_err_sd_avg_less10km),0,data_counts$annual_ind_err_sd_avg_less10km)
  data_counts$annual_ind_err_2sd_avg_less10km <- ifelse(is.na(data_counts$annual_ind_err_2sd_avg_less10km),0,data_counts$annual_ind_err_2sd_avg_less10km)
  data_counts$annual_ind_err_3sd_avg_less10km <- ifelse(is.na(data_counts$annual_ind_err_3sd_avg_less10km),0,data_counts$annual_ind_err_3sd_avg_less10km)
  data_counts$annual_ind_err_sd_avg_10_20km <- ifelse(is.na(data_counts$annual_ind_err_sd_avg_10_20km),0,data_counts$annual_ind_err_sd_avg_10_20km)
  data_counts$annual_ind_err_2sd_avg_10_20km <- ifelse(is.na(data_counts$annual_ind_err_2sd_avg_10_20km),0,data_counts$annual_ind_err_2sd_avg_10_20km)
  data_counts$annual_ind_err_3sd_avg_10_20km <- ifelse(is.na(data_counts$annual_ind_err_3sd_avg_10_20km),0,data_counts$annual_ind_err_3sd_avg_10_20km)
  data_counts$annual_ind_err_sd_avg_20_40km <- ifelse(is.na(data_counts$annual_ind_err_sd_avg_20_40km),0,data_counts$annual_ind_err_sd_avg_20_40km)
  data_counts$annual_ind_err_2sd_avg_20_40km <- ifelse(is.na(data_counts$annual_ind_err_2sd_avg_20_40km),0,data_counts$annual_ind_err_2sd_avg_20_40km)
  data_counts$annual_ind_err_3sd_avg_20_40km <- ifelse(is.na(data_counts$annual_ind_err_3sd_avg_20_40km),0,data_counts$annual_ind_err_3sd_avg_20_40km)
  
  # add spatial correlation to independent error using the monitor spatial correlation coefficients b0, b1, b2, b3
  # new error = old error-0.006525+0.873326*err_less10km+0.033130*err_10_20km+0.046112*err_20_40km
  # if any of the <=10km, 10-20km, or 20-40km average error takes missing values, do not apply the correction formula
  data_counts$annual_sp_err_sd <- ifelse(((data_counts$annual_ind_err_sd_avg_less10km==0)|(data_counts$annual_ind_err_sd_avg_10_20km==0)|(data_counts$annual_ind_err_sd_avg_20_40km==0)), # 0 indicates missing originally 
                                         data_counts$annual_ind_err_sd, 
                                         data_counts$annual_ind_err_sd+b0+b1*data_counts$annual_ind_err_sd_avg_less10km+b2*data_counts$annual_ind_err_sd_avg_10_20km+b3*data_counts$annual_ind_err_sd_avg_20_40km)
  data_counts$annual_sp_err_2sd <- ifelse(((data_counts$annual_ind_err_2sd_avg_less10km==0)|(data_counts$annual_ind_err_2sd_avg_10_20km==0)|(data_counts$annual_ind_err_2sd_avg_20_40km==0)), # 0 indicates missing originally 
                                          data_counts$annual_ind_err_2sd, 
                                          data_counts$annual_ind_err_2sd+b0+b1*data_counts$annual_ind_err_2sd_avg_less10km+b2*data_counts$annual_ind_err_2sd_avg_10_20km+b3*data_counts$annual_ind_err_2sd_avg_20_40km)
  data_counts$annual_sp_err_3sd <- ifelse(((data_counts$annual_ind_err_3sd_avg_less10km==0)|(data_counts$annual_ind_err_3sd_avg_10_20km==0)|(data_counts$annual_ind_err_3sd_avg_20_40km==0)), # 0 indicates missing originally 
                                          data_counts$annual_ind_err_3sd, 
                                          data_counts$annual_ind_err_3sd+b0+b1*data_counts$annual_ind_err_3sd_avg_less10km+b2*data_counts$annual_ind_err_3sd_avg_10_20km+b3*data_counts$annual_ind_err_3sd_avg_20_40km)
  
  data_counts <- data_counts[,c("zip","year","pm25",
                                "annual_ind_err_sd","annual_ind_err_2sd","annual_ind_err_3sd",
                                "annual_sp_err_sd","annual_sp_err_2sd","annual_sp_err_3sd")]
  
  return(data_counts)
}

# generate spatially adjusted errors for each dataset each piece
for (i in 1:20){
  counts_covar_simu_set_50 <- readRDS(paste0(dir_simulation,'counts_covar_simu_set_50_p',i,'.rds'))
  new_data <- data.frame()
  for (j in 1:50) {
    n1 <- 649910*(j-1)+1
    n2 <- 649910*j
    add_data <- spatial_corr_adj_add(data_counts=counts_covar_simu_set_50[n1:n2,],
                                     join_zip_less10km=join_zip_less10km,
                                     join_zip_10_20km=join_zip_10_20km,
                                     join_zip_20_40km=join_zip_20_40km,
                                     b0=-0.006525,
                                     b1=0.873326,
                                     b2=0.033130,
                                     b3=0.046112)
    new_data <- rbind(new_data,add_data)
    print(paste0("set ",i,'-',j," is done"))
    rm(add_data)
    gc()
  }
  saveRDS(new_data,file = paste0(dir_simulation,'counts_covar_simu_set_50_p',i,'.rds'))
  rm(counts_covar_simu_set_50,new_data)
  gc()
}


