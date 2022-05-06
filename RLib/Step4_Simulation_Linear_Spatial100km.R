###############################################################################
# Project: Exposure measurement error                                         #
# Code: Step 4 - simulation - linear dr relation                              #
# Machine: QNAP                                                               #
###############################################################################

######################################### 0. set up ################################################
rm(list=ls())
gc()

library(stats)
library(MASS)
library(nlme)
library(splines)
library(mgcv)
library(dplyr)
library(doParallel)

dir_data <- '/media/qnap4/Yaguang/EME/data/SimulationData/'
dir_data_sensitivity <- '/media/qnap4/Yaguang/EME/data/SimulationData_Sensitivity/'
dir_results <- '/media/qnap4/Yaguang/EME/results/Linear_Spatial100km/'  ### create directory to save results



######################################### 1. function ################################################
simulate_results <- function(key,covar,b0,b1,n.reps){
  ##### the following codes are for testing only, comment out when running real analysis
  # covar <- readRDS(paste0(dir_data,'counts_weight_pred_sd_20211221.rds'))
  # # create binary variables for calendar year
  # covar$year2001 <- ifelse(covar$year==2001,1,0)
  # covar$year2002 <- ifelse(covar$year==2002,1,0)
  # covar$year2003 <- ifelse(covar$year==2003,1,0)
  # covar$year2004 <- ifelse(covar$year==2004,1,0)
  # covar$year2005 <- ifelse(covar$year==2005,1,0)
  # covar$year2006 <- ifelse(covar$year==2006,1,0)
  # covar$year2007 <- ifelse(covar$year==2007,1,0)
  # covar$year2008 <- ifelse(covar$year==2008,1,0)
  # covar$year2009 <- ifelse(covar$year==2009,1,0)
  # covar$year2010 <- ifelse(covar$year==2010,1,0)
  # covar$year2011 <- ifelse(covar$year==2011,1,0)
  # covar$year2012 <- ifelse(covar$year==2012,1,0)
  # covar$year2013 <- ifelse(covar$year==2013,1,0)
  # covar$year2014 <- ifelse(covar$year==2014,1,0)
  # covar$year2015 <- ifelse(covar$year==2015,1,0)
  # covar$year2016 <- ifelse(covar$year==2016,1,0)
  # # create binary variables for region5
  # covar$region5_2 <- ifelse(covar$region5==2,1,0)
  # covar$region5_3 <- ifelse(covar$region5==3,1,0)
  # covar$region5_4 <- ifelse(covar$region5==4,1,0)
  # covar$region5_5 <- ifelse(covar$region5==5,1,0)
  # 
  # key_files <- list.files(path=dir_data,pattern = "^counts_covar_simu_set_50_(.*)rds$")
  # key <- readRDS(paste0(dir_data,key_files[1]))
  # 
  # b0 <- log(8)
  # b1 <- 0.005
  # n.reps <- 50
  ##### stop commenting for testing
  
  # dataframe to store resutls
  results <- data.frame(matrix(NA, nrow=n.reps*6, ncol=5))
  names(results) <- c('type','corr_magnitude','pollutants','b0_est','b1_est')
  results$b0 <- b0
  results$b1 <- b1
  
  # run in loops 
  for (i in 1:n.reps){
    # i=1
    cat(paste0("set ",i," starts running \n"))
    # create start/end index 
    start <- (i-1)*649910+1
    end <- 649910*i
    
    # link covar with measurement error set 
    DATA <- merge(key[start:end,],covar,all.x=TRUE,by.x=c('zip','year'),by.y=c('zip','year'))
    # head(DATA,20)
    DATA <- na.omit(DATA)
    
    # set IDs for each ob
    DATA$id <- 1:nrow(DATA)
    
    ### classical errors with ozone and no2 adjustment ###
    # generate true counts 
    DATA$y_classical_multi <- rpois(nrow(DATA),exp(b0+b1*DATA$pm25.x+4.795e-04*DATA$ozone_summer+7.526e-04*DATA$no2+6.693e-04*DATA$temp+7.395e-04*DATA$rh
                                                   -1.665e-03*DATA$PctEye-1.015e-03*DATA$PctLDL-1.525e-03*DATA$Pctmam+6.110e-01*DATA$LungCancerRate
                                                   +1.920e-01*DATA$poverty-4.245e-06*DATA$popdensity-2.574e-07*DATA$medianhousevalue-4.967e-02*DATA$pct_blk
                                                   -9.385e-07*DATA$medhouseholdincome-2.910e-01*DATA$pct_owner_occ-2.623e-01*DATA$hispanic+1.448e-01*DATA$education
                                                   +9.841e-02*DATA$smoke_rate+6.779e-03*DATA$mean_bmi+5.998e-04*DATA$amb_visit_pct+4.921e-04*DATA$a1c_exm_pct
                                                   -3.336e-03*DATA$nearest_hospital_km-5.609e-03*DATA$year2001-1.341e-03*DATA$year2002-1.192e-02*DATA$year2003
                                                   -4.005e-02*DATA$year2004-2.362e-02*DATA$year2005-1.901e-02*DATA$year2006-2.230e-02*DATA$year2007-1.388e-02*DATA$year2008
                                                   -6.234e-02*DATA$year2009-6.470e-02*DATA$year2010-7.065e-02*DATA$year2011-9.103e-02*DATA$year2012
                                                   -1.029e-01*DATA$year2013-1.207e-01*DATA$year2014-1.134e-01*DATA$year2015-1.329e-01*DATA$year2016
                                                   -5.513e-02*DATA$region5_2-2.685e-02*DATA$region5_3-5.430e-02*DATA$region5_4-3.998e-02*DATA$region5_5))
    
    # run Poisson regressions 
    mod_classical_multi_sp_error_sd <- glm(y_classical_multi ~ I(pm25.x+annual_sp_err_sd) + ozone_summer + no2 + temp + rh + PctEye + PctLDL + Pctmam + LungCancerRate + poverty 
                                           + popdensity + medianhousevalue + pct_blk + medhouseholdincome + pct_owner_occ + hispanic + education + smoke_rate 
                                           + mean_bmi + amb_visit_pct + a1c_exm_pct + nearest_hospital_km + year2001 + year2002 + year2003 + year2004 + year2005 
                                           + year2006 + year2007 + year2008 + year2009 + year2010 + year2011 + year2012 + year2013 + year2014 + year2015 + year2016 
                                           + region5_2 + region5_3 + region5_4 + region5_5, data=DATA, family=quasipoisson) 
    results[((i-1)*6+1),c('type','corr_magnitude','pollutants')] <- c('classical','sp_error_sd','multi')
    results[((i-1)*6+1),c('b0_est','b1_est')] <- mod_classical_multi_sp_error_sd$coefficients[1:2]
    
    mod_classical_multi_sp_error_2sd <- glm(y_classical_multi ~ I(pm25.x+annual_sp_err_2sd) + ozone_summer + no2 + temp + rh + PctEye + PctLDL + Pctmam + LungCancerRate + poverty 
                                            + popdensity + medianhousevalue + pct_blk + medhouseholdincome + pct_owner_occ + hispanic + education + smoke_rate 
                                            + mean_bmi + amb_visit_pct + a1c_exm_pct + nearest_hospital_km + year2001 + year2002 + year2003 + year2004 + year2005 
                                            + year2006 + year2007 + year2008 + year2009 + year2010 + year2011 + year2012 + year2013 + year2014 + year2015 + year2016 
                                            + region5_2 + region5_3 + region5_4 + region5_5, data=DATA, family=quasipoisson) 
    results[((i-1)*6+2),c('type','corr_magnitude','pollutants')] <- c('classical','sp_error_2sd','multi')
    results[((i-1)*6+2),c('b0_est','b1_est')] <- mod_classical_multi_sp_error_2sd$coefficients[1:2]
    
    mod_classical_multi_sp_error_3sd <- glm(y_classical_multi ~ I(pm25.x+annual_sp_err_3sd) + ozone_summer + no2 + temp + rh + PctEye + PctLDL + Pctmam + LungCancerRate + poverty 
                                            + popdensity + medianhousevalue + pct_blk + medhouseholdincome + pct_owner_occ + hispanic + education + smoke_rate 
                                            + mean_bmi + amb_visit_pct + a1c_exm_pct + nearest_hospital_km + year2001 + year2002 + year2003 + year2004 + year2005 
                                            + year2006 + year2007 + year2008 + year2009 + year2010 + year2011 + year2012 + year2013 + year2014 + year2015 + year2016 
                                            + region5_2 + region5_3 + region5_4 + region5_5, data=DATA, family=quasipoisson) 
    results[((i-1)*6+3),c('type','corr_magnitude','pollutants')] <- c('classical','sp_error_3sd','multi')
    results[((i-1)*6+3),c('b0_est','b1_est')] <- mod_classical_multi_sp_error_3sd$coefficients[1:2]
    
    # clear memory
    DATA$y_classical_multi <- NULL
    rm(mod_classical_multi_sp_error_sd,mod_classical_multi_sp_error_2sd,mod_classical_multi_sp_error_3sd)
    gc()
    
    
    ### Berkson errors with ozone and no2 adjustment ###
    DATA$y_Berkson_multi_sp_error_sd <- rpois(nrow(DATA),exp(b0+b1*(DATA$pm25.x+DATA$annual_sp_err_sd)+4.795e-04*DATA$ozone_summer+7.526e-04*DATA$no2+6.693e-04*DATA$temp+7.395e-04*DATA$rh
                                                             -1.665e-03*DATA$PctEye-1.015e-03*DATA$PctLDL-1.525e-03*DATA$Pctmam+6.110e-01*DATA$LungCancerRate
                                                             +1.920e-01*DATA$poverty-4.245e-06*DATA$popdensity-2.574e-07*DATA$medianhousevalue-4.967e-02*DATA$pct_blk
                                                             -9.385e-07*DATA$medhouseholdincome-2.910e-01*DATA$pct_owner_occ-2.623e-01*DATA$hispanic+1.448e-01*DATA$education
                                                             +9.841e-02*DATA$smoke_rate+6.779e-03*DATA$mean_bmi+5.998e-04*DATA$amb_visit_pct+4.921e-04*DATA$a1c_exm_pct
                                                             -3.336e-03*DATA$nearest_hospital_km-5.609e-03*DATA$year2001-1.341e-03*DATA$year2002-1.192e-02*DATA$year2003
                                                             -4.005e-02*DATA$year2004-2.362e-02*DATA$year2005-1.901e-02*DATA$year2006-2.230e-02*DATA$year2007-1.388e-02*DATA$year2008
                                                             -6.234e-02*DATA$year2009-6.470e-02*DATA$year2010-7.065e-02*DATA$year2011-9.103e-02*DATA$year2012
                                                             -1.029e-01*DATA$year2013-1.207e-01*DATA$year2014-1.134e-01*DATA$year2015-1.329e-01*DATA$year2016
                                                             -5.513e-02*DATA$region5_2-2.685e-02*DATA$region5_3-5.430e-02*DATA$region5_4-3.998e-02*DATA$region5_5))
    mod_Berkson_multi_sp_error_sd <- glm(y_Berkson_multi_sp_error_sd ~ pm25.x + ozone_summer + no2 + temp + rh + PctEye + PctLDL + Pctmam + LungCancerRate + poverty 
                                         + popdensity + medianhousevalue + pct_blk + medhouseholdincome + pct_owner_occ + hispanic + education + smoke_rate 
                                         + mean_bmi + amb_visit_pct + a1c_exm_pct + nearest_hospital_km + year2001 + year2002 + year2003 + year2004 + year2005 
                                         + year2006 + year2007 + year2008 + year2009 + year2010 + year2011 + year2012 + year2013 + year2014 + year2015 + year2016 
                                         + region5_2 + region5_3 + region5_4 + region5_5, data=DATA, family=quasipoisson) 
    results[((i-1)*6+4),c('type','corr_magnitude','pollutants')] <- c('Berkson','sp_error_sd','multi')
    results[((i-1)*6+4),c('b0_est','b1_est')] <- mod_Berkson_multi_sp_error_sd$coefficients[1:2]
    
    DATA$y_Berkson_multi_sp_error_2sd <- rpois(nrow(DATA),exp(b0+b1*(DATA$pm25.x+DATA$annual_sp_err_2sd)+4.795e-04*DATA$ozone_summer+7.526e-04*DATA$no2+6.693e-04*DATA$temp+7.395e-04*DATA$rh
                                                              -1.665e-03*DATA$PctEye-1.015e-03*DATA$PctLDL-1.525e-03*DATA$Pctmam+6.110e-01*DATA$LungCancerRate
                                                              +1.920e-01*DATA$poverty-4.245e-06*DATA$popdensity-2.574e-07*DATA$medianhousevalue-4.967e-02*DATA$pct_blk
                                                              -9.385e-07*DATA$medhouseholdincome-2.910e-01*DATA$pct_owner_occ-2.623e-01*DATA$hispanic+1.448e-01*DATA$education
                                                              +9.841e-02*DATA$smoke_rate+6.779e-03*DATA$mean_bmi+5.998e-04*DATA$amb_visit_pct+4.921e-04*DATA$a1c_exm_pct
                                                              -3.336e-03*DATA$nearest_hospital_km-5.609e-03*DATA$year2001-1.341e-03*DATA$year2002-1.192e-02*DATA$year2003
                                                              -4.005e-02*DATA$year2004-2.362e-02*DATA$year2005-1.901e-02*DATA$year2006-2.230e-02*DATA$year2007-1.388e-02*DATA$year2008
                                                              -6.234e-02*DATA$year2009-6.470e-02*DATA$year2010-7.065e-02*DATA$year2011-9.103e-02*DATA$year2012
                                                              -1.029e-01*DATA$year2013-1.207e-01*DATA$year2014-1.134e-01*DATA$year2015-1.329e-01*DATA$year2016
                                                              -5.513e-02*DATA$region5_2-2.685e-02*DATA$region5_3-5.430e-02*DATA$region5_4-3.998e-02*DATA$region5_5))
    mod_Berkson_multi_sp_error_2sd <- glm(y_Berkson_multi_sp_error_2sd ~ pm25.x + ozone_summer + no2 + temp + rh + PctEye + PctLDL + Pctmam + LungCancerRate + poverty 
                                          + popdensity + medianhousevalue + pct_blk + medhouseholdincome + pct_owner_occ + hispanic + education + smoke_rate 
                                          + mean_bmi + amb_visit_pct + a1c_exm_pct + nearest_hospital_km + year2001 + year2002 + year2003 + year2004 + year2005 
                                          + year2006 + year2007 + year2008 + year2009 + year2010 + year2011 + year2012 + year2013 + year2014 + year2015 + year2016 
                                          + region5_2 + region5_3 + region5_4 + region5_5, data=DATA, family=quasipoisson) 
    results[((i-1)*6+5),c('type','corr_magnitude','pollutants')] <- c('Berkson','sp_error_2sd','multi')
    results[((i-1)*6+5),c('b0_est','b1_est')] <- mod_Berkson_multi_sp_error_2sd$coefficients[1:2]
    
    DATA$y_Berkson_multi_sp_error_3sd <- rpois(nrow(DATA),exp(b0+b1*(DATA$pm25.x+DATA$annual_sp_err_3sd)+4.795e-04*DATA$ozone_summer+7.526e-04*DATA$no2+6.693e-04*DATA$temp+7.395e-04*DATA$rh
                                                              -1.665e-03*DATA$PctEye-1.015e-03*DATA$PctLDL-1.525e-03*DATA$Pctmam+6.110e-01*DATA$LungCancerRate
                                                              +1.920e-01*DATA$poverty-4.245e-06*DATA$popdensity-2.574e-07*DATA$medianhousevalue-4.967e-02*DATA$pct_blk
                                                              -9.385e-07*DATA$medhouseholdincome-2.910e-01*DATA$pct_owner_occ-2.623e-01*DATA$hispanic+1.448e-01*DATA$education
                                                              +9.841e-02*DATA$smoke_rate+6.779e-03*DATA$mean_bmi+5.998e-04*DATA$amb_visit_pct+4.921e-04*DATA$a1c_exm_pct
                                                              -3.336e-03*DATA$nearest_hospital_km-5.609e-03*DATA$year2001-1.341e-03*DATA$year2002-1.192e-02*DATA$year2003
                                                              -4.005e-02*DATA$year2004-2.362e-02*DATA$year2005-1.901e-02*DATA$year2006-2.230e-02*DATA$year2007-1.388e-02*DATA$year2008
                                                              -6.234e-02*DATA$year2009-6.470e-02*DATA$year2010-7.065e-02*DATA$year2011-9.103e-02*DATA$year2012
                                                              -1.029e-01*DATA$year2013-1.207e-01*DATA$year2014-1.134e-01*DATA$year2015-1.329e-01*DATA$year2016
                                                              -5.513e-02*DATA$region5_2-2.685e-02*DATA$region5_3-5.430e-02*DATA$region5_4-3.998e-02*DATA$region5_5))
    mod_Berkson_multi_sp_error_3sd <- glm(y_Berkson_multi_sp_error_3sd ~ pm25.x + ozone_summer + no2 + temp + rh + PctEye + PctLDL + Pctmam + LungCancerRate + poverty 
                                          + popdensity + medianhousevalue + pct_blk + medhouseholdincome + pct_owner_occ + hispanic + education + smoke_rate 
                                          + mean_bmi + amb_visit_pct + a1c_exm_pct + nearest_hospital_km + year2001 + year2002 + year2003 + year2004 + year2005 
                                          + year2006 + year2007 + year2008 + year2009 + year2010 + year2011 + year2012 + year2013 + year2014 + year2015 + year2016 
                                          + region5_2 + region5_3 + region5_4 + region5_5, data=DATA, family=quasipoisson) 
    results[((i-1)*6+6),c('type','corr_magnitude','pollutants')] <- c('Berkson','sp_error_3sd','multi')
    results[((i-1)*6+6),c('b0_est','b1_est')] <- mod_Berkson_multi_sp_error_3sd$coefficients[1:2]
    
    # clear memory
    DATA$y_Berkson_multi_sp_error_sd <- NULL
    DATA$y_Berkson_multi_sp_error_2sd <- NULL
    DATA$y_Berkson_multi_sp_error_3sd <- NULL
    rm(mod_Berkson_multi_sp_error_sd,mod_Berkson_multi_sp_error_2sd,mod_Berkson_multi_sp_error_3sd)
    gc()
  }
  cat(paste0("dataset is done \n"))
  
  # return the results 
  return(results)
  gc()
}



######################################### 2. test ################################################
# # read in covariates dataset
# covar <- readRDS(paste0(dir_data,'counts_weight_pred_sd_20211221.rds'))
# # create binary variables for calendar year
# covar$year2001 <- ifelse(covar$year==2001,1,0)
# covar$year2002 <- ifelse(covar$year==2002,1,0)
# covar$year2003 <- ifelse(covar$year==2003,1,0)
# covar$year2004 <- ifelse(covar$year==2004,1,0)
# covar$year2005 <- ifelse(covar$year==2005,1,0)
# covar$year2006 <- ifelse(covar$year==2006,1,0)
# covar$year2007 <- ifelse(covar$year==2007,1,0)
# covar$year2008 <- ifelse(covar$year==2008,1,0)
# covar$year2009 <- ifelse(covar$year==2009,1,0)
# covar$year2010 <- ifelse(covar$year==2010,1,0)
# covar$year2011 <- ifelse(covar$year==2011,1,0)
# covar$year2012 <- ifelse(covar$year==2012,1,0)
# covar$year2013 <- ifelse(covar$year==2013,1,0)
# covar$year2014 <- ifelse(covar$year==2014,1,0)
# covar$year2015 <- ifelse(covar$year==2015,1,0)
# covar$year2016 <- ifelse(covar$year==2016,1,0)
# # create binary variables for region5
# covar$region5_2 <- ifelse(covar$region5==2,1,0)
# covar$region5_3 <- ifelse(covar$region5==3,1,0)
# covar$region5_4 <- ifelse(covar$region5==4,1,0)
# covar$region5_5 <- ifelse(covar$region5==5,1,0)
# 
# # read in measurement error datasets
# key_files <- list.files(path=dir_data_sensitivity,pattern = "^counts_covar_simu_set_50_(.*)rds$")
# key <- readRDS(paste0(dir_data_sensitivity,key_files[1]))
# 
# # parameter setting
# b0 <- log(8)
# b1 <- 0.005
# n.reps <- 50  # 50 or 100, 200 (depending how many pieces we combine into a single loop simulation)
# 
# start.time <- Sys.time()
# test <- simulate_results(key,covar,b0,b1,n.reps)
# end.time <- Sys.time()
# end.time - start.time
# # Time difference of 20.60969 mins



######################################### 3. analysis ################################################
# read in covariates dataset
covar <- readRDS(paste0(dir_data,'counts_weight_pred_sd_20211221.rds'))
# create binary variables for calendar year
covar$year2001 <- ifelse(covar$year==2001,1,0)
covar$year2002 <- ifelse(covar$year==2002,1,0)
covar$year2003 <- ifelse(covar$year==2003,1,0)
covar$year2004 <- ifelse(covar$year==2004,1,0)
covar$year2005 <- ifelse(covar$year==2005,1,0)
covar$year2006 <- ifelse(covar$year==2006,1,0)
covar$year2007 <- ifelse(covar$year==2007,1,0)
covar$year2008 <- ifelse(covar$year==2008,1,0)
covar$year2009 <- ifelse(covar$year==2009,1,0)
covar$year2010 <- ifelse(covar$year==2010,1,0)
covar$year2011 <- ifelse(covar$year==2011,1,0)
covar$year2012 <- ifelse(covar$year==2012,1,0)
covar$year2013 <- ifelse(covar$year==2013,1,0)
covar$year2014 <- ifelse(covar$year==2014,1,0)
covar$year2015 <- ifelse(covar$year==2015,1,0)
covar$year2016 <- ifelse(covar$year==2016,1,0)
# create binary variables for region5
covar$region5_2 <- ifelse(covar$region5==2,1,0)
covar$region5_3 <- ifelse(covar$region5==3,1,0)
covar$region5_4 <- ifelse(covar$region5==4,1,0)
covar$region5_5 <- ifelse(covar$region5==5,1,0)

# read in measurement error datasets
key_files <- list.files(path=dir_data,pattern = "^counts_covar_simu_set_50_(.*)rds$")

# parameter setting
b0_list <- c(log(8),log(12),log(20))
b1_list <- c(0,0.005,0.012,0.019)
n.reps <- 50

cl = makeCluster(10,outfile='')    ### change number of clusters
registerDoParallel(cl)

tmp <- foreach(i=1:length(key_files))%dopar%{
  key <- readRDS(paste0(dir_data,key_files[i]))
  for (b0 in b0_list) {
    for (b1 in b1_list) {
      if (!file.exists(paste0(dir_results,'results_set',i,'_b0',b0,'_b1',b1,'.rds'))){
        results_temp <- simulate_results(key,covar,b0,b1,n.reps)
        saveRDS(results_temp,paste0(dir_results,'results_set',i,'_b0',b0,'_b1',b1,'.rds'))
        rm(results_temp)
        gc()
      }
    }
  }
  print(paste0("Key file ",i," is done"))
  rm(key)
  gc()
}

stopCluster(cl)
