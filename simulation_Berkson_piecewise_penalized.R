###########################################################################################################
# Project: EME                                                                                            #
# Code: simulate mortality with piecewise linear PM and test with penalized splines - Berkson error       #
# Author: Yaguang Wei                                                                                     #
###########################################################################################################

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

dir_data <- '/media/qnap4/Yaguang/EME/data/'
dir_results <- '/media/qnap4/Yaguang/EME/results/Berkson/piecewise_penalized/'  ### create directory to save results



######################################### 1. function ################################################
simulate_results<-function(key,covar,threshold,b1,n.reps){
  ##### the following codes are only for testing, comment out when running real analysis
  # covar<-readRDS(paste0(dir_data,'counts_weight_pred_sd_0127.rds'))
  # #create dummy variables for calendar year
  # covar$year2001=ifelse(covar$year==2001,1,0)
  # covar$year2002=ifelse(covar$year==2002,1,0)
  # covar$year2003=ifelse(covar$year==2003,1,0)
  # covar$year2004=ifelse(covar$year==2004,1,0)
  # covar$year2005=ifelse(covar$year==2005,1,0)
  # covar$year2006=ifelse(covar$year==2006,1,0)
  # covar$year2007=ifelse(covar$year==2007,1,0)
  # covar$year2008=ifelse(covar$year==2008,1,0)
  # covar$year2009=ifelse(covar$year==2009,1,0)
  # covar$year2010=ifelse(covar$year==2010,1,0)
  # covar$year2011=ifelse(covar$year==2011,1,0)
  # covar$year2012=ifelse(covar$year==2012,1,0)
  # covar$year2013=ifelse(covar$year==2013,1,0)
  # covar$year2014=ifelse(covar$year==2014,1,0)
  # covar$year2015=ifelse(covar$year==2015,1,0)
  # covar$year2016=ifelse(covar$year==2016,1,0)
  # 
  # #create dummy variables for region5
  # covar$region5_2=ifelse(covar$region5==2,1,0)
  # covar$region5_3=ifelse(covar$region5==3,1,0)
  # covar$region5_4=ifelse(covar$region5==4,1,0)
  # covar$region5_5=ifelse(covar$region5==5,1,0)
  # 
  # #remove sds
  # covar$sd_m1=covar$sd_m2=covar$sd_m3=covar$sd_m4=covar$sd_m5=covar$sd_m6=
  #   covar$sd_m7=covar$sd_m8=covar$sd_m9=covar$sd_m10=covar$sd_m11=covar$sd_m12=NULL
  # 
  # key_files <- list.files(path=dir_data,pattern = "^counts_covar_simu_set_50_(.*)rds$")
  # key <- readRDS(paste0(dir_data,key_files[1]))
  # threshold <- 8
  # b1=0.005
  # n.reps=50
  #####stop commenting for testing
  
  #dataframe to store resutls
  results <- data.frame(matrix(NA, nrow=n.reps*6, ncol=501))
  names(results) <- seq(0, 50, by=0.1)
  results$threshold <- NA
  results$b1 <- NA
  results$exp <- NA
  
  #run it in loops 
  for (i in 1:n.reps){
    # i=2
    cat(paste0("set ",i," starts running \n"))
    #create start/end index 
    start<-(i-1)*649910+1
    end<-649910*i
    
    #link covar with measurement error set 
    DATA<-merge(key[start:end,],covar,all.x=TRUE,by.x=c('zip','year'),by.y=c('zip','year'))
    #head(DATA,20)
    DATA<-na.omit(DATA) #remove the missing, 639926 obs left
    
    #compute new simulated exposure 
    DATA$exposure1_star<-DATA$pm25.x + DATA$annual_tempcorr_adj_err 
    DATA$exposure2_star<-DATA$pm25.x + DATA$annual_tempcorr_adj_err_2sd 
    DATA$exposure3_star<-DATA$pm25.x + DATA$annual_tempcorr_adj_err_3sd 
    DATA$exposure4_star<-DATA$pm25.x + DATA$annual_tempcorr_adj_err_sp
    DATA$exposure5_star<-DATA$pm25.x + DATA$annual_tempcorr_adj_err_2sd_sp
    DATA$exposure6_star<-DATA$pm25.x + DATA$annual_tempcorr_adj_err_3sd_sp
    
    #set IDs for each ob
    DATA$id <- 1:nrow(DATA)
    
    ###exposure1_star
    #inidcate obs with exposure below the threshold
    DATA$low_exp_1 <- 0
    DATA$low_exp_1[DATA$exposure1_star<=threshold] <- 1
    #generate the true counts for low-exposure obs
    DATA_low_exp_1 <- DATA[DATA$low_exp_1==1,]
    DATA_low_exp_1$y_1 <- rpois(nrow(DATA_low_exp_1),exp(4.795e-04*DATA_low_exp_1$ozone_summer+7.526e-04*DATA_low_exp_1$no2+6.693e-04*DATA_low_exp_1$temp+7.395e-04*DATA_low_exp_1$rh
                                                         -1.665e-03*DATA_low_exp_1$PctEye-1.015e-03*DATA_low_exp_1$PctLDL-1.525e-03*DATA_low_exp_1$Pctmam+6.110e-01*DATA_low_exp_1$LungCancerRate
                                                         +1.920e-01*DATA_low_exp_1$poverty-4.245e-06*DATA_low_exp_1$popdensity-2.574e-07*DATA_low_exp_1$medianhousevalue-4.967e-02*DATA_low_exp_1$pct_blk
                                                         -9.385e-07*DATA_low_exp_1$medhouseholdincome-2.910e-01*DATA_low_exp_1$pct_owner_occ-2.623e-01*DATA_low_exp_1$hispanic+1.448e-01*DATA_low_exp_1$education
                                                         +9.841e-02*DATA_low_exp_1$smoke_rate+6.779e-03*DATA_low_exp_1$mean_bmi+5.998e-04*DATA_low_exp_1$amb_visit_pct+4.921e-04*DATA_low_exp_1$a1c_exm_pct-3.336e-03*DATA_low_exp_1$nearest_hospital_km
                                                         -5.609e-03*DATA_low_exp_1$year2001-1.341e-03*DATA_low_exp_1$year2002-1.192e-02*DATA_low_exp_1$year2003-4.005e-02*DATA_low_exp_1$year2004
                                                         -2.362e-02*DATA_low_exp_1$year2005-1.901e-02*DATA_low_exp_1$year2006-2.230e-02*DATA_low_exp_1$year2007-1.388e-02*DATA_low_exp_1$year2008
                                                         -6.234e-02*DATA_low_exp_1$year2009-6.470e-02*DATA_low_exp_1$year2010-7.065e-02*DATA_low_exp_1$year2011-9.103e-02*DATA_low_exp_1$year2012
                                                         -1.029e-01*DATA_low_exp_1$year2013-1.207e-01*DATA_low_exp_1$year2014-1.134e-01*DATA_low_exp_1$year2015-1.329e-01*DATA_low_exp_1$year2016
                                                         -5.513e-02*DATA_low_exp_1$region5_2-2.685e-02*DATA_low_exp_1$region5_3-5.430e-02*DATA_low_exp_1$region5_4-3.998e-02*DATA_low_exp_1$region5_5)) 
    #generate the true counts for high-exposure obs
    DATA_high_exp_1 <- DATA[DATA$low_exp_1==0,]
    DATA_high_exp_1$y_1 <- rpois(nrow(DATA_high_exp_1),exp(-b1*threshold+b1*(DATA_high_exp_1$exposure1_star-threshold)+4.795e-04*DATA_high_exp_1$ozone_summer+7.526e-04*DATA_high_exp_1$no2+6.693e-04*DATA_high_exp_1$temp+7.395e-04*DATA_high_exp_1$rh
                                                           -1.665e-03*DATA_high_exp_1$PctEye-1.015e-03*DATA_high_exp_1$PctLDL-1.525e-03*DATA_high_exp_1$Pctmam+6.110e-01*DATA_high_exp_1$LungCancerRate
                                                           +1.920e-01*DATA_high_exp_1$poverty-4.245e-06*DATA_high_exp_1$popdensity-2.574e-07*DATA_high_exp_1$medianhousevalue-4.967e-02*DATA_high_exp_1$pct_blk
                                                           -9.385e-07*DATA_high_exp_1$medhouseholdincome-2.910e-01*DATA_high_exp_1$pct_owner_occ-2.623e-01*DATA_high_exp_1$hispanic+1.448e-01*DATA_high_exp_1$education
                                                           +9.841e-02*DATA_high_exp_1$smoke_rate+6.779e-03*DATA_high_exp_1$mean_bmi+5.998e-04*DATA_high_exp_1$amb_visit_pct+4.921e-04*DATA_high_exp_1$a1c_exm_pct-3.336e-03*DATA_high_exp_1$nearest_hospital_km
                                                           -5.609e-03*DATA_high_exp_1$year2001-1.341e-03*DATA_high_exp_1$year2002-1.192e-02*DATA_high_exp_1$year2003-4.005e-02*DATA_high_exp_1$year2004
                                                           -2.362e-02*DATA_high_exp_1$year2005-1.901e-02*DATA_high_exp_1$year2006-2.230e-02*DATA_high_exp_1$year2007-1.388e-02*DATA_high_exp_1$year2008
                                                           -6.234e-02*DATA_high_exp_1$year2009-6.470e-02*DATA_high_exp_1$year2010-7.065e-02*DATA_high_exp_1$year2011-9.103e-02*DATA_high_exp_1$year2012
                                                           -1.029e-01*DATA_high_exp_1$year2013-1.207e-01*DATA_high_exp_1$year2014-1.134e-01*DATA_high_exp_1$year2015-1.329e-01*DATA_high_exp_1$year2016
                                                           -5.513e-02*DATA_high_exp_1$region5_2-2.685e-02*DATA_high_exp_1$region5_3-5.430e-02*DATA_high_exp_1$region5_4-3.998e-02*DATA_high_exp_1$region5_5)) 
    #combine and sort
    DATA <- rbind(DATA_low_exp_1,DATA_high_exp_1)
    DATA <- DATA[order(DATA$id),]
    #clear junk
    rm(DATA_low_exp_1,DATA_high_exp_1)
    
    ###exposure2_star
    #inidcate obs with exposure below the threshold
    DATA$low_exp_2 <- 0
    DATA$low_exp_2[DATA$exposure2_star<=threshold] <- 1
    #generate the true counts for low-exposure obs
    DATA_low_exp_2 <- DATA[DATA$low_exp_2==1,]
    DATA_low_exp_2$y_2 <- rpois(nrow(DATA_low_exp_2),exp(4.795e-04*DATA_low_exp_2$ozone_summer+7.526e-04*DATA_low_exp_2$no2+6.693e-04*DATA_low_exp_2$temp+7.395e-04*DATA_low_exp_2$rh
                                                         -1.665e-03*DATA_low_exp_2$PctEye-1.015e-03*DATA_low_exp_2$PctLDL-1.525e-03*DATA_low_exp_2$Pctmam+6.110e-01*DATA_low_exp_2$LungCancerRate
                                                         +1.920e-01*DATA_low_exp_2$poverty-4.245e-06*DATA_low_exp_2$popdensity-2.574e-07*DATA_low_exp_2$medianhousevalue-4.967e-02*DATA_low_exp_2$pct_blk
                                                         -9.385e-07*DATA_low_exp_2$medhouseholdincome-2.910e-01*DATA_low_exp_2$pct_owner_occ-2.623e-01*DATA_low_exp_2$hispanic+1.448e-01*DATA_low_exp_2$education
                                                         +9.841e-02*DATA_low_exp_2$smoke_rate+6.779e-03*DATA_low_exp_2$mean_bmi+5.998e-04*DATA_low_exp_2$amb_visit_pct+4.921e-04*DATA_low_exp_2$a1c_exm_pct-3.336e-03*DATA_low_exp_2$nearest_hospital_km
                                                         -5.609e-03*DATA_low_exp_2$year2001-1.341e-03*DATA_low_exp_2$year2002-1.192e-02*DATA_low_exp_2$year2003-4.005e-02*DATA_low_exp_2$year2004
                                                         -2.362e-02*DATA_low_exp_2$year2005-1.901e-02*DATA_low_exp_2$year2006-2.230e-02*DATA_low_exp_2$year2007-1.388e-02*DATA_low_exp_2$year2008
                                                         -6.234e-02*DATA_low_exp_2$year2009-6.470e-02*DATA_low_exp_2$year2010-7.065e-02*DATA_low_exp_2$year2011-9.103e-02*DATA_low_exp_2$year2012
                                                         -1.029e-01*DATA_low_exp_2$year2013-1.207e-01*DATA_low_exp_2$year2014-1.134e-01*DATA_low_exp_2$year2015-1.329e-01*DATA_low_exp_2$year2016
                                                         -5.513e-02*DATA_low_exp_2$region5_2-2.685e-02*DATA_low_exp_2$region5_3-5.430e-02*DATA_low_exp_2$region5_4-3.998e-02*DATA_low_exp_2$region5_5)) 
    #generate the true counts for high-exposure obs
    DATA_high_exp_2 <- DATA[DATA$low_exp_2==0,]
    DATA_high_exp_2$y_2 <- rpois(nrow(DATA_high_exp_2),exp(-b1*threshold+b1*(DATA_high_exp_2$exposure2_star-threshold)+4.795e-04*DATA_high_exp_2$ozone_summer+7.526e-04*DATA_high_exp_2$no2+6.693e-04*DATA_high_exp_2$temp+7.395e-04*DATA_high_exp_2$rh
                                                           -1.665e-03*DATA_high_exp_2$PctEye-1.015e-03*DATA_high_exp_2$PctLDL-1.525e-03*DATA_high_exp_2$Pctmam+6.110e-01*DATA_high_exp_2$LungCancerRate
                                                           +1.920e-01*DATA_high_exp_2$poverty-4.245e-06*DATA_high_exp_2$popdensity-2.574e-07*DATA_high_exp_2$medianhousevalue-4.967e-02*DATA_high_exp_2$pct_blk
                                                           -9.385e-07*DATA_high_exp_2$medhouseholdincome-2.910e-01*DATA_high_exp_2$pct_owner_occ-2.623e-01*DATA_high_exp_2$hispanic+1.448e-01*DATA_high_exp_2$education
                                                           +9.841e-02*DATA_high_exp_2$smoke_rate+6.779e-03*DATA_high_exp_2$mean_bmi+5.998e-04*DATA_high_exp_2$amb_visit_pct+4.921e-04*DATA_high_exp_2$a1c_exm_pct-3.336e-03*DATA_high_exp_2$nearest_hospital_km
                                                           -5.609e-03*DATA_high_exp_2$year2001-1.341e-03*DATA_high_exp_2$year2002-1.192e-02*DATA_high_exp_2$year2003-4.005e-02*DATA_high_exp_2$year2004
                                                           -2.362e-02*DATA_high_exp_2$year2005-1.901e-02*DATA_high_exp_2$year2006-2.230e-02*DATA_high_exp_2$year2007-1.388e-02*DATA_high_exp_2$year2008
                                                           -6.234e-02*DATA_high_exp_2$year2009-6.470e-02*DATA_high_exp_2$year2010-7.065e-02*DATA_high_exp_2$year2011-9.103e-02*DATA_high_exp_2$year2012
                                                           -1.029e-01*DATA_high_exp_2$year2013-1.207e-01*DATA_high_exp_2$year2014-1.134e-01*DATA_high_exp_2$year2015-1.329e-01*DATA_high_exp_2$year2016
                                                           -5.513e-02*DATA_high_exp_2$region5_2-2.685e-02*DATA_high_exp_2$region5_3-5.430e-02*DATA_high_exp_2$region5_4-3.998e-02*DATA_high_exp_2$region5_5)) 
    #combine and sort
    DATA <- rbind(DATA_low_exp_2,DATA_high_exp_2)
    DATA <- DATA[order(DATA$id),]
    #clear junk
    rm(DATA_low_exp_2,DATA_high_exp_2)
    
    ###exposure3_star
    #inidcate obs with exposure below the threshold
    DATA$low_exp_3 <- 0
    DATA$low_exp_3[DATA$exposure3_star<=threshold] <- 1
    #generate the true counts for low-exposure obs
    DATA_low_exp_3 <- DATA[DATA$low_exp_3==1,]
    DATA_low_exp_3$y_3 <- rpois(nrow(DATA_low_exp_3),exp(4.795e-04*DATA_low_exp_3$ozone_summer+7.526e-04*DATA_low_exp_3$no2+6.693e-04*DATA_low_exp_3$temp+7.395e-04*DATA_low_exp_3$rh
                                                         -1.665e-03*DATA_low_exp_3$PctEye-1.015e-03*DATA_low_exp_3$PctLDL-1.525e-03*DATA_low_exp_3$Pctmam+6.110e-01*DATA_low_exp_3$LungCancerRate
                                                         +1.920e-01*DATA_low_exp_3$poverty-4.245e-06*DATA_low_exp_3$popdensity-2.574e-07*DATA_low_exp_3$medianhousevalue-4.967e-02*DATA_low_exp_3$pct_blk
                                                         -9.385e-07*DATA_low_exp_3$medhouseholdincome-2.910e-01*DATA_low_exp_3$pct_owner_occ-2.623e-01*DATA_low_exp_3$hispanic+1.448e-01*DATA_low_exp_3$education
                                                         +9.841e-02*DATA_low_exp_3$smoke_rate+6.779e-03*DATA_low_exp_3$mean_bmi+5.998e-04*DATA_low_exp_3$amb_visit_pct+4.921e-04*DATA_low_exp_3$a1c_exm_pct-3.336e-03*DATA_low_exp_3$nearest_hospital_km
                                                         -5.609e-03*DATA_low_exp_3$year2001-1.341e-03*DATA_low_exp_3$year2002-1.192e-02*DATA_low_exp_3$year2003-4.005e-02*DATA_low_exp_3$year2004
                                                         -2.362e-02*DATA_low_exp_3$year2005-1.901e-02*DATA_low_exp_3$year2006-2.230e-02*DATA_low_exp_3$year2007-1.388e-02*DATA_low_exp_3$year2008
                                                         -6.234e-02*DATA_low_exp_3$year2009-6.470e-02*DATA_low_exp_3$year2010-7.065e-02*DATA_low_exp_3$year2011-9.103e-02*DATA_low_exp_3$year2012
                                                         -1.029e-01*DATA_low_exp_3$year2013-1.207e-01*DATA_low_exp_3$year2014-1.134e-01*DATA_low_exp_3$year2015-1.329e-01*DATA_low_exp_3$year2016
                                                         -5.513e-02*DATA_low_exp_3$region5_2-2.685e-02*DATA_low_exp_3$region5_3-5.430e-02*DATA_low_exp_3$region5_4-3.998e-02*DATA_low_exp_3$region5_5)) 
    #generate the true counts for high-exposure obs
    DATA_high_exp_3 <- DATA[DATA$low_exp_3==0,]
    DATA_high_exp_3$y_3 <- rpois(nrow(DATA_high_exp_3),exp(-b1*threshold+b1*(DATA_high_exp_3$exposure3_star-threshold)+4.795e-04*DATA_high_exp_3$ozone_summer+7.526e-04*DATA_high_exp_3$no2+6.693e-04*DATA_high_exp_3$temp+7.395e-04*DATA_high_exp_3$rh
                                                           -1.665e-03*DATA_high_exp_3$PctEye-1.015e-03*DATA_high_exp_3$PctLDL-1.525e-03*DATA_high_exp_3$Pctmam+6.110e-01*DATA_high_exp_3$LungCancerRate
                                                           +1.920e-01*DATA_high_exp_3$poverty-4.245e-06*DATA_high_exp_3$popdensity-2.574e-07*DATA_high_exp_3$medianhousevalue-4.967e-02*DATA_high_exp_3$pct_blk
                                                           -9.385e-07*DATA_high_exp_3$medhouseholdincome-2.910e-01*DATA_high_exp_3$pct_owner_occ-2.623e-01*DATA_high_exp_3$hispanic+1.448e-01*DATA_high_exp_3$education
                                                           +9.841e-02*DATA_high_exp_3$smoke_rate+6.779e-03*DATA_high_exp_3$mean_bmi+5.998e-04*DATA_high_exp_3$amb_visit_pct+4.921e-04*DATA_high_exp_3$a1c_exm_pct-3.336e-03*DATA_high_exp_3$nearest_hospital_km
                                                           -5.609e-03*DATA_high_exp_3$year2001-1.341e-03*DATA_high_exp_3$year2002-1.192e-02*DATA_high_exp_3$year2003-4.005e-02*DATA_high_exp_3$year2004
                                                           -2.362e-02*DATA_high_exp_3$year2005-1.901e-02*DATA_high_exp_3$year2006-2.230e-02*DATA_high_exp_3$year2007-1.388e-02*DATA_high_exp_3$year2008
                                                           -6.234e-02*DATA_high_exp_3$year2009-6.470e-02*DATA_high_exp_3$year2010-7.065e-02*DATA_high_exp_3$year2011-9.103e-02*DATA_high_exp_3$year2012
                                                           -1.029e-01*DATA_high_exp_3$year2013-1.207e-01*DATA_high_exp_3$year2014-1.134e-01*DATA_high_exp_3$year2015-1.329e-01*DATA_high_exp_3$year2016
                                                           -5.513e-02*DATA_high_exp_3$region5_2-2.685e-02*DATA_high_exp_3$region5_3-5.430e-02*DATA_high_exp_3$region5_4-3.998e-02*DATA_high_exp_3$region5_5)) 
    #combine and sort
    DATA <- rbind(DATA_low_exp_3,DATA_high_exp_3)
    DATA <- DATA[order(DATA$id),]
    #clear junk
    rm(DATA_low_exp_3,DATA_high_exp_3)
    
    ###exposure4_star
    #inidcate obs with exposure below the threshold
    DATA$low_exp_4 <- 0
    DATA$low_exp_4[DATA$exposure4_star<=threshold] <- 1
    #generate the true counts for low-exposure obs
    DATA_low_exp_4 <- DATA[DATA$low_exp_4==1,]
    DATA_low_exp_4$y_4 <- rpois(nrow(DATA_low_exp_4),exp(4.795e-04*DATA_low_exp_4$ozone_summer+7.526e-04*DATA_low_exp_4$no2+6.693e-04*DATA_low_exp_4$temp+7.395e-04*DATA_low_exp_4$rh
                                                         -1.665e-03*DATA_low_exp_4$PctEye-1.015e-03*DATA_low_exp_4$PctLDL-1.525e-03*DATA_low_exp_4$Pctmam+6.110e-01*DATA_low_exp_4$LungCancerRate
                                                         +1.920e-01*DATA_low_exp_4$poverty-4.245e-06*DATA_low_exp_4$popdensity-2.574e-07*DATA_low_exp_4$medianhousevalue-4.967e-02*DATA_low_exp_4$pct_blk
                                                         -9.385e-07*DATA_low_exp_4$medhouseholdincome-2.910e-01*DATA_low_exp_4$pct_owner_occ-2.623e-01*DATA_low_exp_4$hispanic+1.448e-01*DATA_low_exp_4$education
                                                         +9.841e-02*DATA_low_exp_4$smoke_rate+6.779e-03*DATA_low_exp_4$mean_bmi+5.998e-04*DATA_low_exp_4$amb_visit_pct+4.921e-04*DATA_low_exp_4$a1c_exm_pct-3.336e-03*DATA_low_exp_4$nearest_hospital_km
                                                         -5.609e-03*DATA_low_exp_4$year2001-1.341e-03*DATA_low_exp_4$year2002-1.192e-02*DATA_low_exp_4$year2003-4.005e-02*DATA_low_exp_4$year2004
                                                         -2.362e-02*DATA_low_exp_4$year2005-1.901e-02*DATA_low_exp_4$year2006-2.230e-02*DATA_low_exp_4$year2007-1.388e-02*DATA_low_exp_4$year2008
                                                         -6.234e-02*DATA_low_exp_4$year2009-6.470e-02*DATA_low_exp_4$year2010-7.065e-02*DATA_low_exp_4$year2011-9.103e-02*DATA_low_exp_4$year2012
                                                         -1.029e-01*DATA_low_exp_4$year2013-1.207e-01*DATA_low_exp_4$year2014-1.134e-01*DATA_low_exp_4$year2015-1.329e-01*DATA_low_exp_4$year2016
                                                         -5.513e-02*DATA_low_exp_4$region5_2-2.685e-02*DATA_low_exp_4$region5_3-5.430e-02*DATA_low_exp_4$region5_4-3.998e-02*DATA_low_exp_4$region5_5)) 
    #generate the true counts for high-exposure obs
    DATA_high_exp_4 <- DATA[DATA$low_exp_4==0,]
    DATA_high_exp_4$y_4 <- rpois(nrow(DATA_high_exp_4),exp(-b1*threshold+b1*(DATA_high_exp_4$exposure4_star-threshold)+4.795e-04*DATA_high_exp_4$ozone_summer+7.526e-04*DATA_high_exp_4$no2+6.693e-04*DATA_high_exp_4$temp+7.395e-04*DATA_high_exp_4$rh
                                                           -1.665e-03*DATA_high_exp_4$PctEye-1.015e-03*DATA_high_exp_4$PctLDL-1.525e-03*DATA_high_exp_4$Pctmam+6.110e-01*DATA_high_exp_4$LungCancerRate
                                                           +1.920e-01*DATA_high_exp_4$poverty-4.245e-06*DATA_high_exp_4$popdensity-2.574e-07*DATA_high_exp_4$medianhousevalue-4.967e-02*DATA_high_exp_4$pct_blk
                                                           -9.385e-07*DATA_high_exp_4$medhouseholdincome-2.910e-01*DATA_high_exp_4$pct_owner_occ-2.623e-01*DATA_high_exp_4$hispanic+1.448e-01*DATA_high_exp_4$education
                                                           +9.841e-02*DATA_high_exp_4$smoke_rate+6.779e-03*DATA_high_exp_4$mean_bmi+5.998e-04*DATA_high_exp_4$amb_visit_pct+4.921e-04*DATA_high_exp_4$a1c_exm_pct-3.336e-03*DATA_high_exp_4$nearest_hospital_km
                                                           -5.609e-03*DATA_high_exp_4$year2001-1.341e-03*DATA_high_exp_4$year2002-1.192e-02*DATA_high_exp_4$year2003-4.005e-02*DATA_high_exp_4$year2004
                                                           -2.362e-02*DATA_high_exp_4$year2005-1.901e-02*DATA_high_exp_4$year2006-2.230e-02*DATA_high_exp_4$year2007-1.388e-02*DATA_high_exp_4$year2008
                                                           -6.234e-02*DATA_high_exp_4$year2009-6.470e-02*DATA_high_exp_4$year2010-7.065e-02*DATA_high_exp_4$year2011-9.103e-02*DATA_high_exp_4$year2012
                                                           -1.029e-01*DATA_high_exp_4$year2013-1.207e-01*DATA_high_exp_4$year2014-1.134e-01*DATA_high_exp_4$year2015-1.329e-01*DATA_high_exp_4$year2016
                                                           -5.513e-02*DATA_high_exp_4$region5_2-2.685e-02*DATA_high_exp_4$region5_3-5.430e-02*DATA_high_exp_4$region5_4-3.998e-02*DATA_high_exp_4$region5_5)) 
    #combine and sort
    DATA <- rbind(DATA_low_exp_4,DATA_high_exp_4)
    DATA <- DATA[order(DATA$id),]
    #clear junk
    rm(DATA_low_exp_4,DATA_high_exp_4)
    
    ###exposure5_star
    #inidcate obs with exposure below the threshold
    DATA$low_exp_5 <- 0
    DATA$low_exp_5[DATA$exposure5_star<=threshold] <- 1
    #generate the true counts for low-exposure obs
    DATA_low_exp_5 <- DATA[DATA$low_exp_5==1,]
    DATA_low_exp_5$y_5 <- rpois(nrow(DATA_low_exp_5),exp(4.795e-04*DATA_low_exp_5$ozone_summer+7.526e-04*DATA_low_exp_5$no2+6.693e-04*DATA_low_exp_5$temp+7.395e-04*DATA_low_exp_5$rh
                                                         -1.665e-03*DATA_low_exp_5$PctEye-1.015e-03*DATA_low_exp_5$PctLDL-1.525e-03*DATA_low_exp_5$Pctmam+6.110e-01*DATA_low_exp_5$LungCancerRate
                                                         +1.920e-01*DATA_low_exp_5$poverty-4.245e-06*DATA_low_exp_5$popdensity-2.574e-07*DATA_low_exp_5$medianhousevalue-4.967e-02*DATA_low_exp_5$pct_blk
                                                         -9.385e-07*DATA_low_exp_5$medhouseholdincome-2.910e-01*DATA_low_exp_5$pct_owner_occ-2.623e-01*DATA_low_exp_5$hispanic+1.448e-01*DATA_low_exp_5$education
                                                         +9.841e-02*DATA_low_exp_5$smoke_rate+6.779e-03*DATA_low_exp_5$mean_bmi+5.998e-04*DATA_low_exp_5$amb_visit_pct+4.921e-04*DATA_low_exp_5$a1c_exm_pct-3.336e-03*DATA_low_exp_5$nearest_hospital_km
                                                         -5.609e-03*DATA_low_exp_5$year2001-1.341e-03*DATA_low_exp_5$year2002-1.192e-02*DATA_low_exp_5$year2003-4.005e-02*DATA_low_exp_5$year2004
                                                         -2.362e-02*DATA_low_exp_5$year2005-1.901e-02*DATA_low_exp_5$year2006-2.230e-02*DATA_low_exp_5$year2007-1.388e-02*DATA_low_exp_5$year2008
                                                         -6.234e-02*DATA_low_exp_5$year2009-6.470e-02*DATA_low_exp_5$year2010-7.065e-02*DATA_low_exp_5$year2011-9.103e-02*DATA_low_exp_5$year2012
                                                         -1.029e-01*DATA_low_exp_5$year2013-1.207e-01*DATA_low_exp_5$year2014-1.134e-01*DATA_low_exp_5$year2015-1.329e-01*DATA_low_exp_5$year2016
                                                         -5.513e-02*DATA_low_exp_5$region5_2-2.685e-02*DATA_low_exp_5$region5_3-5.430e-02*DATA_low_exp_5$region5_4-3.998e-02*DATA_low_exp_5$region5_5)) 
    #generate the true counts for high-exposure obs
    DATA_high_exp_5 <- DATA[DATA$low_exp_5==0,]
    DATA_high_exp_5$y_5 <- rpois(nrow(DATA_high_exp_5),exp(-b1*threshold+b1*(DATA_high_exp_5$exposure5_star-threshold)+4.795e-04*DATA_high_exp_5$ozone_summer+7.526e-04*DATA_high_exp_5$no2+6.693e-04*DATA_high_exp_5$temp+7.395e-04*DATA_high_exp_5$rh
                                                           -1.665e-03*DATA_high_exp_5$PctEye-1.015e-03*DATA_high_exp_5$PctLDL-1.525e-03*DATA_high_exp_5$Pctmam+6.110e-01*DATA_high_exp_5$LungCancerRate
                                                           +1.920e-01*DATA_high_exp_5$poverty-4.245e-06*DATA_high_exp_5$popdensity-2.574e-07*DATA_high_exp_5$medianhousevalue-4.967e-02*DATA_high_exp_5$pct_blk
                                                           -9.385e-07*DATA_high_exp_5$medhouseholdincome-2.910e-01*DATA_high_exp_5$pct_owner_occ-2.623e-01*DATA_high_exp_5$hispanic+1.448e-01*DATA_high_exp_5$education
                                                           +9.841e-02*DATA_high_exp_5$smoke_rate+6.779e-03*DATA_high_exp_5$mean_bmi+5.998e-04*DATA_high_exp_5$amb_visit_pct+4.921e-04*DATA_high_exp_5$a1c_exm_pct-3.336e-03*DATA_high_exp_5$nearest_hospital_km
                                                           -5.609e-03*DATA_high_exp_5$year2001-1.341e-03*DATA_high_exp_5$year2002-1.192e-02*DATA_high_exp_5$year2003-4.005e-02*DATA_high_exp_5$year2004
                                                           -2.362e-02*DATA_high_exp_5$year2005-1.901e-02*DATA_high_exp_5$year2006-2.230e-02*DATA_high_exp_5$year2007-1.388e-02*DATA_high_exp_5$year2008
                                                           -6.234e-02*DATA_high_exp_5$year2009-6.470e-02*DATA_high_exp_5$year2010-7.065e-02*DATA_high_exp_5$year2011-9.103e-02*DATA_high_exp_5$year2012
                                                           -1.029e-01*DATA_high_exp_5$year2013-1.207e-01*DATA_high_exp_5$year2014-1.134e-01*DATA_high_exp_5$year2015-1.329e-01*DATA_high_exp_5$year2016
                                                           -5.513e-02*DATA_high_exp_5$region5_2-2.685e-02*DATA_high_exp_5$region5_3-5.430e-02*DATA_high_exp_5$region5_4-3.998e-02*DATA_high_exp_5$region5_5)) 
    #combine and sort
    DATA <- rbind(DATA_low_exp_5,DATA_high_exp_5)
    DATA <- DATA[order(DATA$id),]
    #clear junk
    rm(DATA_low_exp_5,DATA_high_exp_5)
    
    ###exposure6_star
    #inidcate obs with exposure below the threshold
    DATA$low_exp_6 <- 0
    DATA$low_exp_6[DATA$exposure6_star<=threshold] <- 1
    #generate the true counts for low-exposure obs
    DATA_low_exp_6 <- DATA[DATA$low_exp_6==1,]
    DATA_low_exp_6$y_6 <- rpois(nrow(DATA_low_exp_6),exp(4.795e-04*DATA_low_exp_6$ozone_summer+7.526e-04*DATA_low_exp_6$no2+6.693e-04*DATA_low_exp_6$temp+7.395e-04*DATA_low_exp_6$rh
                                                         -1.665e-03*DATA_low_exp_6$PctEye-1.015e-03*DATA_low_exp_6$PctLDL-1.525e-03*DATA_low_exp_6$Pctmam+6.110e-01*DATA_low_exp_6$LungCancerRate
                                                         +1.920e-01*DATA_low_exp_6$poverty-4.245e-06*DATA_low_exp_6$popdensity-2.574e-07*DATA_low_exp_6$medianhousevalue-4.967e-02*DATA_low_exp_6$pct_blk
                                                         -9.385e-07*DATA_low_exp_6$medhouseholdincome-2.910e-01*DATA_low_exp_6$pct_owner_occ-2.623e-01*DATA_low_exp_6$hispanic+1.448e-01*DATA_low_exp_6$education
                                                         +9.841e-02*DATA_low_exp_6$smoke_rate+6.779e-03*DATA_low_exp_6$mean_bmi+5.998e-04*DATA_low_exp_6$amb_visit_pct+4.921e-04*DATA_low_exp_6$a1c_exm_pct-3.336e-03*DATA_low_exp_6$nearest_hospital_km
                                                         -5.609e-03*DATA_low_exp_6$year2001-1.341e-03*DATA_low_exp_6$year2002-1.192e-02*DATA_low_exp_6$year2003-4.005e-02*DATA_low_exp_6$year2004
                                                         -2.362e-02*DATA_low_exp_6$year2005-1.901e-02*DATA_low_exp_6$year2006-2.230e-02*DATA_low_exp_6$year2007-1.388e-02*DATA_low_exp_6$year2008
                                                         -6.234e-02*DATA_low_exp_6$year2009-6.470e-02*DATA_low_exp_6$year2010-7.065e-02*DATA_low_exp_6$year2011-9.103e-02*DATA_low_exp_6$year2012
                                                         -1.029e-01*DATA_low_exp_6$year2013-1.207e-01*DATA_low_exp_6$year2014-1.134e-01*DATA_low_exp_6$year2015-1.329e-01*DATA_low_exp_6$year2016
                                                         -5.513e-02*DATA_low_exp_6$region5_2-2.685e-02*DATA_low_exp_6$region5_3-5.430e-02*DATA_low_exp_6$region5_4-3.998e-02*DATA_low_exp_6$region5_5)) 
    #generate the true counts for high-exposure obs
    DATA_high_exp_6 <- DATA[DATA$low_exp_6==0,]
    DATA_high_exp_6$y_6 <- rpois(nrow(DATA_high_exp_6),exp(-b1*threshold+b1*(DATA_high_exp_6$exposure6_star-threshold)+4.795e-04*DATA_high_exp_6$ozone_summer+7.526e-04*DATA_high_exp_6$no2+6.693e-04*DATA_high_exp_6$temp+7.395e-04*DATA_high_exp_6$rh
                                                           -1.665e-03*DATA_high_exp_6$PctEye-1.015e-03*DATA_high_exp_6$PctLDL-1.525e-03*DATA_high_exp_6$Pctmam+6.110e-01*DATA_high_exp_6$LungCancerRate
                                                           +1.920e-01*DATA_high_exp_6$poverty-4.245e-06*DATA_high_exp_6$popdensity-2.574e-07*DATA_high_exp_6$medianhousevalue-4.967e-02*DATA_high_exp_6$pct_blk
                                                           -9.385e-07*DATA_high_exp_6$medhouseholdincome-2.910e-01*DATA_high_exp_6$pct_owner_occ-2.623e-01*DATA_high_exp_6$hispanic+1.448e-01*DATA_high_exp_6$education
                                                           +9.841e-02*DATA_high_exp_6$smoke_rate+6.779e-03*DATA_high_exp_6$mean_bmi+5.998e-04*DATA_high_exp_6$amb_visit_pct+4.921e-04*DATA_high_exp_6$a1c_exm_pct-3.336e-03*DATA_high_exp_6$nearest_hospital_km
                                                           -5.609e-03*DATA_high_exp_6$year2001-1.341e-03*DATA_high_exp_6$year2002-1.192e-02*DATA_high_exp_6$year2003-4.005e-02*DATA_high_exp_6$year2004
                                                           -2.362e-02*DATA_high_exp_6$year2005-1.901e-02*DATA_high_exp_6$year2006-2.230e-02*DATA_high_exp_6$year2007-1.388e-02*DATA_high_exp_6$year2008
                                                           -6.234e-02*DATA_high_exp_6$year2009-6.470e-02*DATA_high_exp_6$year2010-7.065e-02*DATA_high_exp_6$year2011-9.103e-02*DATA_high_exp_6$year2012
                                                           -1.029e-01*DATA_high_exp_6$year2013-1.207e-01*DATA_high_exp_6$year2014-1.134e-01*DATA_high_exp_6$year2015-1.329e-01*DATA_high_exp_6$year2016
                                                           -5.513e-02*DATA_high_exp_6$region5_2-2.685e-02*DATA_high_exp_6$region5_3-5.430e-02*DATA_high_exp_6$region5_4-3.998e-02*DATA_high_exp_6$region5_5)) 
    #combine and sort
    DATA <- rbind(DATA_low_exp_6,DATA_high_exp_6)
    DATA <- DATA[order(DATA$id),]
    #clear junk
    rm(DATA_low_exp_6,DATA_high_exp_6)
    
    #assign covar value 
    ozone_summer<-DATA$ozone_summer
    no2<-DATA$no2
    temp<-DATA$temp
    rh<-DATA$rh
    PctEye<-DATA$PctEye
    PctLDL<-DATA$PctLDL
    Pctmam<-DATA$Pctmam
    LungCancerRate<-DATA$LungCancerRate
    poverty<-DATA$poverty
    popdensity<-DATA$popdensity
    medianhousevalue<-DATA$medianhousevalue
    pct_blk<-DATA$pct_blk
    medhouseholdincome<-DATA$medhouseholdincome
    pct_owner_occ<-DATA$pct_owner_occ
    hispanic<-DATA$hispanic
    education<-DATA$education
    smoke_rate<-DATA$smoke_rate
    mean_bmi<-DATA$mean_bmi
    amb_visit_pct<-DATA$amb_visit_pct
    a1c_exm_pct<-DATA$a1c_exm_pct
    nearest_hospital_km<-DATA$nearest_hospital_km
    year2001<-DATA$year2001
    year2002<-DATA$year2002
    year2003<-DATA$year2003
    year2004<-DATA$year2004
    year2005<-DATA$year2005
    year2006<-DATA$year2006
    year2007<-DATA$year2007
    year2008<-DATA$year2008
    year2009<-DATA$year2009
    year2010<-DATA$year2010
    year2011<-DATA$year2011
    year2012<-DATA$year2012
    year2013<-DATA$year2013
    year2014<-DATA$year2014
    year2015<-DATA$year2015
    year2016<-DATA$year2016
    region5_2<-DATA$region5_2
    region5_3<-DATA$region5_3
    region5_4<-DATA$region5_4
    region5_5<-DATA$region5_5
    
    #simulated outcomes
    y_1<-DATA$y_1
    y_2<-DATA$y_2
    y_3<-DATA$y_3
    y_4<-DATA$y_4
    y_5<-DATA$y_5
    y_6<-DATA$y_6
    
    #error-prone exposures 
    exposure<-DATA$pm25.x
    
    #build new dataset to predict rr
    newdata_try_1 <- data.frame(exposure=seq(0,50,by=0.1))
    newdata_try_2 <- data.frame(matrix(0, nrow=length(seq(0,50,by=0.1)), ncol=42))
    names(newdata_try_2) <- c('ozone_summer','no2','temp','rh','PctEye','PctLDL','Pctmam','LungCancerRate','poverty','popdensity','medianhousevalue',
                               'pct_blk','medhouseholdincome','pct_owner_occ','hispanic','education','smoke_rate','mean_bmi','amb_visit_pct',
                               'a1c_exm_pct','nearest_hospital_km','amb_visit_pct','year2001','year2002','year2003','year2004','year2005','year2006',
                               'year2007','year2008','year2009','year2010','year2011','year2012','year2013','year2014','year2015','year2016',
                               'region5_2','region5_3','region5_4','region5_5')
    newdata_try <- cbind(newdata_try_1,newdata_try_2)
    
    #run the Poisson models 
    try1_bam<-bam(y_1 ~ s(exposure,bs='cr') + ozone_summer + no2 + temp + rh + PctEye + PctLDL + Pctmam + LungCancerRate + poverty + popdensity + medianhousevalue 
                  + pct_blk + medhouseholdincome + pct_owner_occ + hispanic + education + smoke_rate + mean_bmi + amb_visit_pct + a1c_exm_pct + nearest_hospital_km 
                  + year2001 + year2002 + year2003 + year2004 + year2005 + year2006 + year2007 + year2008 
                  + year2009 + year2010 + year2011 + year2012 + year2013 + year2014 + year2015 + year2016 
                  + region5_2 + region5_3 + region5_4 + region5_5, family=quasipoisson,method="fREML") 
    predict_try1_bam <- data.frame(predict(try1_bam, newdata = newdata_try,type="terms", se.fit = TRUE))
    
    try2_bam<-bam(y_2 ~ s(exposure,bs='cr') + ozone_summer + no2 + temp + rh + PctEye + PctLDL + Pctmam + LungCancerRate + poverty + popdensity + medianhousevalue 
                  + pct_blk + medhouseholdincome + pct_owner_occ + hispanic + education + smoke_rate + mean_bmi + amb_visit_pct + a1c_exm_pct + nearest_hospital_km 
                  + year2001 + year2002 + year2003 + year2004 + year2005 + year2006 + year2007 + year2008 
                  + year2009 + year2010 + year2011 + year2012 + year2013 + year2014 + year2015 + year2016 
                  + region5_2 + region5_3 + region5_4 + region5_5, family=quasipoisson,method="fREML") 
    predict_try2_bam <- data.frame(predict(try2_bam, newdata = newdata_try,type="terms", se.fit = TRUE))
    
    try3_bam<-bam(y_3 ~ s(exposure,bs='cr') + ozone_summer + no2 + temp + rh + PctEye + PctLDL + Pctmam + LungCancerRate + poverty + popdensity + medianhousevalue 
                  + pct_blk + medhouseholdincome + pct_owner_occ + hispanic + education + smoke_rate + mean_bmi + amb_visit_pct + a1c_exm_pct + nearest_hospital_km 
                  + year2001 + year2002 + year2003 + year2004 + year2005 + year2006 + year2007 + year2008 
                  + year2009 + year2010 + year2011 + year2012 + year2013 + year2014 + year2015 + year2016 
                  + region5_2 + region5_3 + region5_4 + region5_5, family=quasipoisson,method="fREML") 
    predict_try3_bam <- data.frame(predict(try3_bam, newdata = newdata_try,type="terms", se.fit = TRUE))
    
    try4_bam<-bam(y_4 ~ s(exposure,bs='cr') + ozone_summer + no2 + temp + rh + PctEye + PctLDL + Pctmam + LungCancerRate + poverty + popdensity + medianhousevalue 
                  + pct_blk + medhouseholdincome + pct_owner_occ + hispanic + education + smoke_rate + mean_bmi + amb_visit_pct + a1c_exm_pct + nearest_hospital_km 
                  + year2001 + year2002 + year2003 + year2004 + year2005 + year2006 + year2007 + year2008 
                  + year2009 + year2010 + year2011 + year2012 + year2013 + year2014 + year2015 + year2016 
                  + region5_2 + region5_3 + region5_4 + region5_5, family=quasipoisson,method="fREML") 
    predict_try4_bam <- data.frame(predict(try4_bam, newdata = newdata_try,type="terms", se.fit = TRUE))
    
    try5_bam<-bam(y_5 ~ s(exposure,bs='cr') + ozone_summer + no2 + temp + rh + PctEye + PctLDL + Pctmam + LungCancerRate + poverty + popdensity + medianhousevalue 
                  + pct_blk + medhouseholdincome + pct_owner_occ + hispanic + education + smoke_rate + mean_bmi + amb_visit_pct + a1c_exm_pct + nearest_hospital_km 
                  + year2001 + year2002 + year2003 + year2004 + year2005 + year2006 + year2007 + year2008 
                  + year2009 + year2010 + year2011 + year2012 + year2013 + year2014 + year2015 + year2016 
                  + region5_2 + region5_3 + region5_4 + region5_5, family=quasipoisson,method="fREML") 
    predict_try5_bam <- data.frame(predict(try5_bam, newdata = newdata_try,type="terms", se.fit = TRUE))
    
    try6_bam<-bam(y_6 ~ s(exposure,bs='cr') + ozone_summer + no2 + temp + rh + PctEye + PctLDL + Pctmam + LungCancerRate + poverty + popdensity + medianhousevalue 
                  + pct_blk + medhouseholdincome + pct_owner_occ + hispanic + education + smoke_rate + mean_bmi + amb_visit_pct + a1c_exm_pct + nearest_hospital_km 
                  + year2001 + year2002 + year2003 + year2004 + year2005 + year2006 + year2007 + year2008 
                  + year2009 + year2010 + year2011 + year2012 + year2013 + year2014 + year2015 + year2016 
                  + region5_2 + region5_3 + region5_4 + region5_5, family=quasipoisson,method="fREML") 
    predict_try6_bam <- data.frame(predict(try6_bam, newdata = newdata_try,type="terms", se.fit = TRUE))
    
    ###plot
    # x=seq(0,50,by=0.1)
    # plot(0, 0, col = "white", xlab="Annual PM2.5", ylab="Log(RR)", xlim=c(0,50), ylim=c(-0.1,0.8))
    # segments(x0 = 0, y0 = 0, x1 = threshold, y1 = 0)
    # segments(x0 = threshold, y0 = 0, x1 = 50, y1 = b1*(50-threshold))
    # lines(x,predict_try1_bam$fit.s.exposure.-(predict_try1_bam$fit.s.exposure.[1]),type = "l", lty=2, col = "blue")
    # lines(x,predict_try2_bam$fit.s.exposure.-(predict_try2_bam$fit.s.exposure.[1]),type = "l", lty=2, col = "red")
    # lines(x,predict_try3_bam$fit.s.exposure.-(predict_try3_bam$fit.s.exposure.[1]),type = "l", lty=2, col = "green")
    # legend(30, 0.15, legend=c("exp", "exp+1sd", "exp+2sd", "exp+3sd"), cex=0.75,lty=c(1,2,2,2), col = c('black','blue','red','green'),border="white")

    #extract results
    results[((i-1)*6+1),(1:501)] <- predict_try1_bam$fit.s.exposure.
    results[((i-1)*6+2),(1:501)] <- predict_try2_bam$fit.s.exposure.
    results[((i-1)*6+3),(1:501)] <- predict_try3_bam$fit.s.exposure.
    results[((i-1)*6+4),(1:501)] <- predict_try4_bam$fit.s.exposure.
    results[((i-1)*6+5),(1:501)] <- predict_try5_bam$fit.s.exposure.
    results[((i-1)*6+6),(1:501)] <- predict_try6_bam$fit.s.exposure.
    results[(((i-1)*6+1):((i-1)*6+6)),'threshold'] <- threshold
    results[(((i-1)*6+1):((i-1)*6+6)),'b1'] <- b1
    results[((i-1)*6+1),'exp'] <- 'exposure1_star'
    results[((i-1)*6+2),'exp'] <- 'exposure2_star'
    results[((i-1)*6+3),'exp'] <- 'exposure3_star'
    results[((i-1)*6+4),'exp'] <- 'exposure4_star'
    results[((i-1)*6+5),'exp'] <- 'exposure5_star'
    results[((i-1)*6+6),'exp'] <- 'exposure6_star'
    
    rm(newdata_try1,try1_bam,predict_try1_bam)
    rm(newdata_try2,try2_bam,predict_try2_bam)
    rm(newdata_try3,try3_bam,predict_try3_bam)
    rm(newdata_try4,try4_bam,predict_try4_bam)
    rm(newdata_try5,try5_bam,predict_try5_bam)
    rm(newdata_try6,try6_bam,predict_try6_bam)
    
    gc()
  }
  cat(paste0("dataset is done \n"))
  #return the results 
  return(results)
  gc()
}



######################################### 2. test ################################################
# #read in the new covariates dataset
# covar<-readRDS(paste0(dir_data,'counts_weight_pred_sd_0127.rds'))
# #read in the new measurement error datasets 
# key_files <- list.files(path=dir_data,pattern = "^counts_covar_simu_set_50_(.*)rds$")
# key <- readRDS(paste0(dir_data,key_files[1]))
# 
# #create dummy variables for calendar year 
# covar$year2001=ifelse(covar$year==2001,1,0)
# covar$year2002=ifelse(covar$year==2002,1,0)
# covar$year2003=ifelse(covar$year==2003,1,0)
# covar$year2004=ifelse(covar$year==2004,1,0)
# covar$year2005=ifelse(covar$year==2005,1,0)
# covar$year2006=ifelse(covar$year==2006,1,0)
# covar$year2007=ifelse(covar$year==2007,1,0)
# covar$year2008=ifelse(covar$year==2008,1,0)
# covar$year2009=ifelse(covar$year==2009,1,0)
# covar$year2010=ifelse(covar$year==2010,1,0)
# covar$year2011=ifelse(covar$year==2011,1,0)
# covar$year2012=ifelse(covar$year==2012,1,0)
# covar$year2013=ifelse(covar$year==2013,1,0)
# covar$year2014=ifelse(covar$year==2014,1,0)
# covar$year2015=ifelse(covar$year==2015,1,0)
# covar$year2016=ifelse(covar$year==2016,1,0)
# 
# #create dummy variables for region5
# covar$region5_2=ifelse(covar$region5==2,1,0)
# covar$region5_3=ifelse(covar$region5==3,1,0)
# covar$region5_4=ifelse(covar$region5==4,1,0)
# covar$region5_5=ifelse(covar$region5==5,1,0)
# 
# #remove sds
# covar$sd_m1=covar$sd_m2=covar$sd_m3=covar$sd_m4=covar$sd_m5=covar$sd_m6=
#   covar$sd_m7=covar$sd_m8=covar$sd_m9=covar$sd_m10=covar$sd_m11=covar$sd_m12=NULL
# 
# #summary(key[1:649910,]) #9984 obs missing, 1.54% 
# 
# #parameter setting 
# n.reps=50 # 50 or 100, 200 (depending how many pieces we combine into a single loop simulation)
# 
# threshold <- 8
# b1=0.005
# 
# start.time <- Sys.time()
# test=simulate_results(key,covar,threshold,b1,n.reps)
# end.time <- Sys.time()
# end.time - start.time
# # Time difference of 1.282171 hours



######################################### 3. analysis ################################################
#read in the new covariates dataset
covar<-readRDS(paste0(dir_data,'counts_weight_pred_sd_0127.rds'))
#read in the new measurement error datasets 
key_files <- list.files(path=dir_data,pattern = "^counts_covar_simu_set_50_(.*)rds$")

#create dummy variables for calendar year 
covar$year2001=ifelse(covar$year==2001,1,0)
covar$year2002=ifelse(covar$year==2002,1,0)
covar$year2003=ifelse(covar$year==2003,1,0)
covar$year2004=ifelse(covar$year==2004,1,0)
covar$year2005=ifelse(covar$year==2005,1,0)
covar$year2006=ifelse(covar$year==2006,1,0)
covar$year2007=ifelse(covar$year==2007,1,0)
covar$year2008=ifelse(covar$year==2008,1,0)
covar$year2009=ifelse(covar$year==2009,1,0)
covar$year2010=ifelse(covar$year==2010,1,0)
covar$year2011=ifelse(covar$year==2011,1,0)
covar$year2012=ifelse(covar$year==2012,1,0)
covar$year2013=ifelse(covar$year==2013,1,0)
covar$year2014=ifelse(covar$year==2014,1,0)
covar$year2015=ifelse(covar$year==2015,1,0)
covar$year2016=ifelse(covar$year==2016,1,0)
#create dummy variables for region5
covar$region5_2=ifelse(covar$region5==2,1,0)
covar$region5_3=ifelse(covar$region5==3,1,0)
covar$region5_4=ifelse(covar$region5==4,1,0)
covar$region5_5=ifelse(covar$region5==5,1,0)
#remove sds
covar$sd_m1=covar$sd_m2=covar$sd_m3=covar$sd_m4=covar$sd_m5=covar$sd_m6=
  covar$sd_m7=covar$sd_m8=covar$sd_m9=covar$sd_m10=covar$sd_m11=covar$sd_m12=NULL

threshold_list <- c(8,9,10)
b1_list <- c(0.005,0.012,0.019)

n.reps=50 # 50 or 100, 200 (depending how many pieces we combine into a single loop simulation)


cl = makeCluster(20,outfile='')    ### change number of clusters; it equals to the number of "key" files to parallel
registerDoParallel(cl)

tmp = foreach(i=1:length(key_files), .packages=c('mgcv'))%dopar%{     ### you can select which "key" files to parallel
  key <- readRDS(paste0(dir_data,key_files[i]))
  for (threshold in threshold_list) {
    for (b1 in b1_list) {
      if (!file.exists(paste0(dir_results,'results_set',i,'_threshold',threshold,'_b1',b1,'.rds'))){
        results_temp=simulate_results(key,covar,threshold,b1,n.reps)
        saveRDS(results_temp,paste0(dir_results,'results_set',i,'_threshold',threshold,'_b1',b1,'.rds'))
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
