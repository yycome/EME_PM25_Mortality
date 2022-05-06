###############################################################################
# Project: Exposure measurement error                                         #
# Code: Step 4 - simulation - soft-threshold cr relation, penalized epi model #
# Machine: Cannon                                                             #
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

dir_data <- '~/data/'    ### change data directory
dir_results <- '~/results/Softplus/'   ### change directory for saving results



######################################### 1. function ################################################
simulate_results <- function(key,covar,k,b,n.reps){
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
  # k <- 0.2
  # b <- 24
  # n.reps <- 50
  ##### stop commenting for testing
  
  # dataframe to store resutls
  results <- data.frame(matrix(NA, nrow=n.reps*12, ncol=504))
  names(results) <- c('type','corr_magnitude','pollutants',seq(0, 50, by=0.1))
  results$k <- k
  results$b <- b
  
  # run in loops 
  for (i in 1:n.reps){
    # i=2
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
    
    # build new dataset for bam function prediction - to predict rr at each exposure level of pm2.5
    newdata_1 <- data.frame(pm25_star=seq(0,50,by=0.1))
    newdata_2 <- data.frame(pm25.x=seq(0,50,by=0.1))
    newdata_3 <- data.frame(matrix(0, nrow=length(seq(0,50,by=0.1)), ncol=42))
    names(newdata_3) <- c('ozone_summer','no2','temp','rh','PctEye','PctLDL','Pctmam','LungCancerRate','poverty','popdensity','medianhousevalue',
                          'pct_blk','medhouseholdincome','pct_owner_occ','hispanic','education','smoke_rate','mean_bmi','amb_visit_pct',
                          'a1c_exm_pct','nearest_hospital_km','amb_visit_pct','year2001','year2002','year2003','year2004','year2005','year2006',
                          'year2007','year2008','year2009','year2010','year2011','year2012','year2013','year2014','year2015','year2016',
                          'region5_2','region5_3','region5_4','region5_5')
    newdata <- cbind(newdata_1,newdata_2,newdata_3)
    
    
    ### classical errors ###
    # generate true counts 
    DATA$y_classical_multi <- rpois(nrow(DATA),exp(log(1+exp(k*(DATA$pm25.x-b)))+4.795e-04*DATA$ozone_summer+7.526e-04*DATA$no2+6.693e-04*DATA$temp+7.395e-04*DATA$rh
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
    DATA$pm25_star <- DATA$pm25.x+DATA$annual_ind_err_sd
    mod_classical_multi_ind_error_sd <- bam(y_classical_multi ~ s(pm25_star,bs='cr') + ozone_summer + no2 + temp + rh + PctEye + PctLDL + Pctmam + LungCancerRate + poverty
                                            + popdensity + medianhousevalue + pct_blk + medhouseholdincome + pct_owner_occ + hispanic + education + smoke_rate 
                                            + mean_bmi + amb_visit_pct + a1c_exm_pct + nearest_hospital_km + year2001 + year2002 + year2003 + year2004 + year2005 
                                            + year2006 + year2007 + year2008 + year2009 + year2010 + year2011 + year2012 + year2013 + year2014 + year2015 + year2016 
                                            + region5_2 + region5_3 + region5_4 + region5_5, data=DATA, family=quasipoisson, method="fREML", discrete=TRUE) 
    results[((i-1)*12+1),c('type','corr_magnitude','pollutants')] <- c('classical','ind_error_sd','multi')
    results[((i-1)*12+1),4:504] <- as.data.frame(predict(mod_classical_multi_ind_error_sd, newdata = newdata,type="terms", se.fit = TRUE))$fit.s.pm25_star.
    DATA$pm25_star <- NULL
    rm(mod_classical_multi_ind_error_sd)
    gc()
    
    DATA$pm25_star <- DATA$pm25.x+DATA$annual_ind_err_2sd
    mod_classical_multi_ind_error_2sd <- bam(y_classical_multi ~ s(pm25_star,bs='cr') + ozone_summer + no2 + temp + rh + PctEye + PctLDL + Pctmam + LungCancerRate + poverty
                                             + popdensity + medianhousevalue + pct_blk + medhouseholdincome + pct_owner_occ + hispanic + education + smoke_rate 
                                             + mean_bmi + amb_visit_pct + a1c_exm_pct + nearest_hospital_km + year2001 + year2002 + year2003 + year2004 + year2005 
                                             + year2006 + year2007 + year2008 + year2009 + year2010 + year2011 + year2012 + year2013 + year2014 + year2015 + year2016 
                                             + region5_2 + region5_3 + region5_4 + region5_5, data=DATA, family=quasipoisson, method="fREML", discrete=TRUE) 
    results[((i-1)*12+2),c('type','corr_magnitude','pollutants')] <- c('classical','ind_error_2sd','multi')
    results[((i-1)*12+2),4:504] <- as.data.frame(predict(mod_classical_multi_ind_error_2sd, newdata = newdata,type="terms", se.fit = TRUE))$fit.s.pm25_star.
    DATA$pm25_star <- NULL
    rm(mod_classical_multi_ind_error_2sd)
    gc()
    
    DATA$pm25_star <- DATA$pm25.x+DATA$annual_ind_err_3sd
    mod_classical_multi_ind_error_3sd <- bam(y_classical_multi ~ s(pm25_star,bs='cr') + ozone_summer + no2 + temp + rh + PctEye + PctLDL + Pctmam + LungCancerRate + poverty
                                             + popdensity + medianhousevalue + pct_blk + medhouseholdincome + pct_owner_occ + hispanic + education + smoke_rate 
                                             + mean_bmi + amb_visit_pct + a1c_exm_pct + nearest_hospital_km + year2001 + year2002 + year2003 + year2004 + year2005 
                                             + year2006 + year2007 + year2008 + year2009 + year2010 + year2011 + year2012 + year2013 + year2014 + year2015 + year2016 
                                             + region5_2 + region5_3 + region5_4 + region5_5, data=DATA, family=quasipoisson, method="fREML", discrete=TRUE) 
    results[((i-1)*12+3),c('type','corr_magnitude','pollutants')] <- c('classical','ind_error_3sd','multi')
    results[((i-1)*12+3),4:504] <- as.data.frame(predict(mod_classical_multi_ind_error_3sd, newdata = newdata,type="terms", se.fit = TRUE))$fit.s.pm25_star.
    DATA$pm25_star <- NULL
    rm(mod_classical_multi_ind_error_3sd)
    gc()
    
    DATA$pm25_star <- DATA$pm25.x+DATA$annual_sp_err_sd
    mod_classical_multi_sp_error_sd <- bam(y_classical_multi ~ s(pm25_star,bs='cr') + ozone_summer + no2 + temp + rh + PctEye + PctLDL + Pctmam + LungCancerRate + poverty
                                           + popdensity + medianhousevalue + pct_blk + medhouseholdincome + pct_owner_occ + hispanic + education + smoke_rate 
                                           + mean_bmi + amb_visit_pct + a1c_exm_pct + nearest_hospital_km + year2001 + year2002 + year2003 + year2004 + year2005 
                                           + year2006 + year2007 + year2008 + year2009 + year2010 + year2011 + year2012 + year2013 + year2014 + year2015 + year2016 
                                           + region5_2 + region5_3 + region5_4 + region5_5, data=DATA, family=quasipoisson, method="fREML", discrete=TRUE) 
    results[((i-1)*12+4),c('type','corr_magnitude','pollutants')] <- c('classical','sp_error_sd','multi')
    results[((i-1)*12+4),4:504] <- as.data.frame(predict(mod_classical_multi_sp_error_sd, newdata = newdata,type="terms", se.fit = TRUE))$fit.s.pm25_star.
    DATA$pm25_star <- NULL
    rm(mod_classical_multi_sp_error_sd)
    gc()
    
    DATA$pm25_star <- DATA$pm25.x+DATA$annual_sp_err_2sd
    mod_classical_multi_sp_error_2sd <- bam(y_classical_multi ~ s(pm25_star,bs='cr') + ozone_summer + no2 + temp + rh + PctEye + PctLDL + Pctmam + LungCancerRate + poverty
                                            + popdensity + medianhousevalue + pct_blk + medhouseholdincome + pct_owner_occ + hispanic + education + smoke_rate 
                                            + mean_bmi + amb_visit_pct + a1c_exm_pct + nearest_hospital_km + year2001 + year2002 + year2003 + year2004 + year2005 
                                            + year2006 + year2007 + year2008 + year2009 + year2010 + year2011 + year2012 + year2013 + year2014 + year2015 + year2016 
                                            + region5_2 + region5_3 + region5_4 + region5_5, data=DATA, family=quasipoisson, method="fREML", discrete=TRUE) 
    results[((i-1)*12+5),c('type','corr_magnitude','pollutants')] <- c('classical','sp_error_2sd','multi')
    results[((i-1)*12+5),4:504] <- as.data.frame(predict(mod_classical_multi_sp_error_2sd, newdata = newdata,type="terms", se.fit = TRUE))$fit.s.pm25_star.
    DATA$pm25_star <- NULL
    rm(mod_classical_multi_sp_error_2sd)
    gc()
    
    DATA$pm25_star <- DATA$pm25.x+DATA$annual_sp_err_3sd
    mod_classical_multi_sp_error_3sd <- bam(y_classical_multi ~ s(pm25_star,bs='cr') + ozone_summer + no2 + temp + rh + PctEye + PctLDL + Pctmam + LungCancerRate + poverty
                                            + popdensity + medianhousevalue + pct_blk + medhouseholdincome + pct_owner_occ + hispanic + education + smoke_rate 
                                            + mean_bmi + amb_visit_pct + a1c_exm_pct + nearest_hospital_km + year2001 + year2002 + year2003 + year2004 + year2005 
                                            + year2006 + year2007 + year2008 + year2009 + year2010 + year2011 + year2012 + year2013 + year2014 + year2015 + year2016 
                                            + region5_2 + region5_3 + region5_4 + region5_5, data=DATA, family=quasipoisson, method="fREML", discrete=TRUE) 
    results[((i-1)*12+6),c('type','corr_magnitude','pollutants')] <- c('classical','sp_error_3sd','multi')
    results[((i-1)*12+6),4:504] <- as.data.frame(predict(mod_classical_multi_sp_error_3sd, newdata = newdata,type="terms", se.fit = TRUE))$fit.s.pm25_star.
    DATA$pm25_star <- NULL
    rm(mod_classical_multi_sp_error_3sd)
    gc()
    
    
    ### Berkson errors ###
    DATA$y_Berkson_multi_ind_error_sd <- rpois(nrow(DATA),exp(log(1+exp(k*(DATA$pm25.x+DATA$annual_ind_err_sd-b)))+4.795e-04*DATA$ozone_summer+7.526e-04*DATA$no2+6.693e-04*DATA$temp+7.395e-04*DATA$rh
                                                              -1.665e-03*DATA$PctEye-1.015e-03*DATA$PctLDL-1.525e-03*DATA$Pctmam+6.110e-01*DATA$LungCancerRate
                                                              +1.920e-01*DATA$poverty-4.245e-06*DATA$popdensity-2.574e-07*DATA$medianhousevalue-4.967e-02*DATA$pct_blk
                                                              -9.385e-07*DATA$medhouseholdincome-2.910e-01*DATA$pct_owner_occ-2.623e-01*DATA$hispanic+1.448e-01*DATA$education
                                                              +9.841e-02*DATA$smoke_rate+6.779e-03*DATA$mean_bmi+5.998e-04*DATA$amb_visit_pct+4.921e-04*DATA$a1c_exm_pct
                                                              -3.336e-03*DATA$nearest_hospital_km-5.609e-03*DATA$year2001-1.341e-03*DATA$year2002-1.192e-02*DATA$year2003
                                                              -4.005e-02*DATA$year2004-2.362e-02*DATA$year2005-1.901e-02*DATA$year2006-2.230e-02*DATA$year2007-1.388e-02*DATA$year2008
                                                              -6.234e-02*DATA$year2009-6.470e-02*DATA$year2010-7.065e-02*DATA$year2011-9.103e-02*DATA$year2012
                                                              -1.029e-01*DATA$year2013-1.207e-01*DATA$year2014-1.134e-01*DATA$year2015-1.329e-01*DATA$year2016
                                                              -5.513e-02*DATA$region5_2-2.685e-02*DATA$region5_3-5.430e-02*DATA$region5_4-3.998e-02*DATA$region5_5))
    mod_Berkson_multi_ind_error_sd <- bam(y_Berkson_multi_ind_error_sd ~ s(pm25.x,bs='cr') + ozone_summer + no2 + temp + rh + PctEye + PctLDL + Pctmam + LungCancerRate + poverty 
                                          + popdensity + medianhousevalue + pct_blk + medhouseholdincome + pct_owner_occ + hispanic + education + smoke_rate 
                                          + mean_bmi + amb_visit_pct + a1c_exm_pct + nearest_hospital_km + year2001 + year2002 + year2003 + year2004 + year2005 
                                          + year2006 + year2007 + year2008 + year2009 + year2010 + year2011 + year2012 + year2013 + year2014 + year2015 + year2016 
                                          + region5_2 + region5_3 + region5_4 + region5_5, data=DATA, family=quasipoisson, method="fREML", discrete=TRUE) 
    results[((i-1)*12+7),c('type','corr_magnitude','pollutants')] <- c('Berkson','ind_error_sd','multi')
    results[((i-1)*12+7),4:504] <- as.data.frame(predict(mod_Berkson_multi_ind_error_sd, newdata = newdata,type="terms", se.fit = TRUE))$fit.s.pm25.x.
    DATA$y_Berkson_multi_ind_error_sd <- NULL
    rm(mod_Berkson_multi_ind_error_sd)
    gc()
    
    DATA$y_Berkson_multi_ind_error_2sd <- rpois(nrow(DATA),exp(log(1+exp(k*(DATA$pm25.x+DATA$annual_ind_err_2sd-b)))+4.795e-04*DATA$ozone_summer+7.526e-04*DATA$no2+6.693e-04*DATA$temp+7.395e-04*DATA$rh
                                                               -1.665e-03*DATA$PctEye-1.015e-03*DATA$PctLDL-1.525e-03*DATA$Pctmam+6.110e-01*DATA$LungCancerRate
                                                               +1.920e-01*DATA$poverty-4.245e-06*DATA$popdensity-2.574e-07*DATA$medianhousevalue-4.967e-02*DATA$pct_blk
                                                               -9.385e-07*DATA$medhouseholdincome-2.910e-01*DATA$pct_owner_occ-2.623e-01*DATA$hispanic+1.448e-01*DATA$education
                                                               +9.841e-02*DATA$smoke_rate+6.779e-03*DATA$mean_bmi+5.998e-04*DATA$amb_visit_pct+4.921e-04*DATA$a1c_exm_pct
                                                               -3.336e-03*DATA$nearest_hospital_km-5.609e-03*DATA$year2001-1.341e-03*DATA$year2002-1.192e-02*DATA$year2003
                                                               -4.005e-02*DATA$year2004-2.362e-02*DATA$year2005-1.901e-02*DATA$year2006-2.230e-02*DATA$year2007-1.388e-02*DATA$year2008
                                                               -6.234e-02*DATA$year2009-6.470e-02*DATA$year2010-7.065e-02*DATA$year2011-9.103e-02*DATA$year2012
                                                               -1.029e-01*DATA$year2013-1.207e-01*DATA$year2014-1.134e-01*DATA$year2015-1.329e-01*DATA$year2016
                                                               -5.513e-02*DATA$region5_2-2.685e-02*DATA$region5_3-5.430e-02*DATA$region5_4-3.998e-02*DATA$region5_5))
    mod_Berkson_multi_ind_error_2sd <- bam(y_Berkson_multi_ind_error_2sd ~ s(pm25.x,bs='cr') + ozone_summer + no2 + temp + rh + PctEye + PctLDL + Pctmam + LungCancerRate + poverty 
                                           + popdensity + medianhousevalue + pct_blk + medhouseholdincome + pct_owner_occ + hispanic + education + smoke_rate 
                                           + mean_bmi + amb_visit_pct + a1c_exm_pct + nearest_hospital_km + year2001 + year2002 + year2003 + year2004 + year2005 
                                           + year2006 + year2007 + year2008 + year2009 + year2010 + year2011 + year2012 + year2013 + year2014 + year2015 + year2016 
                                           + region5_2 + region5_3 + region5_4 + region5_5, data=DATA, family=quasipoisson, method="fREML", discrete=TRUE) 
    results[((i-1)*12+8),c('type','corr_magnitude','pollutants')] <- c('Berkson','ind_error_2sd','multi')
    results[((i-1)*12+8),4:504] <- as.data.frame(predict(mod_Berkson_multi_ind_error_2sd, newdata = newdata,type="terms", se.fit = TRUE))$fit.s.pm25.x.
    DATA$y_Berkson_multi_ind_error_2sd <- NULL
    rm(mod_Berkson_multi_ind_error_2sd)
    gc()
    
    DATA$y_Berkson_multi_ind_error_3sd <- rpois(nrow(DATA),exp(log(1+exp(k*(DATA$pm25.x+DATA$annual_ind_err_3sd-b)))+4.795e-04*DATA$ozone_summer+7.526e-04*DATA$no2+6.693e-04*DATA$temp+7.395e-04*DATA$rh
                                                               -1.665e-03*DATA$PctEye-1.015e-03*DATA$PctLDL-1.525e-03*DATA$Pctmam+6.110e-01*DATA$LungCancerRate
                                                               +1.920e-01*DATA$poverty-4.245e-06*DATA$popdensity-2.574e-07*DATA$medianhousevalue-4.967e-02*DATA$pct_blk
                                                               -9.385e-07*DATA$medhouseholdincome-2.910e-01*DATA$pct_owner_occ-2.623e-01*DATA$hispanic+1.448e-01*DATA$education
                                                               +9.841e-02*DATA$smoke_rate+6.779e-03*DATA$mean_bmi+5.998e-04*DATA$amb_visit_pct+4.921e-04*DATA$a1c_exm_pct
                                                               -3.336e-03*DATA$nearest_hospital_km-5.609e-03*DATA$year2001-1.341e-03*DATA$year2002-1.192e-02*DATA$year2003
                                                               -4.005e-02*DATA$year2004-2.362e-02*DATA$year2005-1.901e-02*DATA$year2006-2.230e-02*DATA$year2007-1.388e-02*DATA$year2008
                                                               -6.234e-02*DATA$year2009-6.470e-02*DATA$year2010-7.065e-02*DATA$year2011-9.103e-02*DATA$year2012
                                                               -1.029e-01*DATA$year2013-1.207e-01*DATA$year2014-1.134e-01*DATA$year2015-1.329e-01*DATA$year2016
                                                               -5.513e-02*DATA$region5_2-2.685e-02*DATA$region5_3-5.430e-02*DATA$region5_4-3.998e-02*DATA$region5_5))
    mod_Berkson_multi_ind_error_3sd <- bam(y_Berkson_multi_ind_error_3sd ~ s(pm25.x,bs='cr') + ozone_summer + no2 + temp + rh + PctEye + PctLDL + Pctmam + LungCancerRate + poverty 
                                           + popdensity + medianhousevalue + pct_blk + medhouseholdincome + pct_owner_occ + hispanic + education + smoke_rate 
                                           + mean_bmi + amb_visit_pct + a1c_exm_pct + nearest_hospital_km + year2001 + year2002 + year2003 + year2004 + year2005 
                                           + year2006 + year2007 + year2008 + year2009 + year2010 + year2011 + year2012 + year2013 + year2014 + year2015 + year2016 
                                           + region5_2 + region5_3 + region5_4 + region5_5, data=DATA, family=quasipoisson, method="fREML", discrete=TRUE) 
    results[((i-1)*12+9),c('type','corr_magnitude','pollutants')] <- c('Berkson','ind_error_3sd','multi')
    results[((i-1)*12+9),4:504] <- as.data.frame(predict(mod_Berkson_multi_ind_error_3sd, newdata = newdata,type="terms", se.fit = TRUE))$fit.s.pm25.x.
    DATA$y_Berkson_multi_ind_error_3sd <- NULL
    rm(mod_Berkson_multi_ind_error_3sd)
    gc()
    
    DATA$y_Berkson_multi_sp_error_sd <- rpois(nrow(DATA),exp(log(1+exp(k*(DATA$pm25.x+DATA$annual_sp_err_sd-b)))+4.795e-04*DATA$ozone_summer+7.526e-04*DATA$no2+6.693e-04*DATA$temp+7.395e-04*DATA$rh
                                                             -1.665e-03*DATA$PctEye-1.015e-03*DATA$PctLDL-1.525e-03*DATA$Pctmam+6.110e-01*DATA$LungCancerRate
                                                             +1.920e-01*DATA$poverty-4.245e-06*DATA$popdensity-2.574e-07*DATA$medianhousevalue-4.967e-02*DATA$pct_blk
                                                             -9.385e-07*DATA$medhouseholdincome-2.910e-01*DATA$pct_owner_occ-2.623e-01*DATA$hispanic+1.448e-01*DATA$education
                                                             +9.841e-02*DATA$smoke_rate+6.779e-03*DATA$mean_bmi+5.998e-04*DATA$amb_visit_pct+4.921e-04*DATA$a1c_exm_pct
                                                             -3.336e-03*DATA$nearest_hospital_km-5.609e-03*DATA$year2001-1.341e-03*DATA$year2002-1.192e-02*DATA$year2003
                                                             -4.005e-02*DATA$year2004-2.362e-02*DATA$year2005-1.901e-02*DATA$year2006-2.230e-02*DATA$year2007-1.388e-02*DATA$year2008
                                                             -6.234e-02*DATA$year2009-6.470e-02*DATA$year2010-7.065e-02*DATA$year2011-9.103e-02*DATA$year2012
                                                             -1.029e-01*DATA$year2013-1.207e-01*DATA$year2014-1.134e-01*DATA$year2015-1.329e-01*DATA$year2016
                                                             -5.513e-02*DATA$region5_2-2.685e-02*DATA$region5_3-5.430e-02*DATA$region5_4-3.998e-02*DATA$region5_5))
    mod_Berkson_multi_sp_error_sd <- bam(y_Berkson_multi_sp_error_sd ~ s(pm25.x,bs='cr') + ozone_summer + no2 + temp + rh + PctEye + PctLDL + Pctmam + LungCancerRate + poverty 
                                         + popdensity + medianhousevalue + pct_blk + medhouseholdincome + pct_owner_occ + hispanic + education + smoke_rate 
                                         + mean_bmi + amb_visit_pct + a1c_exm_pct + nearest_hospital_km + year2001 + year2002 + year2003 + year2004 + year2005 
                                         + year2006 + year2007 + year2008 + year2009 + year2010 + year2011 + year2012 + year2013 + year2014 + year2015 + year2016 
                                         + region5_2 + region5_3 + region5_4 + region5_5, data=DATA, family=quasipoisson, method="fREML", discrete=TRUE) 
    results[((i-1)*12+10),c('type','corr_magnitude','pollutants')] <- c('Berkson','sp_error_sd','multi')
    results[((i-1)*12+10),4:504] <- as.data.frame(predict(mod_Berkson_multi_sp_error_sd, newdata = newdata,type="terms", se.fit = TRUE))$fit.s.pm25.x.
    DATA$y_Berkson_multi_sp_error_sd <- NULL
    rm(mod_Berkson_multi_sp_error_sd)
    gc()
    
    DATA$y_Berkson_multi_sp_error_2sd <- rpois(nrow(DATA),exp(log(1+exp(k*(DATA$pm25.x+DATA$annual_sp_err_2sd-b)))+4.795e-04*DATA$ozone_summer+7.526e-04*DATA$no2+6.693e-04*DATA$temp+7.395e-04*DATA$rh
                                                              -1.665e-03*DATA$PctEye-1.015e-03*DATA$PctLDL-1.525e-03*DATA$Pctmam+6.110e-01*DATA$LungCancerRate
                                                              +1.920e-01*DATA$poverty-4.245e-06*DATA$popdensity-2.574e-07*DATA$medianhousevalue-4.967e-02*DATA$pct_blk
                                                              -9.385e-07*DATA$medhouseholdincome-2.910e-01*DATA$pct_owner_occ-2.623e-01*DATA$hispanic+1.448e-01*DATA$education
                                                              +9.841e-02*DATA$smoke_rate+6.779e-03*DATA$mean_bmi+5.998e-04*DATA$amb_visit_pct+4.921e-04*DATA$a1c_exm_pct
                                                              -3.336e-03*DATA$nearest_hospital_km-5.609e-03*DATA$year2001-1.341e-03*DATA$year2002-1.192e-02*DATA$year2003
                                                              -4.005e-02*DATA$year2004-2.362e-02*DATA$year2005-1.901e-02*DATA$year2006-2.230e-02*DATA$year2007-1.388e-02*DATA$year2008
                                                              -6.234e-02*DATA$year2009-6.470e-02*DATA$year2010-7.065e-02*DATA$year2011-9.103e-02*DATA$year2012
                                                              -1.029e-01*DATA$year2013-1.207e-01*DATA$year2014-1.134e-01*DATA$year2015-1.329e-01*DATA$year2016
                                                              -5.513e-02*DATA$region5_2-2.685e-02*DATA$region5_3-5.430e-02*DATA$region5_4-3.998e-02*DATA$region5_5))
    mod_Berkson_multi_sp_error_2sd <- bam(y_Berkson_multi_sp_error_2sd ~ s(pm25.x,bs='cr') + ozone_summer + no2 + temp + rh + PctEye + PctLDL + Pctmam + LungCancerRate + poverty 
                                          + popdensity + medianhousevalue + pct_blk + medhouseholdincome + pct_owner_occ + hispanic + education + smoke_rate 
                                          + mean_bmi + amb_visit_pct + a1c_exm_pct + nearest_hospital_km + year2001 + year2002 + year2003 + year2004 + year2005 
                                          + year2006 + year2007 + year2008 + year2009 + year2010 + year2011 + year2012 + year2013 + year2014 + year2015 + year2016 
                                          + region5_2 + region5_3 + region5_4 + region5_5, data=DATA, family=quasipoisson, method="fREML", discrete=TRUE) 
    results[((i-1)*12+11),c('type','corr_magnitude','pollutants')] <- c('Berkson','sp_error_2sd','multi')
    results[((i-1)*12+11),4:504] <- as.data.frame(predict(mod_Berkson_multi_sp_error_2sd, newdata = newdata,type="terms", se.fit = TRUE))$fit.s.pm25.x.
    DATA$y_Berkson_multi_sp_error_2sd <- NULL
    rm(mod_Berkson_multi_sp_error_2sd)
    gc()
    
    DATA$y_Berkson_multi_sp_error_3sd <- rpois(nrow(DATA),exp(log(1+exp(k*(DATA$pm25.x+DATA$annual_sp_err_3sd-b)))+4.795e-04*DATA$ozone_summer+7.526e-04*DATA$no2+6.693e-04*DATA$temp+7.395e-04*DATA$rh
                                                              -1.665e-03*DATA$PctEye-1.015e-03*DATA$PctLDL-1.525e-03*DATA$Pctmam+6.110e-01*DATA$LungCancerRate
                                                              +1.920e-01*DATA$poverty-4.245e-06*DATA$popdensity-2.574e-07*DATA$medianhousevalue-4.967e-02*DATA$pct_blk
                                                              -9.385e-07*DATA$medhouseholdincome-2.910e-01*DATA$pct_owner_occ-2.623e-01*DATA$hispanic+1.448e-01*DATA$education
                                                              +9.841e-02*DATA$smoke_rate+6.779e-03*DATA$mean_bmi+5.998e-04*DATA$amb_visit_pct+4.921e-04*DATA$a1c_exm_pct
                                                              -3.336e-03*DATA$nearest_hospital_km-5.609e-03*DATA$year2001-1.341e-03*DATA$year2002-1.192e-02*DATA$year2003
                                                              -4.005e-02*DATA$year2004-2.362e-02*DATA$year2005-1.901e-02*DATA$year2006-2.230e-02*DATA$year2007-1.388e-02*DATA$year2008
                                                              -6.234e-02*DATA$year2009-6.470e-02*DATA$year2010-7.065e-02*DATA$year2011-9.103e-02*DATA$year2012
                                                              -1.029e-01*DATA$year2013-1.207e-01*DATA$year2014-1.134e-01*DATA$year2015-1.329e-01*DATA$year2016
                                                              -5.513e-02*DATA$region5_2-2.685e-02*DATA$region5_3-5.430e-02*DATA$region5_4-3.998e-02*DATA$region5_5))
    mod_Berkson_multi_sp_error_3sd <- bam(y_Berkson_multi_sp_error_3sd ~ s(pm25.x,bs='cr') + ozone_summer + no2 + temp + rh + PctEye + PctLDL + Pctmam + LungCancerRate + poverty 
                                          + popdensity + medianhousevalue + pct_blk + medhouseholdincome + pct_owner_occ + hispanic + education + smoke_rate 
                                          + mean_bmi + amb_visit_pct + a1c_exm_pct + nearest_hospital_km + year2001 + year2002 + year2003 + year2004 + year2005 
                                          + year2006 + year2007 + year2008 + year2009 + year2010 + year2011 + year2012 + year2013 + year2014 + year2015 + year2016 
                                          + region5_2 + region5_3 + region5_4 + region5_5, data=DATA, family=quasipoisson, method="fREML", discrete=TRUE) 
    results[((i-1)*12+12),c('type','corr_magnitude','pollutants')] <- c('Berkson','sp_error_3sd','multi')
    results[((i-1)*12+12),4:504] <- as.data.frame(predict(mod_Berkson_multi_sp_error_3sd, newdata = newdata,type="terms", se.fit = TRUE))$fit.s.pm25.x.
    DATA$y_Berkson_multi_sp_error_3sd <- NULL
    rm(mod_Berkson_multi_sp_error_3sd)
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
# key_files <- list.files(path=dir_data,pattern = "^counts_covar_simu_set_50_(.*)rds$")
# key <- readRDS(paste0(dir_data,key_files[1]))
# 
# # parameter setting
# k <- 0.2
# b <- 24
# n.reps <- 50  # 50 or 100, 200 (depending how many pieces we combine into a single loop simulation)
# 
# start.time <- Sys.time()
# test <- simulate_results(key,covar,k,b,n.reps)
# end.time <- Sys.time()
# end.time - start.time
# # Time difference of 3.750365 hours



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

# list measurement error datasets
key_files <- list.files(path=dir_data,pattern = "^counts_covar_simu_set_50_(.*)rds$")

# parameter setting
k_list <- c(0.2,0.3,0.4)
b_list <- c(22,23,24,25,26)
n.reps <- 50

cl = makeCluster(10,outfile='')    ### change number of clusters
registerDoParallel(cl)

tmp <- foreach(i=1:length(key_files), .packages=c('mgcv'))%dopar%{
  key <- readRDS(paste0(dir_data,key_files[i]))
  for (k in k_list) {
    for (b in b_list) {
      if (!file.exists(paste0(dir_results,'results_set',i,'_k',k,'_b',b,'.rds'))){
        results_temp <- simulate_results(key,covar,k,b,n.reps)
        saveRDS(results_temp,paste0(dir_results,'results_set',i,'_k',k,'_b',b,'.rds'))
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
