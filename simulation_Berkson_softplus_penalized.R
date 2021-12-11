#################################################################################################
# Project: EME                                                                                  #
# Code: simulate mortality with softplus model and test with penalized splines - Berkson error  #
# Author: Yaguang Wei                                                                           #
#################################################################################################

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

dir_data <- '/n/home07/yaw624/EME/data/'                                  ### change directory
dir_results <- '/n/home07/yaw624/EME/results/Berkson/softplus_penalized/'     ### change directory



######################################### 1. function ################################################
simulate_results<-function(key,covar,k,b,n.reps,n.obs){
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
  # k <- 0.3
  # b <- 26
  # n.obs=(649910-9984)
  # n.reps=50
  #####stop commenting for testing
  
  #dataframe to store resutls
  results <- data.frame(matrix(NA, nrow=n.reps*6, ncol=501))
  names(results) <- seq(0, 50, by=0.1)
  results$k <- NA
  results$b <- NA
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
    
    #assign exposure and error value 
    exposure<-DATA$pm25.x #QD's predicted annual pm2.5 exposure 
    error1sd<-DATA$annual_tempcorr_adj_err
    error2sd<-DATA$annual_tempcorr_adj_err_2sd
    error3sd<-DATA$annual_tempcorr_adj_err_3sd
    error1sd_sp<-DATA$annual_tempcorr_adj_err_sp
    error2sd_sp<-DATA$annual_tempcorr_adj_err_2sd_sp
    error3sd_sp<-DATA$annual_tempcorr_adj_err_3sd_sp
    
    #compute new simulated exposure 
    exposure1_star<-exposure + error1sd 
    exposure2_star<-exposure + error2sd  
    exposure3_star<-exposure + error3sd 
    exposure4_star<-exposure + error1sd_sp 
    exposure5_star<-exposure + error2sd_sp 
    exposure6_star<-exposure + error3sd_sp 
    
    #generate the true counts 
    y1<-rpois(n.obs,exp(log(1+exp(k*(exposure1_star-b)))+4.795e-04*ozone_summer+7.526e-04*no2+6.693e-04*temp+7.395e-04*rh
                        -1.665e-03*PctEye-1.015e-03*PctLDL-1.525e-03*Pctmam+6.110e-01*LungCancerRate
                        +1.920e-01*poverty-4.245e-06*popdensity-2.574e-07*medianhousevalue-4.967e-02*pct_blk
                        -9.385e-07*medhouseholdincome-2.910e-01*pct_owner_occ-2.623e-01*hispanic+1.448e-01*education
                        +9.841e-02*smoke_rate+6.779e-03*mean_bmi+5.998e-04*amb_visit_pct+4.921e-04*a1c_exm_pct-3.336e-03*nearest_hospital_km
                        -5.609e-03*year2001-1.341e-03*year2002-1.192e-02*year2003-4.005e-02*year2004
                        -2.362e-02*year2005-1.901e-02*year2006-2.230e-02*year2007-1.388e-02*year2008
                        -6.234e-02*year2009-6.470e-02*year2010-7.065e-02*year2011-9.103e-02*year2012
                        -1.029e-01*year2013-1.207e-01*year2014-1.134e-01*year2015-1.329e-01*year2016
                        -5.513e-02*region5_2-2.685e-02*region5_3-5.430e-02*region5_4-3.998e-02*region5_5)) #remember the exp(), here we assume the exposure is the truth
    
    y2<-rpois(n.obs,exp(log(1+exp(k*(exposure2_star-b)))+4.795e-04*ozone_summer+7.526e-04*no2+6.693e-04*temp+7.395e-04*rh
                        -1.665e-03*PctEye-1.015e-03*PctLDL-1.525e-03*Pctmam+6.110e-01*LungCancerRate
                        +1.920e-01*poverty-4.245e-06*popdensity-2.574e-07*medianhousevalue-4.967e-02*pct_blk
                        -9.385e-07*medhouseholdincome-2.910e-01*pct_owner_occ-2.623e-01*hispanic+1.448e-01*education
                        +9.841e-02*smoke_rate+6.779e-03*mean_bmi+5.998e-04*amb_visit_pct+4.921e-04*a1c_exm_pct-3.336e-03*nearest_hospital_km
                        -5.609e-03*year2001-1.341e-03*year2002-1.192e-02*year2003-4.005e-02*year2004
                        -2.362e-02*year2005-1.901e-02*year2006-2.230e-02*year2007-1.388e-02*year2008
                        -6.234e-02*year2009-6.470e-02*year2010-7.065e-02*year2011-9.103e-02*year2012
                        -1.029e-01*year2013-1.207e-01*year2014-1.134e-01*year2015-1.329e-01*year2016
                        -5.513e-02*region5_2-2.685e-02*region5_3-5.430e-02*region5_4-3.998e-02*region5_5))
    
    y3<-rpois(n.obs,exp(log(1+exp(k*(exposure3_star-b)))+4.795e-04*ozone_summer+7.526e-04*no2+6.693e-04*temp+7.395e-04*rh
                        -1.665e-03*PctEye-1.015e-03*PctLDL-1.525e-03*Pctmam+6.110e-01*LungCancerRate
                        +1.920e-01*poverty-4.245e-06*popdensity-2.574e-07*medianhousevalue-4.967e-02*pct_blk
                        -9.385e-07*medhouseholdincome-2.910e-01*pct_owner_occ-2.623e-01*hispanic+1.448e-01*education
                        +9.841e-02*smoke_rate+6.779e-03*mean_bmi+5.998e-04*amb_visit_pct+4.921e-04*a1c_exm_pct-3.336e-03*nearest_hospital_km
                        -5.609e-03*year2001-1.341e-03*year2002-1.192e-02*year2003-4.005e-02*year2004
                        -2.362e-02*year2005-1.901e-02*year2006-2.230e-02*year2007-1.388e-02*year2008
                        -6.234e-02*year2009-6.470e-02*year2010-7.065e-02*year2011-9.103e-02*year2012
                        -1.029e-01*year2013-1.207e-01*year2014-1.134e-01*year2015-1.329e-01*year2016
                        -5.513e-02*region5_2-2.685e-02*region5_3-5.430e-02*region5_4-3.998e-02*region5_5))
    
    y4<-rpois(n.obs,exp(log(1+exp(k*(exposure4_star-b)))+4.795e-04*ozone_summer+7.526e-04*no2+6.693e-04*temp+7.395e-04*rh
                        -1.665e-03*PctEye-1.015e-03*PctLDL-1.525e-03*Pctmam+6.110e-01*LungCancerRate
                        +1.920e-01*poverty-4.245e-06*popdensity-2.574e-07*medianhousevalue-4.967e-02*pct_blk
                        -9.385e-07*medhouseholdincome-2.910e-01*pct_owner_occ-2.623e-01*hispanic+1.448e-01*education
                        +9.841e-02*smoke_rate+6.779e-03*mean_bmi+5.998e-04*amb_visit_pct+4.921e-04*a1c_exm_pct-3.336e-03*nearest_hospital_km
                        -5.609e-03*year2001-1.341e-03*year2002-1.192e-02*year2003-4.005e-02*year2004
                        -2.362e-02*year2005-1.901e-02*year2006-2.230e-02*year2007-1.388e-02*year2008
                        -6.234e-02*year2009-6.470e-02*year2010-7.065e-02*year2011-9.103e-02*year2012
                        -1.029e-01*year2013-1.207e-01*year2014-1.134e-01*year2015-1.329e-01*year2016
                        -5.513e-02*region5_2-2.685e-02*region5_3-5.430e-02*region5_4-3.998e-02*region5_5))
    
    y5<-rpois(n.obs,exp(log(1+exp(k*(exposure5_star-b)))+4.795e-04*ozone_summer+7.526e-04*no2+6.693e-04*temp+7.395e-04*rh
                        -1.665e-03*PctEye-1.015e-03*PctLDL-1.525e-03*Pctmam+6.110e-01*LungCancerRate
                        +1.920e-01*poverty-4.245e-06*popdensity-2.574e-07*medianhousevalue-4.967e-02*pct_blk
                        -9.385e-07*medhouseholdincome-2.910e-01*pct_owner_occ-2.623e-01*hispanic+1.448e-01*education
                        +9.841e-02*smoke_rate+6.779e-03*mean_bmi+5.998e-04*amb_visit_pct+4.921e-04*a1c_exm_pct-3.336e-03*nearest_hospital_km
                        -5.609e-03*year2001-1.341e-03*year2002-1.192e-02*year2003-4.005e-02*year2004
                        -2.362e-02*year2005-1.901e-02*year2006-2.230e-02*year2007-1.388e-02*year2008
                        -6.234e-02*year2009-6.470e-02*year2010-7.065e-02*year2011-9.103e-02*year2012
                        -1.029e-01*year2013-1.207e-01*year2014-1.134e-01*year2015-1.329e-01*year2016
                        -5.513e-02*region5_2-2.685e-02*region5_3-5.430e-02*region5_4-3.998e-02*region5_5))
    
    y6<-rpois(n.obs,exp(log(1+exp(k*(exposure6_star-b)))+4.795e-04*ozone_summer+7.526e-04*no2+6.693e-04*temp+7.395e-04*rh
                        -1.665e-03*PctEye-1.015e-03*PctLDL-1.525e-03*Pctmam+6.110e-01*LungCancerRate
                        +1.920e-01*poverty-4.245e-06*popdensity-2.574e-07*medianhousevalue-4.967e-02*pct_blk
                        -9.385e-07*medhouseholdincome-2.910e-01*pct_owner_occ-2.623e-01*hispanic+1.448e-01*education
                        +9.841e-02*smoke_rate+6.779e-03*mean_bmi+5.998e-04*amb_visit_pct+4.921e-04*a1c_exm_pct-3.336e-03*nearest_hospital_km
                        -5.609e-03*year2001-1.341e-03*year2002-1.192e-02*year2003-4.005e-02*year2004
                        -2.362e-02*year2005-1.901e-02*year2006-2.230e-02*year2007-1.388e-02*year2008
                        -6.234e-02*year2009-6.470e-02*year2010-7.065e-02*year2011-9.103e-02*year2012
                        -1.029e-01*year2013-1.207e-01*year2014-1.134e-01*year2015-1.329e-01*year2016
                        -5.513e-02*region5_2-2.685e-02*region5_3-5.430e-02*region5_4-3.998e-02*region5_5))
    
    #build new dataset to predict rr
    newdata_1 <- data.frame(exposure=seq(0,50,by=0.1))
    newdata_2 <- data.frame(matrix(0, nrow=length(seq(0,50,by=0.1)), ncol=42))
    names(newdata_2) <- c('ozone_summer','no2','temp','rh','PctEye','PctLDL','Pctmam','LungCancerRate','poverty','popdensity','medianhousevalue',
                          'pct_blk','medhouseholdincome','pct_owner_occ','hispanic','education','smoke_rate','mean_bmi','amb_visit_pct',
                          'a1c_exm_pct','nearest_hospital_km','amb_visit_pct','year2001','year2002','year2003','year2004','year2005','year2006',
                          'year2007','year2008','year2009','year2010','year2011','year2012','year2013','year2014','year2015','year2016',
                          'region5_2','region5_3','region5_4','region5_5')
    newdata <- cbind(newdata_1,newdata_2)
    
    #run the Poisson models 
    try1_bam<-bam(y1 ~ s(exposure,bs='cr') + ozone_summer + no2 + temp + rh + PctEye + PctLDL + Pctmam + LungCancerRate + poverty + popdensity + medianhousevalue 
                  + pct_blk + medhouseholdincome + pct_owner_occ + hispanic + education + smoke_rate + mean_bmi + amb_visit_pct + a1c_exm_pct + nearest_hospital_km 
                  + year2001 + year2002 + year2003 + year2004 + year2005 + year2006 + year2007 + year2008 
                  + year2009 + year2010 + year2011 + year2012 + year2013 + year2014 + year2015 + year2016 
                  + region5_2 + region5_3 + region5_4 + region5_5, family=quasipoisson,method="fREML") 
    predict_try1_bam <- data.frame(predict(try1_bam, newdata = newdata,type="terms", se.fit = TRUE))
    
    try2_bam<-bam(y2 ~ s(exposure,bs='cr') + ozone_summer + no2 + temp + rh + PctEye + PctLDL + Pctmam + LungCancerRate + poverty + popdensity + medianhousevalue 
                  + pct_blk + medhouseholdincome + pct_owner_occ + hispanic + education + smoke_rate + mean_bmi + amb_visit_pct + a1c_exm_pct + nearest_hospital_km 
                  + year2001 + year2002 + year2003 + year2004 + year2005 + year2006 + year2007 + year2008 
                  + year2009 + year2010 + year2011 + year2012 + year2013 + year2014 + year2015 + year2016 
                  + region5_2 + region5_3 + region5_4 + region5_5, family=quasipoisson,method="fREML") 
    predict_try2_bam <- data.frame(predict(try2_bam, newdata = newdata,type="terms", se.fit = TRUE))
    
    try3_bam<-bam(y3 ~ s(exposure,bs='cr') + ozone_summer + no2 + temp + rh + PctEye + PctLDL + Pctmam + LungCancerRate + poverty + popdensity + medianhousevalue 
                  + pct_blk + medhouseholdincome + pct_owner_occ + hispanic + education + smoke_rate + mean_bmi + amb_visit_pct + a1c_exm_pct + nearest_hospital_km 
                  + year2001 + year2002 + year2003 + year2004 + year2005 + year2006 + year2007 + year2008 
                  + year2009 + year2010 + year2011 + year2012 + year2013 + year2014 + year2015 + year2016 
                  + region5_2 + region5_3 + region5_4 + region5_5, family=quasipoisson,method="fREML") 
    predict_try3_bam <- data.frame(predict(try3_bam, newdata = newdata,type="terms", se.fit = TRUE))
    
    try4_bam<-bam(y4 ~ s(exposure,bs='cr') + ozone_summer + no2 + temp + rh + PctEye + PctLDL + Pctmam + LungCancerRate + poverty + popdensity + medianhousevalue 
                  + pct_blk + medhouseholdincome + pct_owner_occ + hispanic + education + smoke_rate + mean_bmi + amb_visit_pct + a1c_exm_pct + nearest_hospital_km 
                  + year2001 + year2002 + year2003 + year2004 + year2005 + year2006 + year2007 + year2008 
                  + year2009 + year2010 + year2011 + year2012 + year2013 + year2014 + year2015 + year2016 
                  + region5_2 + region5_3 + region5_4 + region5_5, family=quasipoisson,method="fREML") 
    predict_try4_bam <- data.frame(predict(try4_bam, newdata = newdata,type="terms", se.fit = TRUE))
    
    try5_bam<-bam(y5 ~ s(exposure,bs='cr') + ozone_summer + no2 + temp + rh + PctEye + PctLDL + Pctmam + LungCancerRate + poverty + popdensity + medianhousevalue 
                  + pct_blk + medhouseholdincome + pct_owner_occ + hispanic + education + smoke_rate + mean_bmi + amb_visit_pct + a1c_exm_pct + nearest_hospital_km 
                  + year2001 + year2002 + year2003 + year2004 + year2005 + year2006 + year2007 + year2008 
                  + year2009 + year2010 + year2011 + year2012 + year2013 + year2014 + year2015 + year2016 
                  + region5_2 + region5_3 + region5_4 + region5_5, family=quasipoisson,method="fREML") 
    predict_try5_bam <- data.frame(predict(try5_bam, newdata = newdata,type="terms", se.fit = TRUE))
    
    try6_bam<-bam(y6 ~ s(exposure,bs='cr') + ozone_summer + no2 + temp + rh + PctEye + PctLDL + Pctmam + LungCancerRate + poverty + popdensity + medianhousevalue 
                  + pct_blk + medhouseholdincome + pct_owner_occ + hispanic + education + smoke_rate + mean_bmi + amb_visit_pct + a1c_exm_pct + nearest_hospital_km 
                  + year2001 + year2002 + year2003 + year2004 + year2005 + year2006 + year2007 + year2008 
                  + year2009 + year2010 + year2011 + year2012 + year2013 + year2014 + year2015 + year2016 
                  + region5_2 + region5_3 + region5_4 + region5_5, family=quasipoisson,method="fREML") 
    predict_try6_bam <- data.frame(predict(try6_bam, newdata = newdata,type="terms", se.fit = TRUE))
    
    #extract results
    results[((i-1)*6+1),(1:501)] <- predict_try1_bam$fit.s.exposure
    results[((i-1)*6+2),(1:501)] <- predict_try2_bam$fit.s.exposure
    results[((i-1)*6+3),(1:501)] <- predict_try3_bam$fit.s.exposure
    results[((i-1)*6+4),(1:501)] <- predict_try4_bam$fit.s.exposure
    results[((i-1)*6+5),(1:501)] <- predict_try5_bam$fit.s.exposure
    results[((i-1)*6+6),(1:501)] <- predict_try6_bam$fit.s.exposure
    results[(((i-1)*6+1):((i-1)*6+6)),'k'] <- k
    results[(((i-1)*6+1):((i-1)*6+6)),'b'] <- b
    results[((i-1)*6+1),'exp'] <- 'exposure1_star'
    results[((i-1)*6+2),'exp'] <- 'exposure2_star'
    results[((i-1)*6+3),'exp'] <- 'exposure3_star'
    results[((i-1)*6+4),'exp'] <- 'exposure4_star'
    results[((i-1)*6+5),'exp'] <- 'exposure5_star'
    results[((i-1)*6+6),'exp'] <- 'exposure6_star'
    
    rm(newdata,
       try1_bam,predict_try1_bam,
       try2_bam,predict_try2_bam,
       try3_bam,predict_try3_bam,
       try4_bam,predict_try4_bam,
       try5_bam,predict_try5_bam,
       try6_bam,predict_try6_bam)
    
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
# n.obs=(649910-9984) #length of obs in each chunked data piece after removing missing
# n.reps=50 # 50 or 100, 200 (depending how many pieces we combine into a single loop simulation)
# 
# k_list <- c(0.2,0.3,0.4)
# b_list <- c(22,23,24,25,26)
# 
# k=k_list[1]
# b=b_list[1]
# 
# start.time <- Sys.time()
# test=simulate_results(key,covar,k,b,n.reps,n.obs)
# end.time <- Sys.time()
# end.time - start.time
# # Time difference of 1.340672 hours



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

k_list <- c(0.2,0.3,0.4)
b_list <- c(22,23,24,25,26)

n.obs=(649910-9984) #length of obs in each chunked data piece after removing missing
n.reps=50 # 50 or 100, 200 (depending how many pieces we combine into a single loop simulation)


cl = makeCluster(5,outfile='')    ### change number of clusters
registerDoParallel(cl)

tmp = foreach(i=1:length(key_files), .packages=c('mgcv'))%dopar%{
  key <- readRDS(paste0(dir_data,key_files[i]))
  for (k in k_list) {
    for (b in b_list) {
      if (!file.exists(paste0(dir_results,'results_set',i,'_k',k,'_b',b,'.rds'))){
        results_temp=simulate_results(key,covar,k,b,n.reps,n.obs)
        saveRDS(results_temp,paste0(dir_results,'results_set',i,'_k',k,'_b',b,'.rds'))
        rm(results_temp)
        gc()
      }
    }
  }
  print(paste0("key file ",i," is done \n"))
  rm(key)
  gc()
}

stopCluster(cl)
