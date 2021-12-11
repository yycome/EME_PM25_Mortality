###########################################################################################################
# Project: EME                                                                                            #
# Code: simulate mortality with piecewise linear PM and test with penalized splines - classical error     #
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
dir_results <- '/media/qnap4/Yaguang/EME/results/piecewise_penalized/'  ### create directory to save results



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
    
    #inidcate obs with exposure below the threshold
    DATA$low_exp <- 0
    DATA$low_exp[DATA$pm25.x<=threshold] <- 1
    
    #set IDs for each ob
    DATA$id <- 1:nrow(DATA)
    
    #generate the true counts for low-exposure obs
    DATA_low_exp <- DATA[DATA$low_exp==1,]
    DATA_low_exp$y <- rpois(nrow(DATA_low_exp),exp(4.795e-04*DATA_low_exp$ozone_summer+7.526e-04*DATA_low_exp$no2+6.693e-04*DATA_low_exp$temp+7.395e-04*DATA_low_exp$rh
                                                   -1.665e-03*DATA_low_exp$PctEye-1.015e-03*DATA_low_exp$PctLDL-1.525e-03*DATA_low_exp$Pctmam+6.110e-01*DATA_low_exp$LungCancerRate
                                                   +1.920e-01*DATA_low_exp$poverty-4.245e-06*DATA_low_exp$popdensity-2.574e-07*DATA_low_exp$medianhousevalue-4.967e-02*DATA_low_exp$pct_blk
                                                   -9.385e-07*DATA_low_exp$medhouseholdincome-2.910e-01*DATA_low_exp$pct_owner_occ-2.623e-01*DATA_low_exp$hispanic+1.448e-01*DATA_low_exp$education
                                                   +9.841e-02*DATA_low_exp$smoke_rate+6.779e-03*DATA_low_exp$mean_bmi+5.998e-04*DATA_low_exp$amb_visit_pct+4.921e-04*DATA_low_exp$a1c_exm_pct-3.336e-03*DATA_low_exp$nearest_hospital_km
                                                   -5.609e-03*DATA_low_exp$year2001-1.341e-03*DATA_low_exp$year2002-1.192e-02*DATA_low_exp$year2003-4.005e-02*DATA_low_exp$year2004
                                                   -2.362e-02*DATA_low_exp$year2005-1.901e-02*DATA_low_exp$year2006-2.230e-02*DATA_low_exp$year2007-1.388e-02*DATA_low_exp$year2008
                                                   -6.234e-02*DATA_low_exp$year2009-6.470e-02*DATA_low_exp$year2010-7.065e-02*DATA_low_exp$year2011-9.103e-02*DATA_low_exp$year2012
                                                   -1.029e-01*DATA_low_exp$year2013-1.207e-01*DATA_low_exp$year2014-1.134e-01*DATA_low_exp$year2015-1.329e-01*DATA_low_exp$year2016
                                                   -5.513e-02*DATA_low_exp$region5_2-2.685e-02*DATA_low_exp$region5_3-5.430e-02*DATA_low_exp$region5_4-3.998e-02*DATA_low_exp$region5_5)) 
    
    #generate the true counts for high-exposure obs
    DATA_high_exp <- DATA[DATA$low_exp==0,]
    DATA_high_exp$y <- rpois(nrow(DATA_high_exp),exp(-b1*threshold+b1*(DATA_high_exp$pm25.x-threshold)+4.795e-04*DATA_high_exp$ozone_summer+7.526e-04*DATA_high_exp$no2+6.693e-04*DATA_high_exp$temp+7.395e-04*DATA_high_exp$rh
                                                     -1.665e-03*DATA_high_exp$PctEye-1.015e-03*DATA_high_exp$PctLDL-1.525e-03*DATA_high_exp$Pctmam+6.110e-01*DATA_high_exp$LungCancerRate
                                                     +1.920e-01*DATA_high_exp$poverty-4.245e-06*DATA_high_exp$popdensity-2.574e-07*DATA_high_exp$medianhousevalue-4.967e-02*DATA_high_exp$pct_blk
                                                     -9.385e-07*DATA_high_exp$medhouseholdincome-2.910e-01*DATA_high_exp$pct_owner_occ-2.623e-01*DATA_high_exp$hispanic+1.448e-01*DATA_high_exp$education
                                                     +9.841e-02*DATA_high_exp$smoke_rate+6.779e-03*DATA_high_exp$mean_bmi+5.998e-04*DATA_high_exp$amb_visit_pct+4.921e-04*DATA_high_exp$a1c_exm_pct-3.336e-03*DATA_high_exp$nearest_hospital_km
                                                     -5.609e-03*DATA_high_exp$year2001-1.341e-03*DATA_high_exp$year2002-1.192e-02*DATA_high_exp$year2003-4.005e-02*DATA_high_exp$year2004
                                                     -2.362e-02*DATA_high_exp$year2005-1.901e-02*DATA_high_exp$year2006-2.230e-02*DATA_high_exp$year2007-1.388e-02*DATA_high_exp$year2008
                                                     -6.234e-02*DATA_high_exp$year2009-6.470e-02*DATA_high_exp$year2010-7.065e-02*DATA_high_exp$year2011-9.103e-02*DATA_high_exp$year2012
                                                     -1.029e-01*DATA_high_exp$year2013-1.207e-01*DATA_high_exp$year2014-1.134e-01*DATA_high_exp$year2015-1.329e-01*DATA_high_exp$year2016
                                                     -5.513e-02*DATA_high_exp$region5_2-2.685e-02*DATA_high_exp$region5_3-5.430e-02*DATA_high_exp$region5_4-3.998e-02*DATA_high_exp$region5_5)) 
    
    #combine and sort
    DATA <- rbind(DATA_low_exp,DATA_high_exp)
    DATA <- DATA[order(DATA$id),]
    
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
    y<-DATA$y
    
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
    
    #build new dataset to predict rr
    newdata_try1_1 <- data.frame(exposure1_star=seq(0,50,by=0.1))
    newdata_try1_2 <- data.frame(matrix(0, nrow=length(seq(0,50,by=0.1)), ncol=42))
    names(newdata_try1_2) <- c('ozone_summer','no2','temp','rh','PctEye','PctLDL','Pctmam','LungCancerRate','poverty','popdensity','medianhousevalue',
                               'pct_blk','medhouseholdincome','pct_owner_occ','hispanic','education','smoke_rate','mean_bmi','amb_visit_pct',
                               'a1c_exm_pct','nearest_hospital_km','amb_visit_pct','year2001','year2002','year2003','year2004','year2005','year2006',
                               'year2007','year2008','year2009','year2010','year2011','year2012','year2013','year2014','year2015','year2016',
                               'region5_2','region5_3','region5_4','region5_5')
    newdata_try1 <- cbind(newdata_try1_1,newdata_try1_2)
    
    newdata_try2 <- newdata_try1
    names(newdata_try2)[1] <- 'exposure2_star'
    
    newdata_try3 <- newdata_try1
    names(newdata_try3)[1] <- 'exposure3_star'
    
    newdata_try4 <- newdata_try1
    names(newdata_try4)[1] <- 'exposure4_star'
    
    newdata_try5 <- newdata_try1
    names(newdata_try5)[1] <- 'exposure5_star'
    
    newdata_try6 <- newdata_try1
    names(newdata_try6)[1] <- 'exposure6_star'
    
    #run the Poisson models 
    try1_bam<-bam(y ~ s(exposure1_star,bs='cr') + ozone_summer + no2 + temp + rh + PctEye + PctLDL + Pctmam + LungCancerRate + poverty + popdensity + medianhousevalue 
                  + pct_blk + medhouseholdincome + pct_owner_occ + hispanic + education + smoke_rate + mean_bmi + amb_visit_pct + a1c_exm_pct + nearest_hospital_km 
                  + year2001 + year2002 + year2003 + year2004 + year2005 + year2006 + year2007 + year2008 
                  + year2009 + year2010 + year2011 + year2012 + year2013 + year2014 + year2015 + year2016 
                  + region5_2 + region5_3 + region5_4 + region5_5, family=quasipoisson,method="fREML") 
    predict_try1_bam <- data.frame(predict(try1_bam, newdata = newdata_try1,type="terms", se.fit = TRUE))
    
    try2_bam<-bam(y ~ s(exposure2_star,bs='cr') + ozone_summer + no2 + temp + rh + PctEye + PctLDL + Pctmam + LungCancerRate + poverty + popdensity + medianhousevalue 
                  + pct_blk + medhouseholdincome + pct_owner_occ + hispanic + education + smoke_rate + mean_bmi + amb_visit_pct + a1c_exm_pct + nearest_hospital_km 
                  + year2001 + year2002 + year2003 + year2004 + year2005 + year2006 + year2007 + year2008 
                  + year2009 + year2010 + year2011 + year2012 + year2013 + year2014 + year2015 + year2016 
                  + region5_2 + region5_3 + region5_4 + region5_5, family=quasipoisson,method="fREML") 
    predict_try2_bam <- data.frame(predict(try2_bam, newdata = newdata_try2,type="terms", se.fit = TRUE))
    
    try3_bam<-bam(y ~ s(exposure3_star,bs='cr') + ozone_summer + no2 + temp + rh + PctEye + PctLDL + Pctmam + LungCancerRate + poverty + popdensity + medianhousevalue 
                  + pct_blk + medhouseholdincome + pct_owner_occ + hispanic + education + smoke_rate + mean_bmi + amb_visit_pct + a1c_exm_pct + nearest_hospital_km 
                  + year2001 + year2002 + year2003 + year2004 + year2005 + year2006 + year2007 + year2008 
                  + year2009 + year2010 + year2011 + year2012 + year2013 + year2014 + year2015 + year2016 
                  + region5_2 + region5_3 + region5_4 + region5_5, family=quasipoisson,method="fREML") 
    predict_try3_bam <- data.frame(predict(try3_bam, newdata = newdata_try3,type="terms", se.fit = TRUE))
    
    try4_bam<-bam(y ~ s(exposure4_star,bs='cr') + ozone_summer + no2 + temp + rh + PctEye + PctLDL + Pctmam + LungCancerRate + poverty + popdensity + medianhousevalue 
                  + pct_blk + medhouseholdincome + pct_owner_occ + hispanic + education + smoke_rate + mean_bmi + amb_visit_pct + a1c_exm_pct + nearest_hospital_km 
                  + year2001 + year2002 + year2003 + year2004 + year2005 + year2006 + year2007 + year2008 
                  + year2009 + year2010 + year2011 + year2012 + year2013 + year2014 + year2015 + year2016 
                  + region5_2 + region5_3 + region5_4 + region5_5, family=quasipoisson,method="fREML") 
    predict_try4_bam <- data.frame(predict(try4_bam, newdata = newdata_try4,type="terms", se.fit = TRUE))
    
    try5_bam<-bam(y ~ s(exposure5_star,bs='cr') + ozone_summer + no2 + temp + rh + PctEye + PctLDL + Pctmam + LungCancerRate + poverty + popdensity + medianhousevalue 
                  + pct_blk + medhouseholdincome + pct_owner_occ + hispanic + education + smoke_rate + mean_bmi + amb_visit_pct + a1c_exm_pct + nearest_hospital_km 
                  + year2001 + year2002 + year2003 + year2004 + year2005 + year2006 + year2007 + year2008 
                  + year2009 + year2010 + year2011 + year2012 + year2013 + year2014 + year2015 + year2016 
                  + region5_2 + region5_3 + region5_4 + region5_5, family=quasipoisson,method="fREML") 
    predict_try5_bam <- data.frame(predict(try5_bam, newdata = newdata_try5,type="terms", se.fit = TRUE))
    
    try6_bam<-bam(y ~ s(exposure6_star,bs='cr') + ozone_summer + no2 + temp + rh + PctEye + PctLDL + Pctmam + LungCancerRate + poverty + popdensity + medianhousevalue 
                  + pct_blk + medhouseholdincome + pct_owner_occ + hispanic + education + smoke_rate + mean_bmi + amb_visit_pct + a1c_exm_pct + nearest_hospital_km 
                  + year2001 + year2002 + year2003 + year2004 + year2005 + year2006 + year2007 + year2008 
                  + year2009 + year2010 + year2011 + year2012 + year2013 + year2014 + year2015 + year2016 
                  + region5_2 + region5_3 + region5_4 + region5_5, family=quasipoisson,method="fREML") 
    predict_try6_bam <- data.frame(predict(try6_bam, newdata = newdata_try6,type="terms", se.fit = TRUE))
    
    ###plot
    # x=seq(0,100,by=0.1)
    # plot(0, 0, col = "white", xlab="Annual PM2.5", ylab="Log(RR)", xlim=c(0,100), ylim=c(-0.1,0.8))
    # segments(x0 = 0, y0 = 0, x1 = threshold, y1 = 0) 
    # segments(x0 = threshold, y0 = 0, x1 = 100, y1 = b1*(100-threshold)) 
    # lines(x,predict_try1_bam$fit.s.exposure1_star.-(predict_try1_bam$fit.s.exposure1_star.[1]),type = "l", lty=2, col = "blue")
    # lines(x,predict_try2_bam$fit.s.exposure2_star.-(predict_try2_bam$fit.s.exposure2_star.[1]),type = "l", lty=2, col = "red")
    # lines(x,predict_try3_bam$fit.s.exposure3_star.-(predict_try3_bam$fit.s.exposure3_star.[1]),type = "l", lty=2, col = "green")
    # legend(60, 0.15, legend=c("exp", "exp+1sd", "exp+2sd", "exp+3sd"), cex=0.75,lty=c(1,2,2,2), col = c('black','blue','red','green'),border="white")
    
    #extract results
    results[((i-1)*6+1),(1:501)] <- predict_try1_bam$fit.s.exposure1_star.
    results[((i-1)*6+2),(1:501)] <- predict_try2_bam$fit.s.exposure2_star.
    results[((i-1)*6+3),(1:501)] <- predict_try3_bam$fit.s.exposure3_star.
    results[((i-1)*6+4),(1:501)] <- predict_try4_bam$fit.s.exposure4_star.
    results[((i-1)*6+5),(1:501)] <- predict_try5_bam$fit.s.exposure5_star.
    results[((i-1)*6+6),(1:501)] <- predict_try6_bam$fit.s.exposure6_star.
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
#read in the new covariates dataset
covar<-readRDS(paste0(dir_data,'counts_weight_pred_sd_0127.rds'))
#read in the new measurement error datasets 
key_files <- list.files(path=dir_data,pattern = "^counts_covar_simu_set_50_(.*)rds$")
key <- readRDS(paste0(dir_data,key_files[1]))

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

#summary(key[1:649910,]) #9984 obs missing, 1.54% 

#parameter setting 
n.reps=50 # 50 or 100, 200 (depending how many pieces we combine into a single loop simulation)

threshold <- 8
b1=0.005

start.time <- Sys.time()
test=simulate_results(key,covar,threshold,b1,n.reps)
end.time <- Sys.time()
end.time - start.time
# Time difference of 1.282171 hours



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
