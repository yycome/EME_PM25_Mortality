###############################################################################
# Project: Exposure measurement error                                         #
# Code: Step 1 - estimate annual-level uncertainty                            #
# Machine: QNAP                                                               #
# Author: Yaguang Wei                                                         #
###############################################################################

############################# 0. Setup ##############################
rm(list=ls())
gc()

require(magrittr)
require(lubridate)
require(readr)
require(haven)
require(dplyr)
require(data.table)
require(ranger)
library(caret)
library(raster)
library(spatialEco)
library(sf)
library(sp)
library(stringr)
library(nlme)

set.seed(1234)

dir_uncertainty_raw <- '/media/qnap2/assembled_data/prediction/PM25_USGrid_Uncertainty/'
dir_uncertainty_save <- '/media/qnap4/Yaguang/EME/data/uncertainty/'
dir_shp <- '/media/qnap4/Yaguang/ZIPCODE_INFO/polygon/'
dir_pobox <- '/media/qnap4/Yaguang/ZIPCODE_INFO/pobox_csv/'



################ 1. calculate SD for grids with monitoring sites ###############
InputData <- readRDS(paste0(dir_uncertainty_raw,"InputData_Step4_1.rds"))
OutputData <- readRDS(paste0(dir_uncertainty_raw,"OutputData_CV.rds"))
OutputData$Date <- as.Date(OutputData$CalendarDay,origin = as.Date("1970-01-01"))
OutputData$Year <- as.numeric(format(OutputData$Date,"%Y"))
OutputData$Month <- as.numeric(format(OutputData$Date,"%m"))
OutputData$MonitorData <- exp(OutputData$MonitorData)
OutputData$pred_ensemble_2 <- exp(OutputData$pred_ensemble_2)
OutputData$res <- OutputData$MonitorData - OutputData$pred_ensemble_2
OutputData <- OutputData[,c(1:3,11:15)]
OutputData$CalendarDay <- NULL
Data <- cbind(OutputData,InputData) #2156 monitors daily data

# compute SD per month & summarize other covariates
res_month <- aggregate(res ~ Year + Month + SiteCode,data = Data,FUN = mean, na.action = na.omit)
Data_month <- aggregate(cbind(Other_Lat,Other_Lon,pred_ensemble_2,USElevation_med100,USElevation_std100,RoadDensity_roads1000,NLCD_Impervious100, NLCD_Developed10000,PM25_Region,MOD09A1,REANALYSIS_shum_2m_DailyMean,NLCD_canopy100,MOD13A2_Nearest4) ~ Year + Month + SiteCode,
                        data = Data,FUN = mean,na.action = na.omit)
Data_uncertainty_month <- merge(x = Data_month,y = res_month,by.x = c("Year","Month","SiteCode"),by.y = c("Year","Month","SiteCode"),all.x = TRUE)

# Data_uncertainty_month$Month <- as.factor(Data_uncertainty_month$Month)
# Data_uncertainty_month$Year <- as.factor(Data_uncertainty_month$Year)
Data_uncertainty_month$PM25_Region <- as.factor(Data_uncertainty_month$PM25_Region)

Data_uncertainty_month_cov <- Data_uncertainty_month[,c("Year","Month","SiteCode","res")]
Data_uncertainty_month_cov <- Data_uncertainty_month_cov[complete.cases(Data_uncertainty_month_cov),]

# compute annual-level uncertainty for each monitored grid
monitored_grid_id <- as.character(unique(Data_uncertainty_month_cov$SiteCode))
for (i in monitored_grid_id) {
  Data_uncertainty_month_test <- Data_uncertainty_month[Data_uncertainty_month$SiteCode==i,c("Year","Month","SiteCode","res")]
  Data_uncertainty_month_test <- Data_uncertainty_month_test[complete.cases(Data_uncertainty_month_test),]
  cov_test <- lme(res~1,data = Data_uncertainty_month_test,random = reStruct(~1|Year),
                  corr = corAR1(form = ~Month | Year), weight = varIdent(form = ~ 1 | Month))
}


# # test heterogeneous AR1 model
Data_uncertainty_month_test <- Data_uncertainty_month[Data_uncertainty_month$SiteCode=='010030010',c("Year","Month","SiteCode","res")]
Data_uncertainty_month_test <- Data_uncertainty_month_test[complete.cases(Data_uncertainty_month_test),]
cov_test <- lme(res~1,data = Data_uncertainty_month_test,random = reStruct(~1|Year),
                corr = corAR1(form = ~Month | Year), weight = varIdent(form = ~ 1 | Month))
summary(cov_test)
# Linear mixed-effects model fit by REML
# Data: Data_uncertainty_month_test 
# AIC      BIC  logLik
# 422.5601 471.1863 -196.28
# 
# Random effects:
#   Formula: ~1 | Year
# (Intercept)  Residual
# StdDev:    0.266763 0.4829384
# 
# Correlation Structure: AR(1)
# Formula: ~Month | Year 
# Parameter estimate(s):
#   Phi 
# 0.5588569 
# Variance function:
#   Structure: Different standard deviations per stratum
# Formula: ~1 | Month 
# Parameter estimates:
#   1         2         3         4         5         6         7         8         9        10        11 
# 1.0000000 4.4462562 0.8166911 1.4149060 1.8080483 1.5035663 1.4751037 1.4874167 1.7763149 1.8375228 1.6505948 
# 12 
# 1.9193342 
# Fixed effects:  res ~ 1 
# Value  Std.Error  DF  t-value p-value
# (Intercept) 0.3463556 0.09218031 174 3.757371   2e-04
# 
# Standardized Within-Group Residuals:
#   Min          Q1         Med          Q3         Max 
# -3.55788183 -0.56046850 -0.02181684  0.41720198  3.20114235 
# 
# Number of Observations: 190
# Number of Groups: 16 

cov_test$sigma
# [1] 0.4829384

coef(cov_test$modelStruct$varStruct, uncons=FALSE)
#  2         3         4         5         6         7         8         9        10        11        12 
# 4.4462562 0.8166911 1.4149060 1.8080483 1.5035663 1.4751037 1.4874167 1.7763149 1.8375228 1.6505948 1.9193342

cov_test$modelStruct$corStruct
# Correlation structure of class corAR1 representing
# Phi 
# 0.5588569





sig2.vect <- cov_test$sigma^2*c(1,coef(cov_test$modelStruct$varStruct, uncons=FALSE))
rows <- matrix(rep(seq(1:12), 12),nrow=12)
columns <- t(rows)
abs.diffs <- abs(rows - columns)
corr.mat <- cov_test$modelStruct$corStruct^abs.diffs
sigma.mat <- diag(sqrt(sig2.vect))%*%corr.mat%*%diag(sqrt(sig2.vect))
avg.vector <- matrix(rep(1/12,12))
sigma2.yearly <- t(avg.vector)%*%sigma.mat%*%avg.vector

