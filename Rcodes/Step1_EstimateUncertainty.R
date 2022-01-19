###############################################################################
# Project: Exposure measurement error                                         #
# Code: Step 1 - estimate exposure uncertainty                                #
# Machine: QNAP                                                               #
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
library(sf)
library(sp)
library(stringr)
library(tseries)
library(spatialEco)

set.seed(1234)

dir_uncertainty_raw <- '~/PM25_USGrid_Uncertainty/'
dir_uncertainty_save <- '~/data/Uncertainty/'
dir_shp <- '~/ZIPCODE_INFO/polygon/'
dir_pobox <- '~/ZIPCODE_INFO/pobox_csv/'
dir_results_save <- '~/results/'



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

# # check autocorrelation & plot Fig S1
# Data_acf <- Data[,c("SiteCode","Date","res")]
# Data_acf <- Data_acf[complete.cases(Data_acf),]
# Data_acf$doy <- yday(Data_acf$Date)
# Data_acf_group_n <- Data_acf %>% group_by(SiteCode) %>% summarise(n=n())
# summary(Data_acf_group_n$n)
# # Min. 1st Qu.  Median  Mean   3rd Qu.  Max.
# # 4.0   372.5  1088.0  1334.5  1826.5  5728.0
# Data_acf_group_n[Data_acf_group_n$n==max(Data_acf_group_n$n,na.rm=TRUE),]
# # # A tibble: 1 × 2
# #   SiteCode      n
# #   <fct>     <int>
# # 1 220330009  5728
# Data_acf_group <- Data_acf %>% group_by(doy) %>% summarise(mean_res=mean(res))
# png(filename = paste0(dir_results_save,"Fig_S1_new.png"), height = 4, width = 6, units = "in",res = 300)
# acf(Data_acf_group$mean_res,lag.max = 100,main='',ylab='Autocorrelation')
# dev.off()

# compute SD per year for every other day & summarize other covariates
Data_res <- Data[,c("SiteCode","Date","Year","Month","res")]
Data_res <- Data_res[order(Data_res$SiteCode,Data_res$Date),]
Data_res <- Data_res %>% group_by(SiteCode) %>% slice(seq(1, n(), by = 2)) %>% ungroup()
Data_res_sd <- aggregate(res ~ Year + SiteCode,data = Data_res,FUN = sd, na.action = na.omit)
Data_res_sd$res <- Data_res_sd$res/sqrt((365/2))
names(Data_res_sd)[which(names(Data_res_sd) == "res")] <- "res_sd"  # rename

Data_cov <- aggregate(cbind(Other_Lat,Other_Lon,pred_ensemble_2,USElevation_med100,USElevation_std100,RoadDensity_roads1000,NLCD_Impervious100, NLCD_Developed10000,PM25_Region,MOD09A1,REANALYSIS_shum_2m_DailyMean,NLCD_canopy100,MOD13A2_Nearest4) ~ Year + SiteCode,
                      data = Data,FUN = mean,na.action = na.omit)

Data_uncertainty <- merge(x = Data_cov,y = Data_res_sd,by.x = c("Year","SiteCode"),by.y = c("Year","SiteCode"),all.x = TRUE)

Data_uncertainty$Year <- as.factor(Data_uncertainty$Year)
Data_uncertainty$PM25_Region <- as.factor(Data_uncertainty$PM25_Region)
Data_uncertainty <- Data_uncertainty[complete.cases(Data_uncertainty),]

saveRDS(Data_uncertainty,file = paste0(dir_uncertainty_save,'Data_uncertainty.rds'))



######################## 2. variable selection & tuning parameters #######################
Data_uncertainty <- readRDS(paste0(dir_uncertainty_save,'Data_uncertainty.rds'))

### Variable selection
rf_rfe <- rfe(res_sd~.-SiteCode,
              data = Data_uncertainty, 
              method = "ranger",
              rfeControl = rfeControl(functions=rfFuncs, method="cv", number=5))

# list the chosen variables
predictors(rf_rfe)
# [1] "NLCD_canopy100"               "pred_ensemble_2"              "MOD09A1"                     
# [4] "NLCD_Developed10000"          "Other_Lat"                    "MOD13A2_Nearest4"
# [7] "USElevation_med100"           "Year"                         "REANALYSIS_shum_2m_DailyMean"
# [10] "Other_Lon"                    "PM25_Region"                                      


### Tuning parameters
rf_grid <- expand.grid(mtry = c(2, 3, 4, 5),
                       splitrule = "extratrees",
                       min.node.size = c(1, 5, 10, 50, 100))
modellist <- list()
for (ntree in c(200,300,500)){
  fit <- train(res_sd~NLCD_canopy100+pred_ensemble_2+MOD09A1+NLCD_Developed10000+Other_Lat+
                 MOD13A2_Nearest4+USElevation_med100+Year+REANALYSIS_shum_2m_DailyMean+
                 Other_Lon+PM25_Region,
               data = Data_uncertainty, 
               method = "ranger",
               trControl = trainControl(method = "cv",number = 5), # 5-fold CV
               tuneGrid = rf_grid,
               num.trees = ntree,
               importance = "impurity")
  key <- toString(ntree)
  modellist[[key]] <- fit
}
min(modellist$`200`$results$RMSE)
# [1] 1.220879
min(modellist$`300`$results$RMSE)
# [1] 1.220507
min(modellist$`500`$results$RMSE)
# [1] 1.221235

modellist$`300`$bestTune
#    mtry  splitrule min.node.size
# 17    5 extratrees             5



######################## 3. fit rf #######################
Data_uncertainty <- readRDS(paste0(dir_uncertainty_save,'Data_uncertainty.rds'))

# fit rf
rf_uncertainty <- ranger(res_sd~NLCD_canopy100+pred_ensemble_2+MOD09A1+NLCD_Developed10000+Other_Lat+
                           MOD13A2_Nearest4+USElevation_med100+Year+REANALYSIS_shum_2m_DailyMean+
                           Other_Lon+PM25_Region,
                         data=Data_uncertainty, num.trees = 300, mtry = 5, min.node.size = 5)

# compare predicted vs observed uncertainties
Data_uncertainty$pred_res_sd <- predict(rf_uncertainty,Data_uncertainty,type="response")$predictions
summary(Data_uncertainty[,c('res_sd','pred_res_sd')])
#     res_sd          pred_res_sd     
# Min.   : 0.01337   Min.   : 0.4026  
# 1st Qu.: 1.23370   1st Qu.: 1.3782  
# Median : 1.69937   Median : 1.7986  
# Mean   : 2.07659   Mean   : 2.0858  
# 3rd Qu.: 2.44344   3rd Qu.: 2.4413  
# Max.   :23.23642   Max.   :14.5050  

quantile(Data_uncertainty$pred_res_sd,probs=c(0.0001,0.9999))
# 0.01%     99.99% 
# 0.4800794 12.2951850 



###################### 4. predict grid-level uncertainty ###############
years_char <- c('00','01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16')

for (i in years_char) {
  if (!file.exists(paste0(dir_uncertainty_save,"PredUncertainty_PM25_USGrid_20",i,".rds"))){
    files_tmp <- Sys.glob(paste0(dir_uncertainty_raw,"InputUncertainty_PM25_USGrid_20",i,"*.rds"))
    datafiles_tmp <- lapply(files_tmp,readRDS)
    data_tmp <- rbindlist(datafiles_tmp)
    data_annual_tmp <- aggregate(cbind(NLCD_canopy100,pred_ensemble_2,MOD09A1,NLCD_Developed10000,Other_Lat,MOD13A2_Nearest4,USElevation_med100,REANALYSIS_shum_2m_DailyMean,Other_Lon,PM25_Region) ~ Year + SiteCode,
                                 data = data_tmp,FUN = mean,na.action = na.omit)
    data_annual_tmp$Year <- as.factor(data_annual_tmp$Year)
    data_annual_tmp$PM25_Region <- as.factor(data_annual_tmp$PM25_Region)
    data_annual_tmp$pred_res_sd <- predict(rf_uncertainty,data_annual_tmp,type="response")$predictions
    data_annual_tmp <- data_annual_tmp[,c("Other_Lat","Other_Lon","PM25_Region","SiteCode","Year","pred_res_sd")]
    saveRDS(data_annual_tmp,paste0(dir_uncertainty_save,"PredUncertainty_PM25_USGrid_20",i,".rds"))
    rm(files_tmp,datafiles_tmp,data_tmp,data_annual_tmp)
    gc()
  }
  print(paste0("20", i, " done"))
}



###################### 5. aggregate to zip code-level uncertainty ###############
SiteData <- readRDS('/media/qnap3/Exposure modeling/3 National scale/USA/4 PM2.5 v2.2000_16/7 Final predictions 1 km/Annual/USGridSite.rds')
sf_use_s2(FALSE) # make zip code geometry valid without using S2

years_char <- c('00','01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16')

for (i in years_char) {
  if (!file.exists(paste0(dir_uncertainty_save,"PredUncertainty_PM25_USZIP_20",i,".rds"))){
    file_tmp <- readRDS(paste0(dir_uncertainty_save,"PredUncertainty_PM25_USGrid_20",i,".rds"))
    file_tmp <- left_join(file_tmp,SiteData,by='SiteCode')
    file_tmp <- file_tmp[,c("Lon","Lat","pred_res_sd")]
    ### aggregate for standard zipcode
    file_tmp_zipcode <- file_tmp
    coordinates(file_tmp_zipcode)<-~Lon+Lat
    cty <- shapefile(paste0(dir_shp,"ESRI",i,"USZIP5_POLY_WGS84.shp"))
    proj4string(file_tmp_zipcode)<-proj4string(cty)
    file_tmp_zipcode <- point.in.poly(file_tmp_zipcode, cty)
    file_tmp_zipcode <- as.data.frame(file_tmp_zipcode)[,c('ZIP','pred_res_sd')]
    file_tmp_zipcode <- file_tmp_zipcode %>% dplyr::group_by(ZIP) %>% dplyr::summarise(pred_res_sd=mean(pred_res_sd,na.rm=T))
    file_tmp_zipcode <- file_tmp_zipcode[!is.na(file_tmp_zipcode$ZIP),]
    ### aggregate for po box
    file_tmp_pobox <- file_tmp
    pobox <- read.csv(paste0(dir_pobox,'ESRI',i,'USZIP5_POINT_WGS84_POBOX.csv'))
    pobox <- pobox[,c(1,3,4)]
    names(pobox) <- c('ZIP','Lon','Lat')
    link <- nabor::knn(file_tmp_pobox[,c("Lon","Lat")],pobox[,c("Lon","Lat")],k=1,radius=2*sqrt(2)/10*0.1)
    link <- cbind.data.frame(link$nn.idx,link$nn.dists)
    names(link) <- c("id","dis")
    pobox$pobox_id <- 1:dim(pobox)[1]
    pobox <- cbind.data.frame(pobox,link)
    file_tmp_pobox$id <- 1:dim(file_tmp_pobox)[1]
    pobox <- left_join(pobox,file_tmp_pobox,by=c("id"))
    pobox <- pobox[pobox$dis!=Inf,]
    pobox <- pobox[,c('ZIP','pred_res_sd')]
    pobox <- pobox[!is.na(pobox$ZIP),]
    pobox$ZIP <- as.character(pobox$ZIP)
    pobox$ZIP <- str_pad(pobox$ZIP, width=5, side="left", pad="0")
    ### merge and save
    file_tmp_combine <- rbind(file_tmp_zipcode,pobox)
    file_tmp_combine <- file_tmp_combine[order(file_tmp_combine$ZIP),]
    saveRDS(file_tmp_combine,file = paste0(dir_uncertainty_save,'PredUncertainty_PM25_USZIP_20',i,'.rds'))
    ### release memory
    rm(file_tmp,file_tmp_zipcode,cty,file_tmp_pobox,pobox,link,file_tmp_combine)
    gc()
  }
  print(paste0("20", i, " done"))
}



###################### 6. Relative frequency for uncertainty ########################
uncertainty_files <- list.files(path=dir_uncertainty_save,pattern = "^PredUncertainty_PM25_USZIP_20(.*)rds$")
dat_uncertainty <- readRDS(paste0(dir_uncertainty_save,uncertainty_files[1]))
dat_uncertainty$ZIP <- NULL
for (i in 2:length(uncertainty_files)) {
  dat_uncertainty_tmp <- readRDS(paste0(dir_uncertainty_save,uncertainty_files[i]))
  dat_uncertainty_tmp$ZIP <- NULL
  dat_uncertainty <- bind_rows(dat_uncertainty,dat_uncertainty_tmp)
  rm(dat_uncertainty_tmp)
  gc()
  print(paste0("uncertainty file ",i," is done"))
}
hist(dat_uncertainty$pred_res_sd, col=rgb(0,0,1,0.2), xlab='Estimated uncertainty', ylab='Relative frequency',
     main='', freq=FALSE)


