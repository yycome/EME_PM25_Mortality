##################################################################################
# Project: Exposure measurement error                                            #
# Code: Step 2. sensitivity analysis - estimate spatial correlation coefficients # 
#               for monitored grids, link each zip code to those 0-10, 10-20,    # 
#               20-40, 40-100km around                                           #
# Machine: QNAP                                                                  #
##################################################################################

############################# 0. Setup ##############################
rm(list=ls())
gc()

library(utils)
library(data.table)
library(dplyr)
library(nngeo)
library(leaflet)
library(sf)

dir_monitor_correlation <- '/media/qnap4/Yaguang/EME/data/SpatialCorrelation/'



#################### 1. pair monitors and calculate distance ######################
# load annual monitored and predicted estimates
DATA <- readRDS(paste0(dir_monitor_correlation,'pm25_moni_pred_2000_2016.rds'))

DATA$res <- DATA$pm25_annual_moni-DATA$pm25_annual_pred
DATA$YEAR <- as.character(DATA$YEAR)
DATA$STUSPS <- as.character(DATA$STUSPS)

# label the monitor in order
monitor_geo <- DATA[c('SITECODE','LAT','LON','XCoord','YCoord','STUSPS')]
monitor_geo$index <- 1:nrow(monitor_geo)

# pair monitors by creating all combinations of rows
combinations <- t(combn(monitor_geo$index, 2, FUN = NULL))
monitor_geo <- as.data.table(monitor_geo)
combinations_monitor1 <- monitor_geo[combinations[,1],c("SITECODE","LAT","LON","index")]
combinations_monitor2 <- monitor_geo[combinations[,2],c("SITECODE","LAT","LON","index")]
comb_data_monitor <- cbind(combinations_monitor1,combinations_monitor2)
names(comb_data_monitor) <- c("SITECODE1","LAT1","LON1","index1","SITECODE2","LAT2","LON2","index2")

rm(combinations_monitor1,combinations_monitor2)
gc()

# function to calculate distance in kilometers between two points with longitude and latitude
earth.dist <- function (long1, lat1, long2, lat2)
{
  rad <- pi/180
  a1 <- lat1 * rad
  a2 <- long1 * rad
  b1 <- lat2 * rad
  b2 <- long2 * rad
  dlon <- b2 - a2
  dlat <- b1 - a1
  a <- (sin(dlat/2))^2 + cos(a1) * cos(b1) * (sin(dlon/2))^2
  c <- 2 * atan2(sqrt(a), sqrt(1 - a))
  R <- 6378.145
  d <- R * c
  return(d)
}

# calculate paired distance 
comb_data_monitor$dist <- earth.dist(comb_data_monitor$LON1,comb_data_monitor$LAT1,
                                     comb_data_monitor$LON2,comb_data_monitor$LAT2)
# summary(comb_data_monitor$dist)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0   925.6  1638.3  1937.8  2766.4 13761.6 



#################### 2. merge residual information ######################
# residual file
data_res <- DATA[,c('SITECODE','YEAR','res')]

# merge to monitor file for each year
years <- 2000:2016

ls_less10km <- list()
ls_10_20km <- list()
ls_20_40km <- list()
ls_40_100km <- list()

for (i in 1:length(years)) {
  comb_data_monitor_tmp <- comb_data_monitor %>% mutate(year=years[i])
  comb_data_monitor_tmp$year <- as.character(comb_data_monitor_tmp$year)
  comb_data_monitor_tmp <- left_join(comb_data_monitor_tmp,data_res,by=c('SITECODE1'='SITECODE','year'='YEAR'))
  comb_data_monitor_tmp <- left_join(comb_data_monitor_tmp,data_res,by=c('SITECODE2'='SITECODE','year'='YEAR'))
  names(comb_data_monitor_tmp)[which(names(comb_data_monitor_tmp)=='res.x')] <- 'res1'
  names(comb_data_monitor_tmp)[which(names(comb_data_monitor_tmp)=='res.y')] <- 'res2'
  
  # create separate files for dist<=10km, 10-20km, 20-40km, and 40-100km
  comb_data_monitor_tmp_less10km <- comb_data_monitor_tmp %>% filter(dist<=10)
  comb_data_monitor_tmp_less10km$dist_gp <- 'less10km'
  comb_data_monitor_tmp_less10km <- comb_data_monitor_tmp_less10km %>% group_by(SITECODE1,year) %>% mutate(res2_avg=mean(res2,na.rm=TRUE))
  comb_data_monitor_tmp_less10km <- comb_data_monitor_tmp_less10km[,c('SITECODE1','year','res1','res2_avg','dist_gp')]
  comb_data_monitor_tmp_less10km <- unique(comb_data_monitor_tmp_less10km)
  
  comb_data_monitor_tmp_10_20km <- comb_data_monitor_tmp %>% filter(dist>10 & dist<=20)
  comb_data_monitor_tmp_10_20km$dist_gp <- '10_20km'
  comb_data_monitor_tmp_10_20km <- comb_data_monitor_tmp_10_20km %>% group_by(SITECODE1,year) %>% mutate(res2_avg=mean(res2,na.rm=TRUE))
  comb_data_monitor_tmp_10_20km <- comb_data_monitor_tmp_10_20km[,c('SITECODE1','year','res1','res2_avg','dist_gp')]
  comb_data_monitor_tmp_10_20km <- unique(comb_data_monitor_tmp_10_20km)
  
  comb_data_monitor_tmp_20_40km <- comb_data_monitor_tmp %>% filter(dist>20 & dist<=40)
  comb_data_monitor_tmp_20_40km$dist_gp <- '20_40km'
  comb_data_monitor_tmp_20_40km <- comb_data_monitor_tmp_20_40km %>% group_by(SITECODE1,year) %>% mutate(res2_avg=mean(res2,na.rm=TRUE))
  comb_data_monitor_tmp_20_40km <- comb_data_monitor_tmp_20_40km[,c('SITECODE1','year','res1','res2_avg','dist_gp')]
  comb_data_monitor_tmp_20_40km <- unique(comb_data_monitor_tmp_20_40km)
  
  comb_data_monitor_tmp_40_100km <- comb_data_monitor_tmp %>% filter(dist>40 & dist<=100)
  comb_data_monitor_tmp_40_100km$dist_gp <- '40_100km'
  comb_data_monitor_tmp_40_100km <- comb_data_monitor_tmp_40_100km %>% group_by(SITECODE1,year) %>% mutate(res2_avg=mean(res2,na.rm=TRUE))
  comb_data_monitor_tmp_40_100km <- comb_data_monitor_tmp_40_100km[,c('SITECODE1','year','res1','res2_avg','dist_gp')]
  comb_data_monitor_tmp_40_100km <- unique(comb_data_monitor_tmp_40_100km)
  
  ls_less10km[[i]] <- comb_data_monitor_tmp_less10km
  ls_10_20km[[i]] <- comb_data_monitor_tmp_10_20km
  ls_20_40km[[i]] <- comb_data_monitor_tmp_20_40km
  ls_40_100km[[i]] <- comb_data_monitor_tmp_40_100km
  
  rm(comb_data_monitor_tmp,comb_data_monitor_tmp_less10km,comb_data_monitor_tmp_10_20km,comb_data_monitor_tmp_20_40km,comb_data_monitor_tmp_40_100km)
  gc()
  
  print(paste0(years[i], " done"))
}

comb_data_monitor_0016_less10km <- rbindlist(ls_less10km)
comb_data_monitor_0016_10_20km <- rbindlist(ls_10_20km)
comb_data_monitor_0016_20_40km <- rbindlist(ls_20_40km)
comb_data_monitor_0016_40_100km <- rbindlist(ls_40_100km)

comb_data_monitor_0016_dist_all <- full_join(comb_data_monitor_0016_less10km,comb_data_monitor_0016_10_20km,
                                             by=c('SITECODE1'='SITECODE1','year'='year'))
comb_data_monitor_0016_dist_all <- full_join(comb_data_monitor_0016_dist_all,comb_data_monitor_0016_20_40km,
                                             by=c('SITECODE1'='SITECODE1','year'='year'))
comb_data_monitor_0016_dist_all <- full_join(comb_data_monitor_0016_dist_all,comb_data_monitor_0016_40_100km,
                                             by=c('SITECODE1'='SITECODE1','year'='year'))

names(comb_data_monitor_0016_dist_all) <- c("SITECODE1","year","res1_less10km","res2_avg_less10km","dist_gp_less10km",
                                            "res1_10_20km","res2_avg_10_20km","dist_gp_10_20km",
                                            "res1_20_40km","res2_avg_20_40km","dist_gp_20_40km",
                                            "res1_40_100km","res2_avg_40_100km","dist_gp_40_100km")

saveRDS(comb_data_monitor_0016_dist_all,file=paste0(dir_monitor_correlation,'comb_data_monitor_0016_dist_all_sensitivity.rds'))



#################### 3. coefficients for spatial correlation among monitored sites ######################
summary(lm(res1_less10km~res2_avg_less10km+res2_avg_10_20km+res2_avg_20_40km+res2_avg_40_100km, 
           na.action=na.omit,data=comb_data_monitor_0016_dist_all))
# Call:
#   lm(formula = res1_less10km ~ res2_avg_less10km + res2_avg_10_20km + 
#        res2_avg_20_40km + res2_avg_40_100km, data = comb_data_monitor_0016_dist_all, 
#      na.action = na.omit)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -10.9989  -0.1570  -0.0062   0.1284  14.4026 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
#   (Intercept)       -0.008100   0.009895  -0.819  0.41306    
#   res2_avg_less10km  0.867077   0.011773  73.647  < 2e-16 ***
#   res2_avg_10_20km   0.029629   0.012639   2.344  0.01910 *  
#   res2_avg_20_40km   0.038560   0.014477   2.664  0.00775 ** 
#   res2_avg_40_100km  0.017282   0.019733   0.876  0.38118    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.7115 on 5684 degrees of freedom
# (30198 observations deleted due to missingness)
# Multiple R-squared:  0.5249,	Adjusted R-squared:  0.5246 
# F-statistic:  1570 on 4 and 5684 DF,  p-value: < 2.2e-16



#################### 4. pair zip codes ######################
dat_zip <- readRDS(paste0(dir_monitor_correlation,'zip_geocode_new.rds'))
dat_zip <- dat_zip[!is.na(dat_zip$Latitude),]

# convert to sf object
dat_zip <- st_as_sf(dat_zip,coords=c("Longitude","Latitude"),crs=4326) # 4326 Stands for WGS84 geographic coordinate system (crs)
# transform to US National Atlas projection crs 
# (it is recommended to run spatial join using sf objects with a projected coordinate system and not a geographic crs) 
dat_zip <- st_transform(dat_zip,2163)  # 2163 is the EPSG code for US National Atlas

# create paired zip codes for all years
years <- 2000:2016
# <=10km
join_zip_less10km <- st_join(dat_zip,dat_zip,join=st_nn,maxdist=10000,k=1000) #maxdist in meters
join_zip_less10km$geometry <- NULL
join_zip_less10km <- join_zip_less10km[join_zip_less10km$zip.x!=join_zip_less10km$zip.y,]
ls_zip_less10km <- list()
for (i in 1:length(years)) {
  join_zip_less10km_tmp <- join_zip_less10km
  join_zip_less10km_tmp$year <- years[i]
  ls_zip_less10km[[i]] <- join_zip_less10km_tmp
  rm(join_zip_less10km_tmp)
  gc()
}
join_zip_less10km_all <- rbindlist(ls_zip_less10km)
names(join_zip_less10km_all) <- c('zipcode','zipcode_near','year')
join_zip_less10km_all$year <- as.character(join_zip_less10km_all$year)
join_zip_less10km_all <- unique(join_zip_less10km_all)
saveRDS(join_zip_less10km_all,paste0(dir_monitor_correlation,"join_zip_less10km_all.rds"))

# 10-20km
join_zip_10_20km <- st_join(dat_zip,dat_zip,join=st_nn,maxdist=20000,k=1000) #maxdist in meters
join_zip_10_20km$geometry <- NULL
join_zip_10_20km <- join_zip_10_20km[join_zip_10_20km$zip.x!=join_zip_10_20km$zip.y,]
join_zip_10_20km <- join_zip_10_20km %>% anti_join(join_zip_less10km)  # exclude those <10km
ls_zip_10_20km <- list()
for (i in 1:length(years)) {
  join_zip_10_20km_tmp <- join_zip_10_20km
  join_zip_10_20km_tmp$year <- years[i]
  ls_zip_10_20km[[i]] <- join_zip_10_20km_tmp
  rm(join_zip_10_20km_tmp)
  gc()
}
join_zip_10_20km_all <- rbindlist(ls_zip_10_20km)
names(join_zip_10_20km_all) <- c('zipcode','zipcode_near','year')
join_zip_10_20km_all$year <- as.character(join_zip_10_20km_all$year)
join_zip_10_20km_all <- unique(join_zip_10_20km_all)
saveRDS(join_zip_10_20km_all,paste0(dir_monitor_correlation,"join_zip_10_20km_all.rds"))

# 20-40km
join_zip_20_40km <- st_join(dat_zip,dat_zip,join=st_nn,maxdist=40000,k=1000) #maxdist in meters
join_zip_20_40km$geometry <- NULL
join_zip_20_40km <- join_zip_20_40km[join_zip_20_40km$zip.x!=join_zip_20_40km$zip.y,]
join_zip_20_40km <- join_zip_20_40km %>% anti_join(join_zip_less10km)  # exclude those <10km
join_zip_20_40km <- join_zip_20_40km %>% anti_join(join_zip_10_20km)  # exclude those between 10-20km
ls_zip_20_40km <- list()
for (i in 1:length(years)) {
  join_zip_20_40km_tmp <- join_zip_20_40km
  join_zip_20_40km_tmp$year <- years[i]
  ls_zip_20_40km[[i]] <- join_zip_20_40km_tmp
  rm(join_zip_20_40km_tmp)
  gc()
}
join_zip_20_40km_all <- rbindlist(ls_zip_20_40km)
names(join_zip_20_40km_all) <- c('zipcode','zipcode_near','year')
join_zip_20_40km_all$year <- as.character(join_zip_20_40km_all$year)
join_zip_20_40km_all <- unique(join_zip_20_40km_all)
saveRDS(join_zip_20_40km_all,paste0(dir_monitor_correlation,"join_zip_20_40km_all.rds"))

# 40-100km
join_zip_40_100km <- st_join(dat_zip,dat_zip,join=st_nn,maxdist=100000,k=1000) #maxdist in meters
join_zip_40_100km$geometry <- NULL
join_zip_40_100km <- join_zip_40_100km[join_zip_40_100km$zip.x!=join_zip_40_100km$zip.y,]
join_zip_40_100km <- join_zip_40_100km %>% anti_join(join_zip_less10km)  # exclude those <10km
join_zip_40_100km <- join_zip_40_100km %>% anti_join(join_zip_10_20km)  # exclude those between 10-20km
join_zip_40_100km <- join_zip_40_100km %>% anti_join(join_zip_20_40km)  # exclude those between 20-40km
ls_zip_40_100km <- list()
for (i in 1:length(years)) {
  join_zip_40_100km_tmp <- join_zip_40_100km
  join_zip_40_100km_tmp$year <- years[i]
  ls_zip_40_100km[[i]] <- join_zip_40_100km_tmp
  rm(join_zip_40_100km_tmp)
  gc()
}
join_zip_40_100km_all <- rbindlist(ls_zip_40_100km)
names(join_zip_40_100km_all) <- c('zipcode','zipcode_near','year')
join_zip_40_100km_all$year <- as.character(join_zip_40_100km_all$year)
join_zip_40_100km_all <- unique(join_zip_40_100km_all)
saveRDS(join_zip_40_100km_all,paste0(dir_monitor_correlation,"join_zip_40_100km_all.rds"))
