###############################################################################
# Project: Exposure measurement error                                         #
# Code: compile results                                                       #
# Machine: local                                                              #
###############################################################################

############################# 0. Setup ##############################
rm(list=ls())
gc()

library(stats)
library(MASS)
library(dplyr)
library(fastDummies)
library(latex2exp)
library(data.table)
library(stringr)
library(nngeo)
library(ggplot2)
library(fst)

dir_data <- '/medizza/qnap4/Yaguang/EME/data/SimulationData/'
dir_uncertainty <- '/media/qnap4/Yaguang/EME/data/Uncertainty/'
dir_results_linear <- '/media/qnap4/Yaguang/EME/results/Linear/'
dir_results_linear_lowlevel <- '/media/qnap4/Yaguang/EME/results/Linear_LowLevel/'
dir_results_linear_spatial100km <- '/media/qnap4/Yaguang/EME/results/Linear_Spatial100km/'
dir_results_quadratic_quadratic <- '/media/qnap4/Yaguang/EME/results/Quadratic_Quadratic/'
dir_results_quadratic_penalized <- '/media/qnap4/Yaguang/EME/results/Quadratic_Penalized/'
dir_results_softplus <- '/media/qnap4/Yaguang/EME/results/Softplus/'
dir_results_save <- '/media/qnap4/Yaguang/EME/results/'



###################### 1. Descriptive statistics ########################
### summarize PM2.5
dat_original <- readRDS(paste0(dir_data,'counts_weight_pred_sd_20211221.rds'))

dim(dat_original)
# [1] 649910     31
summary(dat_original$pm25)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 0.00783  7.52295  9.60506  9.65323 11.78815 30.92493
sd(dat_original$pm25)
# [1] 3.243062
quantile(dat_original$pm25,c(0.01,0.99))
# 1%       99%
# 2.713217 17.035318

### Calculate ratio of variance of error-prone vs variance of error-free exposure
key_files <- list.files(path=dir_data,pattern = "^counts_covar_simu_set_50_(.*)rds$")
dat_errors <- readRDS(paste0(dir_data,key_files[1]))
dat_errors$zip <- NULL
dat_errors$year <- NULL
for (i in 2:length(key_files)) {
  dat_errors_tmp <- readRDS(paste0(dir_data,key_files[i]))
  dat_errors_tmp$zip <- NULL
  dat_errors_tmp$year <- NULL
  dat_errors <- bind_rows(dat_errors,dat_errors_tmp)
  rm(dat_errors_tmp)
  gc()
  print(paste0("Key file ",i," is done"))
}

dat_errors$pm25_error_sp_sd <- dat_errors$pm25+dat_errors$annual_sp_err_sd
dat_errors$pm25_error_sp_2sd <- dat_errors$pm25+dat_errors$annual_sp_err_2sd
dat_errors$pm25_error_sp_3sd <- dat_errors$pm25+dat_errors$annual_sp_err_3sd

sd(dat_errors$pm25,na.rm=TRUE)
# [1] 3.243059
sd(dat_errors$pm25_error_sp_sd,na.rm=TRUE)
# [1] 3.245377
sd(dat_errors$pm25_error_sp_2sd,na.rm=TRUE)
# [1] 3.255555
sd(dat_errors$pm25_error_sp_3sd,na.rm=TRUE)
# [1] 3.2725

(3.2725/3.243059)^2
# [1] 1.018239
(3.243059/3.2725)^2
# [1] 0.982088



###################### 2. Tables 2, S4, S5. Biases of linear cr model ########################
# main analysis and single-pollutant analysis
results_linear_files_list <- Sys.glob(paste0(dir_results_linear,"results*.rds"))
results_linear_files <- lapply(results_linear_files_list,readRDS)
results_linear_raw <- rbindlist(results_linear_files)
results_linear <- results_linear_raw %>% group_by(type,corr_magnitude,pollutants,b0,b1) %>% 
  summarize(b0_est_avg=mean(b0_est),b1_est_avg=mean(b1_est)) %>% 
  ungroup()
results_linear$relative_bias <- (results_linear$b1_est_avg-results_linear$b1)/results_linear$b1*100


# low-level analysis
results_linear_lowlevel_files_list <- Sys.glob(paste0(dir_results_linear_lowlevel,"results*.rds"))
results_linear_lowlevel_files <- lapply(results_linear_lowlevel_files_list,readRDS)
results_linear_lowlevel_raw <- rbindlist(results_linear_lowlevel_files)
results_linear_lowlevel <- results_linear_lowlevel_raw %>% group_by(type,corr_magnitude,pollutants,b0,b1) %>% 
  summarize(b0_est_avg=mean(b0_est),b1_est_avg=mean(b1_est)) %>% 
  ungroup()
results_linear_lowlevel$relative_bias <- (results_linear_lowlevel$b1_est_avg-results_linear_lowlevel$b1)/results_linear_lowlevel$b1*100

# spatial correlation up to 100km
results_linear_spatial100km_files_list <- Sys.glob(paste0(dir_results_linear_spatial100km,"results*.rds"))
results_linear_spatial100km_files <- lapply(results_linear_spatial100km_files_list,readRDS)
results_linear_spatial100km_raw <- rbindlist(results_linear_spatial100km_files)
results_linear_spatial100km <- results_linear_spatial100km_raw %>% group_by(type,corr_magnitude,pollutants,b0,b1) %>% 
  summarize(b0_est_avg=mean(b0_est),b1_est_avg=mean(b1_est)) %>% 
  ungroup()
results_linear_spatial100km$relative_bias <- (results_linear_spatial100km$b1_est_avg-results_linear_spatial100km$b1)/results_linear_spatial100km$b1*100



##################### 3. Fig 2. Density of PM2.5 ########################
covar<-readRDS(paste0(dir_data,'counts_weight_pred_sd_0127.rds'))
pm25 <- covar$pm25[covar$pm25<=30]
### create plot
pdf(paste0(dir_results_save,'Fig_2.pdf'), height = 4, width = 6)
hist(pm25, breaks=25, col=rgb(0,0,0,0.2), xlab=expression(paste("Annual PM"[2.5]," concentration, ",mu,"g/m"^3)), 
     ylim=c(0,0.16),ylab='Density',main='',freq=FALSE)
lines(density(pm25,na.rm=TRUE),lwd=1.5)
dev.off()



###################### 4. Fig 3. Density of errors ########################
key_files <- list.files(path=dir_data,pattern = "^counts_covar_simu_set_50_(.*)rds$")
dat_errors <- readRDS(paste0(dir_data,key_files[1]))
dat_errors$zip <- NULL
dat_errors$year <- NULL
for (i in 2:length(key_files)) {
  dat_errors_tmp <- readRDS(paste0(dir_data,key_files[i]))
  dat_errors_tmp$zip <- NULL
  dat_errors_tmp$year <- NULL
  dat_errors <- bind_rows(dat_errors,dat_errors_tmp)
  rm(dat_errors_tmp)
  gc()
  print(paste0("Key file ",i," is done"))
}

annual_tempcorr_adj_err <- dat_errors$annual_ind_err_sd[dat_errors$annual_ind_err_sd>=-2&dat_errors$annual_ind_err_sd<=2]
annual_tempcorr_adj_err_sp <- dat_errors$annual_sp_err_sd[dat_errors$annual_sp_err_sd>=-2&dat_errors$annual_sp_err_sd<=2]
annual_tempcorr_adj_err_2sd <- dat_errors$annual_ind_err_2sd[dat_errors$annual_ind_err_2sd>=-2&dat_errors$annual_ind_err_2sd<=2]
annual_tempcorr_adj_err_2sd_sp <- dat_errors$annual_sp_err_2sd[dat_errors$annual_sp_err_2sd>=-2&dat_errors$annual_sp_err_2sd<=2]
annual_tempcorr_adj_err_3sd <- dat_errors$annual_ind_err_3sd[dat_errors$annual_ind_err_3sd>=-2&dat_errors$annual_ind_err_3sd<=2]
annual_tempcorr_adj_err_3sd_sp <- dat_errors$annual_sp_err_3sd[dat_errors$annual_sp_err_3sd>=-2&dat_errors$annual_sp_err_3sd<=2]

### create plot
pdf(paste0(dir_results_save,'Fig_3.pdf'), height = 4, width = 7.5)
layout(matrix(c(1,2,3,4,4,4),nrow = 2,ncol = 3,byrow = TRUE),heights = c(0.4,0.08))

plot(density(annual_tempcorr_adj_err,na.rm=TRUE),lty=5,xlim=c(-2,2),ylim=c(0,5),lwd=1.2,cex.lab=1.2, cex.main=1.2, cex.sub=1.2,
     xlab=expression(paste("Exposure measurement error, ",mu,"g/m"^3)), 
     ylab='Density',
     main=TeX("Magnitude of error: $\\sigma_{ik}=\\sigma^{*}_{ik}$"))
lines(density(annual_tempcorr_adj_err_sp,na.rm=TRUE), lty=1,lwd=1)

plot(density(annual_tempcorr_adj_err_2sd,na.rm=TRUE),lty=5,xlim=c(-2,2),ylim=c(0,5),lwd=1.2,cex.lab=1.2, cex.main=1.2, cex.sub=1.2,
     xlab=expression(paste("Exposure measurement error, ",mu,"g/m"^3)), 
     ylab='Density',
     main=TeX("Magnitude of error: $\\sigma_{ik}=2\\sigma^{*}_{ik}$"))
lines(density(annual_tempcorr_adj_err_2sd_sp,na.rm=TRUE), lty=1,lwd=1)

plot(density(annual_tempcorr_adj_err_3sd,na.rm=TRUE),lty=5,xlim=c(-2,2),ylim=c(0,5),lwd=1.2,cex.lab=1.2, cex.main=1.2, cex.sub=1.2,
     xlab=expression(paste("Exposure measurement error, ",mu,"g/m"^3)), 
     ylab='Density',
     main=TeX("Magnitude of error: $\\sigma_{ik}=3\\sigma^{*}_{ik}$"))
lines(density(annual_tempcorr_adj_err_3sd_sp,na.rm=TRUE), lty=1,lwd=1)

par(mar = c(0,0,0,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend(x = "center",x.intersp=0,
       legend=c(TeX("Spatially correlated error ($\\epsilon_{s,ik}$)"),TeX("Independent error ($\\epsilon_{ik}$)")),
       lty=c(1,5), lwd=c(1,1.2), horiz = TRUE, cex=1.2)
dev.off()



###################### 5. Fig 4. Quadratic cr model ########################
# ### reformat quadratic quadratic results
# files_quadratic_quadratic <- list.files(path=dir_results_quadratic_quadratic)
# for (i in 1:length(files_quadratic_quadratic)) {
#   results_tmp <- readRDS(paste0(dir_results_quadratic_quadratic,files_quadratic_quadratic[i]))
#   b1_tmp <- as.numeric(str_match(files_quadratic_quadratic[i], "b1\\s*(.*?)\\s*_b2")[2])
#   b2_tmp <- as.numeric(str_match(files_quadratic_quadratic[i], "b2\\s*(.*?)\\s*.rds")[2])
#   results_tmp$b1 <- b1_tmp
#   results_tmp$b2 <- b2_tmp
#   saveRDS(results_tmp,file = paste0(dir_results_quadratic_quadratic,files_quadratic_quadratic[i]))
#   rm(results_tmp,b1_tmp,b2_tmp)
#   gc()
# }

### extract quadratic quadratic results
files_quadratic_quadratic <- Sys.glob(paste0(dir_results_quadratic_quadratic,'results_set*.rds'))
files_quadratic_quadratic_ls <- lapply(files_quadratic_quadratic,readRDS)
results_quadratic_quadratic <- rbindlist(files_quadratic_quadratic_ls)
results_quadratic_quadratic$pollutants <- NULL
results_quadratic_quadratic <- results_quadratic_quadratic %>% dplyr::group_by(b0,b1,b2,type,corr_magnitude) %>% dplyr::summarise(b0_est=mean(b0_est,na.rm=TRUE),
                                                                                                                                  b1_est=mean(b1_est,na.rm=TRUE),
                                                                                                                                  b2_est=mean(b2_est,na.rm=TRUE))
results_quadratic_quadratic$corr_magnitude <- factor(results_quadratic_quadratic$corr_magnitude,levels=c('sp_error_sd','sp_error_2sd','sp_error_3sd',
                                                                                                         'ind_error_sd','ind_error_2sd','ind_error_3sd'))
results_quadratic_quadratic <- results_quadratic_quadratic[order(results_quadratic_quadratic$b0,results_quadratic_quadratic$b1,results_quadratic_quadratic$b2,
                                                                 results_quadratic_quadratic$type,results_quadratic_quadratic$corr_magnitude),]

# ### reformat quadratic penalized results
# files_quadratic_penalized <- list.files(path=dir_results_quadratic_penalized)
# for (i in 1:length(files_quadratic_penalized)) {
#   results_tmp <- readRDS(paste0(dir_results_quadratic_penalized,files_quadratic_penalized[i]))
#   b1_tmp <- as.numeric(str_match(files_quadratic_penalized[i], "b1\\s*(.*?)\\s*_b2")[2])
#   b2_tmp <- as.numeric(str_match(files_quadratic_penalized[i], "b2\\s*(.*?)\\s*.rds")[2])
#   results_tmp$b1 <- b1_tmp
#   results_tmp$b2 <- b2_tmp
#   saveRDS(results_tmp,file = paste0(dir_results_quadratic_penalized,files_quadratic_penalized[i]))
#   rm(results_tmp,b1_tmp,b2_tmp)
#   gc()
# }

### extract quadratic penalized results
files_quadratic_penalized <- Sys.glob(paste0(dir_results_quadratic_penalized,'results_set*.rds'))
files_quadratic_penalized_ls <- lapply(files_quadratic_penalized,readRDS)
results_quadratic_penalized <- rbindlist(files_quadratic_penalized_ls)
results_quadratic_penalized$pollutants <- NULL
results_quadratic_penalized <- results_quadratic_penalized %>% dplyr::group_by(b0,b1,b2,type,corr_magnitude) %>% dplyr::summarise_all("mean")
results_quadratic_penalized$corr_magnitude <- factor(results_quadratic_penalized$corr_magnitude,levels=c('sp_error_sd','sp_error_2sd','sp_error_3sd',
                                                                                                         'ind_error_sd','ind_error_2sd','ind_error_3sd'))
results_quadratic_penalized <- results_quadratic_penalized[order(results_quadratic_penalized$b0,results_quadratic_penalized$b1,results_quadratic_penalized$b2,
                                                                 results_quadratic_penalized$type,results_quadratic_penalized$corr_magnitude),]

#### create plot
pdf(paste0(dir_results_save,'Fig_4.pdf'), height = 9, width = 7.5)
layout(matrix(c(1,1,17,8,8,
                2,3,17,9,10,
                4,5,17,11,12,
                6,7,17,13,14,
                15,15,15,15,15,
                16,16,16,16,16),
              nrow = 6,ncol = 5,byrow = TRUE),
       heights = c(0.08,0.4,0.4,0.4,0.08,0.08),
       widths = c(0.4,0.4,0.08,0.4,0.4))

x = seq(0,20, by=0.1)

# panel A
par(mar = c(0,0.5,0,0))
plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
text(0,0.4,paste("A."), pos=4,cex = 2, col = "black", font=2)

par(mar = c(5.1, 4.1, 4.1, 2.1),cex.axis=1.2, cex.lab=1.2, cex.main=1.2, cex.sub=1.2)
row_num <- 79
plot(x, x*results_quadratic_quadratic$b1[row_num]+(x^2)*results_quadratic_quadratic$b2[row_num], type = "l",lwd=1.5,cex = 1.2,
     main=bquote(atop("Classical, "~b[0]=="log(8)",b[1]=="0.012, "~b[2]=="-0.0003")),
     xlab=expression(paste("Annual PM"[2.5],", ",mu,"g/m"^3)),ylab="Log(RR)")
lines(x, x*results_quadratic_quadratic$b1_est[row_num]+(x^2)*results_quadratic_quadratic$b2_est[row_num], type = "l", lty=2, col = "#117733")
lines(x, x*results_quadratic_quadratic$b1_est[row_num+1]+(x^2)*results_quadratic_quadratic$b2_est[row_num+1], type = "l", lty=3, col = "#0072B2")
lines(x, x*results_quadratic_quadratic$b1_est[row_num+2]+(x^2)*results_quadratic_quadratic$b2_est[row_num+2], type = "l", lty=6, col = "#D55E00")

row_num <- 73
plot(x, x*results_quadratic_quadratic$b1[row_num]+(x^2)*results_quadratic_quadratic$b2[row_num], type = "l",lwd=1.5,cex = 1.2,
     main=bquote(atop("Berkson, "~b[0]=="log(8)",b[1]=="0.012, "~b[2]=="-0.0003")),
     xlab=expression(paste("Annual PM"[2.5],", ",mu,"g/m"^3)),ylab="Log(RR)")
lines(x, x*results_quadratic_quadratic$b1_est[row_num]+(x^2)*results_quadratic_quadratic$b2_est[row_num], type = "l", lty=2, col = "#117733")
lines(x, x*results_quadratic_quadratic$b1_est[row_num+1]+(x^2)*results_quadratic_quadratic$b2_est[row_num+1], type = "l", lty=3, col = "#0072B2")
lines(x, x*results_quadratic_quadratic$b1_est[row_num+2]+(x^2)*results_quadratic_quadratic$b2_est[row_num+2], type = "l", lty=6, col = "#D55E00")

row_num <- 367
plot(x, x*results_quadratic_quadratic$b1[row_num]+(x^2)*results_quadratic_quadratic$b2[row_num], type = "l",lwd=1.5,cex = 1.2,
     main=bquote(atop("Classical, "~b[0]=="log(12)",b[1]=="0.019, "~b[2]=="-0.0003")),
     xlab=expression(paste("Annual PM"[2.5],", ",mu,"g/m"^3)),ylab="Log(RR)")
lines(x, x*results_quadratic_quadratic$b1_est[row_num]+(x^2)*results_quadratic_quadratic$b2_est[row_num], type = "l", lty=2, col = "#117733")
lines(x, x*results_quadratic_quadratic$b1_est[row_num+1]+(x^2)*results_quadratic_quadratic$b2_est[row_num+1], type = "l", lty=3, col = "#0072B2")
lines(x, x*results_quadratic_quadratic$b1_est[row_num+2]+(x^2)*results_quadratic_quadratic$b2_est[row_num+2], type = "l", lty=6, col = "#D55E00")

row_num <- 361
plot(x, x*results_quadratic_quadratic$b1[row_num]+(x^2)*results_quadratic_quadratic$b2[row_num], type = "l",lwd=1.5,cex = 1.2,
     main=bquote(atop("Berkson, "~b[0]=="log(12)",b[1]=="0.019, "~b[2]=="-0.0003")),
     xlab=expression(paste("Annual PM"[2.5],", ",mu,"g/m"^3)),ylab="Log(RR)")
lines(x, x*results_quadratic_quadratic$b1_est[row_num]+(x^2)*results_quadratic_quadratic$b2_est[row_num], type = "l", lty=2, col = "#117733")
lines(x, x*results_quadratic_quadratic$b1_est[row_num+1]+(x^2)*results_quadratic_quadratic$b2_est[row_num+1], type = "l", lty=3, col = "#0072B2")
lines(x, x*results_quadratic_quadratic$b1_est[row_num+2]+(x^2)*results_quadratic_quadratic$b2_est[row_num+2], type = "l", lty=6, col = "#D55E00")

row_num <- 571
plot(x, x*results_quadratic_quadratic$b1[row_num]+(x^2)*results_quadratic_quadratic$b2[row_num], type = "l",lwd=1.5,cex = 1.2,
     main=bquote(atop("Classical, "~b[0]=="log(20)",b[1]=="0.012, "~b[2]=="0.0003")),
     xlab=expression(paste("Annual PM"[2.5],", ",mu,"g/m"^3)),ylab="Log(RR)")
lines(x, x*results_quadratic_quadratic$b1_est[row_num]+(x^2)*results_quadratic_quadratic$b2_est[row_num], type = "l", lty=2, col = "#117733")
lines(x, x*results_quadratic_quadratic$b1_est[row_num+1]+(x^2)*results_quadratic_quadratic$b2_est[row_num+1], type = "l", lty=3, col = "#0072B2")
lines(x, x*results_quadratic_quadratic$b1_est[row_num+2]+(x^2)*results_quadratic_quadratic$b2_est[row_num+2], type = "l", lty=6, col = "#D55E00")

row_num <- 565
plot(x, x*results_quadratic_quadratic$b1[row_num]+(x^2)*results_quadratic_quadratic$b2[row_num], type = "l",lwd=1.5,cex = 1.2,
     main=bquote(atop("Berkson, "~b[0]=="log(20)",b[1]=="0.012, "~b[2]=="0.0003")),
     xlab=expression(paste("Annual PM"[2.5],", ",mu,"g/m"^3)),ylab="Log(RR)")
lines(x, x*results_quadratic_quadratic$b1_est[row_num]+(x^2)*results_quadratic_quadratic$b2_est[row_num], type = "l", lty=2, col = "#117733")
lines(x, x*results_quadratic_quadratic$b1_est[row_num+1]+(x^2)*results_quadratic_quadratic$b2_est[row_num+1], type = "l", lty=3, col = "#0072B2")
lines(x, x*results_quadratic_quadratic$b1_est[row_num+2]+(x^2)*results_quadratic_quadratic$b2_est[row_num+2], type = "l", lty=6, col = "#D55E00")

# panel B
par(mar = c(0,1,0,0))
plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
text(0,0.4,paste("B."), pos=4,cex = 2, col = "black", font=2)

par(mar = c(5.1, 4.1, 4.1, 2.1),cex.axis=1.2, cex.lab=1.2, cex.main=1.2, cex.sub=1.2)
row_num <- 79
plot(x, x*results_quadratic_penalized$b1[row_num]+(x^2)*results_quadratic_penalized$b2[row_num], type = "l",lwd=1.5,cex = 1.2,
     main=bquote(atop("Classical, "~b[0]=="log(8)",b[1]=="0.012, "~b[2]=="-0.0003")),
     xlab=expression(paste("Annual PM"[2.5],", ",mu,"g/m"^3)), ylab="Log(RR)")
lines(x, as.vector(as.matrix(results_quadratic_penalized[row_num,6:206]))-as.double(results_quadratic_penalized[row_num,6]), type = "l", lty=2, col = "#117733")
lines(x, as.vector(as.matrix(results_quadratic_penalized[row_num+1,6:206]))-as.double(results_quadratic_penalized[row_num+1,6]), type = "l", lty=3, col = "#0072B2")
lines(x, as.vector(as.matrix(results_quadratic_penalized[row_num+2,6:206]))-as.double(results_quadratic_penalized[row_num+2,6]), type = "l", lty=6, col = "#D55E00")

row_num <- 73
plot(x, x*results_quadratic_penalized$b1[row_num]+(x^2)*results_quadratic_penalized$b2[row_num], type = "l",lwd=1.5,cex = 1.2,
     main=bquote(atop("Berkson, "~b[0]=="log(8)",b[1]=="0.012, "~b[2]=="-0.0003")),
     xlab=expression(paste("Annual PM"[2.5],", ",mu,"g/m"^3)), ylab="Log(RR)")
lines(x, as.vector(as.matrix(results_quadratic_penalized[row_num,6:206]))-as.double(results_quadratic_penalized[row_num,6]), type = "l", lty=2, col = "#117733")
lines(x, as.vector(as.matrix(results_quadratic_penalized[row_num+1,6:206]))-as.double(results_quadratic_penalized[row_num+1,6]), type = "l", lty=3, col = "#0072B2")
lines(x, as.vector(as.matrix(results_quadratic_penalized[row_num+2,6:206]))-as.double(results_quadratic_penalized[row_num+2,6]), type = "l", lty=6, col = "#D55E00")

row_num <- 367
plot(x, x*results_quadratic_penalized$b1[row_num]+(x^2)*results_quadratic_penalized$b2[row_num], type = "l",lwd=1.5,cex = 1.2,
     main=bquote(atop("Classical, "~b[0]=="log(12)",b[1]=="0.019, "~b[2]=="-0.0003")),
     xlab=expression(paste("Annual PM"[2.5],", ",mu,"g/m"^3)), ylab="Log(RR)")
lines(x, as.vector(as.matrix(results_quadratic_penalized[row_num,6:206]))-as.double(results_quadratic_penalized[row_num,6]), type = "l", lty=2, col = "#117733")
lines(x, as.vector(as.matrix(results_quadratic_penalized[row_num+1,6:206]))-as.double(results_quadratic_penalized[row_num+1,6]), type = "l", lty=3, col = "#0072B2")
lines(x, as.vector(as.matrix(results_quadratic_penalized[row_num+2,6:206]))-as.double(results_quadratic_penalized[row_num+2,6]), type = "l", lty=6, col = "#D55E00")

row_num <- 361
plot(x, x*results_quadratic_penalized$b1[row_num]+(x^2)*results_quadratic_penalized$b2[row_num], type = "l",lwd=1.5,cex = 1.2,
     main=bquote(atop("Berksonl, "~b[0]=="log(12)",b[1]=="0.019, "~b[2]=="-0.0003")),
     xlab=expression(paste("Annual PM"[2.5],", ",mu,"g/m"^3)), ylab="Log(RR)")
lines(x, as.vector(as.matrix(results_quadratic_penalized[row_num,6:206]))-as.double(results_quadratic_penalized[row_num,6]), type = "l", lty=2, col = "#117733")
lines(x, as.vector(as.matrix(results_quadratic_penalized[row_num+1,6:206]))-as.double(results_quadratic_penalized[row_num+1,6]), type = "l", lty=3, col = "#0072B2")
lines(x, as.vector(as.matrix(results_quadratic_penalized[row_num+2,6:206]))-as.double(results_quadratic_penalized[row_num+2,6]), type = "l", lty=6, col = "#D55E00")

row_num <- 571
plot(x, x*results_quadratic_penalized$b1[row_num]+(x^2)*results_quadratic_penalized$b2[row_num], type = "l",lwd=1.5,cex = 1.2,
     main=bquote(atop("Classical, "~b[0]=="log(20)",b[1]=="0.012, "~b[2]=="0.0003")),
     xlab=expression(paste("Annual PM"[2.5],", ",mu,"g/m"^3)), ylab="Log(RR)")
lines(x, as.vector(as.matrix(results_quadratic_penalized[row_num,6:206]))-as.double(results_quadratic_penalized[row_num,6]), type = "l", lty=2, col = "#117733")
lines(x, as.vector(as.matrix(results_quadratic_penalized[row_num+1,6:206]))-as.double(results_quadratic_penalized[row_num+1,6]), type = "l", lty=3, col = "#0072B2")
lines(x, as.vector(as.matrix(results_quadratic_penalized[row_num+2,6:206]))-as.double(results_quadratic_penalized[row_num+2,6]), type = "l", lty=6, col = "#D55E00")

row_num <- 565
plot(x, x*results_quadratic_penalized$b1[row_num]+(x^2)*results_quadratic_penalized$b2[row_num], type = "l",lwd=1.5,cex = 1.2,
     main=bquote(atop("Berkson, "~b[0]=="log(20)",b[1]=="0.012, "~b[2]=="0.0003")),
     xlab=expression(paste("Annual PM"[2.5],", ",mu,"g/m"^3)), ylab="Log(RR)")
lines(x, as.vector(as.matrix(results_quadratic_penalized[row_num,6:206]))-as.double(results_quadratic_penalized[row_num,6]), type = "l", lty=2, col = "#117733")
lines(x, as.vector(as.matrix(results_quadratic_penalized[row_num+1,6:206]))-as.double(results_quadratic_penalized[row_num+1,6]), type = "l", lty=3, col = "#0072B2")
lines(x, as.vector(as.matrix(results_quadratic_penalized[row_num+2,6:206]))-as.double(results_quadratic_penalized[row_num+2,6]), type = "l", lty=6, col = "#D55E00")

# legend
par(mar = c(0,0,0,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend(x = "center",inset = 0,x.intersp=0,
       legend=c("Error-free (truth)", TeX("Error-prone with $\\sigma_{ik}=\\sigma^{*}_{ik}$")),
       col=c("black", "#117733"),lty=c(1,2),
       lwd=2, cex=1.2, horiz = TRUE)

par(mar = c(0,0,0,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend(x = "center",inset = 0,x.intersp=0,
       legend=c(TeX("Error-prone with $\\sigma_{ik}=2\\sigma^{*}_{ik}$"), TeX("Error-prone with $\\sigma_{ik}=3\\sigma^{*}_{ik}$")),
       col=c("#0072B2", "#D55E00"),lty=c(3,6),
       lwd=2, cex=1.2, horiz = TRUE)

dev.off()



###################### 6. Fig 5. Soft-threshold cr model ########################
### extract results
files_softplus <- Sys.glob(paste0(dir_results_softplus,"results_set*.rds"))
files_softplus_ls <- lapply(files_softplus,readRDS)
results_softplus <- rbindlist(files_softplus_ls)
results_softplus <- results_softplus %>% group_by(k,b,type,corr_magnitude,pollutants) %>% summarise_all("mean")
results_softplus$corr_magnitude <- factor(results_softplus$corr_magnitude,levels=c('sp_error_sd','sp_error_2sd','sp_error_3sd',
                                                                                   'ind_error_sd','ind_error_2sd','ind_error_3sd'))
results_softplus <- results_softplus[order(results_softplus$k,results_softplus$b,results_softplus$type,results_softplus$corr_magnitude),]

### create plot
pdf(paste0(dir_results_save,'Fig_5.pdf'), height = 7.5, width = 4)
layout(matrix(c(1,2,
                3,4,
                5,6,
                7,7,
                8,8),
              nrow = 5,ncol = 2,byrow = TRUE),
       heights = c(0.4,0.4,0.4,0.08,0.08))

x = seq(0,30, by=0.1)

par(mar = c(5.1, 4.1, 4.1, 2.1))
row_num <- 55
plot(x, log(1+exp(results_softplus$k[row_num]*(x-results_softplus$b[row_num])))-log(1+exp(results_softplus$k[row_num]*(-results_softplus$b[row_num]))), 
     type = "l",lwd=2,cex=1.2,main=bquote(atop("Classical",b[0]=="0.2, "~b[1]=="26")),
     xlab=expression(paste("Annual PM"[2.5],", ",mu,"g/m"^3)), ylab="Log(RR)")
lines(x, as.vector(as.matrix(results_softplus[row_num,6:306]))-as.double(results_softplus[row_num,6]), type = "l", lty=2, col = "#117733")
lines(x, as.vector(as.matrix(results_softplus[row_num+1,6:306]))-as.double(results_softplus[row_num+1,6]), type = "l", lty=3, col = "#0072B2")
lines(x, as.vector(as.matrix(results_softplus[row_num+2,6:306]))-as.double(results_softplus[row_num+2,6]), type = "l", lty=6, col = "#D55E00")

row_num <- 49
plot(x, log(1+exp(results_softplus$k[row_num]*(x-results_softplus$b[row_num])))-log(1+exp(results_softplus$k[row_num]*(-results_softplus$b[row_num]))), 
     type = "l",lwd=2,cex=1.2,main=bquote(atop("Berkson",b[0]=="0.2, "~b[1]=="26")),
     xlab=expression(paste("Annual PM"[2.5],", ",mu,"g/m"^3)), ylab="Log(RR)")
lines(x, as.vector(as.matrix(results_softplus[row_num,6:306]))-as.double(results_softplus[row_num,6]), type = "l", lty=2, col = "#117733")
lines(x, as.vector(as.matrix(results_softplus[row_num+1,6:306]))-as.double(results_softplus[row_num+1,6]), type = "l", lty=3, col = "#0072B2")
lines(x, as.vector(as.matrix(results_softplus[row_num+2,6:306]))-as.double(results_softplus[row_num+2,6]), type = "l", lty=6, col = "#D55E00")

row_num <- 91
plot(x, log(1+exp(results_softplus$k[row_num]*(x-results_softplus$b[row_num])))-log(1+exp(results_softplus$k[row_num]*(-results_softplus$b[row_num]))), 
     type = "l",lwd=2,cex=1.2,main=bquote(atop("Classical",b[0]=="0.3, "~b[1]=="24")),
     xlab=expression(paste("Annual PM"[2.5],", ",mu,"g/m"^3)), ylab="Log(RR)")
lines(x, as.vector(as.matrix(results_softplus[row_num,6:306]))-as.double(results_softplus[row_num,6]), type = "l", lty=2, col = "#117733")
lines(x, as.vector(as.matrix(results_softplus[row_num+1,6:306]))-as.double(results_softplus[row_num+1,6]), type = "l", lty=3, col = "#0072B2")
lines(x, as.vector(as.matrix(results_softplus[row_num+2,6:306]))-as.double(results_softplus[row_num+2,6]), type = "l", lty=6, col = "#D55E00")

row_num <- 85
plot(x, log(1+exp(results_softplus$k[row_num]*(x-results_softplus$b[row_num])))-log(1+exp(results_softplus$k[row_num]*(-results_softplus$b[row_num]))), 
     type = "l",lwd=2,cex=1.2,main=bquote(atop("Berkson",b[0]=="0.3, "~b[1]=="24")),
     xlab=expression(paste("Annual PM"[2.5],", ",mu,"g/m"^3)), ylab="Log(RR)")
lines(x, as.vector(as.matrix(results_softplus[row_num,6:306]))-as.double(results_softplus[row_num,6]), type = "l", lty=2, col = "#117733")
lines(x, as.vector(as.matrix(results_softplus[row_num+1,6:306]))-as.double(results_softplus[row_num+1,6]), type = "l", lty=3, col = "#0072B2")
lines(x, as.vector(as.matrix(results_softplus[row_num+2,6:306]))-as.double(results_softplus[row_num+2,6]), type = "l", lty=6, col = "#D55E00")

row_num <- 139
plot(x, log(1+exp(results_softplus$k[row_num]*(x-results_softplus$b[row_num])))-log(1+exp(results_softplus$k[row_num]*(-results_softplus$b[row_num]))), 
     type = "l",lwd=2,cex=1.2,main=bquote(atop("Classical",b[0]=="0.4, "~b[1]=="23")),
     xlab=expression(paste("Annual PM"[2.5],", ",mu,"g/m"^3)), ylab="Log(RR)")
lines(x, as.vector(as.matrix(results_softplus[row_num,6:306]))-as.double(results_softplus[row_num,6]), type = "l", lty=2, col = "#117733")
lines(x, as.vector(as.matrix(results_softplus[row_num+1,6:306]))-as.double(results_softplus[row_num+1,6]), type = "l", lty=3, col = "#0072B2")
lines(x, as.vector(as.matrix(results_softplus[row_num+2,6:306]))-as.double(results_softplus[row_num+2,6]), type = "l", lty=6, col = "#D55E00")

row_num <- 133
plot(x, log(1+exp(results_softplus$k[row_num]*(x-results_softplus$b[row_num])))-log(1+exp(results_softplus$k[row_num]*(-results_softplus$b[row_num]))), 
     type = "l",lwd=2,cex=1.2,main=bquote(atop("Berkson",b[0]=="0.4, "~b[1]=="23")),
     xlab=expression(paste("Annual PM"[2.5],", ",mu,"g/m"^3)), ylab="Log(RR)")
lines(x, as.vector(as.matrix(results_softplus[row_num,6:306]))-as.double(results_softplus[row_num,6]), type = "l", lty=2, col = "#117733")
lines(x, as.vector(as.matrix(results_softplus[row_num+1,6:306]))-as.double(results_softplus[row_num+1,6]), type = "l", lty=3, col = "#0072B2")
lines(x, as.vector(as.matrix(results_softplus[row_num+2,6:306]))-as.double(results_softplus[row_num+2,6]), type = "l", lty=6, col = "#D55E00")

# legend
par(mar = c(0,0,0,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend(x = "center",inset = 0,x.intersp=0,
       legend=c("Error-free (truth)", TeX("Error-prone with $\\sigma_{ik}=\\sigma^{*}_{ik}$")),
       col=c("black", "#117733"),lty=c(1,2),
       lwd=2, cex=1.2, horiz = TRUE)

# legend
par(mar = c(0,0,0,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend(x = "center",inset = 0,x.intersp=0,
       legend=c(TeX("Error-prone with $\\sigma_{ik}=2\\sigma^{*}_{ik}$"), TeX("Error-prone with $\\sigma_{ik}=3\\sigma^{*}_{ik}$")),
       col=c("#0072B2","#D55E00"),lty=c(3,6),
       lwd=2, cex=1.2, horiz = TRUE)

dev.off()



############## 7. Fig S2. Correlation plot ##############
# load data that merges spatially each monitor with all possible monitors around it
Resid <- read_fst(paste0(dir_uncertainty,"redis.fst"))
Resid <- as.data.table(Resid)

n <- c()
cor <- c()
pval <- c()
for (i in unique(Resid$bin)){
  tmp=Resid[bin==i]
  n <- c(n,nrow(tmp))
  cor <- c(cor,cor(tmp$res.x,tmp$res.y))
  pval <- c(pval,  round(cor.test(tmp$res.x,tmp$res.y)$p.value,3))
}
final <- data.frame(bin= unique(Resid$binv),cor=cor,pval=pval,n=n)

# summary(final)
# bin              cor                pval                n            color          
# Min.   :     0   Min.   :0.009206   Min.   :0.000000   Min.   : 5712   Length:61         
# 1st Qu.: 72500   1st Qu.:0.043482   1st Qu.:0.000000   1st Qu.:12054   Class :character  
# Median :147500   Median :0.064018   Median :0.000000   Median :18706   Mode  :character  
# Mean   :147541   Mean   :0.095332   Mean   :0.003934   Mean   :16840                     
# 3rd Qu.:222500   3rd Qu.:0.103493   3rd Qu.:0.000000   3rd Qu.:20954                     
# Max.   :297500   Max.   :0.999203   Max.   :0.172000   Max.   :27076 

range(final[final$bin>=40000,]$cor)
# 0.009206345 0.178043193

ggplot(data = final,aes(x=bin)) +
  geom_point(aes(y = cor), shape = 21, size = 1)+
  # geom_text(aes(y = cor, label = n),size=2, vjust = -1.5) +
  scale_x_continuous("Distance (m)") +
  scale_y_continuous("Pearson correlation coefficient") +
  scale_fill_identity() +
  theme_bw() +
  geom_vline(xintercept = 40000, linetype="dotted", 
             color = "blue", size=1.5)+
  theme(panel.grid = element_blank(),
        text = element_text(size = 8))

ggsave(paste0(dir_results_save,"Fig_S2.png"),width=6,height=4,dpi=600)



###################### 8. Fig S3. Quadratic cr model ########################
### extract quadratic quadratic results
files_quadratic_quadratic <- Sys.glob(paste0(dir_results_quadratic_quadratic,'results_set*.rds'))
files_quadratic_quadratic_ls <- lapply(files_quadratic_quadratic,readRDS)
results_quadratic_quadratic <- rbindlist(files_quadratic_quadratic_ls)
results_quadratic_quadratic$pollutants <- NULL
results_quadratic_quadratic <- results_quadratic_quadratic %>% dplyr::group_by(b0,b1,b2,type,corr_magnitude) %>% dplyr::summarise(b0_est=mean(b0_est,na.rm=TRUE),
                                                                                                                                  b1_est=mean(b1_est,na.rm=TRUE),
                                                                                                                                  b2_est=mean(b2_est,na.rm=TRUE))
results_quadratic_quadratic$corr_magnitude <- factor(results_quadratic_quadratic$corr_magnitude,levels=c('sp_error_sd','sp_error_2sd','sp_error_3sd',
                                                                                                         'ind_error_sd','ind_error_2sd','ind_error_3sd'))
results_quadratic_quadratic <- results_quadratic_quadratic[order(results_quadratic_quadratic$b0,results_quadratic_quadratic$b1,results_quadratic_quadratic$b2,
                                                                 results_quadratic_quadratic$type,results_quadratic_quadratic$corr_magnitude),]

### extract quadratic penalized results
files_quadratic_penalized <- Sys.glob(paste0(dir_results_quadratic_penalized,'results_set*.rds'))
files_quadratic_penalized_ls <- lapply(files_quadratic_penalized,readRDS)
results_quadratic_penalized <- rbindlist(files_quadratic_penalized_ls)
results_quadratic_penalized$pollutants <- NULL
results_quadratic_penalized <- results_quadratic_penalized %>% dplyr::group_by(b0,b1,b2,type,corr_magnitude) %>% dplyr::summarise_all("mean")
results_quadratic_penalized$corr_magnitude <- factor(results_quadratic_penalized$corr_magnitude,levels=c('sp_error_sd','sp_error_2sd','sp_error_3sd',
                                                                                                         'ind_error_sd','ind_error_2sd','ind_error_3sd'))
results_quadratic_penalized <- results_quadratic_penalized[order(results_quadratic_penalized$b0,results_quadratic_penalized$b1,results_quadratic_penalized$b2,
                                                                 results_quadratic_penalized$type,results_quadratic_penalized$corr_magnitude),]

#### create plot
png(filename = paste0(dir_results_save,"Fig_S3.png"), height = 9, width = 7.5,units = "in",res = 600)
layout(matrix(c(1,1,17,8,8,
                2,3,17,9,10,
                4,5,17,11,12,
                6,7,17,13,14,
                15,15,15,15,15,
                16,16,16,16,16),
              nrow = 6,ncol = 5,byrow = TRUE),
       heights = c(0.08,0.4,0.4,0.4,0.08,0.08),
       widths = c(0.4,0.4,0.08,0.4,0.4))

x = seq(0,20, by=0.1)

# panel A
par(mar = c(0,0.5,0,0))
plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
text(0,0.4,paste("A."), pos=4,cex = 2, col = "black", font=2)

par(mar = c(5.1, 4.1, 4.1, 2.1),cex.axis=1.2, cex.lab=1.2, cex.main=1.2, cex.sub=1.2)
row_num <- 82
plot(x, x*results_quadratic_quadratic$b1[row_num]+(x^2)*results_quadratic_quadratic$b2[row_num], type = "l",lwd=1.5,cex = 1.2,
     main=bquote(atop("Classical, "~b[0]=="log(8)",b[1]=="0.012, "~b[2]=="-0.0003")),
     xlab=expression(paste("Annual PM"[2.5],", ",mu,"g/m"^3)),ylab="Log(RR)")
lines(x, x*results_quadratic_quadratic$b1_est[row_num]+(x^2)*results_quadratic_quadratic$b2_est[row_num], type = "l", lty=2, col = "#117733")
lines(x, x*results_quadratic_quadratic$b1_est[row_num+1]+(x^2)*results_quadratic_quadratic$b2_est[row_num+1], type = "l", lty=3, col = "#0072B2")
lines(x, x*results_quadratic_quadratic$b1_est[row_num+2]+(x^2)*results_quadratic_quadratic$b2_est[row_num+2], type = "l", lty=6, col = "#D55E00")

row_num <- 76
plot(x, x*results_quadratic_quadratic$b1[row_num]+(x^2)*results_quadratic_quadratic$b2[row_num], type = "l",lwd=1.5,cex = 1.2,
     main=bquote(atop("Berkson, "~b[0]=="log(8)",b[1]=="0.012, "~b[2]=="-0.0003")),
     xlab=expression(paste("Annual PM"[2.5],", ",mu,"g/m"^3)),ylab="Log(RR)")
lines(x, x*results_quadratic_quadratic$b1_est[row_num]+(x^2)*results_quadratic_quadratic$b2_est[row_num], type = "l", lty=2, col = "#117733")
lines(x, x*results_quadratic_quadratic$b1_est[row_num+1]+(x^2)*results_quadratic_quadratic$b2_est[row_num+1], type = "l", lty=3, col = "#0072B2")
lines(x, x*results_quadratic_quadratic$b1_est[row_num+2]+(x^2)*results_quadratic_quadratic$b2_est[row_num+2], type = "l", lty=6, col = "#D55E00")

row_num <- 370
plot(x, x*results_quadratic_quadratic$b1[row_num]+(x^2)*results_quadratic_quadratic$b2[row_num], type = "l",lwd=1.5,cex = 1.2,
     main=bquote(atop("Classical, "~b[0]=="log(12)",b[1]=="0.019, "~b[2]=="-0.0003")),
     xlab=expression(paste("Annual PM"[2.5],", ",mu,"g/m"^3)),ylab="Log(RR)")
lines(x, x*results_quadratic_quadratic$b1_est[row_num]+(x^2)*results_quadratic_quadratic$b2_est[row_num], type = "l", lty=2, col = "#117733")
lines(x, x*results_quadratic_quadratic$b1_est[row_num+1]+(x^2)*results_quadratic_quadratic$b2_est[row_num+1], type = "l", lty=3, col = "#0072B2")
lines(x, x*results_quadratic_quadratic$b1_est[row_num+2]+(x^2)*results_quadratic_quadratic$b2_est[row_num+2], type = "l", lty=6, col = "#D55E00")

row_num <- 364
plot(x, x*results_quadratic_quadratic$b1[row_num]+(x^2)*results_quadratic_quadratic$b2[row_num], type = "l",lwd=1.5,cex = 1.2,
     main=bquote(atop("Berkson, "~b[0]=="log(12)",b[1]=="0.019, "~b[2]=="-0.0003")),
     xlab=expression(paste("Annual PM"[2.5],", ",mu,"g/m"^3)),ylab="Log(RR)")
lines(x, x*results_quadratic_quadratic$b1_est[row_num]+(x^2)*results_quadratic_quadratic$b2_est[row_num], type = "l", lty=2, col = "#117733")
lines(x, x*results_quadratic_quadratic$b1_est[row_num+1]+(x^2)*results_quadratic_quadratic$b2_est[row_num+1], type = "l", lty=3, col = "#0072B2")
lines(x, x*results_quadratic_quadratic$b1_est[row_num+2]+(x^2)*results_quadratic_quadratic$b2_est[row_num+2], type = "l", lty=6, col = "#D55E00")

row_num <- 574
plot(x, x*results_quadratic_quadratic$b1[row_num]+(x^2)*results_quadratic_quadratic$b2[row_num], type = "l",lwd=1.5,cex = 1.2,
     main=bquote(atop("Classical, "~b[0]=="log(20)",b[1]=="0.012, "~b[2]=="0.0003")),
     xlab=expression(paste("Annual PM"[2.5],", ",mu,"g/m"^3)),ylab="Log(RR)")
lines(x, x*results_quadratic_quadratic$b1_est[row_num]+(x^2)*results_quadratic_quadratic$b2_est[row_num], type = "l", lty=2, col = "#117733")
lines(x, x*results_quadratic_quadratic$b1_est[row_num+1]+(x^2)*results_quadratic_quadratic$b2_est[row_num+1], type = "l", lty=3, col = "#0072B2")
lines(x, x*results_quadratic_quadratic$b1_est[row_num+2]+(x^2)*results_quadratic_quadratic$b2_est[row_num+2], type = "l", lty=6, col = "#D55E00")

row_num <- 568
plot(x, x*results_quadratic_quadratic$b1[row_num]+(x^2)*results_quadratic_quadratic$b2[row_num], type = "l",lwd=1.5,cex = 1.2,
     main=bquote(atop("Berkson, "~b[0]=="log(20)",b[1]=="0.012, "~b[2]=="0.0003")),
     xlab=expression(paste("Annual PM"[2.5],", ",mu,"g/m"^3)),ylab="Log(RR)")
lines(x, x*results_quadratic_quadratic$b1_est[row_num]+(x^2)*results_quadratic_quadratic$b2_est[row_num], type = "l", lty=2, col = "#117733")
lines(x, x*results_quadratic_quadratic$b1_est[row_num+1]+(x^2)*results_quadratic_quadratic$b2_est[row_num+1], type = "l", lty=3, col = "#0072B2")
lines(x, x*results_quadratic_quadratic$b1_est[row_num+2]+(x^2)*results_quadratic_quadratic$b2_est[row_num+2], type = "l", lty=6, col = "#D55E00")

# panel B
par(mar = c(0,1,0,0))
plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
text(0,0.4,paste("B."), pos=4,cex = 2, col = "black", font=2)

par(mar = c(5.1, 4.1, 4.1, 2.1),cex.axis=1.2, cex.lab=1.2, cex.main=1.2, cex.sub=1.2)
row_num <- 82
plot(x, x*results_quadratic_penalized$b1[row_num]+(x^2)*results_quadratic_penalized$b2[row_num], type = "l",lwd=1.5,cex = 1.2,
     main=bquote(atop("Classical, "~b[0]=="log(8)",b[1]=="0.012, "~b[2]=="-0.0003")),
     xlab=expression(paste("Annual PM"[2.5],", ",mu,"g/m"^3)), ylab="Log(RR)")
lines(x, as.vector(as.matrix(results_quadratic_penalized[row_num,6:206]))-as.double(results_quadratic_penalized[row_num,6]), type = "l", lty=2, col = "#117733")
lines(x, as.vector(as.matrix(results_quadratic_penalized[row_num+1,6:206]))-as.double(results_quadratic_penalized[row_num+1,6]), type = "l", lty=3, col = "#0072B2")
lines(x, as.vector(as.matrix(results_quadratic_penalized[row_num+2,6:206]))-as.double(results_quadratic_penalized[row_num+2,6]), type = "l", lty=6, col = "#D55E00")

row_num <- 76
plot(x, x*results_quadratic_penalized$b1[row_num]+(x^2)*results_quadratic_penalized$b2[row_num], type = "l",lwd=1.5,cex = 1.2,
     main=bquote(atop("Berkson, "~b[0]=="log(8)",b[1]=="0.012, "~b[2]=="-0.0003")),
     xlab=expression(paste("Annual PM"[2.5],", ",mu,"g/m"^3)), ylab="Log(RR)")
lines(x, as.vector(as.matrix(results_quadratic_penalized[row_num,6:206]))-as.double(results_quadratic_penalized[row_num,6]), type = "l", lty=2, col = "#117733")
lines(x, as.vector(as.matrix(results_quadratic_penalized[row_num+1,6:206]))-as.double(results_quadratic_penalized[row_num+1,6]), type = "l", lty=3, col = "#0072B2")
lines(x, as.vector(as.matrix(results_quadratic_penalized[row_num+2,6:206]))-as.double(results_quadratic_penalized[row_num+2,6]), type = "l", lty=6, col = "#D55E00")

row_num <- 370
plot(x, x*results_quadratic_penalized$b1[row_num]+(x^2)*results_quadratic_penalized$b2[row_num], type = "l",lwd=1.5,cex = 1.2,
     main=bquote(atop("Classical, "~b[0]=="log(12)",b[1]=="0.019, "~b[2]=="-0.0003")),
     xlab=expression(paste("Annual PM"[2.5],", ",mu,"g/m"^3)), ylab="Log(RR)")
lines(x, as.vector(as.matrix(results_quadratic_penalized[row_num,6:206]))-as.double(results_quadratic_penalized[row_num,6]), type = "l", lty=2, col = "#117733")
lines(x, as.vector(as.matrix(results_quadratic_penalized[row_num+1,6:206]))-as.double(results_quadratic_penalized[row_num+1,6]), type = "l", lty=3, col = "#0072B2")
lines(x, as.vector(as.matrix(results_quadratic_penalized[row_num+2,6:206]))-as.double(results_quadratic_penalized[row_num+2,6]), type = "l", lty=6, col = "#D55E00")

row_num <- 364
plot(x, x*results_quadratic_penalized$b1[row_num]+(x^2)*results_quadratic_penalized$b2[row_num], type = "l",lwd=1.5,cex = 1.2,
     main=bquote(atop("Berksonl, "~b[0]=="log(12)",b[1]=="0.019, "~b[2]=="-0.0003")),
     xlab=expression(paste("Annual PM"[2.5],", ",mu,"g/m"^3)), ylab="Log(RR)")
lines(x, as.vector(as.matrix(results_quadratic_penalized[row_num,6:206]))-as.double(results_quadratic_penalized[row_num,6]), type = "l", lty=2, col = "#117733")
lines(x, as.vector(as.matrix(results_quadratic_penalized[row_num+1,6:206]))-as.double(results_quadratic_penalized[row_num+1,6]), type = "l", lty=3, col = "#0072B2")
lines(x, as.vector(as.matrix(results_quadratic_penalized[row_num+2,6:206]))-as.double(results_quadratic_penalized[row_num+2,6]), type = "l", lty=6, col = "#D55E00")

row_num <- 574
plot(x, x*results_quadratic_penalized$b1[row_num]+(x^2)*results_quadratic_penalized$b2[row_num], type = "l",lwd=1.5,cex = 1.2,
     main=bquote(atop("Classical, "~b[0]=="log(20)",b[1]=="0.012, "~b[2]=="0.0003")),
     xlab=expression(paste("Annual PM"[2.5],", ",mu,"g/m"^3)), ylab="Log(RR)")
lines(x, as.vector(as.matrix(results_quadratic_penalized[row_num,6:206]))-as.double(results_quadratic_penalized[row_num,6]), type = "l", lty=2, col = "#117733")
lines(x, as.vector(as.matrix(results_quadratic_penalized[row_num+1,6:206]))-as.double(results_quadratic_penalized[row_num+1,6]), type = "l", lty=3, col = "#0072B2")
lines(x, as.vector(as.matrix(results_quadratic_penalized[row_num+2,6:206]))-as.double(results_quadratic_penalized[row_num+2,6]), type = "l", lty=6, col = "#D55E00")

row_num <- 568
plot(x, x*results_quadratic_penalized$b1[row_num]+(x^2)*results_quadratic_penalized$b2[row_num], type = "l",lwd=1.5,cex = 1.2,
     main=bquote(atop("Berkson, "~b[0]=="log(20)",b[1]=="0.012, "~b[2]=="0.0003")),
     xlab=expression(paste("Annual PM"[2.5],", ",mu,"g/m"^3)), ylab="Log(RR)")
lines(x, as.vector(as.matrix(results_quadratic_penalized[row_num,6:206]))-as.double(results_quadratic_penalized[row_num,6]), type = "l", lty=2, col = "#117733")
lines(x, as.vector(as.matrix(results_quadratic_penalized[row_num+1,6:206]))-as.double(results_quadratic_penalized[row_num+1,6]), type = "l", lty=3, col = "#0072B2")
lines(x, as.vector(as.matrix(results_quadratic_penalized[row_num+2,6:206]))-as.double(results_quadratic_penalized[row_num+2,6]), type = "l", lty=6, col = "#D55E00")

# legend
par(mar = c(0,0,0,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend(x = "center",inset = 0,x.intersp=0,
       legend=c("Error-free (truth)", TeX("Error-prone with $\\sigma_{ik}=\\sigma^{*}_{ik}$")),
       col=c("black", "#117733"),lty=c(1,2),
       lwd=2, cex=1.2, horiz = TRUE)

par(mar = c(0,0,0,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend(x = "center",inset = 0,x.intersp=0,
       legend=c(TeX("Error-prone with $\\sigma_{ik}=2\\sigma^{*}_{ik}$"), TeX("Error-prone with $\\sigma_{ik}=3\\sigma^{*}_{ik}$")),
       col=c("#0072B2", "#D55E00"),lty=c(3,6),
       lwd=2, cex=1.2, horiz = TRUE)

dev.off()



###################### 9. Fig S4. Soft-threshold cr model ########################
### extract results
files_softplus <- Sys.glob(paste0(dir_results_softplus,"results_set*.rds"))
files_softplus_ls <- lapply(files_softplus,readRDS)
results_softplus <- rbindlist(files_softplus_ls)
results_softplus <- results_softplus %>% group_by(k,b,type,corr_magnitude,pollutants) %>% summarise_all("mean")
results_softplus$corr_magnitude <- factor(results_softplus$corr_magnitude,levels=c('sp_error_sd','sp_error_2sd','sp_error_3sd',
                                                                                   'ind_error_sd','ind_error_2sd','ind_error_3sd'))
results_softplus <- results_softplus[order(results_softplus$k,results_softplus$b,results_softplus$type,results_softplus$corr_magnitude),]

### create plot
png(filename = paste0(dir_results_save,"Fig_S4.png"), height = 7.5, width = 4, units = "in", res = 600)
layout(matrix(c(1,2,
                3,4,
                5,6,
                7,7,
                8,8),
              nrow = 5,ncol = 2,byrow = TRUE),
       heights = c(0.4,0.4,0.4,0.08,0.08))

x = seq(0,30, by=0.1)

par(mar = c(5.1, 4.1, 4.1, 2.1))
row_num <- 58
plot(x, log(1+exp(results_softplus$k[row_num]*(x-results_softplus$b[row_num])))-log(1+exp(results_softplus$k[row_num]*(-results_softplus$b[row_num]))), 
     type = "l",lwd=2,cex=1.2,main=bquote(atop("Classical",b[0]=="0.2, "~b[1]=="26")),
     xlab=expression(paste("Annual PM"[2.5],", ",mu,"g/m"^3)), ylab="Log(RR)")
lines(x, as.vector(as.matrix(results_softplus[row_num,6:306]))-as.double(results_softplus[row_num,6]), type = "l", lty=2, col = "#117733")
lines(x, as.vector(as.matrix(results_softplus[row_num+1,6:306]))-as.double(results_softplus[row_num+1,6]), type = "l", lty=3, col = "#0072B2")
lines(x, as.vector(as.matrix(results_softplus[row_num+2,6:306]))-as.double(results_softplus[row_num+2,6]), type = "l", lty=6, col = "#D55E00")

row_num <- 52
plot(x, log(1+exp(results_softplus$k[row_num]*(x-results_softplus$b[row_num])))-log(1+exp(results_softplus$k[row_num]*(-results_softplus$b[row_num]))), 
     type = "l",lwd=2,cex=1.2,main=bquote(atop("Berkson",b[0]=="0.2, "~b[1]=="26")),
     xlab=expression(paste("Annual PM"[2.5],", ",mu,"g/m"^3)), ylab="Log(RR)")
lines(x, as.vector(as.matrix(results_softplus[row_num,6:306]))-as.double(results_softplus[row_num,6]), type = "l", lty=2, col = "#117733")
lines(x, as.vector(as.matrix(results_softplus[row_num+1,6:306]))-as.double(results_softplus[row_num+1,6]), type = "l", lty=3, col = "#0072B2")
lines(x, as.vector(as.matrix(results_softplus[row_num+2,6:306]))-as.double(results_softplus[row_num+2,6]), type = "l", lty=6, col = "#D55E00")

row_num <- 94
plot(x, log(1+exp(results_softplus$k[row_num]*(x-results_softplus$b[row_num])))-log(1+exp(results_softplus$k[row_num]*(-results_softplus$b[row_num]))), 
     type = "l",lwd=2,cex=1.2,main=bquote(atop("Classical",b[0]=="0.3, "~b[1]=="24")),
     xlab=expression(paste("Annual PM"[2.5],", ",mu,"g/m"^3)), ylab="Log(RR)")
lines(x, as.vector(as.matrix(results_softplus[row_num,6:306]))-as.double(results_softplus[row_num,6]), type = "l", lty=2, col = "#117733")
lines(x, as.vector(as.matrix(results_softplus[row_num+1,6:306]))-as.double(results_softplus[row_num+1,6]), type = "l", lty=3, col = "#0072B2")
lines(x, as.vector(as.matrix(results_softplus[row_num+2,6:306]))-as.double(results_softplus[row_num+2,6]), type = "l", lty=6, col = "#D55E00")

row_num <- 88
plot(x, log(1+exp(results_softplus$k[row_num]*(x-results_softplus$b[row_num])))-log(1+exp(results_softplus$k[row_num]*(-results_softplus$b[row_num]))), 
     type = "l",lwd=2,cex=1.2,main=bquote(atop("Berkson",b[0]=="0.3, "~b[1]=="24")),
     xlab=expression(paste("Annual PM"[2.5],", ",mu,"g/m"^3)), ylab="Log(RR)")
lines(x, as.vector(as.matrix(results_softplus[row_num,6:306]))-as.double(results_softplus[row_num,6]), type = "l", lty=2, col = "#117733")
lines(x, as.vector(as.matrix(results_softplus[row_num+1,6:306]))-as.double(results_softplus[row_num+1,6]), type = "l", lty=3, col = "#0072B2")
lines(x, as.vector(as.matrix(results_softplus[row_num+2,6:306]))-as.double(results_softplus[row_num+2,6]), type = "l", lty=6, col = "#D55E00")

row_num <- 142
plot(x, log(1+exp(results_softplus$k[row_num]*(x-results_softplus$b[row_num])))-log(1+exp(results_softplus$k[row_num]*(-results_softplus$b[row_num]))), 
     type = "l",lwd=2,cex=1.2,main=bquote(atop("Classical",b[0]=="0.4, "~b[1]=="23")),
     xlab=expression(paste("Annual PM"[2.5],", ",mu,"g/m"^3)), ylab="Log(RR)")
lines(x, as.vector(as.matrix(results_softplus[row_num,6:306]))-as.double(results_softplus[row_num,6]), type = "l", lty=2, col = "#117733")
lines(x, as.vector(as.matrix(results_softplus[row_num+1,6:306]))-as.double(results_softplus[row_num+1,6]), type = "l", lty=3, col = "#0072B2")
lines(x, as.vector(as.matrix(results_softplus[row_num+2,6:306]))-as.double(results_softplus[row_num+2,6]), type = "l", lty=6, col = "#D55E00")

row_num <- 136
plot(x, log(1+exp(results_softplus$k[row_num]*(x-results_softplus$b[row_num])))-log(1+exp(results_softplus$k[row_num]*(-results_softplus$b[row_num]))), 
     type = "l",lwd=2,cex=1.2,main=bquote(atop("Berkson",b[0]=="0.4, "~b[1]=="23")),
     xlab=expression(paste("Annual PM"[2.5],", ",mu,"g/m"^3)), ylab="Log(RR)")
lines(x, as.vector(as.matrix(results_softplus[row_num,6:306]))-as.double(results_softplus[row_num,6]), type = "l", lty=2, col = "#117733")
lines(x, as.vector(as.matrix(results_softplus[row_num+1,6:306]))-as.double(results_softplus[row_num+1,6]), type = "l", lty=3, col = "#0072B2")
lines(x, as.vector(as.matrix(results_softplus[row_num+2,6:306]))-as.double(results_softplus[row_num+2,6]), type = "l", lty=6, col = "#D55E00")

# legend
par(mar = c(0,0,0,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend(x = "center",inset = 0,x.intersp=0,
       legend=c("Error-free (truth)", TeX("Error-prone with $\\sigma_{ik}=\\sigma^{*}_{ik}$")),
       col=c("black", "#117733"),lty=c(1,2),
       lwd=2, cex=1.2, horiz = TRUE)

# legend
par(mar = c(0,0,0,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend(x = "center",inset = 0,x.intersp=0,
       legend=c(TeX("Error-prone with $\\sigma_{ik}=2\\sigma^{*}_{ik}$"), TeX("Error-prone with $\\sigma_{ik}=3\\sigma^{*}_{ik}$")),
       col=c("#0072B2","#D55E00"),lty=c(3,6),
       lwd=2, cex=1.2, horiz = TRUE)

dev.off()



############## 10. Table S6 & S7. coefficients and CIs in linear cr model ##############
# main analysis and single-pollutant analysis
results_linear_files_list <- Sys.glob(paste0(dir_results_linear,"results*.rds"))
results_linear_files <- lapply(results_linear_files_list,readRDS)
results_linear_raw <- rbindlist(results_linear_files)
results_linear <- results_linear_raw %>% group_by(type,corr_magnitude,pollutants,b0,b1) %>% 
  summarize(b1_est_avg=mean(b1_est),b1_est_lb=quantile(b1_est, 0.025),b1_est_ub=quantile(b1_est, 0.975)) %>% 
  ungroup()
results_linear$corr_magnitude <- factor(results_linear$corr_magnitude,levels=c('sp_error_sd','sp_error_2sd','sp_error_3sd',
                                                                               'ind_error_sd','ind_error_2sd','ind_error_3sd'))
results_linear <- results_linear[order(results_linear$b0,results_linear$b1,results_linear$type,results_linear$corr_magnitude),]
results_linear_multi <- results_linear[results_linear$pollutants=='multi',]
results_linear_single <- results_linear[results_linear$pollutants=='single',]

# low-level analysis
results_linear_lowlevel_files_list <- Sys.glob(paste0(dir_results_linear_lowlevel,"results*.rds"))
results_linear_lowlevel_files <- lapply(results_linear_lowlevel_files_list,readRDS)
results_linear_lowlevel_raw <- rbindlist(results_linear_lowlevel_files)
results_linear_lowlevel <- results_linear_lowlevel_raw %>% group_by(type,corr_magnitude,b0,b1) %>% 
  summarize(b1_est_avg=mean(b1_est),b1_est_lb=quantile(b1_est, 0.025),b1_est_ub=quantile(b1_est, 0.975)) %>% 
  ungroup()
results_linear_lowlevel$corr_magnitude <- factor(results_linear_lowlevel$corr_magnitude,levels=c('sp_error_sd','sp_error_2sd','sp_error_3sd',
                                                                                                 'ind_error_sd','ind_error_2sd','ind_error_3sd'))
results_linear_lowlevel <- results_linear_lowlevel[order(results_linear_lowlevel$b0,results_linear_lowlevel$b1,results_linear_lowlevel$type,results_linear_lowlevel$corr_magnitude),]




