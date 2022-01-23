###############################################################################
# Project: Exposure measurement error                                         #
# Code: compile results                                                       #
# Machine: QNAP                                                               #
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

dir_data <- '~/data/SimulationData/'
dir_uncertainty <- '~/data/Uncertainty/'
dir_results_linear <- '~/results/Linear/'
dir_results_linear_lowlevel <- '~/results/Linear_LowLevel/'
dir_results_quadratic_quadratic <- '~/results/Quadratic_Quadratic/'
dir_results_quadratic_penalized <- '~/results/Quadratic_Penalized/'
dir_results_softplus <- '~/results/Softplus/'
dir_results_save <- '~/results/'



###################### 1. descriptive statistics ########################
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



###################### 2. Table 2. Biases of linear dr model ########################
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



##################### 3. Fig 2. Relative frequency for PM2.5 ########################
covar<-readRDS(paste0(dir_data,'counts_weight_pred_sd_0127.rds'))
pm25 <- covar$pm25[covar$pm25<=30]
### create plot
pdf(paste0(dir_results_save,'Fig_2_new.pdf'), height = 4, width = 6)
hist(pm25, breaks=25, col=rgb(0,0,1,0.2), xlab=expression(paste("Annual PM"[2.5]," concentration, ",mu,"g/m"^3)), 
     ylim=c(0,0.16),ylab='Relative frequency',main='',freq=FALSE)
dev.off()



###################### 4. Fig 3. Relative frequency for errors ########################
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
pdf(paste0(dir_results_save,'Fig_3_new.pdf'), height = 4, width = 12)
layout(matrix(c(1,2,3,4,4,4),nrow = 2,ncol = 3,byrow = TRUE),heights = c(0.4,0.08))

hist(annual_tempcorr_adj_err, col=rgb(0,0,1,0.2), xlab=expression(paste("Exposure measurement error, ",mu,"g/m"^3)), ylab='Relative frequency',
     main=TeX("Magnitude of error: $\\sigma_{ik}=\\sigma^{*}_{ik}$"),xlim=c(-2,2),ylim=c(0,4.5),breaks=seq(from=-2, to=2, by=0.1), freq=FALSE)
hist(annual_tempcorr_adj_err_sp, col=rgb(1,0,0,0.2),breaks=seq(from=-2, to=2, by=0.1), freq=FALSE, add=TRUE)

hist(annual_tempcorr_adj_err_2sd, col=rgb(0,0,1,0.2), xlab=expression(paste("Exposure measurement error, ",mu,"g/m"^3)), ylab='Relative frequency',
     main=TeX("Magnitude of error: $\\sigma_{ik}=2\\sigma^{*}_{ik}$"),xlim=c(-2,2),ylim=c(0,4.5),breaks=seq(from=-2, to=2, by=0.1), freq=FALSE)
hist(annual_tempcorr_adj_err_2sd_sp, col=rgb(1,0,0,0.2),breaks=seq(from=-2, to=2, by=0.1), freq=FALSE, add=TRUE)

hist(annual_tempcorr_adj_err_3sd, col=rgb(0,0,1,0.2), xlab=expression(paste("Exposure measurement error, ",mu,"g/m"^3)), ylab='Relative frequency',
     main=TeX("Magnitude of error: $\\sigma_{ik}=3\\sigma^{*}_{ik}$"),xlim=c(-2,2),ylim=c(0,4.5),breaks=seq(from=-2, to=2, by=0.1), freq=FALSE)
hist(annual_tempcorr_adj_err_3sd_sp, col=rgb(1,0,0,0.2),breaks=seq(from=-2, to=2, by=0.1), freq=FALSE, add=TRUE)

par(mar = c(0,0,0,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend(x = "center",inset = 0,
       legend=c(TeX("Independent error ($\\epsilon_{ik}$)"), TeX("Spatially correlated error ($\\epsilon_{s,ik}$)")),fill=c(rgb(0,0,1,0.15), rgb(1,0,0,0.15)),
       cex=1.2, horiz = TRUE)
dev.off()



###################### 5. Fig 4. Quadratic dr model ########################
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
pdf(paste0(dir_results_save,'Fig_4_new.pdf'), height = 12, width = 12)
layout(matrix(c(1,1,16,8,8,
                2,3,16,9,10,
                4,5,16,11,12,
                6,7,16,13,14,
                15,15,15,15,15),
              nrow = 5,ncol = 5,byrow = TRUE),
       heights = c(0.08,0.4,0.4,0.4,0.12),
       widths = c(0.4,0.4,0.08,0.4,0.4))

x = seq(0,20, by=0.1)

# panel A
par(mar = c(0,0.5,0,0))
plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
text(0,0.4,paste("A."), pos=4,cex = 2, col = "black", font=2)

par(mar = c(5.1, 4.1, 4.1, 2.1))
row_num <- 79
plot(x, x*results_quadratic_quadratic$b1[row_num]+(x^2)*results_quadratic_quadratic$b2[row_num], type = "l",lwd=2,
     main=bquote(atop("Classical error",b[0]=="log(8), "~b[1]=="0.012, "~b[2]=="-0.0003")),
     xlab=expression(paste("Annual PM"[2.5]," concentration, ",mu,"g/m"^3)),ylab="Log(RR)")
lines(x, x*results_quadratic_quadratic$b1_est[row_num]+(x^2)*results_quadratic_quadratic$b2_est[row_num], type = "l", lty=5, col = rainbow(3)[1])
lines(x, x*results_quadratic_quadratic$b1_est[row_num+1]+(x^2)*results_quadratic_quadratic$b2_est[row_num+1], type = "l", lty=5, col = rainbow(3)[2])
lines(x, x*results_quadratic_quadratic$b1_est[row_num+2]+(x^2)*results_quadratic_quadratic$b2_est[row_num+2], type = "l", lty=5, col = rainbow(3)[3])

row_num <- 73
plot(x, x*results_quadratic_quadratic$b1[row_num]+(x^2)*results_quadratic_quadratic$b2[row_num], type = "l",lwd=2,
     main=bquote(atop("Berkson error",b[0]=="log(8), "~b[1]=="0.012, "~b[2]=="-0.0003")),
     xlab=expression(paste("Annual PM"[2.5]," concentration, ",mu,"g/m"^3)),ylab="Log(RR)")
lines(x, x*results_quadratic_quadratic$b1_est[row_num]+(x^2)*results_quadratic_quadratic$b2_est[row_num], type = "l", lty=5, col = rainbow(3)[1])
lines(x, x*results_quadratic_quadratic$b1_est[row_num+1]+(x^2)*results_quadratic_quadratic$b2_est[row_num+1], type = "l", lty=5, col = rainbow(3)[2])
lines(x, x*results_quadratic_quadratic$b1_est[row_num+2]+(x^2)*results_quadratic_quadratic$b2_est[row_num+2], type = "l", lty=5, col = rainbow(3)[3])

row_num <- 367
plot(x, x*results_quadratic_quadratic$b1[row_num]+(x^2)*results_quadratic_quadratic$b2[row_num], type = "l",lwd=2,
     main=bquote(atop("Classical error",b[0]=="log(12), "~b[1]=="0.019, "~b[2]=="-0.0003")),
     xlab=expression(paste("Annual PM"[2.5]," concentration, ",mu,"g/m"^3)),ylab="Log(RR)")
lines(x, x*results_quadratic_quadratic$b1_est[row_num]+(x^2)*results_quadratic_quadratic$b2_est[row_num], type = "l", lty=5, col = rainbow(3)[1])
lines(x, x*results_quadratic_quadratic$b1_est[row_num+1]+(x^2)*results_quadratic_quadratic$b2_est[row_num+1], type = "l", lty=5, col = rainbow(3)[2])
lines(x, x*results_quadratic_quadratic$b1_est[row_num+2]+(x^2)*results_quadratic_quadratic$b2_est[row_num+2], type = "l", lty=5, col = rainbow(3)[3])

row_num <- 361
plot(x, x*results_quadratic_quadratic$b1[row_num]+(x^2)*results_quadratic_quadratic$b2[row_num], type = "l",lwd=2,
     main=bquote(atop("Berkson error",b[0]=="log(12), "~b[1]=="0.019, "~b[2]=="-0.0003")),
     xlab=expression(paste("Annual PM"[2.5]," concentration, ",mu,"g/m"^3)),ylab="Log(RR)")
lines(x, x*results_quadratic_quadratic$b1_est[row_num]+(x^2)*results_quadratic_quadratic$b2_est[row_num], type = "l", lty=5, col = rainbow(3)[1])
lines(x, x*results_quadratic_quadratic$b1_est[row_num+1]+(x^2)*results_quadratic_quadratic$b2_est[row_num+1], type = "l", lty=5, col = rainbow(3)[2])
lines(x, x*results_quadratic_quadratic$b1_est[row_num+2]+(x^2)*results_quadratic_quadratic$b2_est[row_num+2], type = "l", lty=5, col = rainbow(3)[3])

row_num <- 571
plot(x, x*results_quadratic_quadratic$b1[row_num]+(x^2)*results_quadratic_quadratic$b2[row_num], type = "l",lwd=2,
     main=bquote(atop("Classical error",b[0]=="log(20), "~b[1]=="0.012, "~b[2]=="0.0003")),
     xlab=expression(paste("Annual PM"[2.5]," concentration, ",mu,"g/m"^3)),ylab="Log(RR)")
lines(x, x*results_quadratic_quadratic$b1_est[row_num]+(x^2)*results_quadratic_quadratic$b2_est[row_num], type = "l", lty=5, col = rainbow(3)[1])
lines(x, x*results_quadratic_quadratic$b1_est[row_num+1]+(x^2)*results_quadratic_quadratic$b2_est[row_num+1], type = "l", lty=5, col = rainbow(3)[2])
lines(x, x*results_quadratic_quadratic$b1_est[row_num+2]+(x^2)*results_quadratic_quadratic$b2_est[row_num+2], type = "l", lty=5, col = rainbow(3)[3])

row_num <- 565
plot(x, x*results_quadratic_quadratic$b1[row_num]+(x^2)*results_quadratic_quadratic$b2[row_num], type = "l",lwd=2,
     main=bquote(atop("Berkson error",b[0]=="log(20), "~b[1]=="0.012, "~b[2]=="0.0003")),
     xlab=expression(paste("Annual PM"[2.5]," concentration, ",mu,"g/m"^3)),ylab="Log(RR)")
lines(x, x*results_quadratic_quadratic$b1_est[row_num]+(x^2)*results_quadratic_quadratic$b2_est[row_num], type = "l", lty=5, col = rainbow(3)[1])
lines(x, x*results_quadratic_quadratic$b1_est[row_num+1]+(x^2)*results_quadratic_quadratic$b2_est[row_num+1], type = "l", lty=5, col = rainbow(3)[2])
lines(x, x*results_quadratic_quadratic$b1_est[row_num+2]+(x^2)*results_quadratic_quadratic$b2_est[row_num+2], type = "l", lty=5, col = rainbow(3)[3])

# panel B
par(mar = c(0,1,0,0))
plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
text(0,0.4,paste("B."), pos=4,cex = 2, col = "black", font=2)

par(mar = c(5.1, 4.1, 4.1, 2.1))
row_num <- 79
plot(x, x*results_quadratic_penalized$b1[row_num]+(x^2)*results_quadratic_penalized$b2[row_num], type = "l",lwd=2,
     main=bquote(atop("Classical error",b[0]=="log(8), "~b[1]=="0.012, "~b[2]=="-0.0003")),
     xlab=expression(paste("Annual PM"[2.5]," concentration, ",mu,"g/m"^3)), ylab="Log(RR)")
lines(x, as.vector(as.matrix(results_quadratic_penalized[row_num,6:206]))-as.double(results_quadratic_penalized[row_num,6]), type = "l", lty=5, col = rainbow(3)[1])
lines(x, as.vector(as.matrix(results_quadratic_penalized[row_num+1,6:206]))-as.double(results_quadratic_penalized[row_num+1,6]), type = "l", lty=5, col = rainbow(3)[2])
lines(x, as.vector(as.matrix(results_quadratic_penalized[row_num+2,6:206]))-as.double(results_quadratic_penalized[row_num+2,6]), type = "l", lty=5, col = rainbow(3)[3])

row_num <- 73
plot(x, x*results_quadratic_penalized$b1[row_num]+(x^2)*results_quadratic_penalized$b2[row_num], type = "l",lwd=2,
     main=bquote(atop("Berkson error",b[0]=="log(8), "~b[1]=="0.012, "~b[2]=="-0.0003")),
     xlab=expression(paste("Annual PM"[2.5]," concentration, ",mu,"g/m"^3)), ylab="Log(RR)")
lines(x, as.vector(as.matrix(results_quadratic_penalized[row_num,6:206]))-as.double(results_quadratic_penalized[row_num,6]), type = "l", lty=5, col = rainbow(3)[1])
lines(x, as.vector(as.matrix(results_quadratic_penalized[row_num+1,6:206]))-as.double(results_quadratic_penalized[row_num+1,6]), type = "l", lty=5, col = rainbow(3)[2])
lines(x, as.vector(as.matrix(results_quadratic_penalized[row_num+2,6:206]))-as.double(results_quadratic_penalized[row_num+2,6]), type = "l", lty=5, col = rainbow(3)[3])

row_num <- 367
plot(x, x*results_quadratic_penalized$b1[row_num]+(x^2)*results_quadratic_penalized$b2[row_num], type = "l",lwd=2,
     main=bquote(atop("Classical error",b[0]=="log(12), "~b[1]=="0.019, "~b[2]=="-0.0003")),
     xlab=expression(paste("Annual PM"[2.5]," concentration, ",mu,"g/m"^3)), ylab="Log(RR)")
lines(x, as.vector(as.matrix(results_quadratic_penalized[row_num,6:206]))-as.double(results_quadratic_penalized[row_num,6]), type = "l", lty=5, col = rainbow(3)[1])
lines(x, as.vector(as.matrix(results_quadratic_penalized[row_num+1,6:206]))-as.double(results_quadratic_penalized[row_num+1,6]), type = "l", lty=5, col = rainbow(3)[2])
lines(x, as.vector(as.matrix(results_quadratic_penalized[row_num+2,6:206]))-as.double(results_quadratic_penalized[row_num+2,6]), type = "l", lty=5, col = rainbow(3)[3])

row_num <- 361
plot(x, x*results_quadratic_penalized$b1[row_num]+(x^2)*results_quadratic_penalized$b2[row_num], type = "l",lwd=2,
     main=bquote(atop("Berkson error",b[0]=="log(12), "~b[1]=="0.019, "~b[2]=="-0.0003")),
     xlab=expression(paste("Annual PM"[2.5]," concentration, ",mu,"g/m"^3)), ylab="Log(RR)")
lines(x, as.vector(as.matrix(results_quadratic_penalized[row_num,6:206]))-as.double(results_quadratic_penalized[row_num,6]), type = "l", lty=5, col = rainbow(3)[1])
lines(x, as.vector(as.matrix(results_quadratic_penalized[row_num+1,6:206]))-as.double(results_quadratic_penalized[row_num+1,6]), type = "l", lty=5, col = rainbow(3)[2])
lines(x, as.vector(as.matrix(results_quadratic_penalized[row_num+2,6:206]))-as.double(results_quadratic_penalized[row_num+2,6]), type = "l", lty=5, col = rainbow(3)[3])

row_num <- 571
plot(x, x*results_quadratic_penalized$b1[row_num]+(x^2)*results_quadratic_penalized$b2[row_num], type = "l",lwd=2,
     main=bquote(atop("Classical error",b[0]=="log(20), "~b[1]=="0.012, "~b[2]=="0.0003")),
     xlab=expression(paste("Annual PM"[2.5]," concentration, ",mu,"g/m"^3)), ylab="Log(RR)")
lines(x, as.vector(as.matrix(results_quadratic_penalized[row_num,6:206]))-as.double(results_quadratic_penalized[row_num,6]), type = "l", lty=5, col = rainbow(3)[1])
lines(x, as.vector(as.matrix(results_quadratic_penalized[row_num+1,6:206]))-as.double(results_quadratic_penalized[row_num+1,6]), type = "l", lty=5, col = rainbow(3)[2])
lines(x, as.vector(as.matrix(results_quadratic_penalized[row_num+2,6:206]))-as.double(results_quadratic_penalized[row_num+2,6]), type = "l", lty=5, col = rainbow(3)[3])

row_num <- 565
plot(x, x*results_quadratic_penalized$b1[row_num]+(x^2)*results_quadratic_penalized$b2[row_num], type = "l",lwd=2,
     main=bquote(atop("Berkson error",b[0]=="log(20), "~b[1]=="0.012, "~b[2]=="0.0003")),
     xlab=expression(paste("Annual PM"[2.5]," concentration, ",mu,"g/m"^3)), ylab="Log(RR)")
lines(x, as.vector(as.matrix(results_quadratic_penalized[row_num,6:206]))-as.double(results_quadratic_penalized[row_num,6]), type = "l", lty=5, col = rainbow(3)[1])
lines(x, as.vector(as.matrix(results_quadratic_penalized[row_num+1,6:206]))-as.double(results_quadratic_penalized[row_num+1,6]), type = "l", lty=5, col = rainbow(3)[2])
lines(x, as.vector(as.matrix(results_quadratic_penalized[row_num+2,6:206]))-as.double(results_quadratic_penalized[row_num+2,6]), type = "l", lty=5, col = rainbow(3)[3])

# legend
par(mar = c(0,0,0,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend(x = "center",inset = 0,
       legend=c("Error-free (truth)", TeX("Error-prone with $\\sigma_{ik}=\\sigma^{*}_{ik}$"), TeX("Error-prone with $\\sigma_{ik}=2\\sigma^{*}_{ik}$"), TeX("Error-prone with $\\sigma_{ik}=3\\sigma^{*}_{ik}$")),
       col=c("black", rainbow(3)[1], rainbow(3)[2],rainbow(3)[3]),lty=c(1,5,5,5),
       lwd=2, cex=1.5, horiz = TRUE)

plot(0,type='n',axes=FALSE,ann=FALSE)

dev.off()



###################### 6. Fig 5. Soft-threshold dr model ########################
### extract results
files_softplus <- Sys.glob(paste0(dir_results_softplus,"results_set*.rds"))
files_softplus_ls <- lapply(files_softplus,readRDS)
results_softplus <- rbindlist(files_softplus_ls)
results_softplus <- results_softplus %>% group_by(k,b,type,corr_magnitude,pollutants) %>% summarise_all("mean")
results_softplus$corr_magnitude <- factor(results_softplus$corr_magnitude,levels=c('sp_error_sd','sp_error_2sd','sp_error_3sd',
                                                                                   'ind_error_sd','ind_error_2sd','ind_error_3sd'))
results_softplus <- results_softplus[order(results_softplus$k,results_softplus$b,results_softplus$type,results_softplus$corr_magnitude),]

### create plot
pdf(paste0(dir_results_save,'Fig_5_new.pdf'), height = 12, width = 6)
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
     type = "l",lwd=2,main=bquote(atop("Classical error",a=="0.2, "~b=="26")),
     xlab=expression(paste("Annual PM"[2.5]," concentration, ",mu,"g/m"^3)), ylab="Log(RR)")
lines(x, as.vector(as.matrix(results_softplus[row_num,6:306]))-as.double(results_softplus[row_num,6]), type = "l", lty=5, col = rainbow(3)[1])
lines(x, as.vector(as.matrix(results_softplus[row_num+1,6:306]))-as.double(results_softplus[row_num+1,6]), type = "l", lty=5, col = rainbow(3)[2])
lines(x, as.vector(as.matrix(results_softplus[row_num+2,6:306]))-as.double(results_softplus[row_num+2,6]), type = "l", lty=5, col = rainbow(3)[3])

row_num <- 49
plot(x, log(1+exp(results_softplus$k[row_num]*(x-results_softplus$b[row_num])))-log(1+exp(results_softplus$k[row_num]*(-results_softplus$b[row_num]))), 
     type = "l",lwd=2,main=bquote(atop("Berkson error",a=="0.2, "~b=="26")),
     xlab=expression(paste("Annual PM"[2.5]," concentration, ",mu,"g/m"^3)), ylab="Log(RR)")
lines(x, as.vector(as.matrix(results_softplus[row_num,6:306]))-as.double(results_softplus[row_num,6]), type = "l", lty=5, col = rainbow(3)[1])
lines(x, as.vector(as.matrix(results_softplus[row_num+1,6:306]))-as.double(results_softplus[row_num+1,6]), type = "l", lty=5, col = rainbow(3)[2])
lines(x, as.vector(as.matrix(results_softplus[row_num+2,6:306]))-as.double(results_softplus[row_num+2,6]), type = "l", lty=5, col = rainbow(3)[3])

row_num <- 91
plot(x, log(1+exp(results_softplus$k[row_num]*(x-results_softplus$b[row_num])))-log(1+exp(results_softplus$k[row_num]*(-results_softplus$b[row_num]))), 
     type = "l",lwd=2,main=bquote(atop("Classical error",a=="0.3, "~b=="24")),
     xlab=expression(paste("Annual PM"[2.5]," concentration, ",mu,"g/m"^3)), ylab="Log(RR)")
lines(x, as.vector(as.matrix(results_softplus[row_num,6:306]))-as.double(results_softplus[row_num,6]), type = "l", lty=5, col = rainbow(3)[1])
lines(x, as.vector(as.matrix(results_softplus[row_num+1,6:306]))-as.double(results_softplus[row_num+1,6]), type = "l", lty=5, col = rainbow(3)[2])
lines(x, as.vector(as.matrix(results_softplus[row_num+2,6:306]))-as.double(results_softplus[row_num+2,6]), type = "l", lty=5, col = rainbow(3)[3])

row_num <- 85
plot(x, log(1+exp(results_softplus$k[row_num]*(x-results_softplus$b[row_num])))-log(1+exp(results_softplus$k[row_num]*(-results_softplus$b[row_num]))), 
     type = "l",lwd=2,main=bquote(atop("Berkson error",a=="0.3, "~b=="24")),
     xlab=expression(paste("Annual PM"[2.5]," concentration, ",mu,"g/m"^3)), ylab="Log(RR)")
lines(x, as.vector(as.matrix(results_softplus[row_num,6:306]))-as.double(results_softplus[row_num,6]), type = "l", lty=5, col = rainbow(3)[1])
lines(x, as.vector(as.matrix(results_softplus[row_num+1,6:306]))-as.double(results_softplus[row_num+1,6]), type = "l", lty=5, col = rainbow(3)[2])
lines(x, as.vector(as.matrix(results_softplus[row_num+2,6:306]))-as.double(results_softplus[row_num+2,6]), type = "l", lty=5, col = rainbow(3)[3])

row_num <- 139
plot(x, log(1+exp(results_softplus$k[row_num]*(x-results_softplus$b[row_num])))-log(1+exp(results_softplus$k[row_num]*(-results_softplus$b[row_num]))), 
     type = "l",lwd=2,main=bquote(atop("Classical error",a=="0.4, "~b=="23")),
     xlab=expression(paste("Annual PM"[2.5]," concentration, ",mu,"g/m"^3)), ylab="Log(RR)")
lines(x, as.vector(as.matrix(results_softplus[row_num,6:306]))-as.double(results_softplus[row_num,6]), type = "l", lty=5, col = rainbow(3)[1])
lines(x, as.vector(as.matrix(results_softplus[row_num+1,6:306]))-as.double(results_softplus[row_num+1,6]), type = "l", lty=5, col = rainbow(3)[2])
lines(x, as.vector(as.matrix(results_softplus[row_num+2,6:306]))-as.double(results_softplus[row_num+2,6]), type = "l", lty=5, col = rainbow(3)[3])

row_num <- 133
plot(x, log(1+exp(results_softplus$k[row_num]*(x-results_softplus$b[row_num])))-log(1+exp(results_softplus$k[row_num]*(-results_softplus$b[row_num]))), 
     type = "l",lwd=2,main=bquote(atop("Berkson error",a=="0.4, "~b=="23")),
     xlab=expression(paste("Annual PM"[2.5]," concentration, ",mu,"g/m"^3)), ylab="Log(RR)")
lines(x, as.vector(as.matrix(results_softplus[row_num,6:306]))-as.double(results_softplus[row_num,6]), type = "l", lty=5, col = rainbow(3)[1])
lines(x, as.vector(as.matrix(results_softplus[row_num+1,6:306]))-as.double(results_softplus[row_num+1,6]), type = "l", lty=5, col = rainbow(3)[2])
lines(x, as.vector(as.matrix(results_softplus[row_num+2,6:306]))-as.double(results_softplus[row_num+2,6]), type = "l", lty=5, col = rainbow(3)[3])

# legend
par(mar = c(0,0,0,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend(x = "center",inset = 0,
       legend=c("Error-free (truth)", TeX("Error-prone with $\\sigma_{ik}=\\sigma^{*}_{ik}$")),
       col=c("black", rainbow(3)[1]),lty=c(1,5),
       lwd=2, cex=1.5, horiz = TRUE)

# legend
par(mar = c(0,0,0,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend(x = "center",inset = 0,
       legend=c(TeX("Error-prone with $\\sigma_{ik}=2\\sigma^{*}_{ik}$"), TeX("Error-prone with $\\sigma_{ik}=3\\sigma^{*}_{ik}$")),
       col=c(rainbow(3)[2],rainbow(3)[3]),lty=c(5,5),
       lwd=2, cex=1.5, horiz = TRUE)

dev.off()



###################### 7. Fig S2. Quadratic dr model ########################
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
png(filename = paste0(dir_results_save,"Fig_S2_new.png"), width = 12, height = 12,units = "in",res = 300)
layout(matrix(c(1,1,16,8,8,
                2,3,16,9,10,
                4,5,16,11,12,
                6,7,16,13,14,
                15,15,15,15,15),
              nrow = 5,ncol = 5,byrow = TRUE),
       heights = c(0.08,0.4,0.4,0.4,0.12),
       widths = c(0.4,0.4,0.08,0.4,0.4))

x = seq(0,20, by=0.1)

# panel A
par(mar = c(0,0.5,0,0))
plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
text(0,0.4,paste("A."), pos=4,cex = 2, col = "black", font=2)

par(mar = c(5.1, 4.1, 4.1, 2.1))
row_num <- 82
plot(x, x*results_quadratic_quadratic$b1[row_num]+(x^2)*results_quadratic_quadratic$b2[row_num], type = "l",lwd=2,
     main=bquote(atop("Classical error",b[0]=="log(8), "~b[1]=="0.012, "~b[2]=="-0.0003")),
     xlab=expression(paste("Annual PM"[2.5]," concentration, ",mu,"g/m"^3)),ylab="Log(RR)")
lines(x, x*results_quadratic_quadratic$b1_est[row_num]+(x^2)*results_quadratic_quadratic$b2_est[row_num], type = "l", lty=5, col = rainbow(3)[1])
lines(x, x*results_quadratic_quadratic$b1_est[row_num+1]+(x^2)*results_quadratic_quadratic$b2_est[row_num+1], type = "l", lty=5, col = rainbow(3)[2])
lines(x, x*results_quadratic_quadratic$b1_est[row_num+2]+(x^2)*results_quadratic_quadratic$b2_est[row_num+2], type = "l", lty=5, col = rainbow(3)[3])

row_num <- 76
plot(x, x*results_quadratic_quadratic$b1[row_num]+(x^2)*results_quadratic_quadratic$b2[row_num], type = "l",lwd=2,
     main=bquote(atop("Berkson error",b[0]=="log(8), "~b[1]=="0.012, "~b[2]=="-0.0003")),
     xlab=expression(paste("Annual PM"[2.5]," concentration, ",mu,"g/m"^3)),ylab="Log(RR)")
lines(x, x*results_quadratic_quadratic$b1_est[row_num]+(x^2)*results_quadratic_quadratic$b2_est[row_num], type = "l", lty=5, col = rainbow(3)[1])
lines(x, x*results_quadratic_quadratic$b1_est[row_num+1]+(x^2)*results_quadratic_quadratic$b2_est[row_num+1], type = "l", lty=5, col = rainbow(3)[2])
lines(x, x*results_quadratic_quadratic$b1_est[row_num+2]+(x^2)*results_quadratic_quadratic$b2_est[row_num+2], type = "l", lty=5, col = rainbow(3)[3])

row_num <- 370
plot(x, x*results_quadratic_quadratic$b1[row_num]+(x^2)*results_quadratic_quadratic$b2[row_num], type = "l",lwd=2,
     main=bquote(atop("Classical error",b[0]=="log(12), "~b[1]=="0.019, "~b[2]=="-0.0003")),
     xlab=expression(paste("Annual PM"[2.5]," concentration, ",mu,"g/m"^3)),ylab="Log(RR)")
lines(x, x*results_quadratic_quadratic$b1_est[row_num]+(x^2)*results_quadratic_quadratic$b2_est[row_num], type = "l", lty=5, col = rainbow(3)[1])
lines(x, x*results_quadratic_quadratic$b1_est[row_num+1]+(x^2)*results_quadratic_quadratic$b2_est[row_num+1], type = "l", lty=5, col = rainbow(3)[2])
lines(x, x*results_quadratic_quadratic$b1_est[row_num+2]+(x^2)*results_quadratic_quadratic$b2_est[row_num+2], type = "l", lty=5, col = rainbow(3)[3])

row_num <- 364
plot(x, x*results_quadratic_quadratic$b1[row_num]+(x^2)*results_quadratic_quadratic$b2[row_num], type = "l",lwd=2,
     main=bquote(atop("Berkson error",b[0]=="log(12), "~b[1]=="0.019, "~b[2]=="-0.0003")),
     xlab=expression(paste("Annual PM"[2.5]," concentration, ",mu,"g/m"^3)),ylab="Log(RR)")
lines(x, x*results_quadratic_quadratic$b1_est[row_num]+(x^2)*results_quadratic_quadratic$b2_est[row_num], type = "l", lty=5, col = rainbow(3)[1])
lines(x, x*results_quadratic_quadratic$b1_est[row_num+1]+(x^2)*results_quadratic_quadratic$b2_est[row_num+1], type = "l", lty=5, col = rainbow(3)[2])
lines(x, x*results_quadratic_quadratic$b1_est[row_num+2]+(x^2)*results_quadratic_quadratic$b2_est[row_num+2], type = "l", lty=5, col = rainbow(3)[3])

row_num <- 574
plot(x, x*results_quadratic_quadratic$b1[row_num]+(x^2)*results_quadratic_quadratic$b2[row_num], type = "l",lwd=2,
     main=bquote(atop("Classical error",b[0]=="log(20), "~b[1]=="0.012, "~b[2]=="0.0003")),
     xlab=expression(paste("Annual PM"[2.5]," concentration, ",mu,"g/m"^3)),ylab="Log(RR)")
lines(x, x*results_quadratic_quadratic$b1_est[row_num]+(x^2)*results_quadratic_quadratic$b2_est[row_num], type = "l", lty=5, col = rainbow(3)[1])
lines(x, x*results_quadratic_quadratic$b1_est[row_num+1]+(x^2)*results_quadratic_quadratic$b2_est[row_num+1], type = "l", lty=5, col = rainbow(3)[2])
lines(x, x*results_quadratic_quadratic$b1_est[row_num+2]+(x^2)*results_quadratic_quadratic$b2_est[row_num+2], type = "l", lty=5, col = rainbow(3)[3])

row_num <- 568
plot(x, x*results_quadratic_quadratic$b1[row_num]+(x^2)*results_quadratic_quadratic$b2[row_num], type = "l",lwd=2,
     main=bquote(atop("Berkson error",b[0]=="log(20), "~b[1]=="0.012, "~b[2]=="0.0003")),
     xlab=expression(paste("Annual PM"[2.5]," concentration, ",mu,"g/m"^3)),ylab="Log(RR)")
lines(x, x*results_quadratic_quadratic$b1_est[row_num]+(x^2)*results_quadratic_quadratic$b2_est[row_num], type = "l", lty=5, col = rainbow(3)[1])
lines(x, x*results_quadratic_quadratic$b1_est[row_num+1]+(x^2)*results_quadratic_quadratic$b2_est[row_num+1], type = "l", lty=5, col = rainbow(3)[2])
lines(x, x*results_quadratic_quadratic$b1_est[row_num+2]+(x^2)*results_quadratic_quadratic$b2_est[row_num+2], type = "l", lty=5, col = rainbow(3)[3])

# panel B
par(mar = c(0,1,0,0))
plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
text(0,0.4,paste("B."), pos=4,cex = 2, col = "black", font=2)

par(mar = c(5.1, 4.1, 4.1, 2.1))
row_num <- 82
plot(x, x*results_quadratic_penalized$b1[row_num]+(x^2)*results_quadratic_penalized$b2[row_num], type = "l",lwd=2,
     main=bquote(atop("Classical error",b[0]=="log(8), "~b[1]=="0.012, "~b[2]=="-0.0003")),
     xlab=expression(paste("Annual PM"[2.5]," concentration, ",mu,"g/m"^3)), ylab="Log(RR)")
lines(x, as.vector(as.matrix(results_quadratic_penalized[row_num,6:206]))-as.double(results_quadratic_penalized[row_num,6]), type = "l", lty=5, col = rainbow(3)[1])
lines(x, as.vector(as.matrix(results_quadratic_penalized[row_num+1,6:206]))-as.double(results_quadratic_penalized[row_num+1,6]), type = "l", lty=5, col = rainbow(3)[2])
lines(x, as.vector(as.matrix(results_quadratic_penalized[row_num+2,6:206]))-as.double(results_quadratic_penalized[row_num+2,6]), type = "l", lty=5, col = rainbow(3)[3])

row_num <- 76
plot(x, x*results_quadratic_penalized$b1[row_num]+(x^2)*results_quadratic_penalized$b2[row_num], type = "l",lwd=2,
     main=bquote(atop("Berkson error",b[0]=="log(8), "~b[1]=="0.012, "~b[2]=="-0.0003")),
     xlab=expression(paste("Annual PM"[2.5]," concentration, ",mu,"g/m"^3)), ylab="Log(RR)")
lines(x, as.vector(as.matrix(results_quadratic_penalized[row_num,6:206]))-as.double(results_quadratic_penalized[row_num,6]), type = "l", lty=5, col = rainbow(3)[1])
lines(x, as.vector(as.matrix(results_quadratic_penalized[row_num+1,6:206]))-as.double(results_quadratic_penalized[row_num+1,6]), type = "l", lty=5, col = rainbow(3)[2])
lines(x, as.vector(as.matrix(results_quadratic_penalized[row_num+2,6:206]))-as.double(results_quadratic_penalized[row_num+2,6]), type = "l", lty=5, col = rainbow(3)[3])

row_num <- 370
plot(x, x*results_quadratic_penalized$b1[row_num]+(x^2)*results_quadratic_penalized$b2[row_num], type = "l",lwd=2,
     main=bquote(atop("Classical error",b[0]=="log(12), "~b[1]=="0.019, "~b[2]=="-0.0003")),
     xlab=expression(paste("Annual PM"[2.5]," concentration, ",mu,"g/m"^3)), ylab="Log(RR)")
lines(x, as.vector(as.matrix(results_quadratic_penalized[row_num,6:206]))-as.double(results_quadratic_penalized[row_num,6]), type = "l", lty=5, col = rainbow(3)[1])
lines(x, as.vector(as.matrix(results_quadratic_penalized[row_num+1,6:206]))-as.double(results_quadratic_penalized[row_num+1,6]), type = "l", lty=5, col = rainbow(3)[2])
lines(x, as.vector(as.matrix(results_quadratic_penalized[row_num+2,6:206]))-as.double(results_quadratic_penalized[row_num+2,6]), type = "l", lty=5, col = rainbow(3)[3])

row_num <- 364
plot(x, x*results_quadratic_penalized$b1[row_num]+(x^2)*results_quadratic_penalized$b2[row_num], type = "l",lwd=2,
     main=bquote(atop("Berkson error",b[0]=="log(12), "~b[1]=="0.019, "~b[2]=="-0.0003")),
     xlab=expression(paste("Annual PM"[2.5]," concentration, ",mu,"g/m"^3)), ylab="Log(RR)")
lines(x, as.vector(as.matrix(results_quadratic_penalized[row_num,6:206]))-as.double(results_quadratic_penalized[row_num,6]), type = "l", lty=5, col = rainbow(3)[1])
lines(x, as.vector(as.matrix(results_quadratic_penalized[row_num+1,6:206]))-as.double(results_quadratic_penalized[row_num+1,6]), type = "l", lty=5, col = rainbow(3)[2])
lines(x, as.vector(as.matrix(results_quadratic_penalized[row_num+2,6:206]))-as.double(results_quadratic_penalized[row_num+2,6]), type = "l", lty=5, col = rainbow(3)[3])

row_num <- 574
plot(x, x*results_quadratic_penalized$b1[row_num]+(x^2)*results_quadratic_penalized$b2[row_num], type = "l",lwd=2,
     main=bquote(atop("Classical error",b[0]=="log(20), "~b[1]=="0.012, "~b[2]=="0.0003")),
     xlab=expression(paste("Annual PM"[2.5]," concentration, ",mu,"g/m"^3)), ylab="Log(RR)")
lines(x, as.vector(as.matrix(results_quadratic_penalized[row_num,6:206]))-as.double(results_quadratic_penalized[row_num,6]), type = "l", lty=5, col = rainbow(3)[1])
lines(x, as.vector(as.matrix(results_quadratic_penalized[row_num+1,6:206]))-as.double(results_quadratic_penalized[row_num+1,6]), type = "l", lty=5, col = rainbow(3)[2])
lines(x, as.vector(as.matrix(results_quadratic_penalized[row_num+2,6:206]))-as.double(results_quadratic_penalized[row_num+2,6]), type = "l", lty=5, col = rainbow(3)[3])

row_num <- 568
plot(x, x*results_quadratic_penalized$b1[row_num]+(x^2)*results_quadratic_penalized$b2[row_num], type = "l",lwd=2,
     main=bquote(atop("Berkson error",b[0]=="log(20), "~b[1]=="0.012, "~b[2]=="0.0003")),
     xlab=expression(paste("Annual PM"[2.5]," concentration, ",mu,"g/m"^3)), ylab="Log(RR)")
lines(x, as.vector(as.matrix(results_quadratic_penalized[row_num,6:206]))-as.double(results_quadratic_penalized[row_num,6]), type = "l", lty=5, col = rainbow(3)[1])
lines(x, as.vector(as.matrix(results_quadratic_penalized[row_num+1,6:206]))-as.double(results_quadratic_penalized[row_num+1,6]), type = "l", lty=5, col = rainbow(3)[2])
lines(x, as.vector(as.matrix(results_quadratic_penalized[row_num+2,6:206]))-as.double(results_quadratic_penalized[row_num+2,6]), type = "l", lty=5, col = rainbow(3)[3])

# legend
par(mar = c(0,0,0,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend(x = "center",inset = 0,
       legend=c("Error-free (truth)", TeX("Error-prone with $\\sigma_{ik}=\\sigma^{*}_{ik}$"), TeX("Error-prone with $\\sigma_{ik}=2\\sigma^{*}_{ik}$"), TeX("Error-prone with $\\sigma_{ik}=3\\sigma^{*}_{ik}$")),
       col=c("black", rainbow(3)[1], rainbow(3)[2],rainbow(3)[3]),lty=c(1,5,5,5),
       lwd=2, cex=1.5, horiz = TRUE)

plot(0,type='n',axes=FALSE,ann=FALSE)

dev.off()



###################### 8. Fig S3. Soft-threshold dr model ########################
### extract results
files_softplus <- Sys.glob(paste0(dir_results_softplus,"results_set*.rds"))
files_softplus_ls <- lapply(files_softplus,readRDS)
results_softplus <- rbindlist(files_softplus_ls)
results_softplus <- results_softplus %>% group_by(k,b,type,corr_magnitude,pollutants) %>% summarise_all("mean")
results_softplus$corr_magnitude <- factor(results_softplus$corr_magnitude,levels=c('sp_error_sd','sp_error_2sd','sp_error_3sd',
                                                                                   'ind_error_sd','ind_error_2sd','ind_error_3sd'))
results_softplus <- results_softplus[order(results_softplus$k,results_softplus$b,results_softplus$type,results_softplus$corr_magnitude),]

### create plot
png(filename = paste0(dir_results_save,"Fig_S3_new.png"), width = 6, height = 12,units = "in",res = 300)
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
     type = "l",lwd=2,main=bquote(atop("Classical error",a=="0.2, "~b=="26")),
     xlab=expression(paste("Annual PM"[2.5]," concentration, ",mu,"g/m"^3)), ylab="Log(RR)")
lines(x, as.vector(as.matrix(results_softplus[row_num,6:306]))-as.double(results_softplus[row_num,6]), type = "l", lty=5, col = rainbow(3)[1])
lines(x, as.vector(as.matrix(results_softplus[row_num+1,6:306]))-as.double(results_softplus[row_num+1,6]), type = "l", lty=5, col = rainbow(3)[2])
lines(x, as.vector(as.matrix(results_softplus[row_num+2,6:306]))-as.double(results_softplus[row_num+2,6]), type = "l", lty=5, col = rainbow(3)[3])

row_num <- 52
plot(x, log(1+exp(results_softplus$k[row_num]*(x-results_softplus$b[row_num])))-log(1+exp(results_softplus$k[row_num]*(-results_softplus$b[row_num]))), 
     type = "l",lwd=2,main=bquote(atop("Berkson error",a=="0.2, "~b=="26")),
     xlab=expression(paste("Annual PM"[2.5]," concentration, ",mu,"g/m"^3)), ylab="Log(RR)")
lines(x, as.vector(as.matrix(results_softplus[row_num,6:306]))-as.double(results_softplus[row_num,6]), type = "l", lty=5, col = rainbow(3)[1])
lines(x, as.vector(as.matrix(results_softplus[row_num+1,6:306]))-as.double(results_softplus[row_num+1,6]), type = "l", lty=5, col = rainbow(3)[2])
lines(x, as.vector(as.matrix(results_softplus[row_num+2,6:306]))-as.double(results_softplus[row_num+2,6]), type = "l", lty=5, col = rainbow(3)[3])

row_num <- 94
plot(x, log(1+exp(results_softplus$k[row_num]*(x-results_softplus$b[row_num])))-log(1+exp(results_softplus$k[row_num]*(-results_softplus$b[row_num]))), 
     type = "l",lwd=2,main=bquote(atop("Classical error",a=="0.3, "~b=="24")),
     xlab=expression(paste("Annual PM"[2.5]," concentration, ",mu,"g/m"^3)), ylab="Log(RR)")
lines(x, as.vector(as.matrix(results_softplus[row_num,6:356]))-as.double(results_softplus[row_num,6]), type = "l", lty=5, col = rainbow(3)[1])
lines(x, as.vector(as.matrix(results_softplus[row_num+1,6:306]))-as.double(results_softplus[row_num+1,6]), type = "l", lty=5, col = rainbow(3)[2])
lines(x, as.vector(as.matrix(results_softplus[row_num+2,6:306]))-as.double(results_softplus[row_num+2,6]), type = "l", lty=5, col = rainbow(3)[3])

row_num <- 88
plot(x, log(1+exp(results_softplus$k[row_num]*(x-results_softplus$b[row_num])))-log(1+exp(results_softplus$k[row_num]*(-results_softplus$b[row_num]))), 
     type = "l",lwd=2,main=bquote(atop("Berkson error",a=="0.3, "~b=="24")),
     xlab=expression(paste("Annual PM"[2.5]," concentration, ",mu,"g/m"^3)), ylab="Log(RR)")
lines(x, as.vector(as.matrix(results_softplus[row_num,6:306]))-as.double(results_softplus[row_num,6]), type = "l", lty=5, col = rainbow(3)[1])
lines(x, as.vector(as.matrix(results_softplus[row_num+1,6:306]))-as.double(results_softplus[row_num+1,6]), type = "l", lty=5, col = rainbow(3)[2])
lines(x, as.vector(as.matrix(results_softplus[row_num+2,6:306]))-as.double(results_softplus[row_num+2,6]), type = "l", lty=5, col = rainbow(3)[3])

row_num <- 142
plot(x, log(1+exp(results_softplus$k[row_num]*(x-results_softplus$b[row_num])))-log(1+exp(results_softplus$k[row_num]*(-results_softplus$b[row_num]))), 
     type = "l",lwd=2,main=bquote(atop("Classical error",a=="0.4, "~b=="23")),
     xlab=expression(paste("Annual PM"[2.5]," concentration, ",mu,"g/m"^3)), ylab="Log(RR)")
lines(x, as.vector(as.matrix(results_softplus[row_num,6:306]))-as.double(results_softplus[row_num,6]), type = "l", lty=5, col = rainbow(3)[1])
lines(x, as.vector(as.matrix(results_softplus[row_num+1,6:306]))-as.double(results_softplus[row_num+1,6]), type = "l", lty=5, col = rainbow(3)[2])
lines(x, as.vector(as.matrix(results_softplus[row_num+2,6:306]))-as.double(results_softplus[row_num+2,6]), type = "l", lty=5, col = rainbow(3)[3])

row_num <- 136
plot(x, log(1+exp(results_softplus$k[row_num]*(x-results_softplus$b[row_num])))-log(1+exp(results_softplus$k[row_num]*(-results_softplus$b[row_num]))), 
     type = "l",lwd=2,main=bquote(atop("Berkson error",a=="0.4, "~b=="23")),
     xlab=expression(paste("Annual PM"[2.5]," concentration, ",mu,"g/m"^3)), ylab="Log(RR)")
lines(x, as.vector(as.matrix(results_softplus[row_num,6:306]))-as.double(results_softplus[row_num,6]), type = "l", lty=5, col = rainbow(3)[1])
lines(x, as.vector(as.matrix(results_softplus[row_num+1,6:306]))-as.double(results_softplus[row_num+1,6]), type = "l", lty=5, col = rainbow(3)[2])
lines(x, as.vector(as.matrix(results_softplus[row_num+2,6:306]))-as.double(results_softplus[row_num+2,6]), type = "l", lty=5, col = rainbow(3)[3])

# legend
par(mar = c(0,0,0,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend(x = "center",inset = 0,
       legend=c("Error-free (truth)", TeX("Error-prone with $\\sigma_{ik}=\\sigma^{*}_{ik}$")),
       col=c("black", rainbow(3)[1]),lty=c(1,5),
       lwd=2, cex=1.5, horiz = TRUE)

# legend
par(mar = c(0,0,0,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend(x = "center",inset = 0,
       legend=c(TeX("Error-prone with $\\sigma_{ik}=2\\sigma^{*}_{ik}$"), TeX("Error-prone with $\\sigma_{ik}=3\\sigma^{*}_{ik}$")),
       col=c(rainbow(3)[2],rainbow(3)[3]),lty=c(5,5),
       lwd=2, cex=1.5, horiz = TRUE)

dev.off()



############## 9. Table S2. coefficients and CIs in linear dr model ##############
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
