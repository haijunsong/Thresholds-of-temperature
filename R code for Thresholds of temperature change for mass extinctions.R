#####R code for "Thresholds of temperature change for mass extinctions"
setwd("set Working Directory")

##########################
### Conodont O isotope ###
##########################
library("MASS")
library("tmvtnorm")
library("mvtnorm")
library("gmm")
library("sandwich")
library("Rcpp")
library("ff")
library("bit")
library("propagate")
COdata = read.csv("conodont data.csv",header = TRUE) #load data
esresult = matrix(data = NA,nrow(COdata),1e4) #create a matrix to record estimated temperature data

##propagate parameter uncertainties to temperature estimation 
for (i in 1:nrow(COdata)){
n = 1e4
EXPR1 <- expression(x - y*(ovalue - Os)) ##based on the equation (2) of Puc¨¦at et al., 2010. EPSL
Os = COdata$Os[i]
ovalue = COdata$O[i]
x <- rnorm(n, 118.7, 4.9)
y <- rnorm(n, 4.22, 0.2)
DF1 <- cbind(x, y,ovalue,Os)
RES1 <- propagate(expr = EXPR1, data = DF1, type = "raw",
                  nsim = n)
esresult[i,] = RES1$resSIM
}

##export quantiles of estimated temperature data
T_result_conodont = matrix(data = NA,nrow(COdata),8)
for (i in 1:nrow(COdata)){
  T_result_conodont[i,] = c(quantile(esresult[i,],0.025),quantile(esresult[i,],0.05),
                            quantile(esresult[i,],0.15),mean(esresult[i,]),
                            quantile(esresult[i,],0.85),quantile(esresult[i,],0.95),
                            quantile(esresult[i,],0.975),sd(esresult[i,]))
}
write.csv(T_result_conodont,file = "temperature_conodont.csv")

###########################
### carbonate O isotope ###
###########################

library("bayfoxr") 

MCMCD = get_draws(foram = NULL, seasonal_seatemp = FALSE)

Odata = read.csv("carbonate data.csv",header = TRUE)

#prior was set using SST = 16.1 - 4.64*(deta_O_c -deta_O_sw)*1000 + 0.09*[(deta_O_c -deta_O_sw)*1000]^2
sst <- predict_seatemp(Odata$d18o, d18osw=-1, prior_mean=19.82, prior_std= 4.13)  
seatpr = sst$ensemble  # extract estimated temperature
nrow(seatpr)
# create a matrix to record quantiles of estimated temperature
# we have 147 carbonate O isotope to estimated
T_result_carbonate = matrix(data = NA, 101, 8)
for(i in 1:101){
  T_result_carbonate[i,] = c(quantile(seatpr[i,],0.025),quantile(seatpr[i,],0.05),quantile(seatpr[i,],0.15),
                   mean(seatpr[i,]),quantile(seatpr[i,],0.85),quantile(seatpr[i,],0.95),
                   quantile(seatpr[i,],0.975),sd(seatpr[i,]))
}

#export estimate temperature data
write.csv(T_result_carbonate, file = "temperature_carbonate.csv")


###########################
###         TEX86       ###
###########################

##estimated using matlab package BAYSPAR

#################################################################################
######## Define a function (TRSIM) to estimate deta temperature and rate ########
#################################################################################
##parameters defined in the function
#mean, mean temperature
#sd, standard deviation of temperature
#t, two options t0 or t1, 
#age, the age of each data

#out_put estimated deta temperature, time and rate

Tdata = read.csv("Temperature data.csv",header = TRUE) #load data¡¢

TRSIM = function(Tdata, b){
  bindata =Tdata[which(Tdata$binno == b),] ##extract studied bin data
  bindatat0 = bindata[which(bindata$t == "t0"),] 
  bindatat1 = bindata[which(bindata$t == "t1"),] 
  
  ##deta temperature
  #t1_T
  sumt1T = rep(0,1e4)
  for(i in 1:nrow(bindatat1)){
    sumt1T = sumt1T + rnorm(1e4,bindatat1$mean[i],bindatat1$sd[i])
  }
  #t0_T
  sumt0T = rep(0,1e4)
  for(j in 1:nrow(bindatat0)){
    sumt0T = sumt0T + rnorm(1e4,bindatat0$mean[j],bindatat0$sd[j])
  }
  detaT = sumt1T/nrow(bindatat1) - sumt0T/nrow(bindatat0)
  
  Tsim = c(quantile(detaT,0.025),quantile(detaT,0.05),
           quantile(detaT,0.15),mean(detaT),
           quantile(detaT,0.85),quantile(detaT,0.95),
           quantile(detaT,0.975),sd(detaT))
  
  ##deta Time
  detaTime = abs(runif(1e4,min(bindatat1$age),max(bindatat1$age)) - 
                   runif(1e4,min(bindatat0$age),max(bindatat0$age)))
  
  timesim = c(quantile(detaTime,0.025),quantile(detaTime,0.05),
              quantile(detaTime,0.15),mean(detaTime),
              quantile(detaTime,0.85),quantile(detaTime,0.95),
              quantile(detaTime,0.975),sd(detaTime))
  
  ##Rate
  Rate = detaT/detaTime
  
  ratesim = c(quantile(Rate,0.025),quantile(Rate,0.05),
              quantile(Rate,0.15),mean(Rate),
              quantile(Rate,0.85),quantile(Rate,0.95),
              quantile(Rate,0.975),sd(Rate))
  
  return(c(Tsim,timesim,ratesim))
}

result_m = matrix(data = NA, 100, 24) 
for (i in 1:4){
  result_m[i,] = TRSIM(COdata,i)
}

write.csv(result_m, file = "result_m.csv")

###calculating mean rate and T
data = read.csv("trdata.csv",header = TRUE) #load data

mean_TR = function(k){
  bin_data = data[which(data$binno == k),]
  n = nrow(bin_data)
  sumT = rep(0,1e4)
  sumR = rep(0,1e4)
  for (i in 1:n){
    a = bin_data$meant[i]
    b = bin_data$sdt[i]
    c = bin_data$meanr[i]
    d = bin_data$sdr[i]
    sumT = sumT + rnorm(1e4,a,b)
    sumR = sumR + rnorm(1e4,c,d)
}
mean_T = sumT/n
mean_R = sumR/n

T_data = c(quantile(mean_T,0.025),quantile(mean_T,0.05),quantile(mean_T,0.15),mean(mean_T),
           quantile(mean_T,0.85),quantile(mean_T,0.95),quantile(mean_T,0.975),sd(mean_T),
           quantile(mean_R,0.025),quantile(mean_R,0.05),quantile(mean_R,0.15),
           mean(mean_R),quantile(mean_R,0.85),quantile(mean_R,0.95),
           quantile(mean_R,0.975),sd(mean_R))
return(T_data)
}

Mean_result = matrix(data = NA, 6, 16)
for (i in 1:6){
  Mean_result[i,] = mean_TR(i)
}

write.csv(Mean_result, file = "Mean_result.csv")

###mean time
T_SIM = function(Tdata, b){
  bindata =Tdata[which(Tdata$binno == b),] ##extract studied bin data
  bindatat0 = bindata[which(bindata$t == "t0"),] 
  bindatat1 = bindata[which(bindata$t == "t1"),] 

  
  ##deta Time
  detaTime = abs(runif(1e4,min(bindatat1$age),max(bindatat1$age)) - 
                   runif(1e4,min(bindatat0$age),max(bindatat0$age)))
  
  timesim = c(quantile(detaTime,0.025),quantile(detaTime,0.05),
              quantile(detaTime,0.15),mean(detaTime),
              quantile(detaTime,0.85),quantile(detaTime,0.95),
              quantile(detaTime,0.975),sd(detaTime))
  
  return(detaTime)
}

#example for calculating a mean time of four samples
mean_deta_Time = matrix(NA, 5, 8)
time_1 = T_SIM(data,1)
time_2 = T_SIM(data,2)
time_3 = T_SIM(data,3)
time_4 = T_SIM(data,4)

meantime = (time_1+time_2+time_3+time_4)/4

mean_deta_Time[1,] = c(quantile(meantime,0.025),quantile(meantime,0.05),
            quantile(meantime,0.15),mean(meantime),
            quantile(meantime,0.85),quantile(meantime,0.95),
            quantile(meantime,0.975),sd(meantime))
write.csv(mean_deta_Time, file = "C:\\Users\\DAI Xu\\Desktop\\time_out.csv")

################################
### Autocorrelation function ###
################################

data = read.csv("acf_data.csv", header = TRUE)
log_T = log10(abs(data$T))
log_R = log10(abs(data$R))
log_GF = log10(abs(data$GF))
acf(log_T,type = "correlation")
acf(log_R,type = "correlation")
acf(log_GF,type = "correlation")
acf(data$GF,type = "correlation")
acf(data$T,type = "correlation")
acf(data$R,type = "correlation")

############ comparison correlations using "cocor" ##############
install.packages("cocor")
library(cocor)

TECO = read.csv("acf_data.csv", header = TRUE)
log_R = log10(abs(TECO$R))
TECO = cbind.data.frame(TECO,log_R)

##deta T and rate correlations comparison
cocor.result1 <- cocor(~GF + T | GF + logR,TECO) # comparing the correlations 
as.htest(cocor.result1)  # showing results

##cooling and warming correlations comparison
cocor.result2 <- cocor(~WGF + WT | CGF + CT,TECO)
as.htest(cocor.result2)

cocor.result3 <- cocor(~WGF + WlogR | CGF + ClogR,TECO)
as.htest(cocor.result3)

## paleozoic-mesozoic-cenozoic correlations
cocor.result4 <- cocor(~CeGF + CelogR | MeGF + MelogR,TECO)
as.htest(cocor.result4)

cocor.result5 <- cocor(~CeGF + CelogR | PaGF + PalogR,TECO)
as.htest(cocor.result5)

cocor.result6 <- cocor(~MeGF + MelogR | PaGF + PalogR,TECO)
as.htest(cocor.result6)

cocor.result7 <- cocor(~CeGF + CeT | MeGF + MeT,TECO)
as.htest(cocor.result7)

cocor.result8 <- cocor(~CeGF + CeT | PaGF + PaT,TECO)
as.htest(cocor.result8)

cocor.result9 <- cocor(~MeGF + MeT | PaGF + PaT,TECO)
as.htest(cocor.result9)
