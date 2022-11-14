################################################################################
#############################Main Analysis###################################### 
##############################Air pollution ####################################
################################################################################
################################################################################
#working paper Cortes et al. 
################################################################################
#load packages 
if(!require(dplyr)) install.packages("dplyr"); library(dplyr)
#if(!require(lubridate)) install.packages("lubridate"); library(lubridate)
if(!require(data.table)) install.packages("data.table"); library(data.table)
if(!require(tidyr)) install.packages("tidyr"); library(tidyr)

if(!require(survival)) install.packages("survival"); library(survival)
if(!require(dlnm)) install.packages("dlnm"); library(dlnm)
library(splines)
library(tidyverse)
################################################################################
# read air pollution data (case crossover)
#5lag
pm10_inf= fread( "exposure_pm10_backtlog_cc_lag0_10_2012_2017_5NOna.csv")
pm10_inf=pm10_inf[,-1]
names(pm10_inf)
pm10_inf=pm10_inf[,-c(10:15)]
names(pm10_inf)=gsub("lag", "pm10.lag", names(pm10_inf))

ozone_inf= fread( "exposure_ozone_backtlog_cc_lag0_10_2012_2017_5NOna.csv")
ozone_inf=ozone_inf[,-1]
names(ozone_inf)
ozone_inf=ozone_inf[,-c(10:15)]
names(ozone_inf)=gsub("lag", "ozone.lag", names(ozone_inf))

###########################################covariate data

temp= fread("exposure_temperature_cc_lag0_10_2012_2017_5NOna.csv")
temp=temp[,-1]
temp=temp[,-c(10:14)]
names(temp)
names(temp)=gsub("lag", "temp.lag", names(temp))

temp=mutate(temp, temp.lag0_5 = rowMeans(select(temp, c(temp.lag0, 
                                                        temp.lag1, 
                                                        temp.lag2,
                                                        temp.lag3,
                                                        temp.lag4,
                                                        temp.lag5))))

temp=mutate(temp, temp.lag0_1 = rowMeans(select(temp, c(temp.lag0, 
                                                        temp.lag1))))

temp=mutate(temp, temp.lag0_2 = rowMeans(select(temp, c(temp.lag0, 
                                                        temp.lag1, 
                                                        temp.lag2))))

temp=mutate(temp, temp.lag0_3 = rowMeans(select(temp, c(temp.lag0, 
                                                        temp.lag1, 
                                                        temp.lag2,
                                                        temp.lag3))))

temp=mutate(temp, temp.lag0_4 = rowMeans(select(temp, c(temp.lag0, 
                                                        temp.lag1, 
                                                        temp.lag2, temp.lag3,
                                                        temp.lag4))))

###humid

humid=fread("exposure_humidity_cc_lag0_10_2012_2017_5NOna.csv")
humid=humid[,-1]
names(humid)
humid=humid[,-c(10:15)]
names(humid)=gsub("lag", "humid.lag", names(humid))


#absolute humidity
source("CalculateAbsHumidity.R") # from Ismael script

humid$abs_humid.lag0=CalculateAbsHumidity(humid$humid.lag0, 
                                          temp$temp.lag0, fahrenheit = FALSE, percent = TRUE)

humid$abs_humid.lag1=CalculateAbsHumidity(humid$humid.lag1, 
                                          temp$temp.lag1, fahrenheit = FALSE, percent = TRUE)

humid$abs_humid.lag2=CalculateAbsHumidity(humid$humid.lag2, 
                                          temp$temp.lag2, fahrenheit = FALSE, percent = TRUE)

humid$abs_humid.lag3=CalculateAbsHumidity(humid$humid.lag3, 
                                          temp$temp.lag3, fahrenheit = FALSE, percent = TRUE)

humid$abs_humid.lag4=CalculateAbsHumidity(humid$humid.lag4, 
                                          temp$temp.lag4, fahrenheit = FALSE, percent = TRUE)

humid$abs_humid.lag5=CalculateAbsHumidity(humid$humid.lag5, 
                                          temp$temp.lag5, fahrenheit = FALSE, percent = TRUE)


humid=mutate(humid, abs_humid.lag0_1 = rowMeans(select(humid, c(abs_humid.lag0, 
                                                                abs_humid.lag1))))

humid=mutate(humid, abs_humid.lag0_2 = rowMeans(select(humid, c(abs_humid.lag0,
                                                                abs_humid.lag1, 
                                                                abs_humid.lag2))))

humid=mutate(humid, abs_humid.lag0_3 = rowMeans(select(humid, c(abs_humid.lag0, 
                                                                abs_humid.lag1,
                                                                abs_humid.lag2, 
                                                                abs_humid.lag3))))

humid=mutate(humid, abs_humid.lag0_4 = rowMeans(select(humid, c(abs_humid.lag0, 
                                                                abs_humid.lag1, 
                                                                abs_humid.lag2,
                                                                abs_humid.lag3,
                                                                abs_humid.lag4))))


humid=mutate(humid, abs_humid.lag0_5 = rowMeans(select(humid, c(abs_humid.lag0, 
                                                                abs_humid.lag1, 
                                                                abs_humid.lag2,
                                                                abs_humid.lag3,
                                                                abs_humid.lag4,
                                                                abs_humid.lag5))))




###############################################################################

data_inf=merge(pm10_inf, ozone_inf, by=c("rec_number", "case", "Lag0_date"), all = T)
data_inf=merge(data_inf, temp, by=c("rec_number", "case", "Lag0_date"), all = T)
data_inf=merge(data_inf, humid, by=c("rec_number", "case", "Lag0_date"), all = T)
covar=fread("covariate_data.csv")
data_inf=merge(data_inf, covar, by="rec_number")
data_inf_cvd=subset(data_inf, data_inf$cause_group=="CVD")
data_inf_resp=subset(data_inf, data_inf$cause_group=="Resp")

################################################################################
# pm10 resp analysis
names(data_inf_resp)
pm10_exp_resp_inf= data_inf_resp[,c(4:9)]
names(pm10_exp_resp_inf)

basis.pm10_resp_inf <- crossbasis(pm10_exp_resp_inf, lag= c(0,5),
                                  argvar=list(fun="lin"),
                                  arglag=list(fun="ns", df=3))


model_pm10_resp_inf=clogit(case~ basis.pm10_resp_inf + ns(abs_humid.lag0_5, 3)+
                                   ns(temp.lag0_5, 3)+
                                   strata(rec_number), data_inf_resp)

cp_ad=crosspred(basis.pm10_resp_inf, model_pm10_resp_inf, cen = 0, cumul = T)
summary(cp_ad)
#plot(cp_ad)
aic=round(AIC(model_pm10_resp_inf), digits = 1)

cp_resp_pm10=cp_ad

#plot

tiff("pm10_resp_inf", width = 980, height = 500,)
plot(cp_ad, var=10, type="b", pch=18, cex=1.9, ci="b", col="blue",
     ylab="OR", ylim=c(0.96,1.06), xlab="Lag (days)", xlim=c(0,5), xaxt='n')
axis(side =1 , at=0:5)
mytitle = "Respiratory diseases mortality"
mtext(side=3, line=1, at=-0.07, adj=0, cex=1.1, mytitle)
legend("topright", paste(aic), col=c('black', lwd=1.2))
dev.off()

OR_cum=cbind(cp_ad$allRRfit["10"], cp_ad$allRRlow["10"],cp_ad$allRRhigh["10"])
OR_lag=cbind(cp_ad$matRRfit["10",], cp_ad$matRRlow["10",], cp_ad$matRRhigh["10",])

table=rbind(OR_cum, OR_lag, aic)
names(table)
colnames(table) <- c("OR","IC_inf","IC_sup")
table=formatC(table, format = "f", digits = 3)
write.csv(table, file="pm10_resp.csv")

table=rbind(OR_cum, OR_lag, aic)
names(table)
colnames(table) <- c("OR","IC_inf","IC_sup")
table=formatC(table, format = "f", digits = 2)
write.csv(table, file="pm10_resp_inf_2dg.csv")

#PM10 CVD
pm10_exp_cvd_inf= data_inf_cvd[,c(4:9)]

basis.pm10_cvd_inf <- crossbasis(pm10_exp_cvd_inf, lag= c(0,5),
                                 argvar=list(fun="lin"),
                                 arglag=list(fun="ns", df=3))

model_pm10_cvd_inf=clogit(case~ basis.pm10_cvd_inf + ns(abs_humid.lag0_5, 3)+
                                  ns(temp.lag0_5, 3)+
                                  strata(rec_number), data_inf_cvd)

cp_ad=crosspred(basis.pm10_cvd_inf, model_pm10_cvd_inf, cen = 0, cumul = T)
#summary(cp_ad)
cp_cvd_pm10=cp_ad

tiff("pm10_cvd_inf", width = 980, height = 500,)
plot(cp_ad, var=10, type="b", pch=18, cex=1.9, ci="b", col="blue",
     ylab="OR", ylim=c(0.96,1.06), xlab="Lag (days)", xlim=c(0,5), xaxt='n')
axis(side =1 , at=0:5)
mytitle = "Cardiovascular diseases mortality"
mtext(side=3, line=1, at=-0.07, adj=0, cex=1.1, mytitle)
legend("topright", paste(aic), col=c('black', lwd=1.2))
dev.off()

OR_cum=cbind(cp_ad$allRRfit["10"], cp_ad$allRRlow["10"],cp_ad$allRRhigh["10"])
OR_lag=cbind(cp_ad$matRRfit["10",], cp_ad$matRRlow["10",], cp_ad$matRRhigh["10",])

table=rbind(OR_cum, OR_lag, aic)
names(table)
colnames(table) <- c("OR","IC_inf","IC_sup")
table=formatC(table, format = "f", digits = 3)
write.csv(table, file="pm10_cvd.csv")

table=rbind(OR_cum, OR_lag, aic)
names(table)
colnames(table) <- c("OR","IC_inf","IC_sup")
table=formatC(table, format = "f", digits = 2)
write.csv(table, file="pm10_cvd_2dg.csv")

###################ozone
ozone_exp_resp_inf= data_inf_resp[,c(10:15)]
names(ozone_exp_resp_inf)

basis.ozone_resp_inf <- crossbasis(ozone_exp_resp_inf, lag= c(0,5),
                                   argvar=list(fun="lin"),
                                   arglag=list(fun="ns", df=3))


model_ozone_resp_inf=clogit(case~ basis.ozone_resp_inf + ns(abs_humid.lag0_5, 3)+
                                    ns(temp.lag0_5, 3)+
                                    strata(rec_number), data_inf_resp)

cp_ad=crosspred(basis.ozone_resp_inf, model_ozone_resp_inf, cen = 0, cumul = T)
summary(cp_ad)
#plot(cp_ad)
aic=round(AIC(model_ozone_resp_inf), digits = 1)

cp_resp_ozone=cp_ad

#plot

tiff("ozone_resp_inf", width = 980, height = 500,)
plot(cp_ad, var=10, type="b", pch=18, cex=1.9, ci="b", col="blue",
     ylab="OR", ylim=c(0.96,1.06), xlab="Lag (days)", xlim=c(0,5), xaxt='n')
axis(side =1 , at=0:5)
mytitle = "Respiratory diseases mortality"
mtext(side=3, line=1, at=-0.07, adj=0, cex=1.1, mytitle)
legend("topright", paste(aic), col=c('black', lwd=1.2))
dev.off()

OR_cum=cbind(cp_ad$allRRfit["10"], cp_ad$allRRlow["10"],cp_ad$allRRhigh["10"])
OR_lag=cbind(cp_ad$matRRfit["10",], cp_ad$matRRlow["10",], cp_ad$matRRhigh["10",])

table=rbind(OR_cum, OR_lag, aic)
names(table)
colnames(table) <- c("OR","IC_inf","IC_sup")
table=formatC(table, format = "f", digits = 3)
write.csv(table, file="ozone_resp.csv")

table=rbind(OR_cum, OR_lag, aic)
names(table)
colnames(table) <- c("OR","IC_inf","IC_sup")
table=formatC(table, format = "f", digits = 2)
write.csv(table, file="ozone_resp_2dg.csv")


#OZONE CVD
ozone_exp_cvd_inf= data_inf_cvd[,c(10:15)]
names(ozone_exp_cvd_inf)

basis.ozone_cvd_inf <- crossbasis(ozone_exp_cvd_inf, lag= c(0,5),
                                  argvar=list(fun="lin"),
                                  arglag=list(fun="ns", df=3))

model_ozone_cvd_inf=clogit(case~ basis.ozone_cvd_inf + ns(abs_humid.lag0_5, 3)+
                                   ns(temp.lag0_5, 3)+
                                   strata(rec_number), data_inf_cvd)

cp_ad=crosspred(basis.ozone_cvd_inf, model_ozone_cvd_inf, cen = 0, cumul = T)
summary(cp_ad)
cp_cvd_ozone=cp_ad

tiff("ozone_cvd", width = 980, height = 500,)
plot(cp_ad, var=10, type="b", pch=18, cex=1.9, ci="b", col="blue",
     ylab="OR", ylim=c(0.96,1.06), xlab="Lag (days)", xlim=c(0,5), xaxt='n')
axis(side =1 , at=0:5)
mytitle = "Cardiovascular diseases mortality"
mtext(side=3, line=1, at=-0.07, adj=0, cex=1.1, mytitle)
legend("topright", paste(aic), col=c('black', lwd=1.2))
dev.off()

OR_cum=cbind(cp_ad$allRRfit["10"], cp_ad$allRRlow["10"],cp_ad$allRRhigh["10"])
OR_lag=cbind(cp_ad$matRRfit["10",], cp_ad$matRRlow["10",], cp_ad$matRRhigh["10",])

table=rbind(OR_cum, OR_lag, aic)
names(table)
colnames(table) <- c("OR","IC_inf","IC_sup")
table=formatC(table, format = "f", digits = 3)
write.csv(table, file="ozone_cvd.csv")

table=rbind(OR_cum, OR_lag, aic)
names(table)
colnames(table) <- c("OR","IC_inf","IC_sup")
table=formatC(table, format = "f", digits = 2)
write.csv(table, file="ozone_cvd_5lags_2dg.csv")


###############################################################################
#final plot

#final all 4 plots

tiff("Cardiovascular & Respiratory_pm10_ozone_main", width = 880, height = 800,  compression = "lzw")

par(mfrow=c(2,2))

plot(cp_resp_pm10, var=10, type="p", pch=15, cex=1.9, ci="bars", col=1,
     ylab="OR", ylim=c(0.98,1.04), xlab="Lag (days)", xlim=c(0,5))
axis(side = 1, at=0:5)

title("Respiratory diseases mortality")

mytitle = "(A) PM10"
mtext(side=3, line=0, at=-0.07, adj=0, cex=1, mytitle)


plot(cp_cvd_pm10, var=10, type="p", pch=15, cex=1.9, ci="bars", col=1,
     ylab="OR", ylim=c(0.98,1.04), xlab="Lag (days)", xlim=c(0,5))
axis(side = 1, at=0:5)

mytitle = "(C) PM10"
title("Cardiovascular diseases mortality")
mtext(side=3, line=0, at=-0.07, adj=0, cex=1, mytitle)


plot(cp_resp_ozone, var=10, type="p", pch=19, cex=1.9, ci="bars", col=1,
     ylab="OR", ylim=c(0.98,1.04), xlab="Lag (days)", xlim=c(0,5))

mytitle = "(B) O3"
mtext(side=3, line=0, at=-0.07, adj=0, cex=1, mytitle)
axis(side = 1, at=0:5)


plot(cp_cvd_ozone, var=10, type="p", pch=19, cex=1.9, ci="bars", col=1,
     ylab="OR", ylim=c(0.98,1.04), xlab="Lag (days)", xlim=c(0,5))

mytitle = "(D) O3"
mtext(side=3, line=0, at=-0.07, adj=0, cex=1, mytitle)
axis(side = 1, at=0:5)



dev.off()


table(data_inf_cvd$case)
table(data_inf_resp$case)
