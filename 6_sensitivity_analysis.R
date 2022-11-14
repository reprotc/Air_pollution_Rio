################################################################################
#############################Exposure lags###################################### 
##############################Air pollution ####################################
################################################################################
################################################################################
#working paper Cortes et al


#load packages 
if(!require(dplyr)) install.packages("dplyr"); library(dplyr)
if(!require(data.table)) install.packages("data.table"); library(data.table)
if(!require(survival)) install.packages("survival"); library(survival)
if(!require(dlnm)) install.packages("dlnm"); library(dlnm)
library(splines)
library(lubridate)
library(survival)
library(gnm)
################################################################################
################################################################################
#read exposure lag dates

cvd=read.csv("geocoded_cvd_final.csv", stringsAsFactors = F)
cvd=cvd[,-1]
names(cvd)

cvd=cvd[,c(1,2,3,4,5,6,7,8,9,10,11,12,17,18,20,21)]
resp=read.csv("geocoded_resp_final.csv", stringsAsFactors = F)
names(resp)
resp=resp[,-1]
resp=resp[,c(1,2,3,4,5,6,7,8,9,10,11,12,17,18,20,21)]
death=rbind.data.frame(cvd, resp)

death$cause_group=death$cause_death
death$cause_group=substr(death$cause_group, 1,1)
death$cause_group= ifelse(death$cause_group=='J', "Resp", "CVD")
table(death$cause_group)   

death$Date=death$date_death
death$Date=ifelse(nchar(death$Date)<8, paste0("0",death$Date), death$Date)
death$Date=dmy(death$Date)
#write.csv(death, "death_all.csv")
death=subset(death, year(death$Date)<2018)

cvd=subset(death,  death$cause_group=="CVD")
resp=subset(death, death$cause_group=="Resp")

resp_t=aggregate(resp$rec_number,by=list(resp$Date),length)
cvd_t=aggregate(cvd$rec_number,by=list(cvd$Date),length)
colnames(cvd_t)=c('Date','n_deaths')
colnames(resp_t)=c('Date','n_deaths')
names(cvd_t)
head(resp_t)
head(cvd_t)
#read exposure data period 
pm10=read.csv("pm10_imputed_new.csv")
date=subset(pm10, select="Date")
stations=c("AV", "CA","CG", "IR" ,"PG", "SC", 'SP')
pm10=subset(pm10, select = stations)
pm10=cbind.data.frame(date, pm10)
head(pm10)
pm10=pm10 %>% mutate(exposure_mean=rowMeans(.[ , c("AV","CA", "CG", "IR", "PG","SC","SP")], na.rm=TRUE))
pm10=pm10[,c(1,9)]
names(pm10)
names(pm10)[2]="mean_pm10"


###ozone
ozone=read.csv("ozone_imputed_5df.csv")
date=subset(ozone, select="Date")
stations=c("AV", "BG", "CA","CG", "IR" ,"PG", "SC", 'SP')
ozone=subset(ozone, select = stations)
ozone=cbind.data.frame(date, ozone)
head(ozone)
ozone=ozone %>% mutate(exposure_mean=rowMeans(.[ , c("AV","BG","CA", "CG", "IR", "PG","SC","SP")], na.rm=TRUE))
ozone=ozone[,c(1,10)]
names(ozone)[2]="mean_ozone"

temp= read.csv(file = "temp_imputed_new.csv")
names(temp)
temp=temp[,-1]
names(temp)[1]="Date"
unique(temp$Date)
summary(temp)
names(temp)
temp=temp %>% mutate(exposure_mean=rowMeans(.[ , c("dcea83054","dcea83115","dcea83746", 
                                                           "inmet83743", "dcea83748",  "dcea83755", 
                                                           "inmetA621","inmetA625","AV",
                                                           "CG","IR","PG","SC","SP")], na.rm=TRUE))
temp=temp[,c(1,16)]
names(temp)[2]="mean_temp"




humid=fread("humid_imputed_new.csv")
humid=humid[,-1]
names(humid)
names(humid)[1]="Date"
unique(humid$Date)
summary(humid)

humid=humid %>% mutate(exposure_mean=rowMeans(.[ , c("dcea83054","dcea83115","dcea83746", 
                                                   "dcea83748",  "dcea83755", 
                                                   "inmetA621","inmetA625","AV",
                                                   "BG","IR","SC","SP")], na.rm=TRUE))



#absolute humidity
source("CalculateAbsHumidity.R")

humid$abs_humid=CalculateAbsHumidity(humid$exposure_mean, 
                                          temp$mean_temp, fahrenheit = FALSE, percent = TRUE)
names(humid)
humid=humid[,c(1,14,15)]
names(humid)[2]="mean_humid"
names(humid)[3]="mean_abs_humid"
humid$Date=as.character(humid$Date)

head(pm10)
head(ozone)
head(temp)
head(humid)

exposure

exposure=merge.data.frame(pm10, ozone, by="Date")
exposure=merge.data.frame(exposure, temp, by="Date")
exposure=merge.data.frame(exposure, humid, by="Date")
head(exposure)
 
exposure$time=rep(1:2192)

exposure$Date=as.Date(exposure$Date)
exposure$dow= as.factor(weekdays(exposure$Date))
exposure$month= as.factor(months(exposure$Date))
exposure$year= as.factor(format(exposure$Date, format="%Y") )


exposure %>% group_by(year) %>% 
  mutate(doy = length(Date))  

  
exposure$Date


cvd_t=merge(exposure,cvd_t,by="Date",all.x=T)
summary(cvd_t)

resp_t=merge(exposure,resp_t,by="Date",all.x=T)
summary(resp_t)

######################PM10 CVD

basis_cvd_pm10 <- crossbasis(cvd_t$mean_pm10, lag= c(0,5),
                              argvar=list(fun="lin"),
                              arglag=list(fun="ns", df=3))

summary(basis_cvd_pm10)

range(temp$mean_temp,na.rm=T)

#tknots=quantile(temp$temp_mean,  probs = c(.25,.50,.75))

basis_cvd_temp <- crossbasis(cvd_t$mean_temp, lag= c(0,5),
                             argvar=list(fun="ns", df=3),
                             arglag=list(fun="ns", df=3))




summary(basis_cvd_pm10)
summary(basis_cvd_temp)

library(gnm)
names(cvd_t)
md_cvd <- glm(n_deaths ~ basis_cvd_pm10 ,  data=cvd_t, family=poisson)

cp_ad=crosspred(basis_cvd_pm10, md_cvd, cen=0,  at=10, cumul = T)
range(exposure$mean_pm10)
a=IQR(exposure$mean_pm10)
plot(cp_ad, var=10, type="p", pch=18, cex=1.9, ci="bars", col="black",
     ylab="OR", ylim=c(0.96,1.04), xlab="Lag (days)", xlim=c(0,5), xaxt='n')
axis(side =1 , at=0:5)
mytitle = "CVD & pm10_TS"
mtext(side=3, line=1, at=-0.07, adj=0, cex=1.1, mytitle)


md_cvd <- glm(n_deaths ~  basis_cvd_pm10+ basis_cvd_temp+ ns(mean_abs_humid, df=3)+ ns(time,5) + dow,
              family=poisson(), data=cvd_t)

cp_ad=crosspred(basis_cvd_pm10, md_cvd, cen=0, at=10, cumul = T)

plot(cp_ad, var=10, type="p", pch=18, cex=1.9, ci="bars", col="black",
     ylab="OR", ylim=c(0.96,1.04), xlab="Lag (days)", xlim=c(0,5), xaxt='n')
axis(side =1 , at=0:5)
mytitle = "CVD & pm10_cc_TS"
mtext(side=3, line=1, at=-0.07, adj=0, cex=1.1, mytitle)

summary(md_cvd)
aic=round(AIC(md_cvd), digits = 1)

#plot

tiff("pm10_cvd_5lags_poisson_times_series", width = 980, height = 500,)
plot(cp_ad, var=10, type="b", pch=18, cex=1.9, ci="b", col="blue",
     ylab="OR", ylim=c(0.96,1.06), xlab="Lag (days)", xlim=c(0,5), xaxt='n')
axis(side =1 , at=0:5)
mytitle = " PM10 Cardiovascular diseases mortality TS"
mtext(side=3, line=1, at=-0.07, adj=0, cex=1.1, mytitle)
legend("topright", paste(aic), col=c('black', lwd=1.2))
dev.off()

OR_cum=cbind(cp_ad$allRRfit["10"], cp_ad$allRRlow["10"],cp_ad$allRRhigh["10"])
OR_lag=cbind(cp_ad$matRRfit["10",], cp_ad$matRRlow["10",], cp_ad$matRRhigh["10",])

table=rbind(OR_cum, OR_lag, aic)
names(table)
colnames(table) <- c("OR","IC_inf","IC_sup")
table=formatC(table, format = "f", digits = 3)
write.csv(table, file="pm10_cvd_5lags_poisson_times_series.csv")

table=rbind(OR_cum, OR_lag, aic)
names(table)
colnames(table) <- c("OR","IC_inf","IC_sup")
table=formatC(table, format = "f", digits = 2)
write.csv(table, file="pm10_cvd_5lags_poisson_times_series_2dg.csv")


######################PM10 resp

basis_resp_pm10 <- crossbasis(resp_t$mean_pm10, lag= c(0,5),
                             argvar=list(fun="lin"),
                             arglag=list(fun="ns", df=3))

summary(basis_resp_pm10)

range(temp$mean_temp,na.rm=T)

#tknots=quantile(temp$temp_mean,  probs = c(.25,.50,.75))

basis_resp_temp <- crossbasis(resp_t$mean_temp, lag= c(0,5),
                             argvar=list(fun="ns", df=3),
                             arglag=list(fun="ns", df=3))




summary(basis_resp_pm10)
summary(basis_resp_temp)

library(gnm)
names(resp_t)
md_resp <- glm(n_deaths ~ basis_resp_pm10 ,  data=resp_t, family=poisson, eliminate=factor(stratum))

cp_ad=crosspred(basis_resp_pm10, md_resp, cen=0,  at=10, cumul = T)
range(exposure$mean_pm10)
a=IQR(exposure$mean_pm10)
plot(cp_ad, var=10, type="p", pch=18, cex=1.9, ci="bars", col="black",
     ylab="OR", ylim=c(0.96,1.04), xlab="Lag (days)", xlim=c(0,5), xaxt='n')
axis(side =1 , at=0:5)
mytitle = "resp & pm10_TS"
mtext(side=3, line=1, at=-0.07, adj=0, cex=1.1, mytitle)


md_resp <- glm(n_deaths ~  basis_resp_pm10+ basis_resp_temp+ ns(mean_abs_humid, df=3)+ ns(time,5) + dow,
              family=poisson(), data=resp_t)

cp_ad=crosspred(basis_resp_pm10, md_resp, cen=0, at=10, cumul = T)

plot(cp_ad, var=10, type="p", pch=18, cex=1.9, ci="bars", col="black",
     ylab="OR", ylim=c(0.96,1.04), xlab="Lag (days)", xlim=c(0,5), xaxt='n')
axis(side =1 , at=0:5)
mytitle = "resp & pm10_cc_TS"
mtext(side=3, line=1, at=-0.07, adj=0, cex=1.1, mytitle)

summary(md_resp)
aic=round(AIC(md_resp), digits = 1)

#plot

tiff("pm10_resp_5lags_poisson_times_series", width = 980, height = 500,)
plot(cp_ad, var=10, type="b", pch=18, cex=1.9, ci="b", col="blue",
     ylab="OR", ylim=c(0.96,1.06), xlab="Lag (days)", xlim=c(0,5), xaxt='n')
axis(side =1 , at=0:5)
mytitle = "PM10 Respiratory diseases mortality TS"
mtext(side=3, line=1, at=-0.07, adj=0, cex=1.1, mytitle)
legend("topright", paste(aic), col=c('black', lwd=1.2))
dev.off()

OR_cum=cbind(cp_ad$allRRfit["10"], cp_ad$allRRlow["10"],cp_ad$allRRhigh["10"])
OR_lag=cbind(cp_ad$matRRfit["10",], cp_ad$matRRlow["10",], cp_ad$matRRhigh["10",])

table=rbind(OR_cum, OR_lag, aic)
names(table)
colnames(table) <- c("OR","IC_inf","IC_sup")
table=formatC(table, format = "f", digits = 3)
write.csv(table, file="pm10_resp_5lags_poisson_times_series.csv")

table=rbind(OR_cum, OR_lag, aic)
names(table)
colnames(table) <- c("OR","IC_inf","IC_sup")
table=formatC(table, format = "f", digits = 2)
write.csv(table, file="pm10_resp_5lags_poisson_times_series_2dg.csv")

######################ozone CVD

basis_cvd_ozone <- crossbasis(cvd_t$mean_ozone, lag= c(0,5),
                             argvar=list(fun="lin"),
                             arglag=list(fun="ns", df=3))

summary(basis_cvd_ozone)

range(temp$mean_temp,na.rm=T)

#tknots=quantile(temp$temp_mean,  probs = c(.25,.50,.75))

basis_cvd_temp <- crossbasis(cvd_t$mean_temp, lag= c(0,5),
                             argvar=list(fun="ns", df=3),
                             arglag=list(fun="ns", df=3))




summary(basis_cvd_ozone)
summary(basis_cvd_temp)

library(gnm)
names(cvd_t)
md_cvd <- gnm(n_deaths ~ basis_cvd_ozone ,  data=cvd_t, family=poisson, eliminate=factor(stratum))

cp_ad=crosspred(basis_cvd_ozone, md_cvd, cen=0,  at=10, cumul = T)
range(exposure$mean_ozone)
a=IQR(exposure$mean_ozone)
plot(cp_ad, var=10, type="p", pch=18, cex=1.9, ci="bars", col="black",
     ylab="OR", ylim=c(0.96,1.04), xlab="Lag (days)", xlim=c(0,5), xaxt='n')
axis(side =1 , at=0:5)
mytitle = "CVD & ozone_TS"
mtext(side=3, line=1, at=-0.07, adj=0, cex=1.1, mytitle)

md_cvd <- glm(n_deaths ~  basis_cvd_ozone+ basis_cvd_temp+ ns(mean_abs_humid, df=3)+ ns(time,5) + dow,
               family=poisson(), data=cvd_t)

cp_ad=crosspred(basis_cvd_ozone, md_cvd, cen=0, at=10, cumul = T)

plot(cp_ad, var=10, type="p", pch=18, cex=1.9, ci="bars", col="black",
     ylab="OR", ylim=c(0.96,1.04), xlab="Lag (days)", xlim=c(0,5), xaxt='n')
axis(side =1 , at=0:5)
mytitle = "CVD & ozone_cc_TS"
mtext(side=3, line=1, at=-0.07, adj=0, cex=1.1, mytitle)

summary(md_cvd)
aic=round(AIC(md_cvd), digits = 1)

#plot

tiff("ozone_cvd_5lags_poisson_times_series", width = 980, height = 500,)
plot(cp_ad, var=10, type="b", pch=18, cex=1.9, ci="b", col="blue",
     ylab="OR", ylim=c(0.96,1.06), xlab="Lag (days)", xlim=c(0,5), xaxt='n')
axis(side =1 , at=0:5)
mytitle = "Ozone Cardiovascular diseases mortality TS"
mtext(side=3, line=1, at=-0.07, adj=0, cex=1.1, mytitle)
legend("topright", paste(aic), col=c('black', lwd=1.2))
dev.off()

OR_cum=cbind(cp_ad$allRRfit["10"], cp_ad$allRRlow["10"],cp_ad$allRRhigh["10"])
OR_lag=cbind(cp_ad$matRRfit["10",], cp_ad$matRRlow["10",], cp_ad$matRRhigh["10",])

table=rbind(OR_cum, OR_lag, aic)
names(table)
colnames(table) <- c("OR","IC_inf","IC_sup")
table=formatC(table, format = "f", digits = 3)
write.csv(table, file="ozone_cvd_5lags_poisson_times_series.csv")

table=rbind(OR_cum, OR_lag, aic)
names(table)
colnames(table) <- c("OR","IC_inf","IC_sup")
table=formatC(table, format = "f", digits = 2)
write.csv(table, file="ozone_cvd_5lags_poisson_times_series_2dg.csv")

######################ozone resp

basis_resp_ozone <- crossbasis(resp_t$mean_ozone, lag= c(0,5),
                              argvar=list(fun="lin"),
                              arglag=list(fun="ns", df=3))

summary(basis_resp_ozone)

range(temp$mean_temp,na.rm=T)

#tknots=quantile(temp$temp_mean,  probs = c(.25,.50,.75))

basis_resp_temp <- crossbasis(resp_t$mean_temp, lag= c(0,5),
                             argvar=list(fun="ns", df=3),
                             arglag=list(fun="ns", df=3))




summary(basis_resp_ozone)
summary(basis_resp_temp)

library(gnm)
names(resp_t)
md_resp <- gnm(n_deaths ~ basis_resp_ozone ,  data=resp_t, family=poisson, eliminate=factor(stratum))

cp_ad=crosspred(basis_resp_ozone, md_resp, cen=0,  at=10, cumul = T)
range(exposure$mean_ozone)
a=IQR(exposure$mean_ozone)
plot(cp_ad, var=10, type="p", pch=18, cex=1.9, ci="bars", col="black",
     ylab="OR", ylim=c(0.96,1.04), xlab="Lag (days)", xlim=c(0,5), xaxt='n')
axis(side =1 , at=0:5)
mytitle = "resp & ozone_TS"
mtext(side=3, line=1, at=-0.07, adj=0, cex=1.1, mytitle)


md_resp <- glm(n_deaths ~  basis_resp_ozone+ basis_resp_temp+ ns(mean_abs_humid, df=3)+ ns(time,5) + dow,
              family=poisson(), data=resp_t)

cp_ad=crosspred(basis_resp_ozone, md_resp, cen=0, at=10, cumul = T)

plot(cp_ad, var=10, type="p", pch=18, cex=1.9, ci="bars", col="black",
     ylab="OR", ylim=c(0.96,1.04), xlab="Lag (days)", xlim=c(0,5), xaxt='n')
axis(side =1 , at=0:5)
mytitle = "resp & ozone_cc_TS"
mtext(side=3, line=1, at=-0.07, adj=0, cex=1.1, mytitle)

summary(md_resp)
aic=round(AIC(md_resp), digits = 1)

#plot

tiff("ozone_resp_5lags_poisson_times_series", width = 980, height = 500,)
plot(cp_ad, var=10, type="b", pch=18, cex=1.9, ci="b", col="blue",
     ylab="OR", ylim=c(0.96,1.06), xlab="Lag (days)", xlim=c(0,5), xaxt='n')
axis(side =1 , at=0:5)
mytitle = "Ozone Respiratory diseases mortality TS"
mtext(side=3, line=1, at=-0.07, adj=0, cex=1.1, mytitle)
legend("topright", paste(aic), col=c('black', lwd=1.2))
dev.off()

OR_cum=cbind(cp_ad$allRRfit["10"], cp_ad$allRRlow["10"],cp_ad$allRRhigh["10"])
OR_lag=cbind(cp_ad$matRRfit["10",], cp_ad$matRRlow["10",], cp_ad$matRRhigh["10",])

table=rbind(OR_cum, OR_lag, aic)
names(table)
colnames(table) <- c("OR","IC_inf","IC_sup")
table=formatC(table, format = "f", digits = 3)
write.csv(table, file="ozone_resp_5lags_poisson_times_series.csv")

table=rbind(OR_cum, OR_lag, aic)
names(table)
colnames(table) <- c("OR","IC_inf","IC_sup")
table=formatC(table, format = "f", digits = 2)
write.csv(table, file="ozone_resp_5lags_poisson_times_series_2dg.csv")


