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
library(tidyverse)

setwd("~/Documentos/NPEA_ATUAL/CaseCrossover dataset_pollutant/Analysis Final")

################################################################################
#pm10 data
################################################################################
#crude data
pm10_all= read.csv(file = "pm10_all_transp_stations.csv")
names(pm10_all)
date=subset(pm10_all, select="Date")
stations=c("AV", "CA","CG", "IR" ,"PG", "SC", 'SP')
pm10_all=subset(pm10_all, select = stations)
pm10_all=cbind.data.frame(date, pm10_all)
head(pm10_all)
pm10_all=pm10_all %>% mutate(exposure_mean=rowMeans(.[ , c("AV","CA", "CG", "IR", "PG","SC","SP")], na.rm=TRUE))
pm10_all=pm10_all[,c(1,9)]

###imputed only 
pm10_imp=read.csv("pm10_imputed_new.csv")
date=subset(pm10_imp, select="Date")
stations=c("AV", "CA","CG", "IR" ,"PG", "SC", 'SP')
pm10_imp=subset(pm10_imp, select = stations)
pm10_imp=cbind.data.frame(date, pm10_imp)
head(pm10_imp)
pm10_imp=pm10_imp %>% mutate(exposure_mean=rowMeans(.[ , c("AV","CA", "CG", "IR", "PG","SC","SP")], na.rm=TRUE))
pm10_imp=pm10_imp[,c(1,9)]

# estimated air pollution data (case crossover)
pm10_inf= fread( "exposure_pm10_backtlog_cc_lag0_10_2012_2017_5NOna.csv")
pm10_inf=pm10_inf[,-1]
names(pm10_inf)
pm10_inf=pm10_inf[,-c(10:15)]
names(pm10_inf)=gsub("lag", "pm10.lag", names(pm10_inf))
names(pm10_inf)
pm10_inf=pm10_inf[,c(1,2,4)]
pm10_inf=subset(pm10_inf, pm10_inf$case==1)

summary(pm10_all$exposure_mean)
sd(pm10_all$exposure_mean, na.rm = T)
IQR(pm10_all$exposure_mean,na.rm = T)

summary(pm10_imp$exposure_mean)
sd(pm10_imp$exposure_mean, na.rm = T)
IQR(pm10_imp$exposure_mean, na.rm = T)

summary(pm10_inf$pm10.lag0)
sd(pm10_inf$pm10.lag0)
IQR(pm10_inf$pm10.lag0)

################################################################################
#ozone data
################################################################################
#crude data
ozone_all= read.csv(file = "ozone_all_INEA_SMAC.csv")
names(ozone_all)
date=subset(ozone_all, select="Date")
stations=c("AV", "BG", "CA","CG", "IR" ,"PG", "SC", 'SP')
ozone_all=subset(ozone_all, select = stations)
ozone_all=cbind.data.frame(date, ozone_all)
head(ozone_all)
ozone_all=ozone_all %>% mutate(exposure_mean=rowMeans(.[ , c("AV","BG","CA", "CG", "IR", "PG","SC","SP")], na.rm=TRUE))
ozone_all=ozone_all[,c(1,10)]

###imputed only 
ozone_imp=read.csv("ozone_imputed_5df.csv")
date=subset(ozone_imp, select="Date")
stations=c("AV", "BG", "CA","CG", "IR" ,"PG", "SC", 'SP')
ozone_imp=subset(ozone_imp, select = stations)
ozone_imp=cbind.data.frame(date, ozone_imp)
head(ozone_imp)
ozone_imp=ozone_imp %>% mutate(exposure_mean=rowMeans(.[ , c("AV","BG","CA", "CG", "IR", "PG","SC","SP")], na.rm=TRUE))
ozone_imp=ozone_imp[,c(1,10)]

# estimated air pollution data (case crossover)
ozone_inf= fread( "exposure_ozone_backtlog_cc_lag0_10_2012_2017_5NOna.csv")
ozone_inf=ozone_inf[,-1]
names(ozone_inf)
ozone_inf=ozone_inf[,-c(10:15)]
names(ozone_inf)=gsub("lag", "ozone.lag", names(ozone_inf))
names(ozone_inf)
ozone_inf=ozone_inf[,c(1,2,4)]

ozone_inf=subset(ozone_inf, ozone_inf$case==1)

summary(ozone_all$exposure_mean)
sd(ozone_all$exposure_mean, na.rm = T)
IQR(ozone_all$exposure_mean,na.rm = T)

summary(ozone_imp$exposure_mean)
sd(ozone_imp$exposure_mean, na.rm = T)
IQR(ozone_imp$exposure_mean, na.rm = T)

summary(ozone_inf$ozone.lag0)
sd(ozone_inf$ozone.lag0)
IQR(ozone_inf$ozone.lag0)

################################################################################
#temp data
################################################################################
#crude data

###imputed only 
temp_imp= read.csv(file = "temp_imputed_new.csv")
names(temp_imp)
temp_imp=temp_imp %>% mutate(exposure_mean=rowMeans(.[ , c("dcea83054","dcea83115","dcea83746", 
                                                           "inmet83743", "dcea83748",  "dcea83755", 
                                                          "inmetA621","inmetA625","AV",
                                                          "CG","IR","PG","SC","SP")], na.rm=TRUE))
temp_imp=temp_imp[,c(2,17)]


# estimated met data (case crossover)
temp_inf= fread( "exposure_temperature_cc_lag0_10_2012_2017_5NOna.csv")
temp_inf=temp_inf[,-1]
names(temp_inf)
temp_inf=temp_inf[,-c(10:15)]
names(temp_inf)=gsub("lag", "temp.lag", names(temp_inf))
names(temp_inf)
temp_inf=temp_inf[,c(1,2,4)]
temp_inf=subset(temp_inf, temp_inf$case==1)

# summary(temp_all$exposure_mean)
# sd(temp_all$exposure_mean, na.rm = T)
# IQR(temp_all$exposure_mean,na.rm = T)

summary(temp_imp$exposure_mean)
sd(temp_imp$exposure_mean, na.rm = T)
IQR(temp_imp$exposure_mean, na.rm = T)

summary(temp_inf$temp.lag0)
sd(temp_inf$temp.lag0)
IQR(temp_inf$temp.lag0)

################################################################################
#humid data
################################################################################
#crude data

###imputed only 
humid_imp= read.csv(file = "humid_imputed_new.csv")
names(humid_imp)
humid_imp=humid_imp %>% mutate(exposure_mean=rowMeans(.[ , c("dcea83054","dcea83115","dcea83746", 
                                                           "dcea83748",  "dcea83755", 
                                                           "inmetA621","inmetA625","AV", "BG",
                                                           "IR","SC","SP")], na.rm=TRUE))
humid_imp=humid_imp[,c(2,15)]


# estimated met data (case crossover)
humid_inf= fread( "exposure_humidity_cc_lag0_10_2012_2017_5NOna.csv")
humid_inf=humid_inf[,-1]
names(humid_inf)
humid_inf=humid_inf[,-c(10:15)]
names(humid_inf)=gsub("lag", "humid.lag", names(humid_inf))
names(humid_inf)
humid_inf=humid_inf[,c(1,2,4)]
humid_inf=subset(humid_inf, humid_inf$case==1)

# summary(humid_all$exposure_mean)
# sd(humid_all$exposure_mean, na.rm = T)
# IQR(humid_all$exposure_mean,na.rm = T)

summary(humid_imp$exposure_mean)
sd(humid_imp$exposure_mean, na.rm = T)
IQR(humid_imp$exposure_mean, na.rm = T)

summary(humid_inf$humid.lag0)
sd(humid_inf$humid.lag0)
IQR(humid_inf$humid.lag0)
