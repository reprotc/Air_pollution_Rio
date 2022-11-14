#################################################################################
##################################Exposure###################################### 
##############################Air pollution ####################################
################################################################################
################################################################################
#working paper Cortes et al


#load packages 
if(!require(dplyr)) install.packages("dplyr"); library(dplyr)
if(!require(lubridate)) install.packages("lubridate"); library(lubridate)
if(!require(gstat)) install.packages("gstat"); library(gstat)
if(!require(sp)) install.packages("sp"); library(sp)
if(!require(raster)) install.packages("raster"); library(raster)
if(!require(rgdal)) install.packages("rgdal"); library(rgdal)
if(!require(sf)) install.packages("sf"); library(sf)
if(!require(tidyr)) install.packages("tidyr");library(tidyr)
if(!require(exactextractr)) install.packages("exactextractr");library(exactextractr)
if(!require(rgeos)) install.packages("rgeos");library(rgeos)
if(!require(data.table)) install.packages("data.table");library(data.table)


############################################################
# Read health data
cvd=read.csv("geocoded_cvd_final.csv", stringsAsFactors = F)
cvd=cvd[,-1]
names(cvd)

cvd=cvd[,c(1,2,3,5,6,7,8,9,10,11,12,17,18,20,21)]
cvd=subset(cvd, !is.na(cvd$latitude))
cvd=subset(cvd, cvd$non_Geocoded_reason=="geocoded")
#plot(cvd$longitude)

resp=read.csv("geocoded_resp_final.csv", stringsAsFactors = F)
names(resp)
resp=resp[,-1]
resp=resp[,c(1,2,3,5,6,7,8,9,10,11,12,17,18,20,21)]
resp=subset(resp, !is.na(resp$latitude))
resp=subset(resp, resp$non_Geocoded_reason=="geocoded")

death=rbind.data.frame(cvd, resp)
rm(cvd, resp)

#flags from Kim/Bellgroup
flags=read.csv("CVD_RESP_DEATH_LOCATION_INDICATOR.csv")
names(flags)
flags=subset(flags, select=c(rec_number, WithinRio))
death=merge(flags, death, by="rec_number", all.y=TRUE)
rm(flags)

table(death$non_Geocoded_reason, death$WithinRio)
death=subset(death, death$non_Geocoded_reason=="geocoded" & death$WithinRio==1)
names(death)


death$Date=death$date_death
death$Date=ifelse(nchar(death$Date)<8, paste0("0",death$Date), death$Date)
death$Date=dmy(death$Date)
#write.csv(death, "death_all.csv")

death=subset(death, year(death$Date)<2018)

##############################################################

path_1 <- "/home/tc/Documentos/NPEA_ATUAL/CaseCrossover dataset_pollutant/Pm10_cc"


tif_2012 <- list.files(path_1,
                       full.names = TRUE,
                       pattern = "^x201201|^x201202|^x201203|^x201204|^x201205|^x201206")

stack2012 <- stack(tif_2012)

#2012
names(stack2012)
z_2012  <- gsub("^x","", names(stack2012))
z_2012  <- gsub("_pm10_5km","", z_2012)
z_2012= ymd(z_2012)
stack2012=setZ(stack2012, z_2012 , name='Date')
getZ(stack2012)


seq_stack=names(stack2012)
stack_date=stack2012

 names(death)
 death_day=death[-c(2,4,4:12,15,16, 18:29)]
 names(death_day)
 death_day=death_day[c(1,5,3,4)]
 names(death_day)
 names(death_day)[2]="Date_death"
 death_day=death_day[,-2]

 coordinates(death_day) <- ~longitude + latitude
 proj4string(death_day) <- CRS("+init=epsg:4326") # WGS 84
 death_day<- spTransform(death_day, CRS("+proj=utm +zone=23 +south +ellps=aust_SA +units=m +no_defs"))

buffer_day <- gBuffer(death_day,width=1000, byid = T) # using gBuffer
names(buffer_day)

rm(stack2012)
rm(death, death_day)
exposure=NULL
gc()
for (i in seq_along(seq_stack)){
  print(i)
  stack_day=stack_date[[i]]

  exposure_i <- exact_extract(stack_day, buffer_day, 'weighted_mean', weights = 'area', force_df = T, append_cols = T)
  exposure_i=cbind.data.frame(stack_day@z, exposure_i)
  exposure=rbind.data.frame(exposure, exposure_i)
  
}

#
exposure=subset(exposure,!is.na(exposure$weighted_mean))
fwrite(exposure,"exposure_pm10_period_2012_01_06_5km.csv")

rm(stack2012, stack_day)
gc()
