#################################################################################
#######################################IDW###################################### 
##############################Air pollution ####################################
###########################Temperature/humid####################################
################################################################################
#load packages 
if(!require(dplyr)) install.packages("dplyr"); library(dplyr)
if(!require(lubridate)) install.packages("lubridate"); library(lubridate)
if(!require(gstat)) install.packages("gstat"); library(gstat)
if(!require(sp)) install.packages("sp"); library(sp)
if(!require(raster)) install.packages("raster"); library(raster)
if(!require(rgdal)) install.packages("rgdal"); library(rgdal)
if(!require(sf)) install.packages("sf"); library(sf)
if(!require(stars)) install.packages("stars"); library(stars)
if(!require(tmap)) install.packages("tmap"); library(tmap)
if(!require(tidyr)) install.packages("tidyr");library(tidyr)
if(!require(ggplot2)) install.packages("ggplo2");library(ggplot2)
if(!require(rasterVis)) install.packages("rasterVis");library(rasterVis)
if(!require(DescTools)) install.packages("DescTools");library(DescTools)


############################################################
##tips and sources:
#https://www.geo.fu-berlin.de/en/v/soga/Geodata-analysis/geostatistics/
#https://github.com/Robinlovelace/geocompr/blob/master/06-reproj.Rmd
#https://gis-ops.com/measuring-distances-and-why-projections-matter-practical-examples/
#https://stackoverflow.com/questions/63727886/proj4-to-proj6-upgrade-and-discarded-datum-warnings
#https://www.nceas.ucsb.edu/sites/default/files/2020-04/OverviewCoordinateReferenceSystems.pdf

#############################################################
# Load RJ boudary map
#source: https://www.data.rio/datasets/limite-do-munic%C3%ADpio-do-rio-de-janeiro

rj=readOGR(dsn ="/home/npea/Documentos/NPEA_ATUAL/IDW_humidity", layer = "LimitemunicRJ")
proj4string(rj) <- CRS("+proj=utm +zone=23 +south +ellps=aust_SA +units=m +no_defs")


# Creating Grids for Rio 
grd_2 <- expand.grid(x=seq(from=bbox(rj)[1,1], to=bbox(rj)[1,2], by=1000), 
                     y=seq(from=bbox(rj)[2,1], to=bbox(rj)[2,2], by=1000))

coordinates(grd_2) <- c("x", "y")
gridded(grd_2) <- TRUE
gridded(grd_2)
proj4string(grd_2) <- CRS("+proj=utm +zone=23 +south +ellps=aust_SA +units=m +no_defs")
gridparameters(grd_2)

#read pollutant data
polmet_data=read.csv(file="var_data_imputed.csv")
names(polmet_data)
polmet_data=polmet_data[,-1]

polmet_data$year=year(polmet_data$Date)
lapply(polmet_data, function(x){ length(which(x==0))}) #SP=8days with 0
polmet_data$variable_log=log(polmet_data$variable+1)
hist(polmet_data$variable_log)
hist(polmet_data$variable)


names(polmet_data)
names(polmet_data)[2]="Date"
polmet_data=polmet_data[,c(1,2,3,4,5)]
polmet_data$Date=gsub("-", "",polmet_data$Date)
polmet_data$Date=paste0("x", polmet_data$Date)
names(polmet_data)
polmet_data_T=spread(polmet_data, Date, variable)


colnames(polmet_data_T)[1:10]
polmet_data_T=polmet_data_T[,-1]
coordinates(polmet_data_T) <- ~longitude + latitude
class(polmet_data_T)
crs(polmet_data_T)
proj4string(polmet_data_T) <- CRS("+init=epsg:4326") # WGS 84
polmet_data_T_r<- spTransform(polmet_data_T, CRS("+proj=utm +zone=23 +south +ellps=aust_SA +units=m +no_defs")) 
crs(polmet_data_T_r)

list_daysperiod=unique(polmet_data$Date)
list_daysperiod=as.data.frame(list_daysperiod)
names(list_daysperiod)[1]="Date"

list_daysperiod$Date_2=gsub("^x", "", list_daysperiod$Date)
list_daysperiod$Date_2=ymd(list_daysperiod$Date_2)
list_daysperiod=unique(list_daysperiod$Date)

# loocv 
seldays=list_daysperiod
seqPower = seq(from = 1, to = 3, 0.5)
seqNeighbors = seq(from = 1, to = 13, 1)
maxdist_1=c(10000,15000,20000, Inf)

spatialDF=polmet_data_T_r
polmet_data_T_r
gridded(grd_2)

#create combinations of parameters
cv.Grid <- expand.grid(Power = seqPower,
                       Neighbors = seqNeighbors,
                       Day=seldays,
                       Search_dist=maxdist_1)
cv.Grid$RMSE <- NA
cv.Grid$MAE <- NA
cv.Grid$RMSE_2 <- NA


#loop cv

cv.Grid=read.csv("cv.Grid_variable_multidist.csv")
names(cv.Grid)
cv.Grid=cv.Grid[,-1]

cv_var=subset(cv_var, !is.na(cv_var$RMSE))

for (i in 1:nrow(cv.Grid)){
  print(i)
  idw <- gstat(formula =   
                 as.formula(paste(cv.Grid[i, 'Day'], 1, sep =" ~ ")),
               data = spatialDF, 
               nmax =as.numeric(cv.Grid[i, 'Neighbors']), 
               set = list(idp = as.numeric(cv.Grid[i, 'Power'])), 
               maxdist=as.numeric(cv.Grid[i, "Search_dist"]))
  
  crossval <- gstat.cv(idw, 
                       nmax =as.numeric(cv.Grid[i, 'Neighbors']),
                       Power = as.numeric(cv.Grid[i, 'Power']),
                       maxdist= as.numeric(cv.Grid[i, "Search_dist"]),
                       debug.level = 0)
  
  cv.Grid[i, 'RMSE'] <- DescTools::RMSE(crossval$var1.pred, crossval$observed, na.rm=T)
  cv.Grid[i, 'MAE'] <- DescTools::MAE(crossval$var1.pred, crossval$observed, na.rm=T)
  cv.Grid[i, 'RMSE_2'] <- DescTools::RMSE(crossval$var1.pred, crossval$observed)

  }


cv.Grid$Date=gsub("^x", "", cv.Grid$Day)
cv.Grid$Date=ymd(cv.Grid$Date)
write.csv(cv.Grid, file="cv.Grid_humidity_multdist_nolog_final_3.csv")

#do it a scatter-plot for cv results!

names(cv.Grid)
cv.Grid=cv.Grid[,-c(9,10)]

cv.Grid0=read.csv("cv.Grid_var_2012_2013_ndist.csv")
cv.Grid1=read.csv("cv.Grid_var_2014_2015_ndist.csv")
cv.Grid2=read.csv("cv.Grid_var_2016_2017_ndist.csv")

names(cv.Grid0)
names(cv.Grid1)
cv.Grid1=cv.Grid1[,-c(1)]
cv.Grid0=cv.Grid0[,-c(1)]
cv.Grid2=cv.Grid2[,-c(1)]

cv.Grid=rbind.data.frame(cv.Grid0, cv.Grid1)
cv.Grid=rbind.data.frame(cv.Grid, cv.Grid2)

rm(cv.Grid0, cv.Grid1)


cv.Grid1=read.csv("cv.Grid_variable_log_multdist_1.csv")
cv.Grid2=read.csv("cv.Grid_variable_log_multdist_2.csv")

cv.Grid_dist=rbind.data.frame(cv.Grid1,cv.Grid2)
names(cv.Grid_dist)
cv.Grid_dist=cv.Grid_dist[,-c(1)]
rm(cv.Grid1, cv.Grid2)

#loocv results
cv.Grid$RMSE_nlog=exp(cv.Grid$RMSE)-1
cv.Grid$MAE_nlog=exp(cv.Grid$MAE)-1

cv.Grid_dist$RMSE_nlog=exp(cv.Grid_dist$RMSE)-1
cv.Grid_dist$MAE_nlog=exp(cv.Grid_dist$MAE)-1

summary(best.cv.Grid$RMSE_nlog)


best.cv.Grid=cv.Grid %>%
  group_by(Day) %>%
  slice(which.min(RMSE))

best.cv.Grid=cv.Grid %>%
  group_by(Day, Search_dist) %>%
  slice(which.min(RMSE))

table(best.cv.Grid$Power, best.cv.Grid$Neighbors)


#best.cv.Grid=subset(best.cv.Grid_2, best.cv.Grid_2$Search_dist==Inf)
best.cv.Grid_15=subset(best.cv.Grid, best.cv.Grid$Search_dist==15000)
best.cv.Grid_20=subset(best.cv.Grid, best.cv.Grid$Search_dist==20000)
best.cv.Grid_10=subset(best.cv.Grid, best.cv.Grid$Search_dist==10000)


################################################################################
#rm(cv.Grid, crossval)

best.cv.Grid=as.data.frame(best.cv.Grid)
names(best.cv.Grid)
idw.best.output=data.frame(Longitude=NA, Latitude=NA, variable=NA, var1.var=NA,  Day=NA)

#run the idw
for (i in 1:nrow(best.cv.Grid)) {
  print(i)
  best.idw <- gstat(
    formula =
      as.formula(paste(best.cv.Grid[i, 'Day'], 1, sep = " ~ ")),
    data = spatialDF,
    nmax = best.cv.Grid[i, 'Neighbors'],
    set = list(idp = best.cv.Grid[i, 'Power']))

  idw.best.predict <- predict(object = best.idw,
                              newdata = grd_2)
  r=raster(idw.best.predict)
  rf=writeRaster(r, filename= paste0(best.cv.Grid[i, 'Day'], "_variable_log_dist.tif"), format="GTiff")
  
  idw.best.data= as.data.frame(idw.best.predict)
  idw.best.data$Day=best.cv.Grid[i, 'Day']
  names(idw.best.data)[1:3] <- c("Longitude", "Latitude", "variable")
  idw.best.output=rbind.data.frame(idw.best.output, idw.best.data)
  
}
idw.best.output=subset(idw.best.output, 
                       !is.na(idw.best.output$Day))
idw.best.output$Date=gsub("^x", "", idw.best.output$Day)
idw.best.output$Date=ymd(idw.best.output$Date)
write.csv(idw.best.output, "idw.best.output_log.csv")
