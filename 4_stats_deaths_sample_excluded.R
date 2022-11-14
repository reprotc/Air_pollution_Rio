#################################################################################
##################################Analyses######################################
##############################Air pollution ####################################
################################################################################
################################################################################
#working paper Cortes et al. 

#load packages 
if(!require(dplyr)) install.packages("dplyr"); library(dplyr)
if(!require(lubridate)) install.packages("lubridate"); library(lubridate)
if(!require(data.table)) install.packages("data.table"); library(data.table)
if(!require(tidyr)) install.packages("tidyr"); library(tidyr)
library(tidyverse)


#stats fo all deaths, sample and excluded
# read air pollution data (case crossover)

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

humid=fread("exposure_humidity_cc_lag0_10_2012_2017_5NOna.csv")
humid=humid[,-1]
names(humid)
humid=humid[,-c(10:15)]
names(humid)=gsub("lag", "humid.lag", names(humid))


#inf km
data_inf=merge(pm10_inf, ozone_inf, by=c("rec_number", "case", "Lag0_date"), all = T)
data_inf=merge(data_inf, temp, by=c("rec_number", "case", "Lag0_date"), all = T)
data_inf=merge(data_inf, humid, by=c("rec_number", "case", "Lag0_date"), all = T)
covar=fread("covariate_data.csv")
data_inf=merge(data_inf, covar, by="rec_number")

data_inf_cvd=subset(data_inf, data_inf$cause_group=="CVD")
data_inf_resp=subset(data_inf, data_inf$cause_group=="Resp")

names(data_inf)

data_inf=data_inf[,c(1,2,3,29:48)]

data_inf=data_inf %>% group_by(rec_number) %>% 
  mutate(count2 = length(rec_number[case==1]))

table(data_inf$count2)

see=subset(data_inf, data_inf$count2==0)      
data_inf=subset(data_inf, !data_inf$count2==0)      

data=data_inf
data=subset(data, data$case==1)

table(data$cause_group)
names(data)

cvd=subset(data, data$cause_group=="CVD")
resp=subset(data, data$cause_group=="Resp")
table(resp$causes_ctg)
table(cvd$causes_ctg)

round(prop.table(table(resp$causes_ctg))*100, digits = 0)
round(prop.table(table(cvd$causes_ctg))*100, digits = 0)


data$age_ctg2=ifelse(data$age_edit<65,  "<65", 
                     ifelse(data$age_edit>64, "elderly (>=65 y)", "other"))

table(data$death_race_ctg)
round(prop.table(table(data$death_race_ctg))*100, digits = 0)



#Tem categoria 0???????????Grau de instrução (Escolaridade em anos. 
#Valores: 1 – Nenhuma; 2 – de 1 a 3 anos; 3 –de 4 a 7 anos;
#4 – de 8 a 11 anos; 5 – 12 anos e mais; 9 – Ignorado)

data$education_years
data$school_ctg=ifelse(data$education_years<4, "primary or none" ,
                       ifelse(data$education_years==9, "primary or none", ">7y"))
                     
table(data$school_ctg)
round(prop.table(table(data$school_ctg))*100, digits = 0)

round(prop.table(table(cvd$causes_ctg))*100, digits = 0)
round(prop.table(table(resp$causes_ctg))*100, digits = 0)

table(data$age_ctg2)
table(is.na(data$age_ctg2))
round(prop.table(table(data$age_ctg2))*100)

round(prop.table(table(resp$causes_ctg))*100, digits = 0)

table(data$gender)
table(is.na(data$gender))
round(prop.table(table(data$gender))*100)

names(data)
table(data$education_years)

data$season=ifelse(month(data$Date)==5,  "5-10  winter", 
                             ifelse(month(data$Date)==6,  "5-10  winter", 
                                    ifelse(month(data$Date)==7,  "5-10  winter", 
                                           ifelse(month(data$Date)==8,  "5-10  winter",
                                                  ifelse(month(data$Date)==9,  "5-10  winter", 
                                                         ifelse(month(data$Date)==10,  "5-10  winter", 
                                                                ifelse(month(data$Date)==12, "12-03 warm", 
                                                                       ifelse(month(data$Date)==1, "12-03 warm", 
                                                                              ifelse(month(data$Date)==2, "12-03 warm", 
                                                                                     ifelse(month(data$Date)==3, "12-03 warm", "other")))) ))))))

table(data$season)
round(prop.table(table(data$season))*100)




data$status="included"
###########################################################
cvd=read.csv("geocoded_cvd_final.csv", stringsAsFactors = F)
cvd=cvd[,-1]
names(cvd)

cvd=cvd[,c(1,2,3,4,5,6,7,8,9,10,11,12,17,18,20,21)]
cvd=subset(cvd, !is.na(cvd$latitude))
cvd=subset(cvd, cvd$non_Geocoded_reason=="geocoded")
#plot(cvd$longitude)

resp=read.csv("geocoded_resp_final.csv", stringsAsFactors = F)
names(resp)
resp=resp[,-1]
resp=resp[,c(1,2,3,4,5,6,7,8,9,10,11,12,17,18,20,21)]
resp=subset(resp, !is.na(resp$latitude))
resp=subset(resp, resp$non_Geocoded_reason=="geocoded")

death_all=rbind.data.frame(cvd, resp)



#############################################################

data_red=data[,c(1,3,26)]

death_all=merge(death_all, data_red, by="rec_number", all.x = T)
head(death_all)

death_all$status=ifelse(is.na(death_all$status), "excluded", death_all$status)

death_excluded=subset(death_all, death_all$status=="excluded")

death_excluded$cause_death


death_excluded$Date=death_excluded$date_death
death_excluded$Date=ifelse(nchar(death_excluded$Date)<8, paste0("0",death_excluded$Date), death_excluded$Date)
death_excluded$Date=dmy(death_excluded$Date)
#write.csv(death, "death_all.csv")


names(death_excluded)
table(death_excluded$gender)
table(death_excluded$race_color)
table(death_excluded$marital)
table(death_excluded$education_years)

death_excluded$marital_ctg=ifelse(death_excluded$marital==2, "married", ifelse(death_excluded$marital==5, "married",
                                                             ifelse(death_excluded$marital==9,NA, "not married")))

table(death_excluded$marital_ctg)

death_excluded$cause_group=death_excluded$cause_death
death_excluded$cause_group=substr(death_excluded$cause_group, 1,1)
death_excluded$cause_group= ifelse(death_excluded$cause_group=='J', "Resp", "CVD")
table(death_excluded$cause_group)                         
cvd=subset(death_excluded, death_excluded$cause_group=="CVD")
resp=subset(death_excluded, death_excluded$cause_group=="Resp")


#J44 DPOC
#j45 ASMA
# J20 acute bronchitis
#J21  Bronquiolite Aguda

death_excluded$causes_ctg=death_excluded$cause_death
death_excluded$causes_ctg

death_excluded$causes_ctg=ifelse(substring(death_excluded$cause_death,1,3)=="I20", "Doencas isquemicas do coracao (I20- I25)", 
                        ifelse(substring(death_excluded$cause_death,1,3)=="I21", "Doencas isquemicas do coracao (I20- I25)", 
                               ifelse(substring(death_excluded$cause_death,1,3)=="I22", "Doencas isquemicas do coracao (I20- I25)", 
                                      ifelse(substring(death_excluded$cause_death,1,3)=="I23","Doencas isquemicas do coracao (I20- I25)", 
                                             ifelse(substring(death_excluded$cause_death,1,3)=="I24","Doencas isquemicas do coracao (I20- I25)", 
                                                    ifelse(substring(death_excluded$cause_death,1,3)=="I25","Doencas isquemicas do coracao (I20- I25)", 
                                                           ifelse(substring(death_excluded$cause_death,1,3)=="I50","Insuficiencia Cardiaca  (I50)", 
                                                                  ifelse(substring(death_excluded$cause_death,1,3)=="I60","Doencas Cerebrovasculares I60-I69)", 
                                                                         ifelse(substring(death_excluded$cause_death,1,3)=="I61","Doencas Cerebrovasculares I60-I69)", 
                                                                                ifelse(substring(death_excluded$cause_death,1,3)=="I62","Doencas Cerebrovasculares I60-I69)", 
                                                                                       ifelse(substring(death_excluded$cause_death,1,3)=="I63","Doencas Cerebrovasculares I60-I69)", 
                                                                                              ifelse(substring(death_excluded$cause_death,1,3)=="I64","Doencas Cerebrovasculares I60-I69)", 
                                                                                                     ifelse(substring(death_excluded$cause_death,1,3)=="I65","Doencas Cerebrovasculares I60-I69)", 
                                                                                                            ifelse(substring(death_excluded$cause_death,1,3)=="I66","Doencas Cerebrovasculares I60-I69)", 
                                                                                                                   ifelse(substring(death_excluded$cause_death,1,3)=="I67","Doencas Cerebrovasculares I60-I69)", 
                                                                                                                          ifelse(substring(death_excluded$cause_death,1,3)=="I68","Doencas Cerebrovasculares I60-I69)", 
                                                                                                                                 ifelse(substring(death_excluded$cause_death,1,3)=="I69","Doencas Cerebrovasculares I60-I69)",
                                                                                                                                        ifelse(substring(death_excluded$cause_death,1,3)=="J44", "DPOC", 
                                                                                                                                               ifelse(substring(death_excluded$cause_death,1,3)=="J45", "Asma", 
                                                                                                                                                      ifelse(substring(death_excluded$cause_death,1,3)=="J20", "Bronquite Aguda", 
                                                                                                                                                             ifelse(substring(death_excluded$cause_death,1,3)=="J21", "Bronquiolite Aguda", "outras")))))))))))))))))))))


table(death_excluded$causes_ctg)


names(death_excluded)

summary(death_excluded$birth_date)

death_excluded$birth_date_ctg=death_excluded$birth_date
death_excluded$birth_date_ctg=ifelse(nchar(death_excluded$birth_date)<8, paste0("0",death_excluded$birth_date), death_excluded$birth_date)
death_excluded$birth_date_ctg=dmy(death_excluded$birth_date_ctg)
death_excluded$birth_date_ctg=as.Date(death_excluded$birth_date_ctg, format="%Y-%m-%d")
death_excluded$birth_date_ctg

death_excluded$Date=as.Date(death_excluded$Date, format="%Y-%m-%d")
death_excluded$Date

death_excluded$age_edit = year(as.period(interval(start = death_excluded$birth_date_ctg, end = death_excluded$Date)))

# infants <5 20–64 years; the elderly: ≥ 65 years
# death_excluded$age_ctg=ifelse(death_excluded$age_edit<5, "infants-preschoolers  (<5y)", 
#                      ifelse(death_excluded$age_edit>19 & 
#                             death_excluded$age_edit<65,  "adults (20-64 y)", 
#                             ifelse(death_excluded$age_edit>64, "elderly (>64 y)", "other")))
#                             

death_excluded$age_ctg=ifelse(death_excluded$age_edit<65,  "adults (20-64 y)", 
                     ifelse(death_excluded$age_edit>64, "elderly (>64 y)", "other"))


table(death_excluded$age_ctg)                         
round(prop.table(table(death_excluded$age_ctg))*100, digits = 0)                         

names(death_excluded)

table(death_excluded$race_color)
death_excluded$death_race_ctg=ifelse(death_excluded$race_color==1, "white","non white")
table(death_excluded$death_race_ctg)

table(death_excluded$death_race_ctg)
round(prop.table(table(death_excluded$death_race_ctg))*100, digits = 0)



#Tem categoria 0???????????Grau de instrução (Escolaridade em anos. 
#Valores: 1 – Nenhuma; 2 – de 1 a 3 anos; 3 –de 4 a 7 anos;
#4 – de 8 a 11 anos; 5 – 12 anos e mais; 9 – Ignorado)

death_excluded$education_years
death_excluded$school_ctg=ifelse(death_excluded$education_years<4, "primary or none" ,
                       ifelse(death_excluded$education_years==9, "primary or none", ">7y"))

table(death_excluded$school_ctg)
round(prop.table(table(death_excluded$school_ctg))*100, digits = 0)


table(cvd$causes_ctg)
table(resp$causes_ctg)

round(prop.table(table(cvd$causes_ctg))*100, digits = 0)
round(prop.table(table(resp$causes_ctg))*100, digits = 0)

table(death_excluded$age_ctg2)
table(is.na(death_excluded$age_ctg2))
round(prop.table(table(death_excluded$age_ctg2))*100)

round(prop.table(table(resp$causes_ctg))*100, digits = 0)

table(death_excluded$gender)
table(is.na(death_excluded$gender))
round(prop.table(table(death_excluded$gender))*100)

names(death_excluded)
table(death_excluded$education_years)

death_excluded$season=ifelse(month(death_excluded$Date)==5,  "5-10  winter", 
                              ifelse(month(death_excluded$Date)==6,  "5-10  winter", 
                                      ifelse(month(death_excluded$Date)==7,  "5-10  winter", 
                                              ifelse(month(death_excluded$Date)==8,  "5-10  winter",
                                                      ifelse(month(death_excluded$Date)==9,  "5-10  winter", 
                                                              ifelse(month(death_excluded$Date)==10,  "5-10  winter", 
                              ifelse(month(death_excluded$Date)==12, "12-03 warm", 
                                      ifelse(month(death_excluded$Date)==1, "12-03 warm", 
                                              ifelse(month(death_excluded$Date)==2, "12-03 warm", 
                                                      ifelse(month(death_excluded$Date)==3, "12-03 warm", "other")))) ))))))

table(death_excluded$season)
round(prop.table(table(death_excluded$season))*100)
