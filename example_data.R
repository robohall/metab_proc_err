
#Script  get data together for two sites


library(streamMetabolizer)
library(tidyverse)
library(lubridate)


bpcalc_atm<- function(bpst, alt) {
  bpst*exp((-9.80665*0.0289644*alt)/(8.31447*(273.15+15)))
}
bpcalc_atm(bpst=bpstd, alt=altitude)
Kcor<-function (temp,K600) {
  K600/(600/(1800.6-(temp*120.1)+(3.7818*temp^2)-(0.047608*temp^3)))^-0.5
}


gallatin_downstream<-read.csv("/Users/bob.hall/Dropbox/met_error_models//metabolism process error/GallatinTestDownstream.csv")
oxy<-gallatin_downstream
head(oxy)  #look ok??
lat<- 45.5
long<- -111.2 #(note - for us folks in western hemisphere)
altitude<- 5400/3.28 # m asl
bpstd <- 1013 #Do not change. ignoring weather relate changes in bp for this minimal code.  1013 mb = 760 mm Hg.  
K600_guess <- 30 #your best guess of K600 in units of 1/d.  Big slow rivers, 0.5 to 2.  Medium slow rivers, 2-10. 
#Fun rivers to float 10-20.  Non bubbly streams 20-40.  Bubbly streams 50-200.
mean_depth <- 0.5 
oxy$time<-as_datetime(oxy$unixtime) #converts unit time to a R time object in UTC
oxy$solar.time<-convert_UTC_to_solartime(oxy$time, longitude= long, time.type="mean solar")
oxy$light<- calc_light(oxy$solar.time, latitude=lat, longitude=long, max.PAR =2326, attach.units = F)
oxy$osat<- calc_DO_sat(temp.water=oxy$temp, pressure.air=bpcalc_atm(bpst=bpstd, alt=altitude))
oxy_g<- oxy[oxy$solar.time > ymd_hms("2024-09-01 04:00:00") & oxy$solar.time < ymd_hms("2024-10-01 04:00:00"),  ]
oxy_g$Kc<-Kcor(oxy_g$temp,1)
ntime<-144
nday<-length(oxy_g$oxy)/ntime
oxy_g$day<- rep(1:nday, each=ntime)




###Data for  Spring Creek Wyoming. 

spring17<- read.csv("/Users/bob.hall/Dropbox/met_error_models/metabolism process error/spring_june17.csv")
oxy<-spring17
lat<- 41.304
long<- -105.5676 #(note - for us folks in western hemisphere)
altitude<- 7200/3.28 # m asl
bpstd <- 1013 #Do not change. ignoring weather relate changes in bp for this minimal code.  1013 mb = 760 mm Hg.  
K600_guess <- 30 #your best guess of K600 in units of 1/d.  Big slow rivers, 0.5 to 2.  Medium slow rivers, 2-10. 
#Fun rivers to float 10-20.  Non bubbly streams 20-40.  Bubbly streams 50-200.
mean_depth <- 0.2 # knowing mean river depth (m) is not important for code testing, but it sure is important if you want to publish meaningful estimates of metabolis,m
oxy$time<-as_datetime(oxy$unixtime) #converts unit time to a R time object in UTC
oxy$solar.time<-convert_UTC_to_solartime(oxy$time, longitude= long, time.type="mean solar")
oxy$light<- calc_light(oxy$solar.time, latitude=lat, longitude=long, max.PAR =2326, attach.units = F)
oxy$osat<- calc_DO_sat(temp.water=oxy$temp, pressure.air=bpcalc_atm(bpst=bpstd, alt=altitude))

oxy_s<-oxy[oxy$time > ymd_hms("2017-06-01 11:01:00 ") & oxy$time < ymd_hms("2017-06-30 11:01:00 "), ]
oxy_s$Kc<-Kcor(oxy_s$temp,1)
nday<-length(oxy_g$oxy)/ntime
oxy_s$day<- rep(1:nday, each=ntime)

head(oxy)

plot(oxy_s$time, oxy_s$oxy)

