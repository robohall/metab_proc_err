#AR process error models for Gallatin River


library(rstan)
options(mc.cores = parallel::detectCores())

library(streamMetabolizer)
library(tidyverse)
library(lubridate)

Kcor<-function (temp,K600) {
  K600/(600/(1800.6-(temp*120.1)+(3.7818*temp^2)-(0.047608*temp^3)))^-0.5
}
bpcalc_atm<- function(bpst, alt) {
  bpst*exp((-9.80665*0.0289644*alt)/(8.31447*(273.15+15)))
}

gall<- read.csv("path/file")
lat<- 45.36
long<-  -111.2
gall$time<-as_datetime(gall$unixtime) #converts unit time to a R time object in UTC

gall$solar.time<-convert_UTC_to_solartime(gall$time, longitude= long, time.type="mean solar")
gall$light<- calc_light(gall$solar.time, latitude=lat, longitude=long, max.PAR =2326, attach.units = F)
gall$Kc <- Kcor(gall$temp,1)


  
  ######Site 1

gall1<- gall[gall$site_id==1,]
gall1<- gall1[gall1$solar.time > mdy_hms() &  gall1$solar.time < mdy_hms() ,  ] #start near 4 AM first day.  End one time step earlier on last day
gall1$day <- #need code here to put integer on each day but starting at the 4 am time slot and going to 3:50 the following day
gall1$osat<- calc_DO_sat(temp.water=gall1$temp, pressure.air=bpcalc_atm(bpst=1013, alt=altitude))


             
gall1_metab_list<- list(T=length(gall1$solar.time), D=length(gall1$solar.time)/144, 
                        y=gall1$oxy, light=gall1$light,oxysat=gall1$osat,day=gall1$day,
                                    Kc=oxy_g$Kc,sumlight=colSums(matrix(gall1$light, nrow=ntime)), z=1, ts=1/144   )
 

gall1arfit<- stan(file= "metab_pi_ar_Q_2.stan", data=gall1_metab_list,  iter = 2000, chains = 4)
print(gall1arfit, pars=c("GPP","ER","K", "Kmean", "Ksd", "sigproc", "phi" ), digits=3 ) 


gpp_est_gall1<- summary(gall1arfit, pars = c("GPP"), probs = 0.5)$summary
ER_est_gall1<- summary(gall1arfit, pars = c("ER"), probs = 0.5)$summary
K_est_gall1<- summary(gall1arfit, pars = c("K"), probs = 0.5)$summary


plot(K_est_gall1[,1], ER_est_gall1[,1])

par(mfrow=c(2,1), mai=c(0.65,0.7,0.1,0.1), mgp=c(2,1,0))
plot(K_est_gall1[,1],  gpp_est_gall1[,1], pch=16, ylim=c(0,8), xlim=c(20,50), xlab="K (1/d)", ylab = (expression(paste('GPP (g O'[2] ~m^-2, ~d^-1,')')))   )
plot(K_est_gall1[,1],  ER_est_gall1[,1], pch=16, ylim=c(-15,-4), xlim=c(20,50), xlab="K (1/d)", ylab = (expression(paste('ER (g O'[2] ~m^-2, ~d^-1,')')))   )



# mu_summary_procar <- summary(procarfit, pars = c("mu"), probs = c(0.5))$summary
pred_mu_procar<- mu_summary_procar[,1]

ypred_summary_procar <- summary(procarfit, pars = c("ypred"), probs = c(0.5))$summary
pred_ypred_procar<- ypred_summary_procar[,1]

eta_summary_procar <- summary(procarfit, pars = c("eta"), probs = c(0.5))$summary
pred_eta_procar<- eta_summary_procar[,1]

plot(oxy_g$solar.time, oxy_g$oxy, pch=16, cex=0.5, col="orange")
lines(oxy_g$solar.time, pred_ypred_procar)


acf(pred_eta_procar)
acf(pred_eta_procar[1:144])

sqrt(mean((oxy_g$oxy - pred_ypred_procar)^2))




