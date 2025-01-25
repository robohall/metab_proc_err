

library(rstan)
options(mc.cores = parallel::detectCores())

library(streamMetabolizer)
library(tidyverse)
library(lubridate)

bpcalc_atm<- function(bpst, alt) {
  bpst*exp((-9.80665*0.0289644*alt)/(8.31447*(273.15+15)))
}
bpcalc_atm(bpst=bpstd, alt=altitude)


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

head(oxy)


plot(oxy$solar.time,oxy$oxy)

model_name <- mm_name(type='bayes', pool_K600='normal', err_obs_iid = T, err_proc_iid =T)
model_specs <- specs(model_name, K600_daily_meanlog_meanlog=log(K600_guess), K600_daily_meanlog_sdlog=0.7, K600_daily_sdlog_sigma=0.05, burnin_steps=1000, 
                     saved_steps=1000) 

### Make your data frame

data_sm<-data.frame(DO.obs=oxy$oxy, DO.sat=oxy$osat, 
                    temp.water=oxy$temp, depth=rep(mean_depth,length(oxy$oxy)), 
                    light=oxy$light, solar.time=oxy$solar.time)

head(data_sm) #always take a look at the dataframe!


fit <- metab(model_specs, data=data_sm, info=c(site='spring creek', source='Hall'))
#fit_oipi<-fit
#fit_pi<-fit

plot_DO_preds(predict_DO(fit))
plot_metab_preds(predict_metab(fit))
params<- get_params(fit , uncertainty='ci')

plot(params$K600.daily,params$ER.daily) #if these covary, then equifinality warning!  Indeed my data have that.

mcmc_oipi<-get_mcmc(fit_oipi)
fit_check<-get_fit(fit_oipi)$overall
get_fit(fit_oipi)$overall %>%
  select(ends_with('Rhat'))

print(mcmc_oipi, pars="K600_daily_sdlog")

rstan::stanmodel(mcmc_oipi)


mean(params$GPP.daily, na.rm=T)
mean(params$ER.daily, na.rm=T)


########
#homemade stan model
oxy_s<-oxy[oxy$time > ymd_hms("2017-06-01 11:01:00 ") & oxy$time < ymd_hms("2017-06-30 11:01:00 "), ]
Kcor<-function (temp,K600) {
  K600/(600/(1800.6-(temp*120.1)+(3.7818*temp^2)-(0.047608*temp^3)))^-0.5
}
oxy_s$Kc<-Kcor(oxy_s$temp,1)
  
  
ntime<-144
nday<-length(oxy_s$oxy)/ntime
spring_metab_list<- list(T=ntime, D=nday, y=matrix(oxy_s$oxy, nrow=ntime), light=matrix(oxy_s$light, nrow=ntime),oxysat=matrix(oxy_s$osat, nrow=ntime),
                         Kc=matrix(oxy_s$Kc, nrow=ntime),sumlight=colSums(matrix(oxy_s$light, nrow=ntime)), z=mean_depth, ts=1/ntime   )



####  State space
oipi_fit<- stan(file= "metab_oipi.stan", data=spring_metab_list,  iter = 2000, chains = 4)


print(oipi_fit, pars=c("GPP","ER","K", "Ksd", "sigproc", "sigobs" ) ) 
  ###
gpp_est<- summary(oipi_fit, pars = c("GPP"), probs = 0.5)$summary
K_est<- summary(oipi_fit, pars = c("K"), probs = 0.5)$summary

plot(params$GPP.daily[2:30], gpp_est[,1])
plot(params$K600.daily[2:30], K_est[,1])


mu_summary <- summary(oipi_fit, pars = c("mu"), probs = c(0.5))$summary
pred_mu<- mu_summary[,1]
oxymat<-matrix(oxy_s$oxy, nrow=ntime)

ypred_summary <- summary(oipi_fit, pars = c("ypred"), probs = c(0.5))$summary
pred_ypred<- ypred_summary[,1]


plot(pred_mu,pred_ypred)


####  State space with AR error on state
oipi_ar_fit<- stan(file= "metab_oipi_ar.stan", data=spring_metab_list,  iter = 2000, chains = 4)


print(oipi_ar_fit, pars=c("GPP","ER","K", "Ksd", "sigproc", "sigobs" ) ) 
###
gpp_est<- summary(oipi_ar_fit, pars = c("GPP"), probs = 0.5)$summary
K_est<- summary(oipi_ar_fit, pars = c("K"), probs = 0.5)$summary

plot(params$GPP.daily[2:30], gpp_est[,1])
plot(params$K600.daily[2:30], K_est[,1])






####process error only.  
procfit<- stan(file= "metab_pi.stan", data=spring_metab_list,  iter = 2000, chains = 4)
print(procfit, pars=c("GPP","ER","K", "Ksd", "Kmean", "sigproc" ) ) 


gpp_est<- summary(procfit, pars = c("GPP"), probs = 0.5)$summary
K_est<- summary(procfit, pars = c("K"), probs = 0.5)$summary

plot(params$GPP.daily[2:30], gpp_est[,1])
plot(params$K600.daily[2:30], K_est[,1])

mu_summary <- summary(procfit, pars = c("mu"), probs = c(0.5))$summary
pred_mu<- mu_summary[,1]
oxymat<-matrix(oxy_s$oxy, nrow=ntime)




###########


 
procarfit<- stan(file= "metab_pi_ar.stan", data=spring_metab_list,  iter = 2000, chains = 4)
print(procarfit, pars=c("GPP","ER","K", "Kmean", "Ksd", "sigproc", "phi" ) ) 


###ARCH tries

alpha_0<- 0.00036
alpha_1<-1.74
y<-0
for (i in 2:150){
  
  y[i]<- rnorm(1,0,(alpha_0+alpha_1*y[i-1]^2)^0.5)
  
}

plot(y, type='o')

procarchfit<- stan(file= "metab_pi_arch.stan", data=spring_metab_list,  iter = 2000, chains = 4)
print(procarchfit, pars=c("GPP","ER","K", "Kmean", "Ksd", "alpha_0", "alpha_1" ) ) 

