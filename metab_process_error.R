

library(rstan)
options(mc.cores = parallel::detectCores())

library(streamMetabolizer)
library(tidyverse)
library(lubridate)




##gallatin, uses "oxy_g" from data script

gall_metab_list<- list(T=ntime*nday, D=nday, y=oxy_g$oxy, light=oxy_g$light,oxysat=oxy_g$osat,day=oxy_g$day,
                       Kc=oxy_g$Kc,sumlight=colSums(matrix(oxy_g$light, nrow=ntime)), z=mean_depth, ts=1/ntime   )




########
####process error models.


procfit<- stan(file= "metab_pi_2.stan", data=gall_metab_list_long,  iter = 2000, chains = 4)
print(procfit, pars=c("GPP","ER","K", "Ksd", "Kmean", "sigproc" ) , digits=3) 

K_est_procfit<- summary(procfit, pars = c("K"), probs = 0.5)$summary
gpp_est_procfit<- summary(procfit, pars = c("GPP"), probs = 0.5)$summary
er_est_procfit<- summary(procfit, pars = c("ER"), probs = 0.5)$summary




###summaries for procfit
mu_summary <- summary(procfit, pars = c("mu"), probs = c(0.5))$summary
pred_mu<- mu_summary[,1]


ypred_summary <- summary(procfit, pars = c("ypred"), probs = c(0.5))$summary
pred_ypred<- ypred_summary[,1]

eta_summary <- summary(procfit, pars = c("eta"), probs = c(0.5))$summary
pred_eta<- eta_summary[,1]

acf(pred_eta)
acf(pred_eta[1:144])

plot(oxy_g$solar.time, oxy_g$oxy, pch=16, cex=0.5, col="orange")
lines(oxy_g$solar.time, pred_ypred)

plot(diff(oxy_g$oxy), pch=16, cex=0.5, col="orange")
lines(diff(pred_mu))
lines(diff(pred_ypred))


###
#process error plots
par(mfrow=c(2,1), mai=c(0.65,0.7,0.1,0.1), mgp=c(2,1,0))
plot(K_est[,1], gpp_est[,1], pch=16, ylim=c(0,8), xlim=c(20,50), xlab="K (1/d)", ylab = (expression(paste('GPP (g O'[2] ~m^-2, ~d^-1,')')))   )
plot(K_est[,1], er_est[,1], pch=16, ylim=c(-15,-4), xlim=c(20,50), xlab="K (1/d)", ylab = (expression(paste('ER (g O'[2] ~m^-2, ~d^-1,')')))   )

sqrt(mean((oxy_g$oxy - pred_ypred)^2))



##process error models working great, stick with them




###########



##Process error with AR term
 
#procarfit<- stan(file= "metab_pi_ar.stan", data=spring_metab_list,  iter = 2000, chains = 4)
procarfit<- stan(file= "metab_pi_ar_2.stan", data=gall_metab_list,  iter = 2000, chains = 4)
 print(procarfit, pars=c("GPP","ER","K", "Kmean", "Ksd", "sigproc", "phi" ), digits=3 ) 


gpp_est_procar<- summary(procarfit, pars = c("GPP"), probs = 0.5)$summary
ER_est_procar<- summary(procarfit, pars = c("ER"), probs = 0.5)$summary
K_est_procar<- summary(procarfit, pars = c("K"), probs = 0.5)$summary


plot(gpp_est[,1], gpp_est_procar[,1])
lines(gpp_est[,1],gpp_est[,1])

plot(K_est[,1], K_est_procar[,1])
lines(K_est[,1],K_est[,1])

plot(K_est_procar[,1], ER_est_procar[,1])

par(mfrow=c(2,1), mai=c(0.65,0.7,0.1,0.1), mgp=c(2,1,0))
plot(K_est_procar[,1],  gpp_est_procar[,1], pch=16, ylim=c(0,8), xlim=c(20,50), xlab="K (1/d)", ylab = (expression(paste('GPP (g O'[2] ~m^-2, ~d^-1,')')))   )
plot(K_est_procar[,1],  ER_est_procar[,1], pch=16, ylim=c(-15,-4), xlim=c(20,50), xlab="K (1/d)", ylab = (expression(paste('ER (g O'[2] ~m^-2, ~d^-1,')')))   )



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





####Light and AR process error models.


procar_light_fit<- stan(file= "metab_pi_ar_light_2.stan", data=gall_metab_list,  iter = 2000, chains = 4)

print(procar_light_fit, pars=c("GPP","ER","K", "Kmean", "Ksd","phi", "alpha", "sigproc" ) ) 


gpp_est_procar_light<- summary(procar_light_fit, pars = c("GPP"), probs = 0.5)$summary
ER_est_procar_light<- summary(procar_light_fit, pars = c("ER"), probs = 0.5)$summary
K_est_procar_light<- summary(procar_light_fit, pars = c("K"), probs = 0.5)$summary



plot(K_est_procar_light[,1], ER_est_procar_light[,1])


mu_summary_procar_light <- summary(procar_light_fit, pars = c("mu"), probs = c(0.5))$summary
pred_mu_procar_light<- mu_summary_procar_light[,1]

ypred_summary_procar_light <- summary(procar_light_fit, pars = c("ypred"), probs = c(0.5))$summary
pred_ypred_procar_light<- ypred_summary_procar_light[,1]

eta_summary_procar_light <- summary(procar_light_fit, pars = c("eta"), probs = c(0.5))$summary
pred_eta_proca_light<- eta_summary_procar_light[,1]


par(mfrow=c(2,1), mai=c(0.65,0.7,0.1,0.1), mgp=c(2,1,0))
plot(K_est_procar_light[,1],  gpp_est_procar_light[,1], pch=16, ylim=c(0,8), xlim=c(20,50), xlab="K (1/d)", ylab = (expression(paste('GPP (g O'[2] ~m^-2, ~d^-1,')')))   )
plot(K_est_procar_light[,1],  ER_est_procar_light[,1], pch=16, ylim=c(-15,-4), xlim=c(20,50), xlab="K (1/d)", ylab = (expression(paste('ER (g O'[2] ~m^-2, ~d^-1,')')))   )

sqrt(mean((oxy_g$oxy - pred_ypred_procar_light)^2))



####Spring Creek

spring_metab_list<- list(T=ntime*nday, D=nday, y=oxy_s$oxy, light=oxy_s$light,oxysat=oxy_s$osat,day=oxy_s$day,
                       Kc=oxy_s$Kc,sumlight=colSums(matrix(oxy_s$light, nrow=ntime)), z=spring_mean_depth, ts=1/ntime   )



procfit<- stan(file= "metab_pi_2.stan", data=spring_metab_list,  iter = 2000, chains = 4)
 print(procfit, pars=c("GPP","ER","K", "Ksd", "Kmean", "sigproc" ) , digits=3) 

K_est_procfit<- summary(procfit, pars = c("K"), probs = 0.5)$summary
gpp_est_procfit<- summary(procfit, pars = c("GPP"), probs = 0.5)$summary
er_est_procfit<- summary(procfit, pars = c("ER"), probs = 0.5)$summary




###summaries for procfit
mu_summary <- summary(procfit, pars = c("mu"), probs = c(0.5))$summary
pred_mu<- mu_summary[,1]


ypred_summary <- summary(procfit, pars = c("ypred"), probs = c(0.5))$summary
pred_ypred<- ypred_summary[,1]

eta_summary <- summary(procfit, pars = c("eta"), probs = c(0.5))$summary
pred_eta<- eta_summary[,1]

acf(pred_eta)
acf(pred_eta[1:144])

plot(oxy_g$solar.time, oxy_g$oxy, pch=16, cex=0.5, col="orange")
lines(oxy_g$solar.time, pred_ypred)

plot(diff(oxy_g$oxy), pch=16, cex=0.5, col="orange")
lines(diff(pred_mu))
lines(diff(pred_ypred))


###
#process error plots
par(mfrow=c(2,1), mai=c(0.65,0.7,0.1,0.1), mgp=c(2,1,0))
plot(K_est[,1], gpp_est[,1], pch=16, ylim=c(0,8), xlim=c(20,50), xlab="K (1/d)", ylab = (expression(paste('GPP (g O'[2] ~m^-2, ~d^-1,')')))   )
plot(K_est[,1], er_est[,1], pch=16, ylim=c(-15,-4), xlim=c(20,50), xlab="K (1/d)", ylab = (expression(paste('ER (g O'[2] ~m^-2, ~d^-1,')')))   )

sqrt(mean((oxy_s$oxy - pred_ypred)^2))


##Process error with AR term

procarfit<- stan(file= "metab_pi_ar_2.stan", data=spring_metab_list,  iter = 2000, chains = 4)

print(procarfit, pars=c("GPP","ER","K", "Kmean", "Ksd", "sigproc", "phi" ), digits=3 ) 


gpp_est_procar<- summary(procarfit, pars = c("GPP"), probs = 0.5)$summary
ER_est_procar<- summary(procarfit, pars = c("ER"), probs = 0.5)$summary
K_est_procar<- summary(procarfit, pars = c("K"), probs = 0.5)$summary


plot(gpp_est[,1], gpp_est_procar[,1])
lines(gpp_est[,1],gpp_est[,1])

plot(K_est[,1], K_est_procar[,1])
lines(K_est[,1],K_est[,1])

plot(K_est_procar[,1], ER_est_procar[,1])

par(mfrow=c(2,1), mai=c(0.65,0.7,0.1,0.1), mgp=c(2,1,0))
plot(K_est_procar[,1],  gpp_est_procar[,1], pch=16, ylim=c(0,8), xlim=c(20,50), xlab="K (1/d)", ylab = (expression(paste('GPP (g O'[2] ~m^-2, ~d^-1,')')))   )
plot(K_est_procar[,1],  ER_est_procar[,1], pch=16, ylim=c(-15,-4), xlim=c(20,50), xlab="K (1/d)", ylab = (expression(paste('ER (g O'[2] ~m^-2, ~d^-1,')')))   )



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






