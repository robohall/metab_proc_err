

####Compare my homegrown Stan models with that from streamMetabolizer. Prep data first
##I am using observation and process erro models here (i.e., state space), because that is the usual
#approach with sM

###run data through streammetabolizer



model_name <- mm_name(type='bayes', pool_K600='normal', err_obs_iid = T, err_proc_iid =T)
model_specs <- specs(model_name, K600_daily_meanlog_meanlog=log(K600_guess), K600_daily_meanlog_sdlog=0.7, K600_daily_sdlog_sigma=0.05, burnin_steps=1000, 
                     saved_steps=1000) 

### Make your data sM frame

#Spring Creek
data_sm<-data.frame(DO.obs=oxy$oxy, DO.sat=oxy$osat, 
                    temp.water=oxy$temp, depth=rep(mean_depth,length(oxy$oxy)), 
                    light=oxy$light, solar.time=oxy$solar.time)

#Gallatin River
data_gall_sm<-data.frame(DO.obs=oxy_g$oxy, DO.sat=oxy_g$osat, 
                         temp.water=oxy_g$temp, depth=rep(mean_depth,length(oxy_g$oxy)), 
                         light=oxy_g$light, solar.time=oxy_g$solar.time)


head(data_gall_sm) #always take a look at the dataframe!


#sm_fit <- metab(model_specs, data=data_sm, info=c(site='spring creek', source='Hall'))
sm_fit <- metab(model_specs, data=data_gall_sm, info=c(site='Gallatin', source='Hall'))
#fit_pi<-fit

plot_DO_preds(predict_DO(sm_fit))
plot_metab_preds(predict_metab(sm_fit))
params<- get_params(sm_fit , uncertainty='ci')
plot(params$GPP.daily,params$ER.daily)

plot(params$K600.daily,params$ER.daily) #if these covary, then equifinality warning!  Indeed my data have that.

mcmc_oipi<-get_mcmc(sm_fit)
fit_check<-get_fit(sm_fit)$overall
get_fit(sm_fit)$overall %>%
  select(ends_with('Rhat'))

print(mcmc_oipi, pars=c("K600_daily_sdlog", "err_obs_iid_sigma", "err_proc_iid_sigma"))
##obs error does not converge, makes me think to skip using state space models, despite what
#Appling et al say.
rstan::stanmodel(mcmc_oipi)


mean(params$GPP.daily, na.rm=T)
mean(params$ER.daily, na.rm=T)

#######
#Homemade state space

gall_metab_list<- list(T=ntime*nday, D=nday, y=oxy_g$oxy, light=oxy_g$light,oxysat=oxy_g$osat,day=oxy_g$day,
                            Kc=oxy_g$Kc,sumlight=colSums(matrix(oxy_g$light, nrow=ntime)), z=mean_depth, ts=1/ntime   )




####  State space
#oipi_fit<- stan(file= "metab_oipi_2.stan", data=spring_metab_list,  iter = 2000, chains = 4)
oipi_fit<- stan(file= "metab_oipi_2.stan", data=gall_metab_list,  iter = 2000, chains = 4)


print(oipi_fit, pars=c("GPP","ER","K", "Ksd", "sigproc", "sigobs" ), digits=3 ) 
  ###

###process and obs error struggling to converge
gpp_est<- summary(oipi_fit, pars = c("GPP"), probs = 0.5)$summary
K_est<- summary(oipi_fit, pars = c("K"), probs = 0.5)$summary
er_est<- summary(oipi_fit, pars = c("ER"), probs = 0.5)$summary
  
plot(K_est[,1],er_est[,1])
plot(K_est[,1],er_est[,1], xlab="K", ylab="ER", ylim=c(-12,0), pch=16)
  
plot(params$GPP.daily[1:30], gpp_est[,1])
lines(params$GPP.daily[1:30],params$GPP.daily[1:30])

plot(params$K600.daily[1:30], K_est[,1])
lines(params$K600.daily[1:30], params$K600.daily[1:30])



###Try with process model only

###sM proc
model_name_proc <- mm_name(type='bayes', pool_K600='normal', err_obs_iid = F, err_proc_iid =T)
model_specs_proc <- specs(model_name_proc, K600_daily_meanlog_meanlog=log(K600_guess), K600_daily_meanlog_sdlog=0.7, K600_daily_sdlog_sigma=0.05, burnin_steps=1000, 
                     saved_steps=1000) 

rstan::stanmodel(sm_fit_proc)

sm_fit_proc <- metab(model_specs_proc, data=data_gall_sm, info=c(site='Gallatin', source='Hall'))


params_proc<- get_params(sm_fit_proc , uncertainty='ci')
plot(params_proc$GPP.daily,params_proc$ER.daily)

plot(params_proc$K600.daily,params_proc$ER.daily)  
###proc error model in sM are not working well at all, no idea what is going on there


##Bob's  proc
procfit<- stan(file= "metab_pi_2.stan", data=gall_metab_list_long,  iter = 2000, chains = 4)
print(procfit, pars=c("GPP","ER","K", "Ksd", "Kmean", "sigproc" ) , digits=3) 

gpp_est<- summary(procfit, pars = c("GPP"), probs = 0.5)$summary
K_est<- summary(procfit, pars = c("K"), probs = 0.5)$summary
er_est<- summary(procfit, pars = c("ER"), probs = 0.5)$summary


plot(params$GPP.daily[1:30], gpp_est[,1])
lines(params$GPP.daily[1:30],params$GPP.daily[1:30])

plot(params$K600.daily[1:30], K_est[,1])
lines(params$K600.daily[1:30], params$K600.daily[1:30])

#My model gets better pooiling than sM



