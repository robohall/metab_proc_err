



## plot residuals
mu_summary <- summary(procfit_l, pars = c("mu"), probs = c(0.5))$summary
pred_mu<- mu_summary[,1]
ypred_summary <- summary(procfit_l, pars = c("ypred"), probs = c(0.5))$summary
eta_summary <-  summary(procfit_l, pars = c("eta"), probs = c(0.5))$summary

pred_proc<-data.frame(oxy=oxy_g$oxy, day=oxy_g$day, time=rep(1:ntime, nday), ypred=ypred_summary[,1], mu=mu_summary[,1], 
                      eta=eta_summary[,1])

ggplot(pred_proc, aes(x = time, y =eta)) + geom_point(col="red", size =0.6) +facet_wrap(~day, ncol=10)+
  theme_classic()+
  ylab("Residuals (mg DO/L)")
 
