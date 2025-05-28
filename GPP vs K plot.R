

### Fake data for a figure showing effects of how  K and GPP affect metabolism.  
#Was in my SFS talk, not consequential for the analyses.

set.seed(40)

nday=1
logK_fake<-rnorm(nday, 3.9,0.05)
K_fake<- exp(logK_fake)
P_fake<- rnorm(nday, 5, 0.06)
R_fake<- -0.8*P_fake + rnorm(nday, -5, 0.08)
sigma_oxy<- 0.027
ts<- 1/ntime


oxy_g_pi_mat_eta<- matrix(nrow=ntime, ncol=nday, rnorm(ntime*nday,0,sigma_oxy))
oxy_gHiK_HiP<- matrix(nrow=ntime, ncol=nday, NA)

for (d in 1:nday)	{
  
  oxy_gHiK_HiP[1,d]=gall_metab_list$y[1,d]
  #eta [1,d] = y[1,d]-mu[1,d] ;
  
  
  for (i in 2:ntime){
    
    
    oxy_gHiK_HiP[i,d]=  oxy_gHiK_HiP[i-1,d] + ((P_fake[d]/gall_metab_list$z)*(gall_metab_list$light[i,d]/gall_metab_list$sumlight[d]))+ 
      (R_fake[d]*ts/gall_metab_list$z)+ 
      (gall_metab_list$Kc[i,d]*K_fake[d] * ts*(gall_metab_list$oxysat[i,d]-oxy_gHiK_HiP[i-1,d] )) + oxy_g_pi_mat_eta[i,d]
    
    
  }
}	




#
logK_fake<-rnorm(nday, 3.9,0.05)
K_fake<- exp(logK_fake)
P_fake<- rnorm(nday, 0.5, 0.06)

oxy_g_pi_mat_eta<- matrix(nrow=ntime, ncol=nday, rnorm(ntime*nday,0,sigma_oxy))
oxy_gHiK_LoP<- matrix(nrow=ntime, ncol=nday, NA)

for (d in 1:nday)	{
  
  oxy_gHiK_LoP[1,d]=gall_metab_list$y[1,d]
  #eta [1,d] = y[1,d]-mu[1,d] ;
  
  
  for (i in 2:ntime){
    
    
    oxy_gHiK_LoP[i,d]=  oxy_gHiK_LoP[i-1,d] + ((P_fake[d]/gall_metab_list$z)*(gall_metab_list$light[i,d]/gall_metab_list$sumlight[d]))+ 
      (R_fake[d]*ts/gall_metab_list$z)+ 
      (gall_metab_list$Kc[i,d]*K_fake[d] * ts*(gall_metab_list$oxysat[i,d]-oxy_gHiK_LoP[i-1,d] )) + oxy_g_pi_mat_eta[i,d]
    
    
  }
}	




logK_fake<-rnorm(nday, 2.0,0.05)
K_fake<- exp(logK_fake)
P_fake<- rnorm(nday, 5, 0.06)

oxy_g_pi_mat_eta<- matrix(nrow=ntime, ncol=nday, rnorm(ntime*nday,0,sigma_oxy))
oxy_gLoK_HiP<- matrix(nrow=ntime, ncol=nday, NA)

for (d in 1:nday)	{
  
  oxy_gLoK_HiP[1,d]=6
  #eta [1,d] = y[1,d]-mu[1,d] ;
  
  
  for (i in 2:ntime){
    
    
    oxy_gLoK_HiP[i,d]=  oxy_gLoK_HiP[i-1,d] + ((P_fake[d]/gall_metab_list$z)*(gall_metab_list$light[i,d]/gall_metab_list$sumlight[d]))+ 
      (R_fake[d]*ts/gall_metab_list$z)+ 
      (gall_metab_list$Kc[i,d]*K_fake[d] * ts*(gall_metab_list$oxysat[i,d]-oxy_gLoK_HiP[i-1,d] )) + oxy_g_pi_mat_eta[i,d]
    
    
  }
}	




logK_fake<-rnorm(nday, 2.0,0.05)
K_fake<- exp(logK_fake)
P_fake<- rnorm(nday, 0.5, 0.06)


oxy_g_pi_mat_eta<- matrix(nrow=ntime, ncol=nday, rnorm(ntime*nday,0,sigma_oxy))
oxy_gLoK_LoP<- matrix(nrow=ntime, ncol=nday, NA)

for (d in 1:nday)	{
  
  oxy_gLoK_LoP[1,d]=6
  #eta [1,d] = y[1,d]-mu[1,d] ;
  
  
  for (i in 2:ntime){
    
    
    oxy_gLoK_LoP[i,d]=  oxy_gLoK_LoP[i-1,d] + ((P_fake[d]/gall_metab_list$z)*(gall_metab_list$light[i,d]/gall_metab_list$sumlight[d]))+ 
      (R_fake[d]*ts/gall_metab_list$z)+ 
      (gall_metab_list$Kc[i,d]*K_fake[d] * ts*(gall_metab_list$oxysat[i,d]-oxy_gLoK_LoP[i-1,d] )) + oxy_g_pi_mat_eta[i,d]
    
    
  }
}	

plot(100*oxy_gLoK_LoP/gall_metab_list$oxysat[,d], ylab="DO", pch=16, ylim=c(60,120))
plot(100*oxy_gLoK_HiP/gall_metab_list$oxysat[,d], ylab="DO", pch=16,ylim=c(60,120))
plot(100*oxy_gHiK_LoP/gall_metab_list$oxysat[,d], ylab="DO", pch=16,ylim=c(60,120))
plot(100*oxy_gHiK_HiP/gall_metab_list$oxysat[,d], ylab="DO", pch=16,ylim=c(60,120))





