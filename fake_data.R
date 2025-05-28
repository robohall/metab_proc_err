

###Gallatin fake data.  These derive from the param estimates, but I removed the covariance with K
set.seed(40)
logK_fake<-rnorm(nday, 3.4,0.05)
K_fake<- exp(logK_fake)
P_fake<- rnorm(nday, 3.82, 0.46)
R_fake<- -0.8*P_fake + rnorm(nday, -5, 0.78)
sigma_oxy<- 0.027
ts<- 1/ntime

gall_metab_list<- list(T=ntime*nday, D=nday, y=oxy_g$oxy, light=oxy_g$light,oxysat=oxy_g$osat,day=oxy_g$day,
                            Kc=oxy_g$Kc,sumlight=colSums(matrix(oxy_g$light, nrow=ntime)), z=mean_depth, ts=1/ntime   )


#######Fake proc error data
oxy_g_pi_eta<- rnorm(ntime*nday,0,sigma_oxy)



  oxy_g_pi=oxy_g$oxy[1]
  
  
  
  for (i in 2:length(oxy_g$oxy)){
    if (gall_metab_list$day[i] != gall_metab_list$day[i-1]){
      oxy_g_pi[i]= oxy_g$oxy[i]
    } else {
    
    oxy_g_pi[i]=  oxy_g_pi[i-1] + ((P_fake[gall_metab_list$day[i]]/gall_metab_list$z)*(gall_metab_list$light[i]/gall_metab_list$sumlight[gall_metab_list$day[i]]))+ 
      (R_fake[gall_metab_list$day[i]]*ts/gall_metab_list$z)+ 
      (gall_metab_list$Kc[i]*K_fake[gall_metab_list$day[i]] * ts*(gall_metab_list$oxysat[i]-oxy_g_pi[i-1] )) + oxy_g_pi_eta[i]
    
    }  
  
}	
#############





########## Fake proc with AR data

sigma_oxy<- 0.025
phi<- 0.5



oxy_g_pi_ar=oxy_g$oxy[1]
eta_ar<- 0


for (i in 2:length(oxy_g$oxy)){
  if (gall_metab_list$day[i] != gall_metab_list$day[i-1]){
    oxy_g_pi_ar[i]= oxy_g$oxy[i]
    eta_ar[i]<- 0
  } else {
    
    eta_ar[i]<- phi*eta_ar[i-1] + rnorm(1,0,sigma_oxy)
    
    oxy_g_pi_ar[i]=  oxy_g_pi_ar[i-1] + ((P_fake[gall_metab_list$day[i]]/gall_metab_list$z)*(gall_metab_list$light[i]/gall_metab_list$sumlight[gall_metab_list$day[i]]))+ 
      (R_fake[gall_metab_list$day[i]]*ts/gall_metab_list$z)+ 
      (gall_metab_list$Kc[i]*K_fake[gall_metab_list$day[i]] * ts*(gall_metab_list$oxysat[i]-oxy_g_pi_ar[i-1] )) + eta_ar[i]
    
  }  
  
}	

acf(diff(oxy_g_pi_ar[1:144]))
acf(eta_ar[1:144])


################
#Proc fit with proc error
#####
gall_metab_fake_proc_list<- list(T=ntime*nday, D=nday, y=  oxy_g_pi, light=oxy_g$light,oxysat=oxy_g$osat,day=oxy_g$day,
                                 Kc=oxy_g$Kc,sumlight=colSums(matrix(oxy_g$light, nrow=ntime)), z=mean_depth, ts=1/ntime   )

gall_fake_procdata_procmodel_fit<- stan(file= "metab_pi_2.stan", data=gall_metab_fake_proc_list,  iter = 2000, chains = 4)

print(gall_fake_procdata_procmodel_fit, pars=c("GPP","ER","K", "Ksd", "Kmean", "sigproc" ) , digits=3) 


gall_fake_procdata_procmodel_out<- list(
gpp_est = summary(gall_fake_procdata_procmodel_fit, pars = c("GPP"), probs = c(0.0275, 0.5, 0.975))$summary,
K_est= summary(gall_fake_procdata_procmodel_fit, pars = c("K"), probs = c(0.0275, 0.5, 0.975))$summary,
er_est= summary(gall_fake_procdata_procmodel_fit, pars = c("ER"), probs = c(0.0275, 0.5, 0.975))$summary
)



plot( gall_fake_procdata_procmodel_out$K_est[,1], gall_fake_procdata_procmodel_out$er_est[,1])


#######
#proc fit with procar data

gall_metab_fake_procar_list<- list(T=ntime*nday, D=nday, y=  oxy_g_pi_ar, light=oxy_g$light,oxysat=oxy_g$osat,day=oxy_g$day,
                                 Kc=oxy_g$Kc,sumlight=colSums(matrix(oxy_g$light, nrow=ntime)), z=mean_depth, ts=1/ntime   )

gall_fake_procardata_procmodel_fit<- stan(file= "metab_pi_2.stan", data=gall_metab_fake_procar_list,  iter = 2000, chains = 4)

print(gall_fake_procardata_procmodel_fit, pars=c("GPP","ER","K", "Ksd", "Kmean", "sigproc") , digits=3) 

gall_fake_procardata_procmodel_out<- list(
  gpp_est = summary(gall_fake_procardata_procmodel_fit, pars = c("GPP"), probs = c(0.0275, 0.5, 0.975))$summary,
  K_est= summary(gall_fake_procardata_procmodel_fit, pars = c("K"), probs = c(0.0275, 0.5, 0.975))$summary,
  er_est= summary(gall_fake_procardata_procmodel_fit, pars = c("ER"), probs = c(0.0275, 0.5, 0.975))$summary
)

###########
###procarfit with procar data
#######

gall_fake_procardata_procarmodel_fit<- stan(file= "metab_pi_ar_2.stan", data=gall_metab_fake_procar_list,  iter = 2000, chains = 4)

print(gall_fake_procardata_procarmodel_fit, pars=c("GPP","ER","K", "Ksd", "Kmean", "sigproc", "phi") , digits=3) 

gall_fake_procardata_procarmodel_out<- list(
  gpp_est = summary(gall_fake_procardata_procarmodel_fit, pars = c("GPP"), probs = c(0.0275, 0.5, 0.975))$summary,
  K_est= summary(gall_fake_procardata_procarmodel_fit, pars = c("K"), probs = c(0.0275, 0.5, 0.975))$summary,
  er_est= summary(gall_fake_procardata_procarmodel_fit, pars = c("ER"), probs = c(0.0275, 0.5, 0.975))$summary
)



###Plots

par(mfrow=c(3,1), mai=c(0.5,0.5,0.05,0.05), mgp=c(2,1,0))

plot (P_fake, gall_fake_procdata_procmodel_out$gpp_est[,1], pch=16, xlim=c(2,6), ylim=c(2,6), col="black", ylab=("estimated GPP"), xlab = "" )
lines(P_fake,P_fake)
arrows(x0=P_fake, y0=gall_fake_procdata_procmodel_out$gpp_est[,1], x1 = P_fake, y1 = gall_fake_procdata_procmodel_out$gpp_est[,4], length=0, col="black" )
arrows(x0=P_fake, y0=gall_fake_procdata_procmodel_out$gpp_est[,1], x1 = P_fake, y1 = gall_fake_procdata_procmodel_out$gpp_est[,6], length=0, col="black" )


plot (P_fake, gall_fake_procardata_procmodel_out$gpp_est[,1], pch=16, xlim=c(2,6), ylim=c(2,6), col="red" ,ylab=("estimated GPP"), xlab = "" )
lines(P_fake,P_fake)
arrows(x0=P_fake, y0=gall_fake_procardata_procmodel_out$gpp_est[,1], x1 = P_fake, y1 = gall_fake_procardata_procmodel_out$gpp_est[,4], length=0, col="red" )
arrows(x0=P_fake, y0=gall_fake_procardata_procmodel_out$gpp_est[,1], x1 = P_fake, y1 = gall_fake_procardata_procmodel_out$gpp_est[,6], length=0, col="red" )

plot (P_fake, gall_fake_procardata_procarmodel_out$gpp_est[,1], pch=16, xlim=c(2,6), ylim=c(2,6), col="blue", ylab=("estimated GPP"), xlab = "True (fake) GPP" )
lines(P_fake,P_fake)
arrows(x0=P_fake, y0=gall_fake_procardata_procarmodel_out$gpp_est[,1], x1 = P_fake, y1 = gall_fake_procardata_procarmodel_out$gpp_est[,4], length=0, col="blue")
arrows(x0=P_fake, y0=gall_fake_procardata_procarmodel_out$gpp_est[,1], x1 = P_fake, y1 = gall_fake_procardata_procarmodel_out$gpp_est[,6], length=0, col="blue" )


######
par(mfrow=c(3,3), mai=c(0.4,0.4,0.05,0.05), mgp=c(2,1,0))

plot (P_fake, gall_fake_procdata_procmodel_out$gpp_est[,1], pch=16, xlim=c(2,6), ylim=c(2,6), col="black", ylab=("estimated GPP"), xlab = "" )
lines(P_fake,P_fake)
arrows(x0=P_fake, y0=gall_fake_procdata_procmodel_out$gpp_est[,1], x1 = P_fake, y1 = gall_fake_procdata_procmodel_out$gpp_est[,4], length=0, col="black" )
arrows(x0=P_fake, y0=gall_fake_procdata_procmodel_out$gpp_est[,1], x1 = P_fake, y1 = gall_fake_procdata_procmodel_out$gpp_est[,6], length=0, col="black" )
text(3.5,6,"Data: proc err, Model: proc err", col="black")

plot (P_fake, gall_fake_procardata_procmodel_out$gpp_est[,1], pch=16, xlim=c(2,6), ylim=c(2,6), col="red" ,ylab=("estimated GPP"), xlab = "True (fake) GPP")
lines(P_fake,P_fake)
arrows(x0=P_fake, y0=gall_fake_procardata_procmodel_out$gpp_est[,1], x1 = P_fake, y1 = gall_fake_procardata_procmodel_out$gpp_est[,4], length=0, col="red" )
arrows(x0=P_fake, y0=gall_fake_procardata_procmodel_out$gpp_est[,1], x1 = P_fake, y1 = gall_fake_procardata_procmodel_out$gpp_est[,6], length=0, col="red" )
text(3.5,6,"Data: procAR err, Model: proc err", col="red")

plot (P_fake, gall_fake_procardata_procarmodel_out$gpp_est[,1], pch=16, xlim=c(2,6), ylim=c(2,6), col="blue", ylab=("estimated GPP"), xlab = ""  )
lines(P_fake,P_fake)
arrows(x0=P_fake, y0=gall_fake_procardata_procarmodel_out$gpp_est[,1], x1 = P_fake, y1 = gall_fake_procardata_procarmodel_out$gpp_est[,4], length=0, col="blue")
arrows(x0=P_fake, y0=gall_fake_procardata_procarmodel_out$gpp_est[,1], x1 = P_fake, y1 = gall_fake_procardata_procarmodel_out$gpp_est[,6], length=0, col="blue" )
text(3.5,6,"Data: procAR err, Model: procAR err", col="blue")

##ER
plot (R_fake, gall_fake_procdata_procmodel_out$er_est[,1], pch=16, xlim=c(-11,-4), ylim=c(-11,-4), col="black", ylab=("estimated ER"), xlab = "" )
lines(R_fake,R_fake)
arrows(x0=R_fake, y0=gall_fake_procdata_procmodel_out$er_est[,1], x1 = R_fake, y1 = gall_fake_procdata_procmodel_out$er_est[,4], length=0, col="black" )
arrows(x0=R_fake, y0=gall_fake_procdata_procmodel_out$er_est[,1], x1 = R_fake, y1 = gall_fake_procdata_procmodel_out$er_est[,6], length=0, col="black" )


plot (R_fake, gall_fake_procardata_procmodel_out$er_est[,1], pch=16, xlim=c(-11,-4), ylim=c(-11,-4), col="red" ,ylab="estimated ER",  xlab = "True (fake) ER" )
lines(R_fake,R_fake)
arrows(x0=R_fake, y0=gall_fake_procardata_procmodel_out$er_est[,1], x1 = R_fake, y1 = gall_fake_procardata_procmodel_out$er_est[,4], length=0, col="red" )
arrows(x0=R_fake, y0=gall_fake_procardata_procmodel_out$er_est[,1], x1 = R_fake, y1 = gall_fake_procardata_procmodel_out$er_est[,6], length=0, col="red" )

plot (R_fake, gall_fake_procardata_procarmodel_out$er_est[,1], pch=16, xlim=c(-11,-4), ylim=c(-11,-4), col="blue", ylab=("estimated ER"),xlab = "" )
lines(R_fake,R_fake)
arrows(x0=R_fake, y0=gall_fake_procardata_procarmodel_out$er_est[,1], x1 = R_fake, y1 = gall_fake_procardata_procarmodel_out$er_est[,4], length=0, col="blue")
arrows(x0=R_fake, y0=gall_fake_procardata_procarmodel_out$er_est[,1], x1 = R_fake, y1 = gall_fake_procardata_procarmodel_out$er_est[,6], length=0, col="blue" )

##K
plot (K_fake, gall_fake_procdata_procmodel_out$K_est[,1], pch=16, xlim=c(20,35), ylim=c(20,35), col="black", ylab=("estimated K"), xlab = "" )
lines(K_fake,K_fake)
arrows(x0=K_fake, y0=gall_fake_procdata_procmodel_out$K_est[,1], x1 = K_fake, y1 = gall_fake_procdata_procmodel_out$K_est[,4], length=0, col="black" )
arrows(x0=K_fake, y0=gall_fake_procdata_procmodel_out$K_est[,1], x1 = K_fake, y1 = gall_fake_procdata_procmodel_out$K_est[,6], length=0, col="black" )


plot (K_fake, gall_fake_procardata_procmodel_out$K_est[,1], pch=16,  xlim=c(20,35), ylim=c(20,35), col="red" ,ylab=("estimated K"), xlab = "True (fake) K" )
lines(K_fake,K_fake)
arrows(x0=K_fake, y0=gall_fake_procardata_procmodel_out$K_est[,1], x1 = K_fake, y1 = gall_fake_procardata_procmodel_out$K_est[,4], length=0, col="red" )
arrows(x0=K_fake, y0=gall_fake_procardata_procmodel_out$K_est[,1], x1 = K_fake, y1 = gall_fake_procardata_procmodel_out$K_est[,6], length=0, col="red" )

plot (K_fake, gall_fake_procardata_procarmodel_out$K_est[,1], pch=16,  xlim=c(20,35), ylim=c(20,35), col="blue", ylab=("estimated K"), xlab = ""  )
lines(K_fake,K_fake)
arrows(x0=K_fake, y0=gall_fake_procardata_procarmodel_out$K_est[,1], x1 = K_fake, y1 = gall_fake_procardata_procarmodel_out$K_est[,4], length=0, col="blue")
arrows(x0=K_fake, y0=gall_fake_procardata_procarmodel_out$K_est[,1], x1 = K_fake, y1 = gall_fake_procardata_procarmodel_out$K_est[,6], length=0, col="blue" )





 plot( gall_fake_procardata_procarmodel_out$K_est[,1], gall_fake_procardata_procarmodel_out$er_est[,1], pch=16, col="blue", ylim=c(-12,-5), xlim=c(20,33),
       xlab="Estimated K", ylab = "Estimated ER")
 points( gall_fake_procardata_procmodel_out$K_est[,1], gall_fake_procardata_procmodel_out$er_est[,1], pch=16, col="red")
 points( gall_fake_procdata_procmodel_out$K_est[,1], gall_fake_procdata_procmodel_out$er_est[,1], pch=16, col="black")




