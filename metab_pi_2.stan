
data {
 
 int <lower=0> T;
int <lower=0> D;

  vector [T] y; // oxygen
  vector [T] oxysat; 
  vector [T] light;
  real z;
  vector [T] Kc;
  int <lower=0> day[T];
  vector [D] sumlight;
  real ts;
  
}



parameters {
  
  vector [D] GPP ; 
  vector [D] ER ;
  vector [D] logK ;

  //real<lower=0> sigobs;
real <lower=0> sigproc;
//matrix [T,D]  eta;
real Kmean;
real <lower=0> Ksd;
}


transformed parameters{
  vector [D] K ;
  vector [T]  mu;
  	 vector [T] eta;	
  	
  K=exp(logK);
  



	mu[1]=y[1];
	 eta [1] = y[1]-mu[1] ;
	
	
    for (i in 2:T){
      
      if (day[i] != day[i-1]){
  mu[i]= y[i];
} else { 
     
         mu[i]=  y[i-1] + ((GPP[day[i]]/z)*(light[i]/sumlight[day[i]]))+ (ER[day[i]]*ts/z)+ (Kc[i]*K[day[i]] * ts*(oxysat[i]-y[i-1] ));
}
         eta [i] = y[i]-mu[i] ;
        }
        
    
	
}
model {

	
	
		GPP ~ normal(0, 10);
		ER ~ normal(-2, 8);
		logK ~ normal(Kmean, Ksd);
		
		Kmean ~ normal(log(30),0.7);
		Ksd ~ normal (0,0.05);
  	//sigobs~ normal(0,0.1);
    sigproc ~normal(0,1);
	

		 	   
  //likelihood  
             
        for ( i in 2:T){
          eta[i] ~ normal (0,sigproc);
        }

}


 generated quantities{
   vector [T]  ypred;
   


	ypred[1]=y[1];
	
	
    for (i in 2:T){
      
        if (day[i]>day[i-1]){
        ypred[i]= y[i];
        } else {
          
     ypred[i]=  ypred[i-1] + ((GPP[day[i]]/z)*(light[i]/sumlight[day[i]]))+ (ER[day[i]]*ts/z)+ (Kc[i]*K[day[i]] * ts*(oxysat[i]-ypred[i-1] ));
    }
    }
  
}


