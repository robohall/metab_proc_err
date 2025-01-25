
data {
 
 int <lower=0> T;
int <lower=0> D;

  matrix [T,D] y; // oxygen
  matrix [T,D] oxysat; 
  matrix [T,D] light;
  real z;
  matrix [T,D] Kc;
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
  matrix [T,D]  mu;
  	
  	
  K=exp(logK);
  


for (d in 1:D)	{

	mu[1,d]=y[1,d];
	
	
    for (i in 2:T){
      
        //mu[i,d]= eta[i,d]+ y[i-1,d] + ((GPP[d]/z)*(light[i,d]/sumlight[d]))+ (ER[d]*ts/z)+ (Kc[i,d]*K[d] * ts*(oxysat[i,d]-y[i-1,d] ));
        
         mu[i,d]=  y[i-1,d] + ((GPP[d]/z)*(light[i,d]/sumlight[d]))+ (ER[d]*ts/z)+ (Kc[i,d]*K[d] * ts*(oxysat[i,d]-y[i-1,d] ));
       
        
	}
	}	

	
  
}
model {
	 matrix [T,D] eta;
	
	
		GPP ~ normal(0, 10);
		ER ~ normal(-2, 8);
		logK ~ normal(Kmean, Ksd);
		
		Kmean ~ normal(log(30),0.7);
		Ksd ~ normal (0,0.05);
  	//sigobs~ normal(0,0.1);
    sigproc ~normal(0,1);
	

		 	   
  //likelihood  
  
 
  for (d in 1:D)	{
  	
     for (i in 2:T)  {
     	
          //y[i,d] ~ normal(mu[i,d],sigproc);      
          eta [i,d] = y[i,d]-mu[i,d] ;
          eta[i,d] ~ normal (0,sigproc);
}
}
 
}


