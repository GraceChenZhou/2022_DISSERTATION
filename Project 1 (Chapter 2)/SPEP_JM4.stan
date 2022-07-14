//Fitted: Structure 4: HDSr with spep link + expCor
//Author: Grace C. Zhou
//Date: Apr.7, 2021

functions {

  //Define exponential power CDF expression
  
  real  spep_cdf(real x, real pow_r){
  
    if(pow_r <= 1){
		
		  return(pow(double_exponential_cdf(x/pow_r, 0, 1), pow_r)); //return probability 
		
	} else { 
			return(1-pow(double_exponential_cdf(-pow_r * x, 0, 1), (1/pow_r))); //return probability
		}	
	}
}

data {
  
  int<lower=1> Nobs; //Total # of obs. 
  int<lower=1> N; //Total # of patients
  
  int<lower=1> start_pos[N+1]; //starting point index of each patient
  int<lower=1> T[N]; //number of obs. per patient
  vector[Nobs] t; //time since first PE
   
  int<lower=1> NpredsX; //# of predictors for LME
  int<lower=1> NpredsV; // # of predictors for GLM
  row_vector[NpredsX] X[Nobs]; // design matrix for LME
  row_vector[NpredsV] V[Nobs]; // design matrix for GLM
    
  int<lower=1> Nlev1; // # of centers (level1)
  int<lower=1,upper=Nlev1> levind1[Nobs]; //id index for centers (level1)
  int<lower=1,upper=N> levind2[Nobs];//id index for patients (level2)
  vector[Nobs] y; //Continuous outcome
  int r[Nobs]; //Binary outcome
  real<lower=0> sdscal; //Overall residual
   
}

parameters {
	
  vector[NpredsX] alpha; //fixed coefficients for LME
	vector[NpredsV] beta; //fixed coefficients for GLM
	real<lower=0> sigmalev1; // scale parameter for center
	real<lower=0> sigmalev2; // scale parameter for patient
	real<lower=0> sigmaeps; // scale parameter for epsilon
	real<lower=0> tau; // scale parameter for exponential correlation
	real<lower=0> rho; // 1/range for exponential correlation
	vector[Nlev1] eta1; //latent for site
	vector[N] eta2; //latent for patient
	vector[Nobs] eta; //latent for exponential correlation
	real<lower = -1, upper = 1> rho1; // association parameter for center
	real<lower = -1, upper = 1> rho2; // association parameter for patient
	real<lower=0> pow_r[Nlev1]; //power parameter for flexible link 
	
	
}


transformed parameters {
		
	vector[Nobs] w;
	vector[Nobs] yhat;
	vector[Nlev1] ran1;
	vector[N] ran2;
	vector[Nobs] Fsp;	
	
	for (n in 1:N){
   
		{
			matrix[T[n],T[n]] Sigma;
			matrix[T[n],T[n]] L_Sigma;
			vector[T[n]] t_sub;
			
  			t_sub=t[start_pos[n]:(start_pos[n+1]-1)];
			
			//off-diagonal elements
			for(i in 1:(T[n]-1)){
					for (j in (i+1):T[n]){
						Sigma[i,j] = pow(tau,2) * exp(- rho * fabs(t_sub[i] - t_sub[j]));
						Sigma[j,i] = Sigma[i,j];
					}
			}

			// diagonal elements
			for (k in 1:T[n]){
				Sigma[k,k] = pow(tau,2)+0.000001; // + jitter
			}
				
			L_Sigma=cholesky_decompose(Sigma);
			
			w[start_pos[n]:(start_pos[n+1]-1)]= L_Sigma  * eta[start_pos[n]:(start_pos[n+1]-1)];
		}
  }
  
		ran1 = sigmalev1 * eta1;
		ran2 = sigmalev2 * eta2;
	
		for(i in 1:Nobs){
	
			yhat[i] = X[i] * alpha + ran1[levind1[i]]  + ran2[levind2[i]] + w[i];
			Fsp[i]=spep_cdf(V[i] * beta + rho1 * ran1[levind1[i]] + rho2 * ran2[levind2[i]], pow_r[levind1[i]]);			
		
		}
}


model {

  //priors
  alpha ~ normal(0, 100); 
  beta ~ normal(0, 100);
  eta1 ~ normal(0,1);
  eta2 ~ normal(0,1);
  eta ~ normal(0,1);
  sigmalev1 ~ cauchy(0, 5);
  sigmalev2 ~ cauchy(0, 5);
  sigmaeps ~ cauchy(0, 2.5*sdscal);
  rho1 ~ uniform(-1, 1);
  rho2 ~ uniform(-1, 1);
  pow_r ~ exponential(1);
  tau ~ normal(0, 5);
  rho ~ inv_gamma(2,1);
  
  //likelihood
  
  y ~ normal(yhat,sigmaeps);
  r ~ bernoulli(Fsp);
  
  
}

generated quantities{

  vector[Nobs] log_lik_y; 
  vector[Nobs] log_lik_r; 
  
  for (i in 1:Nobs){
      
	  log_lik_y[i] = normal_lpdf(y[i] | yhat[i], sigmaeps);
	  log_lik_r[i] = bernoulli_lpmf(r[i] | Fsp[i]);
	
  }
  
}


