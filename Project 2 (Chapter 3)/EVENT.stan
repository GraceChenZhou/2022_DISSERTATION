// JM: Multicenter + Longitudinal submodel + Event submodel (recurrent time-to-event)

// Two-stage method: Stage Two-Event submodel

// Grace C. Zhou

// Created based on JM.stan
/********************************************************/


data {

  int<lower=0> e_K; // num. of predictors in event submodel

  int<lower=0> Npat; // num. individuals

  int<lower=0> Nevents; // num. events (ie. not censored)

  int<lower=0> qnodes; // num. of nodes for GK quadrature

  int<lower=0> Nobs_times_qnodes;

  int<lower=0> nrow_e_Xq; // num. rows in event submodel predictor matrix

  matrix[nrow_e_Xq, e_K] e_Xq; // predictor matrix (event submodel)

  vector[e_K] e_Xbar; // predictor means (event submodel)

  real norm_const;

  int<lower=0> basehaz_df; // df for baseline hazard

  vector[nrow_e_Xq] basehaz_X; // design matrix (basis terms) for baseline hazard

  vector[Nobs_times_qnodes] qwts; // GK quadrature weights with (b-a)/2 scaling

  array[nrow_e_Xq] int<lower=0> y1_Cq_id; // center index 

  array[nrow_e_Xq] int<lower=0> y1_Uq_id; // patient index

  vector<lower=0>[e_K] e_prior_scale;

  real<lower=0> e_prior_scale_for_intercept;

  vector<lower=0>[basehaz_df] e_prior_scale_for_aux;
  
  //assoc

  int<lower=0> a_K; // num. of association parameters

  vector<lower=0>[a_K] a_prior_scale;

  vector[nrow_e_Xq] y1_eta_q;

  //frailty

  int<lower=0> u_K;

  vector<lower=0>[u_K] u_prior_scale;

}

parameters {

  real e_gamma; // intercept in event submodel

  vector[e_K] e_z_beta; // primitive coefs in event submodel (log hazard ratios)

  vector<lower=0>[basehaz_df] e_aux_unscaled; // unscaled coefs for baseline hazard

  vector[a_K] a_z_beta; // primitive assoc params (log hazard ratios)

  real<lower=0> v_sd; // frailty  

  vector[Npat] z_v_vec; // unscaled frailty

}

transformed parameters {

  vector[e_K] e_beta; // coefs in event submodel (log hazard ratios)

  vector[a_K] a_beta; // assoc params in event submodel (log hazard ratios)

  vector<lower=0>[basehaz_df] e_aux; // coefs for baseline hazard

  vector[Npat] v_vec; //individual level params for frailty(log(v)~gamma)
 
  e_beta = e_z_beta .* e_prior_scale;

  a_beta = a_z_beta .* a_prior_scale;

  e_aux = e_aux_unscaled .* e_prior_scale_for_aux;

  v_vec = v_sd * z_v_vec;
  
}

model {

  //----- Log likehood (Gauss-Kronrod quadrature)

  vector[nrow_e_Xq] e_eta_q; // linear predictor for both submodels
  
  vector[nrow_e_Xq] log_basehaz; // log baseline hazard AT event time and quadrature points

  vector[nrow_e_Xq] log_haz_q; // log hazard AT event time and quadrature points

  vector[Nevents] log_haz_etimes; // log hazard AT the event time only

  vector[Nobs_times_qnodes] log_haz_qtimes; // log hazard AT the quadrature points

  for (n in 1 : nrow_e_Xq) {
    
    e_eta_q[n] = norm_const + e_gamma + e_Xq[n] * e_beta + y1_eta_q[n] * a_beta[y1_Cq_id[n]]

                   + v_vec[y1_Uq_id[n]];
                   
    log_basehaz[n] = log(e_aux[y1_Cq_id[n]]) + basehaz_X[n] * (e_aux[y1_Cq_id[n]] - 1);

  }

   log_haz_q = log_basehaz + e_eta_q;
    
   log_haz_etimes = head(log_haz_q, Nevents);

   log_haz_qtimes = tail(log_haz_q, Nobs_times_qnodes);

   target += sum(log_haz_etimes) - dot_product(qwts, exp(log_haz_qtimes));

  

  

  //----- Log-priors

  target += normal_lpdf(e_gamma | 0, e_prior_scale_for_intercept);

  target += normal_lpdf(e_z_beta | 0, 1);

  target += normal_lpdf(a_z_beta | 0, 1);

  target += normal_lpdf(e_aux_unscaled | 0, 1);
  
  target += normal_lpdf(z_v_vec | 0, 1);

  target += normal_lpdf(v_sd | 0, u_prior_scale[1]);

}

generated quantities {

  real e_alpha; // transformed intercept for event submodel

  // Transformed intercept for event submodel 

  e_alpha = e_gamma + norm_const - dot_product(e_Xbar, e_beta);

}




