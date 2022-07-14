// JM: Multicenter + Longitudinal submodel + Event submodel (recurrent event)

// Association structure: slope/value

// Risk scale: calendar/gap

// Assumption: random intercept-slope model, weibull baseline hazard

// For: Simulation purpose

// Grace C. Zhou

// Credits to Copyright (C) 2016, 2017 Sam Brilleman; 2015, 2016, 2017 Trustees of Columbia University

/********************************************************/



functions {

  vector evaluate_eta(matrix X, array[] vector Z_u, int Dev_index,

                      array[] int U_id, array[] int C_id, real gamma,

                      vector beta, vector bVec, matrix uMat) {

    int N = rows(X); // num rows in design matrix

    int K = rows(beta); // num predictors

    vector[N] eta;

    if (K > 0) {

      eta = X * beta + gamma * Dev_index;

    } else {

      eta = rep_vector(0.0, N) + gamma * Dev_index;

    }

    for (n in 1 : N) {

      eta[n] = eta[n] + bVec[C_id[n]] * Dev_index

               + uMat[U_id[n], 1] * Z_u[1, n] + uMat[U_id[n], 2] * Z_u[2, n];

    }

    return eta;

  }

  
  /** 

  * Get the indices corresponding to the lower tri of a square matrix

  * @param dim The number of rows in the square matrix

  * @return A vector of indices

  */

  array[] int lower_tri_indices(int dim) {

    array[dim + choose(dim, 2)] int indices;

    int mark = 1;

    for (r in 1 : dim) {

      for (c in r : dim) {

        indices[mark] = (r - 1) * dim + c;

        mark = mark + 1;

      }

    }

    return indices;

  }

}

data {

  //----- Longitudinal submodels

  int<lower=0> y_N; // num observations
  
  int<lower=0> y_K; // num predictors

  vector[y_N] y1; // response vectors

  matrix[y_N, y_K] y1_X; // fix effect design matrix

  vector[y_K] y1_Xbar; // predictor means

  int<lower=0> b_N; // num center

  int<lower=0> b_K; // num center parameter

  int<lower=0> u_N; // num patient

  int<lower=0> u_K; // num patient parameter

  array[y_N] int<lower=0> y1_C_id; // center ids

  array[y_N] int<lower=0> y1_U_id; // patients id 

  array[u_K] vector[y_N] y1_Z; // random effect design matrix


  //----- Event submodel

  // data for calculating event submodel linear predictor in GK quadrature

  // NB these design matrices are evaluated AT the event time and

  // the (unstandardised) quadrature points

  int<lower=0> e_K; // num predictors

  int<lower=0> a_K; // num association parameter

  int<lower=0> Npat; // num patient

  int<lower=0> Nevents; // num events (ie. not censored)

  int<lower=0> qnodes; // num nodes for GK quadrature

  int<lower=0> Nobs_times_qnodes; // obs times nodes

  int<lower=0> nrow_e_Xq; // num rows in predictor matrix

  matrix[nrow_e_Xq, e_K] e_Xq; // predictor matrix
  
  vector[e_K] e_Xbar; // predictor means

  real norm_const; // normal constant

  int<lower=0> basehaz_df; // degree of freedom

  vector[nrow_e_Xq] basehaz_X; // design matrix (basis terms)

  vector[Nobs_times_qnodes] qwts; // GK quadrature weights with (b-a)/2 scaling

  matrix[nrow_e_Xq, y_K] y1_Xq; // fix effect design matrix at quadpoints

  array[u_K] vector[nrow_e_Xq] y1_Zq; // random effect design matrix at quadpoints

  array[nrow_e_Xq] int<lower=0> y1_Cq_id; // center index

  array[nrow_e_Xq] int<lower=0> y1_Uq_id; // patient index

  
  //----- Hyperparameters for prior distributions

  // scale for priors

  vector<lower=0>[y_K] y1_prior_scale;

  vector<lower=0>[e_K] e_prior_scale;

  vector<lower=0>[a_K] a_prior_scale;

  real<lower=0> y_prior_scale_for_intercept;

  real<lower=0> e_prior_scale_for_intercept;

  real<lower=0> y_prior_scale_for_aux;

  vector<lower=0>[basehaz_df] e_prior_scale_for_aux;

  
  // lkj prior stuff

  real<lower=0> b_prior_scale;

  vector<lower=0>[u_K] u_prior_scale;

  vector<lower=0>[u_K] u_prior_df;

  real<lower=0> u_prior_regularization;

  int<lower=0> Dev_index; //derivative indicator

}

transformed data {

  // indexing used to extract lower tri of RE covariance matrix

  array[u_K + choose(u_K, 2)] int u_cov_idx;

  if (u_K > 0) {

    u_cov_idx = lower_tri_indices(u_K);

  }

}

parameters {

  real y1_gamma; // intercept in longitudinal submodel

  vector[y_K] y1_z_beta; // primitive coefs in longitudinal submodel

  real<lower=0> y1_aux_unscaled; // unscaled residual error 

  real<lower=0> b_sd; // center sd  

  vector[b_N] z_b_vec; // unscaled center parameter

  vector<lower=0>[u_K] u_sd; // patient sd 

  matrix[u_K, u_N] z_u_mat; // unscaled patient parameter  

  cholesky_factor_corr[u_K > 1 ? u_K : 0] u_cholesky; // cholesky factor of correlation matrix

  real e_gamma; // intercept in event submodel

  vector[e_K] e_z_beta; // primitive coefs in event submodel (log hazard ratios)

  vector<lower=0>[basehaz_df] e_aux_unscaled; // unscaled coefs for baseline hazard

  vector[a_K] a_z_beta; // primitive association parameter 

  real<lower=0> v_sd; // frailty sd

  vector[Npat] z_v_vec; // unscaled frailty parameter

}

transformed parameters {

  vector[y_K] y1_beta; // coefs in longitudinal submodel

  real<lower=0> y1_aux; // residual error 

  matrix[u_N, u_K] u_mat; // patient parameter

  vector[b_N] b_vec; // center parameter

  vector[e_K] e_beta; // coefs in event submodel (log hazard ratios)

  vector[a_K] a_beta; // association parameter in event submodel

  vector<lower=0>[basehaz_df] e_aux; // coefs for baseline hazard

  vector[Npat] v_vec; //frailty term

  y1_beta = y1_z_beta .* y1_prior_scale;

  y1_aux = y1_aux_unscaled * y_prior_scale_for_aux;

  b_vec = b_sd * z_b_vec;
  
  if (u_K == 1) {

    u_mat = (u_sd[1] * z_u_mat)';

  } else if (u_K > 1) {

    u_mat = (diag_pre_multiply(u_sd, u_cholesky) * z_u_mat)';

  }

  e_beta = e_z_beta .* e_prior_scale;

  a_beta = a_z_beta .* a_prior_scale;

  e_aux = e_aux_unscaled .* e_prior_scale_for_aux;

  v_vec = v_sd * z_v_vec;

}

model {

  //---- Log likelihood for longitudinal submodel
  
  vector[y_N] y1_eta;
  
  y1_eta = evaluate_eta(y1_X, y1_Z, 1, y1_U_id, y1_C_id, y1_gamma, y1_beta,

                        b_vec, u_mat);

  target += normal_lpdf(y1 | y1_eta, y1_aux); // increment the target function

  //----- Log likelihood for event submodel (Gauss-Kronrod quadrature)
  
  {
    vector[nrow_e_Xq] y1_eta_q; // linear predictor for longitudinal submodel
  
    vector[nrow_e_Xq] e_eta_q; // linear predictor for event submodel

    vector[nrow_e_Xq] log_basehaz; // log baseline hazard at event time and quadrature points

    vector[nrow_e_Xq] log_haz_q; // log hazard at event time and quadrature points

    vector[Nevents] log_haz_etimes; // log hazard at the event time only

    vector[Nobs_times_qnodes] log_haz_qtimes; // log hazard at the quadrature points only

    
    // Step 1: linear predictor for longitudinal submodel at event time and quadrature points
	
    y1_eta_q = evaluate_eta(y1_Xq, y1_Zq, Dev_index, y1_Uq_id, y1_Cq_id,

                          y1_gamma, y1_beta, b_vec, u_mat);
  

    for (n in 1 : nrow_e_Xq) {

    // Step 2: linear predictor for event submodel at event time and quadrature points
      
      e_eta_q[n] = e_Xq[n] * e_beta + y1_eta_q[n] * a_beta[y1_Cq_id[n]]

                   + v_vec[y1_Uq_id[n]];

      // Step 3: log baseline hazard (Weibull) at event time and quadrature points

      log_basehaz[n] = e_gamma + norm_const + log(e_aux[y1_Cq_id[n]])

                       + basehaz_X[n] * (e_aux[y1_Cq_id[n]] - 1);

    }

    
    // Step 4: log hazard at event time and quadrature points

    log_haz_q = log_basehaz + e_eta_q;

    
    // Step 5: log hazard at event times only

    log_haz_etimes = head(log_haz_q, Nevents);

    
    // Step 6: log hazard at quadrature points only

    log_haz_qtimes = tail(log_haz_q, Nobs_times_qnodes);

    
    // Step 7: log likelihood for event submodel

    target += sum(log_haz_etimes) - dot_product(qwts, exp(log_haz_qtimes));

  }

  
  //----- Log priors

  // intercept for longitudinal submodel

  target += normal_lpdf(y1_gamma | 0, y_prior_scale_for_intercept);

  // coefs for longitudinal submodel   

  target += normal_lpdf(y1_z_beta | 0, 1);

  // residual error

  target += normal_lpdf(y1_aux_unscaled | 0, 1);

  // standard deviations

  target += normal_lpdf(z_b_vec | 0, 1);

  target += normal_lpdf(b_sd | 0, b_prior_scale); //prior for the intercept following Gelman 2008

  target += student_t_lpdf(u_sd | u_prior_df, 0, u_prior_scale); 

  target += normal_lpdf(to_vector(z_u_mat) | 0, 1);

  // corr matrix

  if (u_K > 1) {

    target += lkj_corr_cholesky_lpdf(u_cholesky | u_prior_regularization);

  }

  
  // intercept for event submodel

  target += normal_lpdf(e_gamma | 0, e_prior_scale_for_intercept);

  // coefs for event submodel

  target += normal_lpdf(e_z_beta | 0, 1);

  target += normal_lpdf(a_z_beta | 0, 1);

  // coefs for baseline hazard

  target += normal_lpdf(e_aux_unscaled | 0, 1);

  //frailty term

  target += normal_lpdf(z_v_vec | 0, 1);

  target += normal_lpdf(v_sd | 0, 10);

}

generated quantities {

  real y1_alpha; // transformed intercepts for longitudinal submodel

  vector[size(u_cov_idx)] u_cov; // transformed variance-covariance matrix for patient

  real rho; // correlation parameter

  real e_alpha; // transformed intercept for event submodel

  y1_alpha = y1_gamma - dot_product(y1_Xbar, y1_beta);

  if (u_K == 1) {

    u_cov[1] = u_sd[1] * u_sd[1];

  } else {

    u_cov = to_vector(quad_form_diag(multiply_lower_tri_self_transpose(

                                     u_cholesky), u_sd))[u_cov_idx];

  }

  rho = u_cov[2] / (u_sd[1] * u_sd[2]);
 
  e_alpha = e_gamma + norm_const - dot_product(e_Xbar, e_beta);

}




