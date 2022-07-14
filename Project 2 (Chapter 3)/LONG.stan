// JM: Multicenter + Longitudinal submodel + Event submodel (recurrent time-to-event)

// Two-stage method: Stage One-longitudinal submodel

// Grace C. Zhou

// Created based on JM.stan
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

  int<lower=0> y_N; // num observations

  int<lower=0> y_K; // num predictors

  vector[y_N] y1; // response vectors

  matrix[y_N, y_K] y1_X; // fix effect design matrix

  vector[y_K] y1_Xbar; // predictor means

  int<lower=0> b_N; // num center

  int<lower=0> b_K; // total num params

  int<lower=0> u_N; // num patients

  int<lower=0> u_K; // total num params

  array[y_N] int<lower=0> y1_C_id; // center ids 

  array[y_N] int<lower=0> y1_U_id; // patients ids 

  array[u_K] vector[y_N] y1_Z; // random effect design matrix

  // scale for priors

  vector<lower=0>[y_K] y1_prior_scale;

  real<lower=0> y_prior_scale_for_intercept;

  real<lower=0> y_prior_scale_for_aux;

  // lkj prior stuff

  real<lower=0> b_prior_scale;

  vector<lower=0>[u_K] u_prior_scale;

  vector<lower=0>[u_K] u_prior_df;

  real<lower=0> u_prior_regularization;

  // for stage two

  int<lower=0> nrow_e_Xq; // num. rows in long. predictor matrix at quadpoints

  matrix[nrow_e_Xq, y_K] y1_Xq; // fix effect design matrix at quadpoints

  array[u_K] vector[nrow_e_Xq] y1_Zq; // random effect design matrix 

  int<lower=0> Dev_index;

  array[nrow_e_Xq] int<lower=0> y1_Cq_id; // center index

  array[nrow_e_Xq] int<lower=0> y1_Uq_id; // patient index

}


transformed data {

  array[u_K + choose(u_K, 2)] int u_cov_idx;

  if (u_K > 0) {

    u_cov_idx = lower_tri_indices(u_K);

  }

}

parameters {

  real y1_gamma; // intercept in long. submodel

  vector[y_K] y1_z_beta; // primitive coefs in long. submodels

  real<lower=0> y1_aux_unscaled; // unscaled residual error SDs 

  real<lower=0> b_sd; // center level sds  

  vector[b_N] z_b_vec; // unscaled patient level params

  matrix[u_K, u_N] z_u_mat; // unscaled patient level params  

  vector<lower=0>[u_K] u_sd; // patient level sds  

  cholesky_factor_corr[u_K > 1 ? u_K : 0] u_cholesky; // cholesky factor of corr matrix

}

transformed parameters {

  vector[y_K] y1_beta; // population level params for long. submodels

  real<lower=0> y1_aux; // residual error SDs for long. submodels

  matrix[u_N, u_K] u_mat; // patient level params

  vector[b_N] b_vec; //center level params

  
  // coefs for long. submodels

  y1_beta = y1_z_beta .* y1_prior_scale;

  
  // residual error SDs for long. submodels

  y1_aux = y1_aux_unscaled * y_prior_scale_for_aux;


  // center level params

  b_vec = b_sd * z_b_vec;

  // patient level params
  
  if (u_K == 1) {

    u_mat = (u_sd[1] * z_u_mat)';

  } else if (u_K > 1) {

    u_mat = (diag_pre_multiply(u_sd, u_cholesky) * z_u_mat)';

  }

}

model {

  //----- Log likelihood 

  vector[y_N] y1_eta;
    
  y1_eta = evaluate_eta(y1_X, y1_Z, 1, y1_U_id, y1_C_id, y1_gamma, y1_beta,
                        b_vec, u_mat);
  
  // increment the target with the log-lik

  target += normal_lpdf(y1 | y1_eta, y1_aux);

  //----- Log-priors
  
  target += normal_lpdf(y1_gamma | 0, y_prior_scale_for_intercept);

  target += normal_lpdf(y1_z_beta | 0, 1);
  
  target += normal_lpdf(y1_aux_unscaled | 0, 1);

  target += normal_lpdf(z_b_vec | 0, 1);

  target += normal_lpdf(b_sd | 0, b_prior_scale);

  target += student_t_lpdf(u_sd | u_prior_df, 0, u_prior_scale);

  target += normal_lpdf(to_vector(z_u_mat) | 0, 1);

  // corr matrix

  if (u_K > 1) {

    target += lkj_corr_cholesky_lpdf(u_cholesky | u_prior_regularization);

  }

}


generated quantities {

  real y1_alpha; // transformed intercepts for long. submodels

  vector[size(u_cov_idx)] u_cov; // var-cov for random effect

  real rho;

  vector[nrow_e_Xq] y1_eta_q; // linear predictor for stage two


  // Transformed intercepts for long. submodels

  y1_alpha = y1_gamma - dot_product(y1_Xbar, y1_beta);

  // Transform variance-covariance matrix for REs

  if (u_K == 1) {

    u_cov[1] = u_sd[1] * u_sd[1];

  } else {

    u_cov = to_vector(quad_form_diag(multiply_lower_tri_self_transpose(

                                     u_cholesky), u_sd))[u_cov_idx];

  }

  
  rho = u_cov[2] / (u_sd[1] * u_sd[2]);

  //For stage two

  y1_eta_q = evaluate_eta(y1_Xq, y1_Zq, Dev_index, y1_Uq_id, y1_Cq_id,

                          y1_gamma, y1_beta, b_vec, u_mat);

  
}




