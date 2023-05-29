data {
  int<lower = 1> N; // number of observations
  int<lower = 1> P; // number of observed indicators
  int<lower = 1> Q; // number of latent variables
  vector[P] y[N]; // observed responses
  int<lower = 0> nC; // number of cross loadings
  real<lower=0> scale_global; // scale for the half-t prior for tau
  real<lower=1> nu_global; // df for the half-t prior for tau
  real<lower=1> nu_local; // df for the half-t priors for lambdas
  real<lower=0> slab_scale; // slab scale for the regularized horseshoe
  real<lower=0> slab_df; // slab df for the regularized horseshoe 
}
parameters {
  corr_matrix[Q] phi; // latent variable correlations (variances fixed to 1 for identification) 
  vector<lower = 0>[P] psi; // residual variances (no correlations assumed for now) 
  vector[P] L_main; 
  // auxiliary parameters/hyperparameters prior
  vector[nC] z;
  real<lower=0> aux1_global;
  real<lower=0> aux2_global;
  vector<lower=0>[nC] aux1_local;
  vector<lower=0>[nC] aux2_local;
  real<lower=0> caux;
}
transformed parameters {
  vector[P] mu; // means measurement model
  matrix[P, P] sigma; // covariance matrix measurement model 
  matrix[P, Q] L; // loading matrix
  vector[nC] L_cross; 
  real<lower=0> tau; // global shrinkage parameter
  vector<lower=0>[nC] lambda; // local shrinkage parameter
  vector<lower=0>[nC] lambda_tilde; // 'truncated' local shrinkage parameter
  real<lower=0> c; // slab scale
  mu = rep_vector(0, P); // no mean structure modelled (yet)
  lambda = aux1_local .* sqrt(aux2_local);
  tau = aux1_global * sqrt(aux2_global) * scale_global; 
  c = slab_scale * sqrt(caux);
  lambda_tilde = sqrt(c^2 * square(lambda) ./ (c^2 + tau^2 * square(lambda)) );
  L_cross = z .* lambda_tilde*tau;
  // create loading matrix: needs automating
  L[1:18, 1] = L_main[1:18];
  L[19:24, 2] = L_main[19:24];
  L[25:36, 3] = L_main[25:36];
  L[19:36, 1] = L_cross[1:18];
  L[1:18, 2] = L_cross[19:36];
  L[25:36, 2] = L_cross[37:48];
  L[1:24, 3] = L_cross[49:72];
  sigma = L * phi * L' + diag_matrix(psi); 
}
model {
  // prior cross-loadings: regularized horseshoe TODO: add scaling to error variance
  // half-t priors for lambdas and tau; inverse-gamma for c^2
  z ~ normal(0, 1);
  aux1_local ~ normal(0, 1);
  aux2_local ~ inv_gamma(0.5*nu_local, 0.5*nu_local);
  aux1_global ~ normal(0, 1);
  aux2_global ~ inv_gamma(0.5*nu_global, 0.5*nu_global);
  caux ~ inv_gamma(0.5*slab_df, 0.5*slab_df);
  // priors nuisance parameters
  target += normal_lpdf(L_main | 0, 5);
  target += cauchy_lpdf(psi | 0, 5);
  // automatic uniform prior latent variable correlations
  // likelihood
  target += multi_normal_lpdf(y | mu, sigma);
}
generated quantities {
  // avoid sign switching by multiplying all loadings of a factor with -1 if a marker item is negative TODO: automate
  // based on this post: https://discourse.mc-stan.org/t/non-convergence-of-latent-variable-model/12450/13
  vector[P] L_main_C;
  vector[nC] L_cross_C;
  corr_matrix[Q] phi_C;
  L_main_C = L_main;
  L_cross_C = L_cross;
  phi_C = phi;
  if(L_main[1] < 0){ // marker item factor 1
    L_main_C[1:18] = -1*L_main[1:18];
    L_cross_C[1:18] = -1*L_cross[1:18];
    if(L_main[19] > 0){
      // the corresponding factor correlation should also be multiplied by -1
      // but only if the loadings on one of the corresponding factors are multiplied by -1
      // otherwise the negatives cancel out
      phi_C[1, 2] = -1*phi[1, 2];
      phi_C[2, 1] = -1*phi[2, 1];
    }
    if(L_main[25] > 0){
      phi_C[1, 3] = -1*phi[1, 3];
      phi_C[3, 1] = -1*phi[3, 1];
    }
  }
  if(L_main[19] < 0){ // marker item factor 2
    L_main_C[19:24] = -1*L_main[19:24];
    L_cross_C[19:48] = -1*L_cross[19:48];
    if(L_main[1] > 0){
      phi_C[1, 2] = -1*phi[1, 2];
      phi_C[2, 1] = -1*phi[2, 1];
    }
    if(L_main[25] > 0){
      phi_C[3, 2] = -1*phi[3, 2];
      phi_C[2, 3] = -1*phi[2, 3];
    }
  }
  if(L_main[25] < 0){ // marker item factor 3
    L_main_C[25:36] = -1*L_main[25:36];
    L_cross_C[49:72] = -1*L_cross[49:72];
    if(L_main[1] > 0){
      phi_C[1, 3] = -1*phi[1, 3];
      phi_C[3, 1] = -1*phi[3, 1];
    }
    if(L_main[19] > 0){
      phi_C[3, 2] = -1*phi[3, 2];
      phi_C[2, 3] = -1*phi[2, 3];
    }
  }
} 
