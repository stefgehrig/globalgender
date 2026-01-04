// fitting model 1 (p. 147) from:
// 'Modeling fertility in modern populations' (Peristera and Kostaki, 2007)
// to data on number of children by age
data {
  int<lower=1> N; // number of data points
  int<lower=1> K; // number of clusters (countries)
  int<lower=1> J; // number of predictor indices
  // data points (counts aggregated within country, age, grb)
  array[N] int<lower=1> nobs;
  array[N] int<lower=0> ys;
  array[N] int<lower=18, upper=49> ages;
  array[N] int<lower=1, upper=J> pred_idcs;
  // predictor values (country and grb for each index)
  array[J] real pred_grb;
  array[J] int<lower=1, upper=K> pred_cty;
}

parameters {
  // scaled country effects (one for each asfr curve model parameter)
  array[K] vector[8] z_b;
  // population average parameters slope
  vector[8] beta;
  // Cholesky factor for the correlation matrix
  cholesky_factor_corr[8] L_Omega;
  // SD among countries of each parameter
  vector<lower=0>[8] L_std;

}

transformed parameters{
  array[K] vector[8] b;
  for(k in 1:K){
      b[k] = beta + diag_pre_multiply(L_std, L_Omega) * z_b[k];
  }
  
  // parameters describing the asfr curves
  array[J,50] real<lower=0> lambda; // 50 lambda values (0 to 49, although we only have data for 18 to 49)
  array[J] real<lower=0> c1, mu, sigma1, sigma2;

  // constants for scaling and centering on log-scale, such that
  // implies indirectly priors on population average paremeters informed by prior literature
  real log_c1_center = log(0.125); // EQT 95%: 0.04691389 0.33305703; exp(log(0.125) + c(-1.96, 1.96) * 0.5)
  real log_mu_center = log(25); // EQT 95%: 16.89260 36.99844; exp(log(25) + c(-1.96, 1.96) * 0.2)
  real log_sigma_center = log(5); // EQT 95%: 1.149627 21.746176; exp(log(5) + c(-1.96, 1.96) * 0.75)
  real log_c1_scale = 0.5;
  real log_mu_scale = 0.2;  
  real log_sigma_scale = 0.75; 

  // obtain the demographic parameters from underlying components on log scale
  for(j in 1:J){
    int k = pred_cty[j];
    if(pred_grb[j] == 0) {
      c1[j] = exp(b[k][1] * log_c1_scale + log_c1_center);
      mu[j] = exp(b[k][2] * log_mu_scale + log_mu_center);
      sigma1[j] = exp(b[k][3] * log_sigma_scale + log_sigma_center);
      sigma2[j] = exp(b[k][4] * log_sigma_scale + log_sigma_center);
    } else {
      c1[j] = exp(b[k][5] * log_c1_scale + log_c1_center);
      mu[j] = exp(b[k][6] * log_mu_scale + log_mu_center);
      sigma1[j] = exp(b[k][7] * log_sigma_scale + log_sigma_center);
      sigma2[j] = exp(b[k][8] * log_sigma_scale + log_sigma_center);
    }
    
    real cum_sum = 0;
    for(m in 0:49){
      // asfr value for each age added to previous
      cum_sum += c1[j] * exp(-pow((m - mu[j]) / ((m <= mu[j]) ? sigma1[j] : sigma2[j]), 2)); // ternary conditional operator
      lambda[j, m+1] = cum_sum;  // shift index by 1
    }
  }
}

model {
  // likelihood
  for(i in 1:N){
target += poisson_lpmf(ys[i] | lambda[pred_idcs[i], ages[i] + 1] * nobs[i]); // shift index by 1 (e.g., index 19 = age 18 etc.)
  }
  // prior on population average paremeters on log-scale
  beta ~ std_normal();
  // prior on cholesky decomposition of correlation matrix
  L_Omega ~ lkj_corr_cholesky(1); // uniform over correlation matrices
  // prior on standard deviations for each parameter
  L_std ~ std_normal(); // for reference category (half-normal, mean: sqrt(2 / pi) = 0.80)
  // prior on multivariate normal, before scaling and centering
  for(k in 1:K){
     z_b[k] ~ std_normal();
  }
  
}

generated quantities {
  // correlation matrix
  matrix[8,8] Omega;
  Omega = multiply_lower_tri_self_transpose(L_Omega);
}
