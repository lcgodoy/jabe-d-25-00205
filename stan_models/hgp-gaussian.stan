functions {
#include utils/pexp_corr.stan
#include utils/gp_lpdf.stan
  real pcp_prec_lpdf(real x, real alpha, real u) {
    real lambda = -log(alpha) / u;
    return log(lambda) - log(2) -
      lambda * exp(-0.5 * x) - 0.5 * x;
  }
}
data {
  int<lower=1> N;
  int<lower=1> k;
  int<lower=1> k2;
  matrix[N, N] dists;
  real<lower=0> rho_0;
  real<lower=0> tau_0;
  real<lower=0, upper=1> pt_u;
  real<lower=0, upper=1> prob_rho;
  vector[N] y;
  matrix[N, k] X; 
  matrix[N, k2] A;
  real<lower=0, upper=1> nu;
}
transformed data {
  real par_exp = - log(prob_rho) / rho_0;
}
parameters {
  real<lower=0> rho;
  real log_tau;
  vector[k] beta;
  vector[k2] alpha;
}
transformed parameters {
  // this notation might be confusing
  // tau is the sqrt of the nugget variance
  real tau = exp(log_tau);
}
model {
  target += gp_heter_xb_lpdf(y | dists, X, A, alpha,
                             rho, nu, tau, beta);
  target += normal_lpdf(beta | 0, 100);
  target += normal_lpdf(alpha | 0, 1);
  target += exponential_lpdf(rho | par_exp);
  target += pcp_prec_lpdf(log_tau | pt_u, tau_0);
}
generated quantities {
}
