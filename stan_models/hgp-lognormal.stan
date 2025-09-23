functions {
#include utils/pexp_corr.stan
#include utils/gp_lpdf.stan
#include utils/lpdfs.stan
}
data {
  int<lower = 1> N;
  int<lower = 1> k;
  int<lower = 1> k2;
  vector[N] y;
  matrix<lower = 0>[N, N] dists;
  matrix[N, k] X;
  real<lower=0, upper=1> prob_rho;
  real<lower=0> rho_0;
  matrix[N, k2] A;
  real<lower = 0, upper = 1> nu;
}
transformed data {
  real par_exp = - log(prob_rho) / rho_0;
}
parameters {
  real<lower = 0> rho;
  real<lower = 0> tau;
  vector[N] z; // Spatial Random Effect
  vector[k] beta;
  vector[k2] alpha;
  // real beta;
}
transformed parameters {
}
model {
  {
    vector[N] lmu;
    lmu = exp(X * beta + z);
    for (i in 1:N)
      target += ln_mu_lpdf(y[i] | lmu[i], tau);
  }
  target += gp_heter_lpdf(z | dists, A, alpha, rho, nu, 0.0);
  target += exponential_lpdf(rho | par_exp);
  if (k > 1) {
    target += normal_lpdf(beta | 0, 10);
  }
  target += normal_lpdf(alpha | 0, 1);
  target += student_t_lpdf(tau | 3, 0, 1) -
    student_t_lccdf(0 | 3, 0, 1);
}
generated quantities {
}
