functions {
#include utils/pois_rng_safe.stan
#include utils/pexp_corr.stan
#include utils/gp_lpdf.stan
}
data {
  int<lower = 1> N;
  int<lower = 1> k;
  int y[N];//observed counts
  vector[N] offset;//expected counts
  matrix<lower = 0>[N, N] dists;//expected counts
  matrix[N, k] X;
  real<lower=0, upper=1> prob_rho;
  real<lower=0> rho_0;
  real<lower = 0, upper = 1> nu;
}
transformed data {
  vector[N] logoff;
  logoff = log(offset);
  real par_exp = - log(prob_rho) / rho_0;
}
parameters {
  real<lower = 0> rho;
  real<lower = 0> sigma;
  vector[N] z; // Spatial Random Effect
  vector[k] beta;
  real beta_0;
}
transformed parameters {
}
model {
  vector[N] intercept;
  intercept = beta_0 + logoff + z;
  target += poisson_log_glm_lpmf(y | X, intercept, beta);
  target += gp_lpdf(z | dists, sigma, rho, nu);
  target += exponential_lpdf(rho | par_exp);
  target += normal_lpdf(beta | 0, 10);
  // target += normal_lpdf(beta_0 | 0, 10);
  target += student_t_lpdf(sigma | 3, 0, 1) -
    student_t_lccdf(0 | 3, 0, 1);
}
generated quantities {
  vector[N] log_lik;
  vector[N] y_rep;
  vector[N] SIR = exp(X * beta + z);
  vector[N] linpred = logoff + beta_0 + X * beta + z;
  for (i in 1:N) {
    if (linpred[i] > 20.7944) {
      y_rep[i] = -1;
    }
    else {
      y_rep[i] =
        poisson_log_rng(linpred[i]);
    }
    log_lik[i] =
      poisson_log_lpmf(y[i] | linpred[i]);
  }
}
