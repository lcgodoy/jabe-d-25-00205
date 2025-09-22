functions{
#include utils/pexp_corr.stan
#include utils/lpdfs.stan
}
data{
  int<lower=1> N;
  int<lower=1> k;
  int<lower=1> k2;
  matrix[N, N] dists;
  vector[N] y;
  matrix[N, k] X; 
  matrix[N, k2] A;
  real<lower=0, upper=1> nu;
  int<lower = 0, upper = 1> ln;
}
parameters {
  real rho;
  real tau;
  vector[N] z;
  vector[k] beta;
  // real beta;
  vector[k2] alpha;
}
generated quantities {// log likelihood & fitted
  array[N] real log_lik;
  array[N] real y_rep;
  {
    vector[N] lmu = exp(X * beta + z);
    for (i in 1:N) {
      if (ln) {
        real mu_ln;
        real sigma_ln;
        sigma_ln = sqrt(log1p(tau * inv_square(lmu[i])));
        mu_ln = log(square(lmu[i]) * inv_sqrt(square(lmu[i]) + tau));
        y_rep[i] = lognormal_rng(mu_ln, sigma_ln);
        log_lik[i] =
          ln_mu_lpdf(y[i] | lmu[i], tau);
      } else {
        real gamma_beta;
        gamma_beta = tau / lmu[i];
        y_rep[i] = gamma_rng(tau, gamma_beta);
        log_lik[i] =
          gamma_mu_lpdf(y[i] | lmu[i], tau);
      }
    }
  }
}
