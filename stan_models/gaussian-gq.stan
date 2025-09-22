functions{
#include utils/pexp_corr.stan
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
}
parameters {
  real rho;
  real tau;
  vector[k] beta;
  vector[k2] alpha;
}
generated quantities {// log likelihood & fitted
  array[N] real log_lik;
  array[N] real y_rep;
  {
      matrix[N, N] Sigma;
      matrix[N, N] R;
      matrix[N, N] chol_sig;
      vector[N] sigmas = exp(A * alpha);
      R = pexp_corr(dists, rho, nu);
      Sigma = add_diag(quad_form_diag(R, sigmas), square(tau));
      chol_sig =
        cholesky_decompose(Sigma);
      vector[N] z;
      z = y - X * beta;
      matrix[N, N] prec = chol2inv(chol_sig);
      vector[N] gi  = prec * z;
      for (i in 1:N) {
        real si  = 1 / prec[i, i];
        real mui = y[i] - gi[i] * si;
        log_lik[i] = normal_lpdf(y[i] | mui, si);
        y_rep[i] = normal_rng(mui, si);
      }
  }
}
