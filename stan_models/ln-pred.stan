functions{
#include utils/pexp_corr.stan
#include utils/gp_pred.stan
}
data{
  int<lower=1> N;
  int<lower=1> Np;
  int<lower=1> k;
  int<lower=1> k2;
  matrix[N, N] dists;
  matrix[Np, Np] distsp;
  matrix[N, Np] cdists;
  matrix[Np, k] Xp;
  matrix[N, k2] A;
  matrix[Np, k2] Ap;
  real<lower=0, upper=1> nu;
}
parameters {
  vector[k] beta;
  vector[k2] alpha;
  real rho;
  real tau;
  vector[N] z;
}
generated quantities {// log likelihood & fitted
  array[Np] real y_rep;
  {
    vector[N] sigmas = exp(A * alpha);
    vector[Np] sigmas_p = exp(Ap * alpha);
    matrix[N, N] Sigma =
      quad_form_diag(pexp_corr(dists, rho, nu), sigmas);
    matrix[Np, Np] Sigma_p =
      quad_form_diag(pexp_corr(distsp, rho, nu), sigmas_p);
    matrix[N, Np] Sigma_c =
      diag_post_multiply(diag_pre_multiply(sigmas,
                                           pexp_corr(cdists, rho, nu)),
                         sigmas_p);
    vector[Np] z_p;
    vector[N] mu = rep_vector(0.0, N);
    vector[Np] mu_p = rep_vector(0.0, Np);
    z_p = gp_pred_rng(z, // observed
                      Sigma,
                      Sigma_p,
                      Sigma_c,
                      mu, // mean observed
                      mu_p);
    vector[Np] lmu = exp(Xp * beta + z_p);
    for (i in 1:Np) {
      real mu_ln;
      real sigma_ln;
      sigma_ln = sqrt(log1p(tau * inv_square(lmu[i])));
      mu_ln = log(square(lmu[i]) * inv_sqrt(square(lmu[i]) + tau));
      y_rep[i] = lognormal_rng(mu_ln, sigma_ln);
    }
  }
}
