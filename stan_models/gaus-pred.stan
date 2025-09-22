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
  matrix[N, k] X;
  matrix[Np, k] Xp;
  matrix[N, k2] A;
  matrix[Np, k2] Ap;
  vector[N] y;
  real<lower=0, upper=1> nu;
}
parameters {
  vector[k] beta;
  vector[k2] alpha;
  real rho;
  real tau;
}
generated quantities {// log likelihood & fitted
  vector[Np] y_rep;
  {
    vector[N] sigmas = exp(A * alpha);
    vector[Np] sigmas_p = exp(Ap * alpha);
    matrix[N, N] Sigma =
      add_diag(quad_form_diag(pexp_corr(dists, rho, nu), sigmas),
               square(tau));
    matrix[Np, Np] Sigma_p =
      add_diag(quad_form_diag(pexp_corr(distsp, rho, nu), sigmas_p),
               square(tau));
    matrix[N, Np] Sigma_c =
      diag_post_multiply(diag_pre_multiply(sigmas,
                                           pexp_corr(cdists, rho, nu)),
                         sigmas_p);
    vector[N] mu = X * beta;
    vector[Np] mu_p = Xp * beta;
    y_rep = gp_pred_rng(y, // observed
                        Sigma,
                        Sigma_p,
                        Sigma_c,
                        mu, // mean observed
                        mu_p);
  }
}
