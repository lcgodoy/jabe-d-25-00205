// https://discourse.mc-stan.org/t/numerical-stability-of-gps-with-negative-binomial-likelihood/19343/
int neg_binomial_2_log_safe_rng(real eta, real phi) {
  real gamma_rate = gamma_rng(phi, phi / exp(eta));
  if (gamma_rate > exp(20.7)) gamma_rate = exp(20.7);
  return poisson_rng(gamma_rate);
}
