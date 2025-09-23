/**
 * Log-Normal lpdf reparametrized (taken from drmr)
 *
 * Reference: https://doi.org/10.1111/jtsa.12242
 * 
 * @param x non-negative random variable
 * @param mu theoretical mean of X
 * @param phi weird parameter parameter.
 * 
 * @return a log-pdf
 */
real ln_mu_lpdf(real x, real mu, real phi) {
  real f_mu_phi;
  f_mu_phi = log1p(phi * inv_square(mu));
  real output = 0;
  output += log2() + log(pi()) + log(f_mu_phi) + log(x) +
    square(log(x) + log2() - 0.5 * f_mu_phi) * inv(f_mu_phi);
  return - 0.5 * output;
}
