real pcp_prec_lpdf(real x, real alpha, real u) {
  real lambda = - log(alpha) / u;
  return log(lambda) - log(2) -
    lambda * exp(-0.5 * x) - 0.5 * x;
}
