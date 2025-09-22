int pois_log_safe_rng(real eta) {
  int out = eta > 20.7944 ? 0 : poisson_log_rng(eta);
  return out;
}
