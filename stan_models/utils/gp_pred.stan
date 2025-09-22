/* Steps to use: build sigma, sigma_p, sigma_c, mu and mu_p beforehand! */
vector gp_pred_rng(vector z, // observed
                   matrix sigma,
                   matrix sigma_p,
                   matrix sigma_c,
                   vector mu, // mean observed
                   vector mu_p) { // mean pred
  int N = rows(sigma);
  int Np = rows(sigma_p);
  matrix[N, N] L = cholesky_decompose(sigma);
  matrix[N, Np] v_pred;
  matrix[Np, Np] sigma_cond;
  vector[Np] mu_cond;
  vector[N] z_c = z - mu;
  v_pred = mdivide_left_tri_low(L, sigma_c);
  mu_cond = mu_p + mdivide_right_tri_low(v_pred', L) * z_c;
  sigma_cond = sigma_p - crossprod(v_pred);
  return multi_normal_rng(mu_cond, sigma_cond);
}
