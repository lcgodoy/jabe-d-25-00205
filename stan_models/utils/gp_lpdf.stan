real gp_lpdf(vector z, data matrix dist,
             real sigma, real rho, real nu) {
  int N;
  N = rows(dist);
  vector[N] mu;
  mu = rep_vector(0.0, N);
  matrix[N, N] Sigma  = square(sigma) *
    pexp_corr(dist, rho, nu);
  return multi_normal_lpdf(z | mu, Sigma);
}

real gp_heter_lpdf(vector z,
                   data matrix dist,
                   data matrix W,
                   vector alpha,
                   real rho, real nu,
                   real b0) {
  int N, K;
  N = rows(dist);
  K = cols(dist);
  vector[N] mu = rep_vector(0.0, N);
  vector[N] y = z - b0;
  vector[N] sigmas;
  sigmas = exp(W * alpha);
  matrix[N, N] R = pexp_corr(dist, rho, nu);
  matrix[N, N] Sigma;
  Sigma = quad_form_diag(R, sigmas);
  return multi_normal_lpdf(y | mu, Sigma);
}

real gp_xb_lpdf(data vector z,
                data matrix dist, data matrix X,
                real sigma, real rho, real nu,
                vector beta) {
  int N, K;
  N = rows(dist);
  K = cols(dist);
  vector[N] mu;
  mu = X * beta;
  matrix[N, N] Sigma  = square(sigma) *
    pexp_corr(dist, rho, nu);
  return multi_normal_lpdf(z | mu, Sigma);
}

real gp_heter_xb_lpdf(vector z,
                      data matrix dist,
                      data matrix X,
                      data matrix W,
                      vector alpha,
                      real rho, real nu,
                      real tau,
                      vector beta) {
  int N, K;
  N = rows(dist);
  K = cols(dist);
  vector[N] mu;
  vector[N] sigmas;
  mu = X * beta;
  sigmas = exp(W * alpha);
  matrix[N, N] R = pexp_corr(dist, rho, nu);
  matrix[N, N] Sigma;
  Sigma = add_diag(quad_form_diag(R, sigmas), square(tau));
  return multi_normal_lpdf(z | mu, Sigma);
}
