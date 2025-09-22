matrix pexp_corr(data matrix dists, real rho, real nu) {
  int N = rows(dists);
  int M = cols(dists);
  real phi = rho / pow(log(10), 1 / nu);
  matrix[N, M] out;
  if (N == M) {
    for (i in 1:N) {
      for (j in 1:i) {
        out[i, j] = exp( - pow(dists[i, j] / phi, nu));
        out[j, i] = out[i, j];
      }
    }
  } else {
    for (i in 1:N) {
      for (j in 1:M) {
        out[i, j] = exp( - pow(dists[i, j] / phi, nu));
      }
    }
  }
  return out;
}
real pexp_corr_s(real dists, real rho, real nu) {
  real phi = rho / pow(log(10), 1 / nu);
  real out = exp( - pow(dists / phi, nu));
  return out;
}
vector pexp_corr_v(vector dists, real rho, real nu) {
  real phi = rho / pow(log(10), 1 / nu);
  vector[num_elements(dists)] out;
  out = exp( - pow(dists / phi, nu));
  return out;
}

matrix pexp_hetcov(data matrix dists,
                   vector sigmas_i,
                   vector sigmas_j,
                   real rho, real nu) {
  int N = rows(dists);
  int M = cols(dists);
  real phi = rho / pow(log(10), 1 / nu);
  matrix[N, M] out;
  if (N == M) {
    for (i in 1:N) {
      for (j in 1:i) {
        out[i, j] = exp( - pow(dists[i, j] / phi, nu)) *
          sigmas_i[i] * sigmas_j[j];
        out[j, i] = out[i, j];
      }
    }
  } else {
    for (i in 1:N) {
      for (j in 1:M) {
        out[i, j] = exp( - pow(dists[i, j] / phi, nu)) *
          sigmas_i[i] * sigmas_j[j];
      }
    }
  }
  return out;
}

real log_exp_lpdf(real x, real beta) {
  return log(beta) + x - beta * exp(x);
}

real beta_4p_lpdf(real y, real alpha, real beta,
                  real lower, real upper) {
    // Scale 4-parameter Beta RV to 2-parameter Beta RV
    real x = (y - lower) / (upper - lower);
    // Return scaled 2-parameter beta lpdf
    return beta_lpdf(x | alpha, beta) - log(upper - lower);
  }
