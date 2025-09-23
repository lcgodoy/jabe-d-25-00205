functions {
  real icar_normal_lpdf(vector phi, int N, int[] node1,
			int[] node2) {
    return -0.5 * dot_self(phi[node1] - phi[node2]) +
      normal_lpdf(sum(phi) | 0, 0.001 * N);
  }
}
data {
  int<lower = 0> N;
  int<lower = 0> N_edges;
  int<lower = 1, upper = N> node1[N_edges];  // node1[i] adjacent to node2[i]
  int<lower = 1, upper = N> node2[N_edges];  // and node1[i] < node2[i]
  int <lower = 1> k;
  int y[N];//observed counts
  vector[N] offset;//expected counts
  matrix[N, k] X; // covariates
  real<lower = 0> sc_f; // scaling factor
}
transformed data {
  vector[N] logoff;
  logoff = log(offset);
}
parameters {
  vector[k] beta;
  real beta_0;
  real<lower = 0> sigma; // precision of heterogeneous effects
  real logit_rho;
  vector[N] phi; // spatial effects
  vector[N] eps; // iid effect
}
transformed parameters {
  // real<lower=0> sigma = inv(sqrt(tau));  // convert precision to sigma
  real xi = inv_logit(logit_rho);
  vector[N] z = sqrt(xi / sc_f) * phi + sqrt(1 - xi) * eps;
  z = sigma * z;
}
model {
  vector[N] intercept;
  intercept = beta_0 + logoff + z;
  target += poisson_log_glm_lpmf(y | X, intercept, beta);
  target += icar_normal_lpdf(phi | N, node1, node2);
  target += std_normal_lpdf(eps);
  target += std_normal_lpdf(logit_rho);
  target += std_normal_lpdf(sigma) - normal_lccdf(0 | 0, 1);
  target += std_normal_lpdf(beta);
  target += std_normal_lpdf(beta_0);
}
generated quantities {
  vector[N] log_lik;
  vector[N] y_rep;
  vector[N] SIR = exp(beta_0 + X * beta + z);
  vector[N] linpred = beta_0 + logoff + X * beta + z;
  for (i in 1:N) {
    if (linpred[i] > 20.7944) {
      y_rep[i] = -1;
    }
    else {
      y_rep[i] =
        poisson_log_rng(linpred[i]);
    }
    log_lik[i] =
      poisson_log_lpmf(y[i] | linpred[i]);
  }
}
