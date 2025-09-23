data {
  int <lower=1> N;
  int <lower=1> k;
  vector[N] N_nei; //number of directed neighbors for each node
  int <lower=1> N_edges;  // total numbere of edges
  int <lower=1, upper=N-1> nei[N_edges]; //vector stacking up the direected neighbors of each node
  int <lower=0, upper=N_edges> adjacency_ends[N]; //Where adjacency of each node ends in the nei vector
  int y[N];//observed counts
  vector[N] offset;//expected counts
  matrix[N, k] X; // single covariate (can be easily extended to multiple)
}
transformed data {
  vector[N] logoff = log(offset);
}
parameters {
  real <lower=0, upper=1> psi; // Spatial correlation parameter
  real <lower=0> tau; // sd sp ref
  vector[N] z; // Spatial Random Effect
  vector[k] beta;
  real beta_0;
}
transformed parameters {
  real <lower=0> sigma = inv_sqrt(tau);
}
model {
  vector[N] b; //
  vector[N] vec_var; //
  vector[N] t_rowsum; // only the rowsum of t is used
  vector[N] std_dev; // Rescaled std_dev by std_dev_s
  // Construct s
  vec_var = (1 - psi * psi) ./ (1 + (N_nei - rep_vector(1,N)) * psi * psi);
  b = psi ./ (1 + (N_nei - rep_vector(1,N)) * psi * psi );
  // Linear in number of edges
  t_rowsum[1] = 0;
  for(i in 2:N){
    t_rowsum[i] =
      sum(z[nei[(adjacency_ends[i-1]+1):adjacency_ends[i]]]) * b[i];
  }
  std_dev = sqrt(vec_var / tau);
  vector[N] intercept;
  intercept = beta_0 + logoff + z;
  target += normal_lpdf(z | t_rowsum, std_dev);
  target += poisson_log_glm_lpmf(y | X, intercept, beta);
  target += normal_lpdf(beta | 0, 1000);
  target += beta_lpdf(psi | 1, 1);
  target += gamma_lpdf(tau | 2, 1);
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
