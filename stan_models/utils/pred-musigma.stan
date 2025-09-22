matrix kronecker_prod(matrix A, matrix B) {
  matrix[rows(A) * rows(B), cols(A) * cols(B)] C;
  int m;
  int n;
  int p;
  int q;
  m = rows(A);
  n = cols(A);
  p = rows(B);
  q = cols(B);
  for (i in 1:m) {
    for (j in 1:n) {
      int row_start;
      int row_end;
      int col_start;
      int col_end;
      row_start = (i - 1) * p + 1;
      row_end = (i - 1) * p + p;
      col_start = (j - 1) * q + 1;
      col_end = (j - 1) * q + q;
      C[row_start:row_end, col_start:col_end] = A[i, j] * B;
    }
  }
  return C;
}

// faster kronecker product for lowertri matrices
// https://discourse.mc-stan.org/t/stan-efficient/26329
matrix chol_kronecker_prod(matrix A, matrix B) {
  int m = rows(A);
  int p = rows(B);
  int N = m * p;
  matrix[N, N] C = rep_matrix(0., N, N);

  for (i in 1:m) {
    for (j in 1:i) {
      if (fabs(A[i, j]) > 1e-12) {
        int row_start = (i - 1) * p + 1;
        int row_end = (i - 1) * p + p;
        int col_start = (j - 1) * p + 1;
        int col_end = (j - 1) * p + p;
        C[row_start:row_end, col_start:col_end] = A[i, j] * B;
      }
    }
  }
  return C;
}

// tuple(vector, matrix) musig(vector mu, vector z,
vector pred_rng(vector mu,
                vector z,
                vector marg_mu_pred,
                matrix marg_Sig_pred,
                matrix L_s,
                matrix L_t,
                matrix Sigmat_cross) {
  int N_t = rows(L_t);
  int N_s = rows(L_s);
  int N_tnew = rows(Sigmat_cross);
  int N_new = N_tnew * N_s;
  vector[N_new] mu_pred;
  matrix[N_new, N_new] L_pred;
  // Sig[1:tst, 1:t] * Sig[1:t, 1:t]^{-1}
  matrix[N_tnew, N_t] v_pred;
  // Sig[1:tst, 1:t] * Sig[1:t, 1:t]^{-1} * Sig[1:t, 1:tst]
  matrix[N_tnew, N_tnew] aux_pred;
  // Sig[t + 1, 1:t]  L_t^{-t} 
  v_pred = mdivide_left_tri_low(L_t, Sigmat_cross')';
  // Sig[t + 1, 1:t] * L_t^{-t} * L_t^{-1} 
  v_pred = mdivide_right_tri_low(v_pred, L_t);
  aux_pred = marg_Sig_pred - v_pred * Sigmat_cross';
  mu_pred = marg_mu_pred +
    kronecker_prod(diag_matrix(rep_vector(1.0, N_s)), v_pred) *
    (z - mu);
  L_pred = chol_kronecker_prod(L_s, cholesky_decompose(aux_pred));
  return multi_normal_cholesky_rng(mu_pred, L_pred);
}
