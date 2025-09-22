real nngp_lpdf(vector w, real sigma, real rho, matrix NN_dist,
               matrix NN_distM, array[,] int NN_ind, int N, int M,
               real nu) {
  vector[N] V;
  vector[N] I_Aw = w;
  int dim;
  int h;
  for (i in 2:N) {
    matrix[ i < (M + 1)? (i - 1) : M, i < (M + 1)? (i - 1): M]
      iNNdistM;
    matrix[ i < (M + 1)? (i - 1) : M, i < (M + 1)? (i - 1): M]
      iNNCholL;
    vector[ i < (M + 1)? (i - 1) : M] iNNcorr;
    vector[ i < (M + 1)? (i - 1) : M] v;
    row_vector[i < (M + 1)? (i - 1) : M] v2;
    dim = (i < (M + 1))? (i - 1) : M;
    if(dim == 1){iNNdistM[1, 1] = 1;}
    else{
      h = 0;
      for (j in 1:(dim - 1)){
        for (k in (j + 1):dim){
          h = h + 1;
          iNNdistM[j, k] = pexp_corr_s(NN_distM[(i - 1), h], rho,
                                       nu);
          iNNdistM[k, j] = iNNdistM[j, k];
        }
      }
      for(j in 1:dim){
        iNNdistM[j, j] = 1;
      }
    }
    iNNCholL = cholesky_decompose(iNNdistM);
    iNNcorr = pexp_corr_v(to_vector(NN_dist[(i - 1), 1:dim]),
                          rho, nu);
    v = mdivide_left_tri_low(iNNCholL, iNNcorr);
    V[i] = 1 - dot_self(v);
    v2 = mdivide_right_tri_low(v', iNNCholL);
    I_Aw[i] = I_Aw[i] - v2 * w[NN_ind[(i - 1), 1:dim]];
  }
  V[1] = 1;
  return - 0.5 * ( 1 / square(sigma) * dot_product(I_Aw, (I_Aw ./ V)) +
                   sum(log(V)) + N * 2 * log(sigma));
}

real nngp_w_lpdf(vector w_b1, real sigma, real rho, matrix NN_dist,
                 matrix NN_distM, array[,] int NN_ind, int N, int M,
                 real nu, real intercept) {
  vector[N] V;
  vector[N] w = w_b1 - intercept;
  vector[N] I_Aw = w;
  int dim;
  int h;
  for (i in 2:N) {
    matrix[ i < (M + 1)? (i - 1) : M, i < (M + 1)? (i - 1): M]
      iNNdistM;
    matrix[ i < (M + 1)? (i - 1) : M, i < (M + 1)? (i - 1): M]
      iNNCholL;
    vector[ i < (M + 1)? (i - 1) : M] iNNcorr;
    vector[ i < (M + 1)? (i - 1) : M] v;
    row_vector[i < (M + 1)? (i - 1) : M] v2;
    dim = (i < (M + 1))? (i - 1) : M;
    if(dim == 1){iNNdistM[1, 1] = 1;}
    else{
      h = 0;
      for (j in 1:(dim - 1)){
        for (k in (j + 1):dim){
          h = h + 1;
          iNNdistM[j, k] = pexp_corr_s(NN_distM[(i - 1), h], rho,
                                       nu);
          iNNdistM[k, j] = iNNdistM[j, k];
        }
      }
      for(j in 1:dim){
        iNNdistM[j, j] = 1;
      }
    }
    iNNCholL = cholesky_decompose(iNNdistM);
    iNNcorr = pexp_corr_v(to_vector(NN_dist[(i - 1), 1:dim]),
                          rho, nu);
    v = mdivide_left_tri_low(iNNCholL, iNNcorr);
    V[i] = 1 - dot_self(v);
    v2 = mdivide_right_tri_low(v', iNNCholL);
    I_Aw[i] = I_Aw[i] - v2 * w[NN_ind[(i - 1), 1:dim]];
  }
  V[1] = 1;
  return - 0.5 * ( 1 / square(sigma) * dot_product(I_Aw, (I_Aw ./ V)) +
                   sum(log(V)) + N * 2 * log(sigma));
}

real nngp_lpdf(vector Y, vector X_beta,
               real sigmasq, real tausq,
               real nu, real rho,
               matrix NN_dist, matrix NN_distM, int[,] NN_ind,
               int N, int M) {

  vector[N] V;
  vector[N] YXb = Y - X_beta;
  vector[N] U = YXb;
  real kappa_p_1 = tausq / sigmasq + 1;
  int dim;
  int h;

  for (i in 2:N) {
    matrix[ i < (M + 1) ? (i - 1) : M, i < (M + 1) ? (i - 1): M]
      iNNdistM;
    matrix[ i < (M + 1) ? (i - 1) : M, i < (M + 1) ? (i - 1): M]
      iNNCholL;
    vector[ i < (M + 1) ? (i - 1) : M] iNNcorr;
    vector[ i < (M + 1) ? (i - 1) : M] v;
    row_vector[i < (M + 1) ? (i - 1) : M] v2;
    dim = (i < (M + 1))? (i - 1) : M;

    if(dim == 1){iNNdistM[1, 1] = kappa_p_1;}
    else{
      h = 0;
      for (j in 1:(dim - 1)){
        for (k in (j + 1):dim){
          h = h + 1;
          iNNdistM[j, k] = pexp_corr_s(NN_distM[(i - 1), h], rho,
                                       nu);
          iNNdistM[k, j] = iNNdistM[j, k];
        }
      }
      for(j in 1:dim){
        iNNdistM[j, j] = kappa_p_1;
      }
    }

    iNNCholL = cholesky_decompose(iNNdistM);
    iNNcorr = pexp_corr_v(to_vector(NN_dist[(i - 1), 1:dim]),
                          rho, nu);

    v = mdivide_left_tri_low(iNNCholL, iNNcorr);

    V[i] = kappa_p_1 - dot_self(v);

    v2 = mdivide_right_tri_low(v', iNNCholL);

    U[i] = U[i] - v2 * YXb[NN_ind[(i - 1), 1:dim]];
  }
  V[1] = kappa_p_1;
  return - 0.5 * ( 1 / sigmasq * dot_product(U, (U ./ V)) +
                   sum(log(V)) + N * log(sigmasq));
}
