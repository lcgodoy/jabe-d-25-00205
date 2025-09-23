### produces vector with list of neighbors
neighbors <- function(Minc) {
    n <- nrow(Minc)
    unlist(lapply(2:n, function(i) which(Minc[i, 1:(i - 1)] == 1)))
}

### produces number of directed neighbors and a vector where each the adjacency
### of each node ends in the nei vector
adj_index <- function(Minc) {
    Minc.low <- Minc
    Minc.low[upper.tri(Minc)] <- 0
    list(N_nei = rowSums(Minc.low), adj_index = cumsum(rowSums(Minc.low)))
}

comp_corr_dagar <- function(N, Minc, rho) {
    a_i <- adj_index(Minc)
    N_nei <- a_i$N_nei
    B <- matrix(0, nrow = N, ncol = N)
    tau_i <- vector(mode = "numeric", length = N)
    b_aux <-  rho / (1 + (N_nei - rep(1, N)) * rho * rho)
    for (i in seq_len(N)) {
        for (j in 1:i) {
            B[i, j] <- Minc[i, j] * b_aux[i]
        }
        tau_i[i] <- rho / (b_aux[i] * (1 - rho * rho))
    }
    L <- diag(N) - B
    Q <- crossprod(L, diag(tau_i)) %*% L
    Sigma <- chol2inv(chol(Q))
    corr_dagar <- matrix(0, nrow = N, ncol = N)
    for (i in 1:N) {
        for (j in 1:i) {
            corr_dagar[i, j] <-
                Sigma[i, j] / sqrt(Sigma[i, i] * Sigma[j, j])
            corr_dagar[j, i] <- corr_dagar[i, j]
        }
    }
    return(corr_dagar)
}
