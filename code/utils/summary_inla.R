aux_inla <- function(x, cl) {
  alpha <- 1 - cl
  x <- INLA::inla.smarginal(x)
  x <- c(mean(x[[1]]),
         as.numeric(quantile(x, probs = c(alpha / 2, 1 - alpha / 2))))
  x
}

inla_summ_stat <- function(x, cl) {
  alpha <- 1 - cl
    c(mean(x),
      quantile(x, probs = c(alpha / 2, 1 - alpha / 2)) |>
      as.numeric())
}

summary_parest_spde <- function(res, cl = .95) {
  svar <-
    INLA::inla.rmarginal(n = 2000,
                         res$marginals.variance.nominal$variance.nominal.1)
  svar <- sqrt(svar) |> inla_summ_stat(cl)
  range <- res$marginals.range.nominal$range.nominal.1 |>
    inla.rmarginal(n = 2000, marginal = _)
  range <- range / 10000 ## 10s of km
  range <- inla_summ_stat(range, cl)
  return(data.frame(parameter = c("sigma", "rho"),
                    mean      = c(svar[1], range[1]),
                    low       = c(svar[2], range[2]),
                    upp       = c(svar[3], range[3])))
}

summary_parest_inla <- function(res, cl = .95) {
  marg1 <- res$marginals.hyperpar[[1]] |>
    INLA::inla.rmarginal(n = 2000,marginal = _)
  marg1 <- inla_summ_stat(1 / marg1, cl)
  fixed <- res$marginals.fixed[[1]] |>
    aux_inla(cl = cl)
  data.frame(parameter = c("beta", "tau"),
             mean = c(fixed[1], marg1[1]),
             low = c(fixed[2], marg1[2]),
             upp = c(fixed[3], marg1[3]))
}
