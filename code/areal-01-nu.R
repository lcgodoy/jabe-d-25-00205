library(sf)
library(rstan)

##--- auxiliary functions ----

rstan_options(auto_write = TRUE)

source("code/utils/aux_icar.R")
source("code/utils/aux_dagar.R")
source("code/utils/process-chains.R")

## beta > 2 * x
alpha_given_beta_mode <- function(.beta, x) {
  (.beta - 2 * x) / (1 - x)
}

##--- loading data ----

data(GGHB.IZ, package = "CARBayesdata")
data(respiratorydata, package = "CARBayesdata")

my_dt <- merge(x = GGHB.IZ,
               y = respiratorydata, by = "IZ",
               all.x = FALSE)

##---+ aux quantities +----

N <- nrow(my_dt)
my_x <- matrix(scale(my_dt$incomedep), ncol = 1)

dists <- st_distance(x = st_geometry(my_dt),
                     by_element = FALSE,
                     which      = "Hausdorff")  |>
  units::set_units(value = "km") |>
  units::set_units(value = NULL) |>
  as.matrix()

range_dist <- range(dists[upper.tri(dists)])

my_x <- matrix(scale(my_dt$incomedep), ncol = 1)

areas <- st_area(my_dt) |>
  units::set_units("km^2") |>
  units::set_units(NULL) |>
  scale() |>
  as.numeric()

nu <- .7
rho_0 <- 0.8 * max(dists[upper.tri(dists)])

hgp_dat <- list(
    N = N,
    k = ncol(my_x),
    y = my_dt$observed,
    offset = my_dt$expected,
    dists = dists,
    X = my_x,
    prob_rho = .05,
    rho_0 = rho_0,
    aux_res = sqrt(.Machine$double.eps),
    nu = nu
)

## hgp_dat$prob_rho <- .01

##---+ Compiling models +----

hgp <- stan_model(
    file = "stan_models/hgp-poisson.stan",
    allow_optimizations = TRUE
)

##---+ Running MCMCs for different \nu values +----

n_iter <- 1000L
burn <- 1000L
thin <- 1

nu <- c(.6, .7, .8, .9)

all_hgps <-
  lapply(nu, \(x) {
    hgp_dat$nu <- x
    sampling(
        object = hgp,
        data = hgp_dat,
        iter = n_iter + burn,
        warmup = burn,
        thin = thin,
        pars = "chol_sig",
        include = FALSE,
        save_warmup = FALSE,
        chains = 4,
        cores = 4,
        seed = 202310
    )
  })

params <- c("beta_0", "beta", "rho", "sigma")

lapply(all_hgps, \(x) print(x, pars = params))

##---+ Goodness-of-fit +---

loo_all <- lapply(all_hgps,
                  \(x) {
                    loo::loo(x,
                             ## moment_match = TRUE,
                             save_psis = TRUE,
                             cores = 8)
                  })
ll_all <- lapply(all_hgps, loo::extract_log_lik)
waic_all <- lapply(ll_all, loo::waic)

loo::loo_compare(loo_all)

## Table S.4
tbl_s4 <-
  rbind(sapply(loo_all, \(x) x$estimates[3, 1]),
        sapply(waic_all, \(x) x$estimates[3, 1]))

colnames(tbl_s4) <- round(nu, 1)
