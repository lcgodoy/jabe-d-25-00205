if (!interactive())
  .libPaths("~/rlibs")

## Load necessary packages
library(parallel)
library(rstan)

source("utils/cv-task-funs.R")

if (interactive()) {
  kf <- 10
  samples <- 50
  warmup <- 100
  chains <- 2
  cores <- 1
  .debug <- 1
} else {
  kf <- Sys.getenv("K") |>
    as.integer()
  samples <- Sys.getenv("SAMPLES") |>
    as.integer()
  warmup <- Sys.getenv("WARMUP") |>
    as.integer()
  chains <- Sys.getenv("CHAINS") |>
    as.integer()
  cores <- Sys.getenv("CORES") |>
    as.integer()
  .debug <- Sys.getenv("DEBUG") |>
    as.integer()
}

## cross-validation and model parameters
models <- if (.debug == 1) c(1, 6) else 1:8

ntasks <- Sys.getenv("SLURM_NTASKS")
ntasks <- ifelse(ntasks == "", 2, as.integer(ntasks))

##--- load data ----

xy 	<- readRDS("data/pm25/cluster/xy.rds")
W  	<- readRDS("data/pm25/cluster/W.rds")
dists_ct <- readRDS("data/pm25/cluster/ct_mat.rds") / 10 ## in 10s of km
dists_h  <- readRDS("data/pm25/cluster/hmat.rds") / 10  ## in 10s of km

##--- compile stan models ----
stan_models <- list(
  hgp_gaus = stan_model("stan_models/hgp-gaussian.stan"),
  hgp_ln = stan_model("stan_models/hgp-lognormal.stan"),
  pred_gaus = stan_model("stan_models/gaus-pred.stan"),
  pred_ln = stan_model("stan_models/ln-pred.stan")
)

##--- execution ----

seq_k <- if(.debug) 1:2 else 1:kf
task_grid <- expand.grid(ki = 1:kf, mod = models)
tasks <- split(task_grid, seq(nrow(task_grid)))

results <- mclapply(tasks,
                    FUN = run_cv_hgp,
                    mc.cores = ntasks)

saveRDS(do.call(rbind, results),
        file = "data/pm25/cluster/cv-hgp.rds")
