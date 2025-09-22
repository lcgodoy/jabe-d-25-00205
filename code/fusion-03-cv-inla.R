## Set library path
if (!interactive())
  .libPaths("~/rlibs")

## Load necessary packages
library(parallel)
library(INLA)
library(sp)

source("code/utils/cv-task-funs.R")

##---+ Configuration ----
## Define all parameters here instead of in the bash script

if (interactive()) {
  kf <- 10
} else {
  kf <- Sys.getenv("K") |>
    as.integer()
}

## Cross-validation and Model parameters
models <- 9:10

ntasks <- Sys.getenv("SLURM_NTASKS")
ntasks <- ifelse(ntasks == "", detectCores() - 1, as.integer(ntasks))

##---+ Data Loading ----

my_dt <- readRDS("data/pm25/dataset.rds")

##---+ Execution ----

## Create a grid of all parameter combinations
seq_k <- 1:kf
task_grid <- expand.grid(ki = 1:kf, mod = models)

## Convert grid to a list of parameters for mclapply
tasks <- split(task_grid, seq(nrow(task_grid)))

results <- mclapply(tasks,
                    FUN = run_cv_inla,
                    mc.cores = ntasks)

saveRDS(do.call(rbind, results),
        file = "data/pm25/cluster/cv-inla.rds")
