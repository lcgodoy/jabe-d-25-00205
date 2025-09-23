library(dplyr)
library(ggplot2)

int_score <- function(y, l, u, alpha) {
  alpha <- 2 / alpha
  ind_l <- as.numeric(y < l)
  ind_u <- as.numeric(y > u)
  (u - l) + alpha * ((l - y) * ind_l + (y - u) * ind_u)
}

cv_results <- rbind(readRDS("data/pm25/cluster/cv-hgp.rds"),
                    readRDS("data/pm25/cluster/cv-inla.rds"))

model_lookup <-
  data.frame(model_id = 1:10) |>
  transform(vart = ifelse(model_id %% 2 == 0, "Het", "Hom"),
            family = c(rep("Gaussian-centroid", 2),
                       rep("Gaussian-HGP", 2),
                       rep("Gamma", 2),
                       rep("Log-Normal", 2),
                       "AGP1", "AGP2"))

cv_results |>
  left_join(model_lookup, by = "model_id") |>
  ## filter(model_id %% 2 == 0) |>
  mutate(bias = median_prd - obs,
         is   = int_score(y = obs, l = lower, u = upper,
                          alpha = .05),
         cvg = between(obs, lower, upper)) |>
  group_by(type, vart, family, k) |>
  summarise(bias = mean(bias),
            rmse = sqrt(mean(bias * bias)),
            is   = mean(is),
            cvg  = 100 * mean(cvg)) |>
  ungroup() |>
  group_by(type, vart, family) |>
  summarise(across(where(is.numeric), mean))

cv_results |>
  left_join(model_lookup, by = "model_id") |>
  ## filter(model_id %% 2 == 0) |>
  mutate(bias = median_prd - obs,
         is   = int_score(y = obs, l = lower, u = upper,
                          alpha = .05),
         cvg = between(obs, lower, upper)) |>
  group_by(vart, family, k) |>
  summarise(bias = mean(bias),
            rmse = sqrt(mean(bias * bias)),
            is   = mean(is),
            cvg  = 100 * mean(cvg)) |>
  ungroup() |>
  group_by(vart, family) |>
  summarise(across(where(is.numeric), mean)) |>
  ungroup() |>
  arrange(rmse) |>
  mutate(rmse = 10 * rmse) |>
  select(- k, - bias) |>
    knitr::kable(format = "latex",
                 booktabs = TRUE,
                 linesep = "",
                 digits = 3)
