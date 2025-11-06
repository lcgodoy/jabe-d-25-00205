library(sf)
library(rstan)
library(ggplot2)
library(bayesplot)
library(INLA)
library(sp)
library(dplyr)

set.seed(2024)

theme_map <- theme_minimal() +
  theme(axis.text = element_blank(),
        strip.background = element_rect(color = 1))

theme_set(theme_bw() +
          theme(strip.background = element_rect(fill = "white")))

gr <- (1 + sqrt(5)) / 2

fig_dir <- "~/git-projects/hausdorff-GP/img"

##--- auxiliary functions ----

int_score <- function(y, l, u, alpha) {
  alpha <- 2 / alpha
  ind_l <- as.numeric(y < l)
  ind_u <- as.numeric(y > u)
  (u - l) + alpha * ((l - y) * ind_l + (y - u) * ind_u)
}

##--- reading and pre-processing data ----

crs_epsg <- 26911

my_dt <- readRDS("data/pm25/dataset.rds")
border_poly <- readRDS("data/pm25/border_counties.rds")

ggplot(data = dplyr::arrange(my_dt, id),
       aes(fill = pm25)) +
  geom_sf(size = 2, pch = 21,
          lwd = .1) +
  ## scale_fill_distiller(palette = "Spectral",
  scale_fill_viridis_c(option = "H",
                       guide = guide_colorbar(expression(PM[2.5]),
                                              title.vjust = .75,
                                              ticks = FALSE)) +
  theme_map +
  theme(legend.position = "inside",
        legend.position.inside = c(.2, .125),
        legend.direction = "horizontal",
        panel.grid = element_blank(),
        legend.key.size = unit(.5, "cm"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10))

ggsave(filename = file.path(fig_dir, "pm25_data.pdf"),
       width = 3.5,
       height = 3.5)

dists <- my_dt |>
  ## st_geometry() |>
  st_distance(by_element = FALSE,
              which = "Hausdorff") |>
  units::set_units(value = "km") |>
  units::set_units(value = NULL) |>
  as.matrix()

hd_lim <- range(dists[upper.tri(dists)])

N <- NROW(my_dt)

x <- matrix(1, nrow = N)

areas <- st_area(my_dt) |>
  units::set_units(value = "km^2") |>
  units::set_units(value = NULL) |>
  as.numeric()

areas_bin <- ifelse(areas == 0, 0, 1)

A <- cbind(1, areas_bin)

##--- regions for predictions ----

pred_regions <- readRDS("data/pm25/pred-locations.rds")

pred_dists <- pred_regions |>
  st_geometry() |>
  st_distance(by_element = FALSE,
              which = "Hausdorff") |>
  units::set_units(value = "km") |>
  units::set_units(value = NULL) |>
  as.matrix()

cross_dists <- pred_regions |>
  st_geometry() |>
  st_distance(x = my_dt, y = _,
              which = "Hausdorff") |>
  units::set_units(value = "km") |>
  units::set_units(value = NULL) |>
  as.matrix()

pred_n <- NROW(pred_regions)

x_new <- matrix(1, nrow = pred_n)

areas_new <- st_area(pred_regions) |>
  units::set_units(value = "km^2") |>
  units::set_units(value = NULL) |>
  as.numeric()

areas_bin_new <- ifelse(areas_new == 0, 0, 1)

A_new <- cbind(1, areas_bin_new)

##---+ HGP +----

nu <- .7

hgp_dat <- list(
    N = N,
    k = NCOL(x),
    k2 = ncol(A),
    dists = dists,
    rho_0 = hd_lim[2] * .8,
    tau_0 = sd(my_dt$pm25) / 2,
    pt_u = .05,
    prob_rho = .05,
    y = my_dt$pm25,
    X = x,
    A = A
)

##----+ Compiling models +---

code_hgp <- stan_model(file = "stan_models/hgp-gaussian.stan")

##----+ Sampling +---

n_iter <- 1000L
burn <- 1000L
thin <- 1

nus <- seq(from = .6, to = .9, by = .1)

pars_hgp <- c("beta", "rho", "alpha", "tau")

all_hgps <-
  lapply(nus, \(x) {
    hgp_dat$nu <- x
    sampling(
        object = code_hgp,
        data = hgp_dat,
        iter = n_iter + burn,
        warmup = burn,
        thin = thin,
        chains = 4,
        cores = 4,
        save_warmup = FALSE,
        seed = 2024
    )})

##---+ Goodness-of-fit +----

stan_ll <-
  stan_model(file = "stan_models/gaussian-gq.stan")

all_gof <-
  Map(\(fit, .nu) {
    hgp_dat$nu <- .nu
    gqs(stan_ll,
        data = hgp_dat,
        draws = as.matrix(fit,
                          pars = pars_hgp))
  }, fit = all_hgps, .nu = nus)

gofs <-
  Map(\(fit, .nu) {
    loo_hgp <- loo::loo(fit,
                        save_psis = TRUE,
                        cores = 8)
    ll_hgp <- loo::extract_log_lik(fit)
    waic_hgp <- loo::waic(ll_hgp)
    c("nu" = .nu,
      "LOOIC" = loo_hgp$estimates[3, 1],
      "WAIC" = waic_hgp$estimates[3, 1])
  }, fit = all_gof,
  .nu = nus)

##--- Table S.5 ----

do.call(rbind, gofs) |>
  t() |>
  knitr::kable(format = "latex",
               booktabs = TRUE,
               linesep = "",
               digits = 1)

best_model <- which.min(sapply(gofs, \(x) x[2]))
fit_hgp <- all_hgps[[best_model]]

print(fit_hgp, pars = pars_hgp)

pairs(fit_hgp, pars = pars_hgp)

rhat(fit_hgp, pars = pars_hgp) |>
  max()

color_scheme_set("mix-blue-red")

mcmc_trace(x = fit_hgp,
           regex_pars = pars_hgp,
           facet_args = list(labeller = label_parsed)) +
  legend_none()
ggsave(filename = file.path(fig_dir, "la-trace-hgp.pdf"),
       width = 6,
       height = 5)

mcmc_dens_overlay(x = fit_hgp,
                  regex_pars = pars_hgp,
                  facet_args = list(labeller = label_parsed)) +
  legend_none()
ggsave(filename = file.path(fig_dir, "/la-dens-hgp.pdf"),
       width = 6,
       height = 5)

##---+ Goodness-of-fit +----

gof_hgp <- gofs[[best_model]]

##---+ Competing method +----

.center <- function(x)
  as.numeric(scale(x, scale = FALSE))

.center_scale <- function(x)
  as.numeric(scale(x))

## inla.setOption(pardiso.license = "/opt/licenses/pardiso.lic")

stations <- my_dt[grepl("^st_", my_dt$id), ]
dt_estim <- my_dt[!grepl("^st_", my_dt$id), ]

coop <- st_coordinates(stations)

yp <- stations$pm25

spol <- dt_estim |>
  st_geometry() |>
  as("Spatial")

spol@proj4string <- CRS()

ya <- dt_estim$pm25

coopred <- pred_regions |>
  st_centroid() |>
  st_coordinates()

sppr <- pred_regions |>
  st_geometry() |>
  as("Spatial")

sppr@proj4string <- CRS()

ypred <- rep(NA_real_, nrow(coopred))

## creating meshes
grid <- "sparse"
max_edge_ct <- ifelse(grid == "fine", 5000, 10000)
max_edge <- c(max_edge_ct,
              10 * max_edge_ct)
cutoff <- max_edge_ct / 100
offset <- c(max_edge_ct,
            5 * max_edge_ct)

mesh_sparse <- inla.mesh.2d(rbind(coopred, coop),
                            ## boundary = borders_xy,
                            max.edge =
                              max_edge,
                            cutoff = cutoff,
                            offset = offset)

spde_sparse <- inla.spde2.matern(mesh = mesh_sparse,
                                 alpha = 2)

## Point observations
Ap <- inla.spde.make.A(mesh = mesh_sparse, loc = coop)
stk_p <- inla.stack(tag = "point",
                    data = list(y = yp),
                    A = list(Ap, 1),
                    effects = list(s = 1:spde_sparse$n.spde,
                                   data.frame(b0 = rep(1, nrow(coop)))))

## Areal observations
locin <- mesh_sparse$loc[as.vector(which(!is.na(over(SpatialPoints(mesh_sparse$loc),
                                                 spol)))), ]
block <- rep(0, nrow(locin))

for (i in seq_len(length(spol))) {
  block[as.vector(which(!is.na(over(SpatialPoints(locin), spol[i]))))] <- i
}

Aa <- inla.spde.make.A(mesh = mesh_sparse,
                       loc = locin,
                       block = block,
                       block.rescale = "sum")
stk_a <- inla.stack(tag = "areal",
                    data = list(y = ya),
                    A = list(Aa, 1),
                    effects = list(s = 1:spde_sparse$n.spde,
                                   data.frame(b0 = rep(1, NROW(ya)))))

## Prediction points
Apred <- inla.spde.make.A(mesh = mesh_sparse, loc = coopred)
stk_pred <- inla.stack(tag = "pred",
                       data = list(y = ypred),
                       A = list(Apred, 1),
                       effects = list(s = 1:spde_sparse$n.spde,
                                      data.frame(b0 = rep(1, NROW(ypred)))))

## it may be dropped
## source("code/utils/moraga_aux.R")

sparse_full <- inla.stack(stk_p, stk_a, stk_pred)

formula <- y ~ -1 + b0 + f(s, model = spde_sparse)

## fitting model with sparse grid

res_sparse <- inla(formula,
                   data = inla.stack.data(sparse_full),
                   control.predictor = list(compute = TRUE,
                                            A = inla.stack.A(sparse_full)),
                   control.compute = list(waic = TRUE,
                                          cpo = TRUE,
                                          dic = TRUE,
                                          return.marginals.predictor = TRUE),
                   verbose = TRUE)

grid <- "fine"
max_edge_ct <- ifelse(grid == "fine", 5000, 10000)
max_edge <- c(max_edge_ct,
              10 * max_edge_ct)
cutoff <- max_edge_ct / 100
offset <- c(max_edge_ct,
            5 * max_edge_ct)

mesh_fine <- inla.mesh.2d(rbind(coopred, coop),
                        ## boundary = borders_xy,
                        max.edge =
                          max_edge,
                        cutoff = cutoff,
                        offset = offset)

spde_fine <- inla.spde2.matern(mesh = mesh_fine,
                               alpha = 2)

## Point observations
Ap <- inla.spde.make.A(mesh = mesh_fine, loc = coop)
stk_p <- inla.stack(tag = "point",
                    data = list(y = yp),
                    A = list(Ap, 1),
                    effects = list(s = 1:spde_fine$n.spde,
                                   data.frame(b0 = rep(1, nrow(coop)))))

## Areal observations
locin <- mesh_fine$loc[as.vector(which(!is.na(over(SpatialPoints(mesh_fine$loc),
                                                 spol)))), ]
block <- rep(0, nrow(locin))

for (i in seq_len(length(spol))) {
  block[as.vector(which(!is.na(over(SpatialPoints(locin), spol[i]))))] <- i
}

Aa <- inla.spde.make.A(mesh = mesh_fine,
                       loc = locin,
                       block = block,
                       block.rescale = "sum")
stk_a <- inla.stack(tag = "areal",
                    data = list(y = ya),
                    A = list(Aa, 1),
                    effects = list(s = 1:spde_fine$n.spde,
                                   data.frame(b0 = rep(1, NROW(ya)))))

## Prediction points
Apred <- inla.spde.make.A(mesh = mesh_fine, loc = coopred)
stk_pred <- inla.stack(tag = "pred",
                       data = list(y = ypred),
                       A = list(Apred, 1),
                       effects = list(s = 1:spde_fine$n.spde,
                                      data.frame(b0 = rep(1, NROW(ypred)))))

fine_full <- inla.stack(stk_p, stk_a, stk_pred)

formula <- y ~ -1 + b0 + f(s, model = spde_fine)

res_fine <- inla(formula,
                 data = inla.stack.data(fine_full),
                 control.predictor = list(compute = TRUE,
                                          A = inla.stack.A(fine_full)),
                 control.compute = list(waic = TRUE,
                                        cpo = TRUE,
                                        dic = TRUE,
                                        return.marginals.predictor = TRUE),
                 verbose = TRUE)

## saveRDS(list(sparse = mesh_sparse,
##              fine   = mesh_fine),
##         file = "~/phd/thesis/talk/meshes-la.rds")

## plotting meshes
pdf(file.path(fig_dir, "meshes.pdf"),
    width = 4,
    height = 4 / 1.618033)
par(mfrow = c(1, 2), mar = c(0, 0, 1, 0))
plot(mesh_sparse, lwd = .25)
## title("Sparse mesh")
plot(mesh_fine, lwd = .25)
## title("Fine mesh")
dev.off()

##--- Table 5 ----

lpml_sparse <- res_sparse$cpo$cpo |>
  log() |>
  sum(na.rm = TRUE)
gof_sparse <-
  c("LOOIC" = -2 * lpml_sparse,
    "WAIC" = res_sparse$waic$waic)

lpml_fine <- res_fine$cpo$cpo |>
  log() |>
  sum(na.rm = TRUE)
gof_fine <-
  c("LOOIC" = -2 * lpml_fine,
    "WAIC" = res_fine$waic$waic)

## something ODD happening to AGP{2}
gof_tbl <-
  cbind("AGP\\textsubscript{1}" = gof_sparse,
        "AGP\\textsubscript{2}" = gof_fine,
        "HGP" = gof_hgp[-1])

gof_tbl |>
  print() |> 
  knitr::kable(format = "latex", booktabs = TRUE,
               digits = 1,
               linesep = "",
               escape = FALSE)

##---+ out-of-sample predictions +----

id_sparse <- inla.stack.index(
    sparse_full,
    "pred"
)$data

id_fine <- inla.stack.index(
    fine_full,
    "pred"
)$data

pred_dat <- c(hgp_dat,
              list(Xp = x_new,
                   Ap = A_new,
                   Np = NROW(pred_regions),
                   distsp = pred_dists,
                   cdists = cross_dists,
                   nu = nus[best_model]))
  
pars_hgp <- c("beta", "alpha", "rho", "tau")

pred_gaus <- stan_model("stan_models/gaus-pred.stan")

pred_samples <- gqs(
    pred_gaus,
    data = pred_dat,
    draws = as.matrix(fit_hgp, pars = pars_hgp)
)

pred_hgp <- as.matrix(pred_samples, pars = "y_rep")

predictions <-
    rbind(transform(pred_regions,
                    model = "AGP[1]",
                    median_prd = res_sparse$summary.fitted.values[id_sparse,
                                                              "0.5quant"],
                    sd_prd = res_sparse$summary.fitted.values$sd[id_sparse],
                    lower = res_sparse$summary.fitted.values[id_sparse,
                                                           "0.025quant"],
                    upper = res_sparse$summary.fitted.values[id_sparse,
                                                           "0.975quant"]),
          transform(pred_regions,
                    model = "AGP[2]",
                    median_prd = res_fine$summary.fitted.values[id_fine,
                                                                "0.5quant"],
                    sd_prd = res_fine$summary.fitted.values$sd[id_fine],
                    lower = res_fine$summary.fitted.values[id_fine,
                                                           "0.025quant"],
                    upper = res_fine$summary.fitted.values[id_fine,
                                                           "0.975quant"]),
          transform(pred_regions,
                    model = "HGP",
                    median_prd = apply(pred_hgp, 2, median),
                    sd_prd = apply(pred_hgp, 2, sd),
                    lower = apply(pred_hgp, 2, quantile,
                                  probs = .025),
                    upper = apply(pred_hgp, 2, quantile,
                                  probs = .975)))

predictions <- predictions |>
  transform(length_ci = upper - lower)

## predictions <- predictions |>
##   dplyr::mutate(model = gsub("AGG", "AGP", model)) |>
##   dplyr::mutate(model = gsub("1", "[1]", model)) |>
##   dplyr::mutate(model = gsub("2", "[2]", model))

## saveRDS(predictions, file = "~/phd/thesis/talk/df-pred.rds")

predictions |>
  dplyr::select(id, model, median_prd) |>
  tidyr::pivot_wider(names_from = "model",
                     values_from = "median_prd") |>
  dplyr::select(-id) |>
  plot(pal = viridis::turbo,
       main = expression("Predicted " * PM[2.5]))

predictions |>
  dplyr::select(id, model, length_ci) |>
  tidyr::pivot_wider(names_from = "model",
                     values_from = "length_ci") |>
  dplyr::select(-id) |>
  plot(pal = viridis::turbo,
       main = expression("CI Length " * PM[2.5]))

##--- Figure 8 ----

p1 <-
  predictions |>
  dplyr::mutate(model =
                  factor(model,
                         levels = c("HGP", "AGP[1]", "AGP[2]"))) |>
  ggplot(data = _,
         aes(fill = median_prd,
             color = median_prd)) +
  geom_sf() +
  scale_fill_viridis_c(option = "H",
                       ## limits = range(my_dt$pm25),
                       guide = guide_colorbar(expression(hat(y)),
                                              title.vjust = 1,
                                              ticks = FALSE)) +
  scale_color_viridis_c(option = "H",
                        ## limits = range(my_dt$pm25),
                        guide = guide_colorbar(expression(hat(y)),
                                               title.vjust = 1,
                                               ticks = FALSE)) +
  facet_wrap(~ model, labeller = label_parsed) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white",
                                        color = 1),
        legend.position = "bottom",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank())


ggsave(filename = file.path(fig_dir, "predictions_pm25.pdf"),
       plot = p1,
       ## width = 6.5,
       ## height = 6.5 / 1.618033
       width = 5,
       height = 2.3)

## ggsave(filename = "img/lengths_pred_pm25.pdf",
##        plot = p2,
##        width = 6.5,
##        height = 6.5 / 1.618033)

##--- cov functions ----

rho <- as.matrix(fit_hgp_v1, pars = "rho")

##--- parameter estimates ----

samples_hgp <- as.matrix(fit_hgp, pars = pars_hgp)
samples_hgp <-
  cbind(samples_hgp,
        "sigma" = exp(samples_hgp[, "alpha[1]"]),
        "sigma_ar" = exp(samples_hgp[, "alpha[1]"] + samples_hgp[, "alpha[2]"]))

samples_hgp[, "rho"] <- samples_hgp[, "rho"] / 10

parest_hgp <- posterior::summarise_draws(samples_hgp,
                                         ~quantile(.x, probs = c(0.5, 0.025,
                                                                 0.975)))

colnames(parest_hgp) <- c("parameter", "median", "low", "upp")
parest_hgp <-
  transform(parest_hgp, ci = sprintf("[%.3f; %.3f]",
                                     low, upp))

## inla sparse

source("code/utils/summary_inla.R")

parest_sparse <-
  res_sparse |>
  inla.spde2.result(name = "s", spde = spde_sparse) |>
  summary_parest_spde() |>
  transform(ci = sprintf("[%.3f; %.3f]", low, upp)) |>
  dplyr::select(-low, -upp)

beta_sparse <-
  res_sparse$summary.fixed[, c("mean", "0.025quant", "0.975quant")] |>
  as.data.frame() |>
  dplyr::mutate(parameter = "beta", .before = mean) |>
  dplyr::mutate(ci = sprintf("[%.2f; %.2f]",
                             `0.025quant`,
                             `0.975quant`)) |>
  dplyr::select(- `0.025quant`,
                - `0.975quant`)

tau_sparse <-
  res_sparse$marginals.hyperpar[[1]] |>
  inla.rmarginal(2000, marginal = _)

tau_sparse <-  1 / sqrt(tau_sparse)

tau_sparse <- data.frame(
    parameter = "tau",
    mean      = mean(tau_sparse),
    ci        = sprintf("[%.3f; %.3f]",
                        quantile(tau_sparse,
                                 probs = .025),
                        quantile(tau_sparse,
                                 probs = .975))
)

parest_sparse <- rbind(beta_sparse, parest_sparse,
                       tau_sparse) |>
  dplyr::rename("expected" = "mean")

## inla fine

parest_fine <-
  res_fine |>
  inla.spde2.result(name = "s", spde = spde_fine) |>
  summary_parest_spde() |>
  transform(ci = sprintf("[%.3f; %.3f]", low, upp)) |>
  dplyr::select(-low, -upp)

beta_fine <-
  res_fine$summary.fixed[, c("mean", "0.025quant", "0.975quant")] |>
  as.data.frame() |>
  dplyr::mutate(parameter = "beta", .before = mean) |>
  dplyr::mutate(ci = sprintf("[%.3f; %.3f]",
                              `0.025quant`,
                              `0.975quant`)) |>
  dplyr::select(- `0.025quant`,
                - `0.975quant`)

tau_fine <-
  res_fine$marginals.hyperpar[[1]] |>
  inla.rmarginal(2000, marginal = _)

tau_fine <-  1 / sqrt(tau_fine)

tau_fine <- data.frame(
    parameter = "tau",
    mean      = mean(tau_fine, na.rm = TRUE),
    ci        = sprintf("[%.3f; %.3f]",
                        quantile(tau_fine,
                                 na.rm = TRUE,
                                 probs = .025),
                        quantile(tau_fine,
                                 na.rm = TRUE,
                                 probs = .975))
)

parest_fine <- rbind(beta_fine, parest_fine,
                     tau_fine) |>
  dplyr::rename("expected" = "mean")

## gathering results

parest_fine <- parest_fine |>
  dplyr::mutate(ci = gsub("\\[", "\\(", ci)) |>
  dplyr::mutate(ci = gsub("\\]", "\\)", ci))
parest_sparse <- parest_sparse |>
  dplyr::mutate(ci = gsub("\\[", "\\(", ci)) |>
  dplyr::mutate(ci = gsub("\\]", "\\)", ci))
final <-
  parest_hgp |>
  dplyr::mutate(ci = gsub("\\[", "\\(", ci)) |>
  dplyr::mutate(ci = gsub("\\]", "\\)", ci)) |>
  dplyr::filter(!grepl("alpha", parameter)) |>
  dplyr::mutate(parameter = stringr::str_extract(parameter,
                                                 "\\w*"))

##--- Table 7 ----

final <-
  final |>
  dplyr::select(-low, -upp) |>
  merge(parest_sparse, by = "parameter",
        all.x = TRUE) |>
  merge(parest_fine, by = "parameter",
        all.x = TRUE) |>
  dplyr::mutate(parameter = sprintf("$\\%s$", parameter))

##--- table 3 ----

final |>
  print(digits = 3) |>
  knitr::kable(format = "latex",
               booktabs = TRUE,
               escape = FALSE,
               digits = 3,
               linesep = "")

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
  filter(grepl("HGP|AGP", family)) |>
  mutate(vart = ifelse(grepl("AGP", family), "Het", vart)) |>
  filter(vart == "Het") |>
  mutate(family = ifelse(grepl("HGP", family), "HGP", family)) |>
  ## filter(model_id %% 2 == 0) |>
  mutate(bias = median_prd - obs,
         is   = int_score(y = obs, l = lower, u = upper,
                          alpha = .05),
         cvg = between(obs, lower, upper)) |>
  group_by(family, k) |>
  summarise(bias = mean(bias),
            rmse = 10 * sqrt(mean(bias * bias)),
            is   = mean(is),
            cvg  = 100 * mean(cvg)) |>
  ungroup() |>
  group_by(family) |>
  summarise(across(where(is.numeric), mean)) |>
  ungroup() |>
  arrange(family) |>
  print(digis = 3) 
