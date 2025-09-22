library(sf)
library(rstan)
library(ggplot2)
library(bayesplot)
library(INLA)

set.seed(2024)

theme_map <- theme_minimal() +
  theme(axis.text = element_blank(),
        strip.background = element_rect(color = 1))
gr <- (1 + sqrt(5)) / 2

##--- auxiliary functions ----

source("code/utils/process-chains.R")

int_score <- function(y, l, u, alpha) {
  alpha <- 2 / alpha
  ind_l <- as.numeric(y < l)
  ind_u <- as.numeric(y > u)
  (u - l) + alpha * ((l - y) * ind_l + (y - u) * ind_u)
}

##--- reading and pre-processing data ----

crs_epsg <- c(4326, 3310, 4269,
              2229, 2874, 3498,
              26911)[7]#[5]

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

ggsave(filename = "img/pm25_data.pdf",
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

## shoud be greater than 1 to predict at a finer grid
## and smaller than 1 to predicet at at more coarse grid
## mf <- 4 ## .5

## pred_regions <- st_sample(x = border_poly,
##                           size = sum(areas_bin) * mf,
##                           type = "regular") |>
##   st_set_crs(crs_epsg)

## pred_regions <- st_sf(pred_regions) |>
##   transform(id = seq_len(NROW(pred_regions)))

## pred_regions <- stars::st_rasterize(pred_regions)[border_poly, ]

## pred_regions <- st_as_sf(pred_regions)
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

## smoothness
nu <- .8

hgp_dat <- list(
    N = N,
    k = NCOL(x),
    k2 = ncol(A),
    dists = dists,
    bp = hd_lim[2] * .8,
    tau_0 = sd(my_dt$pm25) / 2,
    pt_u = .05,
    prob_rho = .01,
    y = my_dt$pm25,
    X = x,
    A = A,
    nu = nu,
    N_new = NROW(A_new),
    X_new = matrix(rep(1, pred_n), ncol = 1),
    A_new = A_new,
    dists_new = pred_dists,
    dists_cross = cross_dists
)

##----+ Compiling models +---

code_hgp <- stan_model(
    file =
      "code/fusion/stan_models/hgp-final-efficient.stan",
    allow_optimizations = TRUE
)

##----+ Sampling +---

n_iter <- 4000L
burn <- 2000L
thin <- 1

pars_hgp <- c("beta", "rho", "alpha", "tau")

fit_hgp <- sampling(
    object = code_hgp,
    data = hgp_dat,
    iter = n_iter,
    warmup = burn,
    thin = thin,
    chains = 4,
    cores = 4,
    pars = c("chol_sig", "mu_new", "sigma_new",
             "sigma_op", "z", "Sigma",
             "sigmas", "sigmas_new",
             "prec", "gi"),
    include = FALSE,
    save_warmup = FALSE,
    ## control = list(adapt_delta = .6),
    ## control = list(max_treedepth = 15),
    seed = 2024
)

print(fit_hgp, pars = pars_hgp)

pairs(fit_hgp, pars = pars_hgp)

stan_trace(fit_hgp, pars = pars_hgp)

rhat(fit_hgp, pars = pars_hgp) |>
  max()

color_scheme_set("mix-blue-red")

extract(fit_hgp, pars = pars_hgp,
        permuted = FALSE) |>
  mcmc_trace(x = _,
             regex_pars = pars_hgp,
             facet_args = list(labeller = label_parsed)) +
  legend_none()
ggsave(filename = "img/la-trace-hgp.pdf",
       width = 6,
       height = 5)

extract(fit_hgp, pars = pars_hgp,
        permuted = FALSE) |>
  mcmc_dens_overlay(x = _,
                    regex_pars = pars_hgp,
                    facet_args = list(labeller = label_parsed)) +
  legend_none()
ggsave(filename = "img/la-dens-hgp.pdf",
       width = 6,
       height = 5)

xx <- as.matrix(fit_hgp, pars = pars_hgp)

##---+ Goodness-of-fit +----

loo_hgp <- loo::loo(fit_hgp,
                    ## moment_match = TRUE,
                    save_psis = TRUE,
                    cores = 8)
ll_hgp <- loo::extract_log_lik(fit_hgp)
waic_hgp <- loo::waic(ll_hgp)
dic_hgp <- marg_dic(stan_fit = fit_hgp,
                    y = hgp_dat$y,
                    x = hgp_dat$X,
                    dists = hgp_dat$dists,
                    a = hgp_dat$A,
                    nu = hgp_dat$nu)

gof_hgp <- c("LOOIC" = loo_hgp$estimates[3, 1],
             "WAIC" = waic_hgp$estimates[3, 1],
             "DIC"  = dic_hgp)

##---+ Competing method +----

source("code/utils/inla-summary.R")

.center <- function(x)
  as.numeric(scale(x, scale = FALSE))

.center_scale <- function(x)
  as.numeric(scale(x))

## inla.setOption(pardiso.license = "/opt/licenses/pardiso.lic")

stations <- my_dt[grepl("^st_", my_dt$id), ] |>
  st_transform(inlabru::fm_crs_set_lengthunit(st_crs(crs_epsg), "km"))
dt_estim <- my_dt[!grepl("^st_", my_dt$id), ] |>
  st_transform(inlabru::fm_crs_set_lengthunit(st_crs(crs_epsg), "km"))

coop <- st_coordinates(stations)

yp <- stations$pm25

spol <- dt_estim |>
  st_geometry() |>
  as("Spatial")

spol@proj4string <- CRS()

ya <- dt_estim$pm25

coopred <- pred_regions |>
  st_transform(inlabru::fm_crs_set_lengthunit(st_crs(crs_epsg), "km")) |>
  st_centroid() |>
  st_coordinates()

sppr <- pred_regions |>
  st_geometry() |>
  as("Spatial")

sppr@proj4string <- CRS()

ypred <- rep(NA_real_, nrow(coopred))

## creating meshes
grid <- "sparse"
max_edge_ct <- ifelse(grid == "fine", 5, 10)
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
max_edge_ct <- ifelse(grid == "fine", 5, 10)
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
pdf("img/meshes.pdf",
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
    "WAIC" = res_sparse$waic$waic,
    "DIC"  = res_sparse$dic$dic)

lpml_fine <- res_fine$cpo$cpo |>
  log() |>
  sum(na.rm = TRUE)
gof_fine <-
  c("LOOIC" = -2 * lpml_fine,
    "WAIC" = res_fine$waic$waic,
    "DIC"  = res_fine$dic$dic)

gof_tbl5 <-
  cbind("AGP\\textsubscript{1}" = gof_sparse,
        "AGP\\textsubscript{2}" = gof_fine,
        "HGP" = gof_hgp)

knitr::kable(gof_tbl5,
             format = "latex", booktabs = TRUE,
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

pred_hgp <- as.matrix(fit_hgp, pars = "y_new")

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
  ## scale_fill_distiller(palette = "Spectral",
  ##                      limits = range(my_dt$pm25),
  ##                      guide = guide_colorbar(expression(hat(y)),
  ##                                             title.vjust = 1,
  ##                                             ticks = FALSE)) +
  ## scale_color_distiller(palette = "Spectral",
  ##                       limits = range(my_dt$pm25),
  ##                       guide = guide_colorbar(expression(hat(y)),
  ##                                              title.vjust = 1,
  ##                                              ticks = FALSE)) +
  facet_wrap(~ model, labeller = label_parsed) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white",
                                        color = 1),
        legend.position = "bottom",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank())


ggsave(filename = "img/predictions_pm25.pdf",
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

phi <- as.matrix(fit_hgp_v1, pars = "phi")

evals <- 400
dists_ev <- seq(from = 1e-16, to = max(dists) * 1.2,
                length.out = evals)
covs <- matrix(nrow = NROW(phi), ncol = length(dists_ev))

for (i in seq_along(phi)) {
  covs[i, ] <-
    smile:::pexp_cov(matrix(dists_ev, nrow = 1),
                     1, phi = phi[i],
                     nu = hgp_dat$nu)
}

plot(x = 0, y = 0, type = "n",
     ylim = c(0, 1),
     xlim = c(0, max(dists_ev)),
     ylab = expression(r(d)),
     xlab = expression(d))
apply(covs, 1,
      \(x) {
        lines(dists_ev, x, type = "l",
              col = scales::alpha(1, .01))
      })
lines(dists_ev, apply(covs, 2, mean),
      col = 2, lwd = 1.2)
rug(dists[upper.tri(dists)],
    col = scales::alpha(1, .1))
abline(h = .1, col = 4, lty = 2)

print(fit_hgp,
      pars = c("beta", "alpha", "rho", "tau"))

##--- parameter estimates ----

samples_hgp <- as.matrix(fit_hgp, pars = c("beta",
                                           "alpha",
                                           "rho",
                                           "tau"))
samples_hgp <-
  cbind(samples_hgp,
        "sigma" = exp(samples_hgp[, "alpha[1]"]),
        "sigma_ar" = exp(samples_hgp[, "alpha[1]"] + samples_hgp[, "alpha[2]"]))

parest_hgp <- my_summary(samples_hgp)

## inla sparse

parest_sparse <-
  res_sparse |>
  inla.spde2.result(name = "s", spde = spde_sparse) |>
  summary_parest_spde() |>
  transform(hpd = sprintf("[%.3f; %.3f]", low, upp)) |>
  dplyr::select(-low, -upp)

beta_sparse <-
  res_sparse$summary.fixed[, c("mean", "0.025quant", "0.975quant")] |>
  as.data.frame() |>
  dplyr::mutate(parameter = "beta", .before = mean) |>
  dplyr::mutate(hpd = sprintf("[%.2f; %.2f]",
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
    hpd       = sprintf("[%.3f; %.3f]",
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
  transform(hpd = sprintf("[%.3f; %.3f]", low, upp)) |>
  dplyr::select(-low, -upp)

beta_fine <-
  res_fine$summary.fixed[, c("mean", "0.025quant", "0.975quant")] |>
  as.data.frame() |>
  dplyr::mutate(parameter = "beta", .before = mean) |>
  dplyr::mutate(hpd = sprintf("[%.3f; %.3f]",
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
    mean      = mean(tau_fine),
    hpd       = sprintf("[%.3f; %.3f]",
                        quantile(tau_fine,
                                 probs = .025),
                        quantile(tau_fine,
                                 probs = .975))
)

parest_fine <- rbind(beta_fine, parest_fine,
                       tau_fine) |>
  dplyr::rename("expected" = "mean")

## gathering results

parest_fine <- parest_fine |>
  dplyr::mutate(hpd = gsub("\\[", "\\(", hpd)) |>
  dplyr::mutate(hpd = gsub("\\]", "\\)", hpd))
parest_sparse <- parest_sparse |>
  dplyr::mutate(hpd = gsub("\\[", "\\(", hpd)) |>
  dplyr::mutate(hpd = gsub("\\]", "\\)", hpd))
final <-
  parest_hgp |>
  dplyr::mutate(hpd = gsub("\\[", "\\(", hpd)) |>
  dplyr::mutate(hpd = gsub("\\]", "\\)", hpd)) |>
  dplyr::filter(!grepl("alpha", parameter)) |>
  dplyr::mutate(parameter = stringr::str_extract(parameter,
                                                 "\\w*"))

##--- Table 7 ----

final <-
  merge(parest_sparse,
        parest_fine, by = "parameter") |>
  merge(final, by = "parameter", all = TRUE) |>
  dplyr::mutate(parameter = sprintf("$\\%s$", parameter))

knitr::kable(final,
             format = "latex",
             booktabs = TRUE,
             escape = FALSE,
             digits = 3,
             linesep = "")

##--- prior vs posterior: rho ----

## bins_sturges <-
##   function(x) diff(range(x)) / nclass.Sturges(x)

## ns <- 2000

## rho <- as.matrix(fit_hgp, pars = "rho") |>
##   as.numeric()

## rho_sparse <-
##   res_sparse |>
##   inla.spde2.result(name = "s", spde = spde_sparse)
## rho_sparse <-
##   as.data.frame(rho_sparse$marginals.range.nominal)
## colnames(rho_sparse) <- c("x", "y")

## rho_fine <-
##   res_fine |>
##   inla.spde2.result(name = "s", spde = spde_fine)
## rho_fine <-
##   as.data.frame(rho_fine$marginals.range.nominal)
## colnames(rho_fine) <- c("x", "y")

## colors <- c("#ffd166", "#ef476f",
##             "#26547c")

## ggplot() +
##   geom_area(data = rho_fine,
##             mapping = aes(x = x,
##                           y = y * 3800),
##             fill = colors[1],
##             color = "transparent",
##             alpha = .5,
##             lwd = 1.1) +
##   geom_area(data = rho_sparse,
##             mapping = aes(x = x,
##                           y = y * 100000),
##             fill = colors[2],
##             color = "transparent",
##             alpha = .5,
##             lwd = 1.1) +
##   geom_density(mapping = aes(x = rho,
##                              y = after_stat(density) * 100000),
##             fill = colors[3],
##             color = "transparent",
##             alpha = .5,
##             adjust = 2,
##             lwd = 1.1) +
##   geom_line(data = rho_fine,
##             mapping = aes(x = x,
##                           y = y * 3800,
##                           color = "AGP[2]"),
##             ## fill = colors[1],
##             ## color = colors[1],
##             ## alpha = .5,
##             lwd = 1.1) +
##   geom_line(data = rho_sparse,
##             mapping = aes(x = x,
##                           y = y * 100000,
##                           color = "AGP[1]"),
##             ## fill = colors[2],
##             ## color = colors[2],
##             ## alpha = .5,
##             lwd = 1.1) +
##   geom_line(mapping = aes(x = rho,
##                           y = after_stat(density) * 100000,
##                           color = "HGP"),
##             stat = "density",
##             ## fill = colors[3],
##             ## color = colors[3],
##             ## alpha = .5,
##             adjust = 2,
##             lwd = 1.1) +
##   scale_y_continuous(sec.axis = sec_axis(trans = ~ .* 1.e-8)) +
##   theme_bw() +
##   scale_color_manual(name = "Model",
##                      breaks = c("AGP[1]", "AGP[2]", "HGP"),
##                      values = c("AGP[1]" = colors[2],
##                                 "AGP[2]" = colors[1],
##                                 "HGP" = colors[3])) +
##   guides(fill = "none") +
##   theme(strip.background = element_rect(fill = "white",
##                                         color = 1),
##         legend.position = c(.85, .8),
##         axis.text.y = element_blank(),
##         axis.ticks.y = element_blank()) +
##   labs(x = expression(rho),
##        y = NULL)

## ggsave(filename = "img/pm25_rho_post.pdf",
##        width = 6,
##        height = 6 / 1.618033)
## ggsave(filename = "img/pm25_rho_post.png",
##        width = 6,
##        height = 6 / 1.618033)

## ##--- priors inla ----

## param2.matern.orig(mesh = mesh_sparse)
## param2.matern.orig(mesh = mesh_fine)
