library(spdep)
library(sf)
library(ggplot2)
library(rstan)
library(bayesplot)
library(dplyr)

theme_set(theme_bw() +
          theme(strip.background = element_rect(fill = "white")))

rstan_options(auto_write = TRUE)

set.seed(2024)

## replace the following with the path where you'd like to have the figures
fig_dir <- "~/git-projects/hausdorff-GP/img"

##--- auxiliary functions ----

source("code/utils/aux_icar.R")
source("code/utils/aux_dagar.R")

##--- loading data ----

data(GGHB.IZ, package = "CARBayesdata")
data(respiratorydata, package = "CARBayesdata")

my_dt <- merge(x = GGHB.IZ,
               y = respiratorydata, by = "IZ",
               all.x = FALSE)

##---+ aux quantities +----

N <- nrow(my_dt)
my_x <- matrix(scale(my_dt$incomedep), ncol = 1)

##---++ HGP ++----

hdists <- st_distance(x = st_geometry(my_dt),
                     by_element = FALSE,
                     which      = "Hausdorff")  |>
  units::set_units(value = "km") |>
  units::set_units(value = NULL) |>
  as.matrix()

my_x <- matrix(scale(my_dt$incomedep), ncol = 1)

nu <- .7
rho_0 <- 0.8 * max(hdists[upper.tri(hdists)])

hgp_dat <- list(
    N = N,
    k = ncol(my_x),
    y = my_dt$observed,
    offset = my_dt$expected,
    dists = hdists,
    X = my_x,
    prob_rho = .05,
    rho_0 = rho_0,
    nu = nu
)

##---++ BYM ++----

## defining neigborhood
neigh <- poly2nb(pl = my_dt)
nb_info <- nb2WB(neigh)
W <- nb2mat(neigh, zero.policy = TRUE,
            style = "B")

plot(st_geometry(my_dt), border = "white",
     col = 1)
plot(neigh,
     st_coordinates(st_centroid(my_dt)),
     add = TRUE, pch = 19,
     col = 2,
     cex = 0.6)
points(st_coordinates(st_centroid(my_dt)),
       pch = 19, col = 2)

adj <- unlist(apply(W == 1, MARGIN = 1,
                    which))
L <- length(adj)
num <- colSums(W)

## ICAR precision matrix
Q <- matrix(0, nrow = nrow(W), ncol = ncol(W))
diag(Q) <- apply(W, 1, sum)
Q <- Q - W

inv_q <- MASS::ginv(Q)
sig_ref <- exp(- mean(log(diag(inv_q))))

aux_dat <- get_nodes(adj = adj, num = num)

bym_dat <- list(N = NROW(my_dt),
                k = NCOL(my_x),
                X = my_x,
                y = my_dt$observed,
                offset = my_dt$expected,
                sc_f = sig_ref) |>
  c(aux_dat)

## BYM2 and ICAR use exactly the same "data" object

##---++ DAGAR ++----

Minc <- nb2mat(neigh, style = "B", zero.policy = TRUE)
N_edges <- sum(Minc) / 2 ### number of edges
nei <- neighbors(Minc)
adj.ends <- adj_index(Minc)$adj_index
N_nei <- adj_index(Minc)$N_nei

dagar_dat <- list(
    N = N,
    N_edges = N_edges,
    nei = nei,
    adjacency_ends = adj.ends,
    k = ncol(my_x),
    y = my_dt$observed,
    offset = my_dt$expected,
    X = my_x,
    N_nei = N_nei
)

##---+ Compiling models +----

hgp <- stan_model(file = "stan_models/hgp-poisson.stan")

bym2 <- stan_model(file = "stan_models/areal/bym2.stan")

dagar <- stan_model(file = "stan_models/areal/dagar.stan")

##---+ Running MCMC +----

##---++ HGP ++----

n_iter <- 5000L
burn <- 1000L
thin <- 10

fit_hgp <- sampling(
    object = hgp,
    data = hgp_dat,
    iter = n_iter + burn,
    warmup = burn,
    thin = thin,
    save_warmup = FALSE,
    chains = 4,
    cores = 4,
    seed = 202506,
    refresh = 100
)

params <- c("beta_0", "beta", "rho", "sigma")

rhat(fit_hgp, pars = params) |>
  max()

print(fit_hgp, pars = params)
pairs(fit_hgp, pars = params)

color_scheme_set("mix-blue-red")

mcmc_trace(x = fit_hgp,
           regex_pars = params,
           facet_args = list(labeller = label_parsed)) +
  legend_none()
ggsave(filename = file.path(fig_dir, "scot-trace-hgp.pdf"),
       width = 6,
       height = 5)

mcmc_dens_overlay(x = fit_hgp,
                  regex_pars = params,
                  facet_args = list(labeller = label_parsed)) +
  legend_none()
ggsave(filename = file.path(fig_dir, "scot-dens-hgp.pdf"),
       width = 6,
       height = 5)

set.seed(2024)
refs_id <- sample(seq_len(nrow(my_dt)), size = 4)
refs <- sprintf("z[%d]", sort(refs_id))

rhat(fit_hgp, pars = refs) |>
  max()

mcmc_trace(x = fit_hgp,
           regex_pars = sprintf("z\\[%d\\]", sort(refs_id)),
           facet_args = list(labeller = label_parsed)) +
  legend_none()
ggsave(filename = file.path(fig_dir, "scot-trace-ref-hgp.pdf"),
       width = 6,
       height = 5)

mcmc_dens_overlay(x = fit_hgp,
                  regex_pars = sprintf("z\\[%d\\]", sort(refs_id)),
                  facet_args = list(labeller = label_parsed)) +
  legend_none()
ggsave(filename = file.path(fig_dir, "scot-dens-ref-hgp.pdf"),
       width = 6,
       height = 5)

##---++ BYM2 ++----

n_iter <- 24000
burn   <- 4000
thin   <- 20

fit_bym2 <- sampling(
    object = bym2,
    data = bym_dat,
    iter = n_iter,
    warmup = burn,
    thin = thin,
    chains = 4,
    cores = 4,
    seed = 202505,
    refresh = 1000
)

print(fit_bym2, pars = c("beta_0", "beta", "sigma", "xi"))
pairs(fit_bym2, pars = c("beta_0", "beta", "sigma", "xi"))

rhat(fit_bym2,
     pars = c("beta_0", "beta", "sigma", "xi")) |>
  max()

mcmc_trace(x = fit_bym2,
           regex_pars = c("beta_0", "beta", "sigma", "xi"),
           facet_args = list(labeller = label_parsed)) +
  legend_none()
ggsave(filename = file.path(fig_dir, "scot-trace-bym.pdf"),
       width = 6,
       height = 5)

mcmc_dens_overlay(x = fit_bym2,
                  regex_pars = c("beta_0", "beta", "sigma", "xi"),
                  facet_args = list(labeller = label_parsed)) +
  legend_none()
ggsave(filename = file.path(fig_dir, "scot-dens-bym.pdf"),
       width = 6,
       height = 5)

rhat(fit_bym2, pars = refs) |>
  max()

mcmc_trace(x = fit_bym2,
           regex_pars = sprintf("z\\[%d\\]", sort(refs_id)),
           facet_args = list(labeller = label_parsed)) +
  legend_none()
ggsave(filename = file.path(fig_dir, "scot-trace-ref-bym.pdf"),
       width = 6,
       height = 5)

mcmc_dens_overlay(x = fit_bym2,
                  regex_pars = sprintf("z\\[%d\\]", sort(refs_id)),
                  facet_args = list(labeller = label_parsed)) +
  legend_none()
ggsave(filename = file.path(fig_dir, "scot-dens-ref-bym.pdf"),
       width = 6,
       height = 5)

##---++ DAGAR ++----

fit_dagar <- sampling(
    object = dagar,
    data = dagar_dat,
    iter = n_iter,
    warmup = burn,
    thin = thin,
    chains = 4,
    cores = 4,
    seed = 202505,
    refresh = 1000
)

print(fit_dagar, pars = c("beta_0", "beta", "sigma", "psi"))
pairs(fit_dagar, pars = c("beta_0", "beta", "sigma",
                          "psi"))

rhat(fit_dagar,
     pars = c("beta_0", "beta", "sigma", "psi")) |>
  max()

mcmc_trace(x = fit_dagar,
           regex_pars = c("beta_0", "beta", "sigma", "psi"),
           facet_args = list(labeller = label_parsed)) +
  legend_none()
ggsave(filename = file.path(fig_dir, "scot-trace-dagar.pdf"),
       width = 6,
       height = 5)

mcmc_dens_overlay(x = fit_dagar,
                  regex_pars = c("beta_0", "beta", "sigma", "psi"),
                  facet_args = list(labeller = label_parsed)) +
  legend_none()
ggsave(filename = file.path(fig_dir, "scot-dens-dagar.pdf"),
       width = 6,
       height = 5)

rhat(fit_dagar, pars = refs) |>
  max()

mcmc_trace(x = fit_dagar,
           regex_pars = sprintf("z\\[%d\\]", sort(refs_id)),
           facet_args = list(labeller = label_parsed)) +
  legend_none()
ggsave(filename = file.path(fig_dir, "scot-trace-ref-dagar.pdf"),
       width = 6,
       height = 5)

mcmc_dens_overlay(x = fit_dagar,
                  regex_pars = sprintf("z\\[%d\\]", sort(refs_id)),
                  facet_args = list(labeller = label_parsed)) +
  legend_none()
ggsave(filename = file.path(fig_dir, "scot-dens-ref-dagar.pdf"),
       width = 6,
       height = 5)

##---+ Goodness-of-fit +---

##---++ HGP ++---

loo_hgp <- loo::loo(fit_hgp,
                    ## moment_match = TRUE,
                    save_psis = TRUE,
                    cores = 4)
ll_hgp <- loo::extract_log_lik(fit_hgp)
waic_hgp <- loo::waic(ll_hgp)

gof_hgp <- c("LPML" = loo_hgp$estimates[3, 1],
             "WAIC" = waic_hgp$estimates[3, 1])

##---++ ICAR and BYM ++----

loo_bym2 <- loo::loo(fit_bym2,
                     ## moment_match = TRUE,
                     save_psis = TRUE,
                     cores = 4)
ll_bym2 <- loo::extract_log_lik(fit_bym2)
waic_bym2 <- loo::waic(ll_bym2)

gof_bym2 <- c("LPML" = loo_bym2$estimates[3, 1],
              "WAIC" = waic_bym2$estimates[3, 1])

##---++ DAGAR ++----

loo_dagar <- loo::loo(fit_dagar,
                      ## moment_match = TRUE,
                      save_psis = TRUE,
                      cores = 4)
ll_dagar <- loo::extract_log_lik(fit_dagar)
waic_dagar <- loo::waic(ll_dagar)

gof_dagar <- c("LPML" = loo_dagar$estimates[3, 1],
               "WAIC" = waic_dagar$estimates[3, 1])

gof_scot <-
  cbind.data.frame("Criteria" = c("LOOIC", "WAIC"),
                   "HGP"  = gof_hgp,
                   "BYM"   = gof_bym2,
                   "DAGAR" = gof_dagar)

rownames(gof_scot) <- NULL

loos <- list("HGP" = loo::loo(fit_hgp),
             "BYM" = loo::loo(fit_bym2),
             "DAG" = loo::loo(fit_dagar))

loo::loo_compare(loos)

##--- ** Table 2 (Bottom) ----

knitr::kable(gof_scot,
             format = "latex",
             digits = c(0, rep(1, 3)),
             booktabs = TRUE,
             linesep = "")

##--- Estimation ----

parest <-
  list(
      fit_hgp |>
      as.matrix(pars = c("beta_0", "beta", "rho", "sigma")) |>
      my_summary(),
      fit_dagar |>
      as.matrix(pars = c("beta_0", "beta", "psi", "sigma")) |>
      my_summary(),
      fit_bym2 |>
      as.matrix(pars = c("beta_0", "beta", "xi", "sigma")) |>
      my_summary()
  )

parest[[2]] <- parest[[2]] |>
  mutate(parameter = ifelse(parameter == "rho", "psi", parameter))

parest[[3]] <- parest[[3]] |>
  mutate(parameter = ifelse(parameter == "rho", "zeta", parameter))

##--- ** Table 2 (top) ----

tbl_4 <- lapply(parest, \(x) {
  mutate(x, hpd = gsub("\\[", "\\(", hpd)) |>
    mutate(hpd = gsub("\\]", "\\)", hpd)) |>
    mutate(hpd = gsub(";", ",", hpd)) |>
    mutate(expected = sprintf("%.2f %s",
                              expected, hpd)) |>
    select(- hpd)
}) |>
  Reduce(f = \(x, y) full_join(x, y, by = "parameter"), x = _) |>
  mutate(parameter = sprintf("$\\%s$", parameter))

tbl_4 |>
  knitr::kable(format = "latex", booktabs = TRUE,
               escape = FALSE,
               digits = 3,
               linesep = "")

##--- ** Figure 6: rho prior and posterior compared to dists ----

rho <- as.matrix(fit_hgp, pars = "rho") |>
  as.numeric()

aux_dists <- hdists[upper.tri(hdists)]

ggplot(mapping = aes(x = hdists[upper.tri(hdists)])) +
  geom_histogram(aes(y = after_stat(density)),
                 binwidth = bins_sturges,
                 color = 1,
                 fill = "gray70") +
  stat_function(# aes(y = after_stat(y) * 5000),
      geom = "area",
      color = "transparent",
      fill = 2,
      alpha = .5,
      lwd = 1.2,
      fun = dexp,
      args = list(rate = - log(hgp_dat$prob_rho) / hgp_dat$rho_0),
      inherit.aes = FALSE) +
  stat_function(# aes(y = after_stat(y) * 5000),
      geom = "line",
      color = 2,
      lwd = 1.2,
      fun = dexp,
      args = list(rate = - log(hgp_dat$prob_rho) / hgp_dat$rho_0),
      inherit.aes = FALSE) +
  geom_density(mapping = aes(x = rho),
               fill = "steelblue",
               color = "steelblue",
               alpha = .5,
               adjust = 2,
               lwd = 1.1) +
  theme(legend.position = c(.85, .8),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  labs(x = "Hausdorff distance (in km)",
       y = "Density")

ggsave(filename = file.path(fig_dir, "scotland_rho_post.pdf"),
       width = .7 * 6.7,
       height = .7 * 6.7 / 1.618033)

##--- ** Figure S.2 ----

cr_dist <- dists[upper.tri(dists)] |>
  sort() |>
  unique()

ndists <- 200

eval_dists <- seq(from = min(cr_dist),
                  ## to = median(cr_dist),
                  to = max(cr_dist),
                  length.out = ndists)

## * HGP

phi <- as.matrix(fit_hgp, pars = "phi")[, 1]

m_pe <- matrix(nrow = length(rho), ncol = ndists)
for (i in seq_len(NROW(m_pe))) {
  m_pe[i, ] <- smile:::pexp_cov(dists = matrix(eval_dists, nrow = 1),
                                1,
                                phi = phi[i],
                                nu = hgp_dat$nu)
}

colnames(m_pe) <- eval_dists

funcs_df <-
  mat_to_df(m_pe, fname = "HGP", nd = ndists)

funcs_df <- funcs_df |>
  group_by(dist, func) |>
  summarise(med = median(corr),
            lower1 = quantile(corr, probs = .025),
            upper1 = quantile(corr, probs = .975),
            lower2 = quantile(corr, probs = .1),
            upper2 = quantile(corr, probs = .9)) |>
  ungroup()

## * DAGAR

upper_dist <- dists[upper.tri(dists)] |>
  cut(breaks = 100) |>
  sapply(get_midpoint)

post_rho <- as.matrix(fit_dagar, pars = "rho")[, 1]

## defining neigborhood
neigh    <- spdep::poly2nb(pl = st_transform(my_dt, 4326))
Minc     <- spdep::nb2mat(neigh, style = "B")
sec_neig <- spdep::nblag(neigh, maxlag = 2) |>
  lapply(spdep::nb2mat)

aux_dagar <- vector(mode = "list", length = NROW(post_rho))
corr_dagar <- array(dim = c(NROW(dists), NROW(dists),
                            NROW(post_rho)))
for (i in seq_along(aux_dagar)) {
  corr_dagar[, , i] <- comp_corr_dagar(N, Minc,
                                       post_rho[i])
  aux_dagar[[i]] <- data.frame(id = i,
                               dist = upper_dist,
                               corrs = corr_dagar[, , i][upper.tri(dists)])
}
corr_dagar <- apply(corr_dagar, c(1, 2), mean)
gc(full = TRUE)
aux_dagar <- aux_dagar |>
  bind_rows() |>
  group_by(id, dist) |>
  summarise(corr = mean(corrs))
aux_dagar <- ungroup(aux_dagar)
aux_dagar <- aux_dagar |>
  group_by(dist) |>
  mutate(transp = ecdf(corr)(corr)) |>
  ungroup() |>
  mutate(transp = if_else(transp > .5, 1 - transp,
                          transp))
funcs_df <- funcs_df |>
  bind_rows(aux_dagar |>
            group_by(dist) |>
            summarise(med = median(corr),
                      lower1 = quantile(corr, probs = .025),
                      upper1 = quantile(corr, probs = .975),
                      lower2 = quantile(corr, probs = .1),
                      upper2 = quantile(corr, probs = .9)) |>
            ungroup() |>
            mutate(func = "DAGAR"))

saveRDS(funcs_df, file = "data/scotland/corr-funs.rds")

ggplot(data = funcs_df,
       aes(x = dist)) +
  geom_ribbon(aes(ymin = lower2, ymax = upper2),
              fill = "lightsteelblue") +
  geom_ribbon(aes(ymin = lower1, ymax = upper1),
              fill = "lightsteelblue",
              alpha = .6) +
  geom_line(aes(y = med),
            color = 1,
            lwd = 1.2) +
  geom_rug(data = data.frame(dst = dists[upper.tri(dists)]),
           aes(x = dst), alpha = .01,
           inherit.aes = FALSE) +
  facet_wrap(~ func, ncol = 1) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white",
                                        color = 1),
        legend.position = "bottom") +
  labs(y = "Spatial correlation",
       x = "Hausdorff distance (in km)")

ggsave(filename = file.path(fig_dir, "scotland_corr_func.pdf"),
       width = 6,
       height = 6 / 1.618033)

##--- Smoothed map ----

SIR_hgp <-
  as.matrix(fit_hgp, pars = "SIR")

SIR_bym2 <-
  as.matrix(fit_bym2, pars = "SIR")

SIR_dagar <-
  as.matrix(fit_dagar, pars = "SIR")

my_dt2 <-
  my_dt |>
  mutate(SIR_hgp = apply(SIR_hgp, 2, median),
         SIR_bym2 = apply(SIR_bym2, 2, median),
         SIR_dagar = apply(SIR_dagar, 2, median),
         SIR = observed / expected)

##--- * Figure 5 ----

my_dt2 |>
  select(starts_with("SIR_"), SIR) |>
  tidyr::pivot_longer(1:4) |>
  mutate(name = factor(name,
                       levels = c("SIR",
                                  sprintf("SIR_%s", c("hgp",
                                                       "dagar",
                                                       "bym2"))),
                       labels = c("Observed", "HGP",
                                  "DAGAR", "BYM"))) |>
  ggplot(data = _,
         aes(fill = value)) +
  geom_sf(col = "gray70", lwd = .0005) +
  ## scale_fill_viridis_c(option = "H",
  ##                      guide = guide_colorbar("SMR",
  ##                                             title.vjust = .75,
  ##                                             ticks = FALSE,
  ##                                             barheight = .5)) +
  scale_fill_distiller(palette = "Spectral",
                       guide = guide_colorbar("SIR",
                                              title.vjust = .75,
                                              ticks = FALSE,
                                              barheight = .5)) +
  facet_wrap(~ name, labeller = label_parsed,
             ncol = 2) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white",
                                        color = 1),
        legend.position = "bottom") +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "bottom")

ggsave(filename = file.path(fig_dir, "scot-smr.pdf"),
       width = 6,
       height = 6)

my_dt2 |>
  select(starts_with("SIR_"), SIR) |>
  tidyr::pivot_longer(1:4) |>
  ggplot(data = _,
         aes(fill = value)) +
  geom_sf() +
  scale_fill_viridis_c(option = "H") +
  facet_wrap(~ name) +
  theme_linedraw()
