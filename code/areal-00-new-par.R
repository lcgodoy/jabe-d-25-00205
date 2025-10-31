library(dplyr)
library(sf)
library(ggplot2)

data(GGHB.IZ, package = "CARBayesdata")
data(respiratorydata, package = "CARBayesdata")

my_dt <- merge(x = GGHB.IZ,
               y = respiratorydata, by = "IZ",
               all.x = FALSE)
N <- nrow(my_dt)

## sf_use_s2(FALSE)

h_areal <-
  sf::st_distance(x = st_transform(sf::st_geometry(my_dt),
                                   27700),
                  by_element = FALSE,
                  which      = "Hausdorff")  |>
  units::set_units(value = "km") |>
  units::set_units(value = NULL) |>
  as.matrix()

h_fusion <- "data/pm25/cluster/hmat.rds" |>
  readRDS()

h_areal <- h_areal / max(h_areal)
h_fusion <- h_fusion / max(h_fusion)

##--- assessing psd ----

range_nu <- c(0.01, 1)
range_rho <- c(0.0001, 1.5)

psd_areal <-
  smile::sev_pexp(range_nu, range_rho,
                  grid_len = 50,
                  dmat = h_areal) |>
  mutate(type = "Original")

##--- Assessing "rho" strategy suggested by R3 ----

rho <- c(seq(from = .15, to = .95, by = .20),
         1 - .Machine$double.eps)

psd_rho <-
  lapply(rho, \(r, dt) {
    mutate(psd_areal,
           lambda = lambda * r + (1.0 - r),
           type   = sprintf("rho==%.2f", r))
  }, dt = psd_areal) |>
  bind_rows()

ggplot(data = psd_rho,
       aes(x = rho,
           y = nu,
           alpha = lambda > 0,
           fill = lambda)) +
  geom_raster() +
  scale_alpha_manual(values = c(0.2, 1)) +
  scale_fill_viridis_c(expression(lambda[1]),
                       option = "H",
                       direction = 1) +
  guides(fill = guide_colorbar(ticks = FALSE,
                               reverse = TRUE,
                               barwidth = .5),
         alpha = "none") +
  facet_wrap(type ~ ., labeller = label_parsed) +
  theme_bw() +
  theme(
      strip.background = element_rect(fill = "white",
                                      color = 1)
  ) +
  labs(x = expression(rho),
       y = expression(nu))
