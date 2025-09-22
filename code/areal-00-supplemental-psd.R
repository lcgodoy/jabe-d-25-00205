library(sf)
library(dplyr)
library(ggplot2)

epsg_proj <- 3857

## The number of sudden infant deaths in each of the counties of North Carolina,
## USA, in 1974 -- refs: `sf` and
## https://www.paulamoraga.com/book-geospatial/sec-spatialdataandCRS.html
nc <- st_read(system.file("shape/nc.shp", package = "sf"),
  quiet = TRUE
  ) |>
  st_transform(st_crs(epsg_proj))

nc_h <- st_distance(nc, which = "Hausdorff")

## Pennsylvania counties:
## https://www.paulamoraga.com/book-geospatial/sec-arealdatatheory.html
data(pennLC, package = "SpatialEpi")
penn <- st_as_sf(pennLC$spatial.polygon) |>
  st_transform(st_crs(epsg_proj))
rm(pennLC)
gc()

penn_h <- st_distance(penn, which = "Hausdorff")

## data from smile package
data(liv_msoa, package = "smile")
data(nl_ct, package = "smile")

liv_msoa <- st_transform(liv_msoa, st_crs(epsg_proj))
nl_ct    <- st_transform(nl_ct, st_crs(epsg_proj))

livm_h <- st_distance(liv_msoa, which = "Hausdorff")
nl_h   <- st_distance(nl_ct , which = "Hausdorff")

##--- bringing distances to the same scale ----

hdists <- lapply(ls(pattern = "_h$"), get)

names(hdists) <- gsub("_h$", "", ls(pattern = "_h$"))

max_dists <- lapply(hdists, max)

std_dists <-
  Map(\(dm, sc) units::drop_units(dm / sc),
      dm = hdists,
      sc = max_dists)

names(std_dists) <- names(hdists)

par(mfrow = c(2, 2))
lapply(std_dists, \(x) hist(x[upper.tri(x)]))

range_nu  <- c(0.01, 1)
range_rho <- c(0.0001, 1.5)

empir_val <-
  purrr::map(.x = std_dists,
             .f = \(x) smile::sev_pexp(range_nu, range_rho,
                                       grid_len = 50,
                                       dmat = x)) |>
  purrr::list_rbind(names_to = "Map")

ssizes <- sapply(std_dists, NROW)

##--- check maps ----

empir_val <- mutate(empir_val,
                    Mapf = factor(Map, levels = names(ssizes)[order(ssizes)]))

ggplot(data = empir_val,
       aes(x = rho,
           y = nu,
           alpha = lambda >= 0,
           fill = lambda)) +
  geom_raster() +
  geom_hline(yintercept = 1) +
  ## scale_color_manual(values = c(2, "transparent")) +
  scale_alpha_manual(values = c(0, 1)) +
  scale_fill_viridis_c(expression(lambda[1]),
                       option = "H",
                       direction = 1) +
  guides(fill = guide_colorbar(ticks = FALSE,
                               reverse = TRUE,
                               ## title.vjust = .75,
                               ## barheight = .2
                               barwidth = .5),
         alpha = "none") +
  facet_wrap(Mapf ~ ., labeller = as_labeller(toupper)) +
  theme_bw() +
  theme(
      ## legend.position = "bottom",
      strip.background = element_rect(fill = "white",
                                      color = 1)
  ) +
  labs(x = expression(rho),
       y = expression(nu))

ggsave(file = "img/psd-areal-supplementary.pdf", width = 6, height = 5)
