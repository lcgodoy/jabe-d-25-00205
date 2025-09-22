library(sf)
library(s2)
library(ggplot2)
library(dplyr)

## functions for calculating Hausdorff distance on the sphere
Rcpp::sourceCpp("src/haus-dist.cpp")

g <- st_as_sfc("POLYGON FULL", crs = 'EPSG:4326')
g

options(s2_oriented = TRUE) # don't change orientation from here on

## for grid
nx <- 10
ny <- 10

set.seed(2025)
my_pts <- st_sample(st_bbox(g),
                    size = .15 * nx * ny, ## less points than polys
                    type = "random",
                    great_circles = TRUE)

my_grid <- st_make_grid(g, n = c(nx, ny))

b <- st_buffer(st_as_sfc("POINT(-30 52)", crs = 'EPSG:4326'),
               9800000)

## visualizing
i <- st_intersection(b, my_grid)
i2 <- st_intersection(b, my_pts)
plot(st_transform(i, "+proj=ortho +lat_0=52 +lon_0=-30"),
     col = "steelblue")
plot(st_transform(i2, "+proj=ortho +lat_0=52 +lon_0=-30"),
     col = "red",
     add = TRUE, pch = 19)

##--- grid to list (easier to work on Rcpp) ----

grid_list <- lapply(my_grid,
                    \(x) {
                      out <- st_coordinates(x)[, 1:2]
                      out[-NROW(out), ]
                    })

##--- building distance matrices ----

haus_grid <- dist_haus(grid_list)
coords_pts <- st_coordinates(my_pts)
haus_point <- crossdist(st_coordinates(my_pts),
                        st_coordinates(my_pts))
haus_cross <- cdist_haus_mp(grid_list, coords_pts)

full_haus <- rbind(
    cbind(haus_grid, haus_cross),
    cbind(t(haus_cross), haus_point)
)

##--- Empirical validation ----

## max at 1 (easier to validade)
full_haus <- full_haus / max(full_haus)
grid_haus <- haus_grid / max(haus_grid)

range_nu  <- c(0.01, 1)
range_rho <- c(0.0001, 1.5)

empir_val <-
  smile::sev_pexp(range_nu, range_rho,
                  grid_len = 50,
                  dmat = full_haus)

empir_val_grid <-
  smile::sev_pexp(range_nu, range_rho,
                  grid_len = 50,
                  dmat = grid_haus)

empir_val <- empir_val |>
  mutate(type = "Fused") |>
  bind_rows(mutate(empir_val_grid, type = "Grid only"))

ggplot(data = empir_val,
       aes(x = rho,
           y = nu,
           alpha = lambda >= 0,
           fill = lambda)) +
  geom_raster() +
  scale_alpha_manual(values = c(0.5, 1)) +
  scale_fill_viridis_c(expression(lambda[1]),
                       option = "H",
                       direction = 1) +
  facet_wrap(. ~ type) +
  guides(fill = guide_colorbar(ticks = FALSE,
                               reverse = TRUE,
                               barwidth = .5),
         alpha = "none") +
  theme_bw() +
  labs(x = expression(rho),
       y = expression(nu)) +
  theme(strip.background = element_rect(fill = "white"))

ggsave(file = "img/psd-sphere.pdf", width = 6, height = 3)
