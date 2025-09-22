library(sf)

set.seed(2024)

##--- reading and pre-processing data ----

## crs_epsg <- 26911

my_dt <- readRDS("data/pm25/dataset.rds")

dists <- my_dt |>
  st_distance(by_element = FALSE,
              which = "Hausdorff") |>
  units::set_units(value = "km") |>
  units::set_units(value = NULL) |>
  as.matrix()

N <- NROW(dists)

x <- matrix(1, nrow = N)

areas <- st_area(my_dt) |>
  units::set_units(value = "km^2") |>
  units::set_units(value = NULL) |>
  as.numeric()

areas_bin <- ifelse(areas == 0, 1, 0)

A <- cbind(1, areas_bin)

ct_dists <- my_dt |>
  st_centroid() |>
  ## st_geometry() |>
  st_distance(which = "Hausdorff") |>
  units::set_units(value = "km") |>
  units::set_units(value = NULL) |>
  as.matrix()

##--- saving quantities ----

saveRDS(cbind(my_dt$pm25, x),
        file = "data/pm25/cluster/xy.rds")

saveRDS(A, file = "data/pm25/cluster/W.rds")

saveRDS(dists, file = "data/pm25/cluster/hmat.rds")

saveRDS(ct_dists, file = "data/pm25/cluster/ct_mat.rds")
