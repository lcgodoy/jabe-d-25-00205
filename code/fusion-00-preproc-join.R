library(dplyr)
library(sf)

crs_epsg <- c(4326, 3310, 4269,
              2229, 2874, 3498,
              26911)[7]#[5]

counties <- readRDS("data/pm25/counties.rds") |>
  st_transform(crs_epsg)

## polygon used to calculate distance to ocean (as a covariate)
water_1 <-
  st_sym_difference(st_union(counties),
                    st_as_sfc(st_bbox(counties)))

shore <- st_union(counties) |>
  st_cast("POLYGON") |>
  st_cast("LINESTRING")

la <- counties |>
  subset(counties$NAME %in% c("Los Angeles",
                              "Ventura")) |>
  ## subset(counties$NAME %in% c("Ventura")) |>
  st_geometry()

border_poly <- la |>
  st_cast("POLYGON")

## ids <- border_poly |>
##   st_touches() |>
##   sapply(\ (x) length(x) > 0)

## border_poly <- border_poly |>
border_poly <- border_poly[c(1, 4)] |>
  st_union()

border_bb <- st_bbox(border_poly)

water <-
  st_intersection(water_1, st_as_sfc(border_bb))
shore <-
  st_intersection(shore, st_as_sfc(border_bb))

water |>
  plot(col = 4, border = "transparent")
plot(shore, col = 2, lwd = 2, add = TRUE)

saveRDS(border_poly, file = "data/pm25/border_counties.rds")
saveRDS(water, file = "data/pm25/ocean.rds")
## not really
saveRDS(shore, file = "data/pm25/shore.rds")

stations <- readRDS("data/pm25/stations_processed.rds") |>
  st_transform(crs_epsg)
dt_sat <- stars::read_stars("data/pm25/satellite.tif") |>
  st_as_sf() |>
  st_transform(crs_epsg)

names(dt_sat)[1] <- "pm25"

stations <- stations[st_within(stations$geometry,
                               ## st_as_sfc(border_bb),
                               border_poly,
                               sparse = FALSE), ]

dt_sat <- dt_sat[st_intersects(dt_sat$geometry,
                               ## st_as_sfc(border_bb),
                               border_poly,
                               sparse = FALSE)[, 1], ]

n_est <- nrow(dt_sat)
n_sts <- nrow(stations)

stations <- stations |>
  mutate(id = sprintf("st_%d",
                      row_number()),
         .before = "site_id") |>
  select(-site_id)

dt_sat <- dt_sat |>
  mutate(id = sprintf("es_%d",
                      row_number()),
         .before = "pm25")

plot(st_geometry(dt_sat))
plot(border_poly, add = TRUE, col = "transparent",
     border = 2, lwd = 2)
plot(water, add = TRUE, col = 4,
     border = "transparent")
plot(st_geometry(stations),
     add = TRUE, col = 1,
     pch = 19)

## creating covariates
stat_wat <- st_distance(stations, shore,
                        which = "Hausdorff") |>
  units::set_units("km") |>
  units::set_units(NULL) |>
  apply(1, min)

sat_wat <- st_distance(dt_sat, shore,
                       which = "Hausdorff") |>
  units::set_units("km") |>
  units::set_units(NULL) |>
  apply(1, min)

stations <- stations |>
  mutate(dist_wat = stat_wat)

dt_sat <- dt_sat |>
  mutate(dist_wat = sat_wat)

## elevation

elev <-
  elevatr::get_elev_raster(locations = st_as_sfc(border_bb),
                           z = 8) |>
  stars::st_as_stars()

elev_sat <-
  aggregate(elev,
            by = st_geometry(dt_sat),
            FUN = median) |>
  st_as_sf()

elev_pts <- stars::st_extract(elev, at = st_geometry(stations))

stars::write_stars(elev, "data/pm25/elevation.tif")

dt_sat <- dt_sat |>
  transform(water = sat_wat,
            elevation = elev_sat[[1]])
stations <- stations |>
  transform(water = stat_wat,
            elevation = elev_pts[[1]])

my_dt <- rbind(dt_sat, stations)

with(my_dt, hist(elevation))
with(my_dt, hist(water))

##--- creating data for predictions ---

mf <- 4 ## .5

pred_regions <- st_sample(x = border_poly,
                          size = n_est * mf,
                          type = "regular") |>
  st_set_crs(crs_epsg)

pred_regions <- st_sf(pred_regions) |>
  transform(id = seq_len(NROW(pred_regions)))

pred_regions <- stars::st_rasterize(pred_regions)
pred_regions <- st_as_sf(pred_regions)

elev_pred <- aggregate(elev, by = st_geometry(pred_regions),
                       FUN = median) |>
  st_as_sf()

pred_wat <- st_distance(elev_pred, shore,
                        which = "Hausdorff") |>
  units::set_units("km") |>
  units::set_units(NULL) |>
  apply(1, min)

pred_within_wat <- st_touches(st_geometry(pred_regions), water,
                              sparse = FALSE) |>
  as.numeric()

pred_wat <- pred_wat * (1 - pred_within_wat)

pred_regions <- pred_regions |>
  transform(water = pred_wat,
            elevation = elev_pred[[1]])

pred_regions <- pred_regions |>
  filter(id > 0)

saveRDS(my_dt, file = "data/pm25/dataset.rds")
saveRDS(pred_regions, file = "data/pm25/pred-locations.rds")
