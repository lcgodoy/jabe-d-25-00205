library(ggplot2)
library(sf)
library(stars)
library(tigris)
library(dplyr)

xx <- read_stars("data/pm25/geotiff/PM25_201001_201212.tif")

us_all <- states()

us_all <- us_all |>
    filter(STATEFP < "60") |> # territories out
    filter(!STATEFP %in% c("02", "15"))

us_all <- us_all |>
    st_geometry() |>
    st_union()

## raster projection
my_crs <- 4326
us_all <- st_transform(us_all, my_crs)

xx <- st_transform(xx, my_crs)

xx <- xx[st_as_sfc(st_bbox(us_all))]

write_stars(xx, "data/pm25/us.tif")

all_counties <- counties(year = 2015, class = "sf") |>
    filter(STATEFP < "60") |> # territories out
    filter(!STATEFP %in% c("02", "15")) |> # AK & HI out
    select(GEOID, STATEFP, geometry) |>
    mutate(GEOID = as.numeric(GEOID)) |>
    st_transform(my_crs)

##--- getting raster data ----

cal_rst <- xx[st_as_sfc(st_bbox(all_counties))]

write_stars(cal_rst, "data/pm25/satellite.tif")
