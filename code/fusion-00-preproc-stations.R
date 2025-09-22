library(ggplot2)
library(data.table)

## cmds <-
##     sprintf("unzip -cq %s",
##             list.files("data/pm25/yearly", full.names = TRUE))

fls <-
    list.files("data/pm25/yearly", full.names = TRUE)

my_dt <- lapply(fls, fread)

my_dt <- rbindlist(my_dt)[STATE == "California"]

setnames(my_dt, old = names(my_dt),
         new = janitor::make_clean_names(names(my_dt)))

my_dt[, year := substr(date, nchar(date) - 3, nchar(date))]

out <- my_dt[, .(pm25 = mean(daily_mean_pm2_5_concentration,
                             na.rm = TRUE)),
             by = c("site_id", "year",
                    "county",
                    "site_latitude",
                    "site_longitude")]

fwrite(x = out, file = "data/pm25/stations2010_2012.csv")
