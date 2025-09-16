library(tidyverse)
library(data.table)

### add precipitation and temperature data to metadata file
### usage: Rscript scripts/climate_data.R <metadata file> <output file path>
### (only works for Africa right now)
### output columns: longitude, latitude, (longlat of sampling data)
### time, (time variable, e.g. year or month)
### precip_min, precip_max, precip_mean, precip_var, (precipitation statistics)
### climatology_min, climatology_max, climatology_mean, climatology_var, (temperature statistics)
### longitude_data, latitude_data (longlat of climate data)


args = commandArgs(trailingOnly=TRUE)
metadata <- args[1] # metadata file with longitude and latitude columns
outpath <- args[2] # output file path
precdata <- 'data/climate/gpcc-precipitation-1992-2022.csv'
landdata <- 'data/climate/land_cover.csv'
tempdata <- 'data/climate/tavg-climatology.csv'

find_nearest <- function(array, value) {
    array <- as.matrix(array)
    idx <- (abs(array - value))[,1] %>% which.min()
    return(array[idx])
}

summarize_loc <- function(location, df, summvar, timevar) {
    lon_value = find_nearest(df$longitude, loc[1])
    lat_value = find_nearest(df$latitude, loc[2])

    subdf <- df %>% filter((longitude == lon_value) & (latitude == lat_value)) %>%
        group_by(!!sym(timevar)) %>% summarise(min_val = min(!!sym(summvar)), max_val = max(!!sym(summvar)),
                                                mean_val = mean(!!sym(summvar)), var_val = var(!!sym(summvar)))
    subdf$longitude <- loc[1]
    subdf$latitude <- loc[2]
    subdf$longitude_data <- lon_value
    subdf$latitude_data <- lat_value
    subdf <- subdf %>% rename(!!paste0(summvar, '_min') := min_val, !!paste0(summvar, '_max') := max_val,
                              !!paste0(summvar, '_mean') := mean_val, !!paste0(summvar, '_var') := var_val,
                              !!sym(timevar) := !!sym(timevar))
    return(subdf)
}

metadata <- fread(metadata, sep = '\t')
sampleLOCs <- as.matrix(metadata[, .(longitude, latitude)])
locs <- unique(sampleLOCs, margin=1)

precipdf <- fread(precdata)
precipdf <- as.data.frame(precipdf)
landdf <- fread(landdata)
landdf <- as.data.frame(landdf)
tempdf <- fread(tempdata)
tempdf <- as.data.frame(tempdf)

loc = locs[1,]
pdf = summarize_loc(loc, precipdf, 'precip', 'time')
for (i in 2:nrow(locs)) {
    loc = locs[i,]
    subdf = summarize_loc(loc, precipdf, 'precip', 'time')
    pdf = rbind(pdf, subdf)
}

tdf = summarize_loc(locs[1,], tempdf, 'climatology', 'time')
for (i in 2:nrow(locs)) {
    loc = locs[i,]
    subdf = summarize_loc(loc, tempdf, 'climatology', 'time')
    tdf = rbind(tdf, subdf)
}

df = merge(pdf, tdf, by = c('longitude', 'latitude', 'time', 'longitude_data', 'latitude_data'))
write.csv(df, outpath, row.names = FALSE)