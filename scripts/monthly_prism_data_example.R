# Script to demonstrate loading and summarizing daily data layers from gridmet

library(tidyverse)
library(terra)
library(sf)
library(downloader)
library(stringr)

#### Constants

# Path to the directory containing the monthly climate data layers
CLIM_DIR = "~/Library/CloudStorage/Box-Box/climate-data-workshop/gridmet_data_daily"

PLSS_DIR = "~/Library/CloudStorage/Box-Box/climate-data-workshop/plss_data"


#### Workflow


##### Load the climate data layers

# Many climate layers, including daily gridmet data, are normally stored in compact nbinary format called netCDF 

# Download netcdf file from the THREDDs climate data server
# see http://thredds.northwestknowledge.net:8080/thredds/reacch_climate_MET_catalog.html 

# First we might need to set the "timeout option" to a longer time to make sure that the download function doesn't just give up before getting the whole file 
options(timeout = max(300, getOption("timeout"))) 

# Download a netCDF file for daily temperature variables for one year of data
variable_to_get <- "tmmx" # daily maximum temperature
year_to_get <- 2023
downloader::download(url = str_c("http://www.northwestknowledge.net/metdata/data/", variable_to_get, "_", year_to_get, ".nc"), destfile = paste0(CLIM_DIR, "/", variable_to_get, "_", year_to_get, ".nc"), mode = "wb")

# Just for reference/comparison, here is an example of getting one netcdf file of nCLIMGRID data. These come in monthly blocks with temp and precip variables combined, so we'll need a different workflow to handle those. 
#downloader::download(url = "https://www.ncei.noaa.gov/pub/data/daily-grids/v1-0-0/grids/2023/ncdd-202301-grd-scaled.nc", destfile ="./ncdd-202301-grd-scaled.nc", mode = "wb")

variable_to_get <- "tmmn" # daily minimum temperature
year_to_get <- 2023
downloader::download(url = str_c("http://www.northwestknowledge.net/metdata/data/", variable_to_get, "_", year_to_get, ".nc"), destfile = paste0(CLIM_DIR, "/", variable_to_get, "_", year_to_get, ".nc"), mode = "wb")

# Load a netCDF file of daily data for one year as one big, multilayer raster
list.files(CLIM_DIR)
r_tmax <- rast(paste0(CLIM_DIR, "/tmmx_2023.nc"))
dim(r_tmax)

# Plot the first few layers of the raster
plot(r_tmax[[1:9]])

# What is the name of the layers?
names(r_tmax)[1:9]

# Let's get the daily mean temperature (approximated as the average of the min and max daily temperatures)
r_tmin <- rast(paste0(CLIM_DIR, "/tmmn_2023.nc"))
r <- (r_tmax + r_tmin) / 2

# Get the average temperature for time period covered by the data set 
r_normal <- mean(r)
plot(r_normal)

# Calculate the daily anomalies
r_anomaly <- r - r_normal 
plot(r_anomaly[[180]]) # plot an early summer day 

names(r)[1:9]
# Hmm noticing some unfamiliar units for days and degrees here. We'll deal with that later!


#### Now we can extract data from the raster using spatial features (points, polygons)

# One source of data might be a field GPS, like Derek's example last week. 

# Another source is polygons from a map such as fire perimeters or PLSS data. Using PLSS data as an example, let's try using centroids of spatial feature polygons to extract data.   

# Load PLSS grid for CA 
plss_grid <- st_read(paste0(PLSS_DIR, "/plss_grid.gpkg"))

# Check out some random set of grid cells 
ggplot(plss_grid[1:2000,]) + geom_sf() + theme_minimal()

head(plss_grid)

# Scenario 1: we have a few specific TRS grid cells in mind to look at. 
# Subset this grid -- for now using some arbitrary locations 
plss_grid_subset = plss_grid |> dplyr::filter(MTRS %in% c("MDM-T15N-R10E-15", "MDM-T18N-R12E-15", "MDM-T20N-R10E-15", "MDM-T22N-R15E-15"))

ggplot(plss_grid_subset) + geom_sf() + theme_minimal()

# Scenario 2: We could also subset the grid based on a region like a single seed zone. 

# Load the USFS seed zones layer
seed_zones <- st_read(paste(PLSS_DIR, "seed_zones.gpkg", sep = "/"))
plot(seed_zones)

# What if we don't know the projection that the seed zones data are in? We don't actually need to know to reproject it. But if we want to check: 
st_crs(seed_zones)

# If we had a spatial dataset that is missing its CRS we can try to find it out 
library(crsuggest)
seed_zone_783 <- seed_zones[seed_zones$SEED_ZONE == 783,]
st_crs(seed_zone_783) = NA
guess_crs(seed_zone_783, target_location = c(-118.990280, 37.592576))
view_crs(3309)

# Transform the seed zones to the same projection as the plss grid 
seed_zones = st_transform(seed_zones, st_crs(plss_grid))

# Take a look at one or more seed zones 
seed_zone_531 <- seed_zones |> 
  filter(SEED_ZONE == 531)
ggplot(seed_zone_531) + geom_sf() + theme_minimal()

# Now grab the PLSS grid cells that intersect with this seed zone
plss_grid_subset = st_intersection(plss_grid, seed_zone_531)
ggplot(plss_grid_subset) + geom_sf() + theme_minimal()

# Get the centroids of our subset of TRS grid squares
seed_zone_531_plss_centroids = plss_grid_subset |> st_centroid()

# Reproject the centroids to the CRS of the climate data 
seed_zone_531_plss_centroids = st_transform(seed_zone_531_plss_centroids, st_crs(r))

# Extract the value of the climate raster at the centroid of the PLSS grid
extracted = terra::extract(r, seed_zone_531_plss_centroids, method = "bilinear")

# Add the site/grid information to the extracted data to keep it all together
d = cbind(seed_zone_531_plss_centroids, extracted)

# Now we have the climate data for our plots in a tabular format! Now let's do some tabular data wrangling.

names(d) # data are very wide because there's a column for each day! 

# Convert the data from wide to long data format 
d_long = pivot_longer(d, cols = starts_with("air_"), names_to = "climate_layer", values_to = "value")
names(d_long)

# Check out the "climate_layer" variable 
d_long$climate_layer[1:10]

# Split this column into info about the variable (air_temperature_day) and the date (days since 1/1/1900)
d_long = d_long  |>
  separate_wider_delim(climate_layer, 
    names = c("clim_var", "days_since_1900"), delim = ".")

# check
d_long |> select(clim_var, days_since_1900) |> head()

# Note that because this data set has a row for each day for each location, it's pretty big!
dim(d_long)

# Convert the "days_since_1900" date info and degrees K data to more understandable formats
d_long = d_long  |>
  mutate(date = lubridate::ymd("1900-01-01") + as.numeric(days_since_1900),
        tmax_c = value - 273.15) |> 
  mutate(year = year(date), month = month(date), day = day(date), 
         julian_date = lubridate::yday(date))

# check result
d_long |> select(days_since_1900, date, year, month, day, julian_date) |> head()

# Plot the time series for a few locations
ggplot(d_long |> filter(ID < 21), aes(x = date, y = tmax_c, group = ID, color = ID)) +
  geom_line() + labs(title = "Time series of tmax", x = "Date", y = "Temperature (C)") +
  theme_minimal()


#### Using daily time series data to make annual weather summaries or indices 

# We can define any kind of function we want!

# Example: growing degree days
calculate_gdd <- function(temp, base_temp = 10) {
  if (is.null(temp)) stop("No temperature data provided")
  gdd = sum(pmax(temp - base_temp, 0))
  return(gdd)
}

# Use the climate index function to summarize by year and point location
d_annual <- d_long |> 
  group_by(year, ID) |> 
  summarize(gdd = calculate_gdd(tmax_c, 10)) |>
  cbind(seed_zone_531_plss_centroids)
  
# To visualize the results, plot growing degree days for the extracted points as a map
ggplot(d_annual) + 
  geom_sf(aes(geometry = geom, color = gdd)) + 
  scale_color_viridis_c(option = "plasma") + 
  labs(title = "Growing degree days", 
       fill = "GDD", 
       x = "Longitude (decimal degrees E)", 
       y = "Latitude (decimal degrees N)") + 
  theme_minimal()


#### Next steps: we could automate downloading and extracting from netCDF files for many different years at once

GRIDMET_URL = "http://www.northwestknowledge.net/metdata/data/"

# Function to get one netCDF file for a given variable and year
get_gridmet_file <- function(variable_to_get, year_to_get, dest_dir) {
  downloader::download(url = str_c("http://www.northwestknowledge.net/metdata/data/", variable_to_get, "_", year_to_get, ".nc"), destfile = paste0(dest_dir, "/", variable_to_get, "_", year_to_get, ".nc"), mode = "wb")
}

# Use the function to get 10 years of data for ppt 
for (year in 2015:2024) {
  get_gridmet_file("pr", year, CLIM_DIR)
}

# Put the 10 years of data together 
raster_files_list = list.files(path = CLIM_DIR, pattern = "pr_.*.nc", full.names = TRUE) # pull out names of all files containing the pattern corresponding to our focal climate variable, here "pr" for precip
raster_files_list

ppt_rasters = rast(raster_files_list) # create the stacked raster 

# Crop to the extent of CA to save space 
ppt_rasters_ca = crop(ppt_rasters, ext(-124.5, -114.1, 32.3, 42.1))
plot(ppt_rasters_ca[[1]])
dim(ppt_rasters_ca)

# Extract the data for the same points as before
extracted_ppt = terra::extract(ppt_rasters_ca, focal_centroids, method = "bilinear")
names(extracted_ppt)[1:10]

# Attach the location information to the extracted data 
d_ppt <- cbind(focal_centroids, extracted_ppt)

# Do some wrangling to put data in long format and convert dates to interpretable units
d_ppt_long = pivot_longer(d_ppt, cols = starts_with("precipitation_"), names_to = "climate_layer", values_to = "value") |> 
  separate_wider_delim(climate_layer, names = c("clim_var", "days_since_1900"), delim = ".") |> 
  mutate(date = lubridate::ymd("1900-01-01") + as.numeric(days_since_1900),
         ppt_mm = value) |> 
  mutate(year = year(date), month = month(date), day = day(date), 
         julian_date = lubridate::yday(date)) 

# Plot an example time series for one location
ggplot(d_ppt_long |> filter(ID == 1000), aes(x = date, y = ppt_mm)) + geom_line(color = "slateblue") + labs(title = "Time series of ppt", x = "Date", y = "Precipitation (mm)") + theme_minimal()
