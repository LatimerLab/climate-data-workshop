### This script draft, which Andrew made by adapting, modifying, 
#     testing, and commenting code generated with Claude Sonnet 4.6. 
#     This demos several ways to explore spatial structure in the residuals
#     of a model including: 
#     - plot empirical variogram
#     - visualize pattern in residuals 
#     - get and plot Moran's I correlogram
#     - use Gaussian Process in brms to characterize spatial pattern 

library(sp)
library(gstat)
library(mgcv)
library(spdep)
library(tidyverse)
library(ncf)
library(brms)

#### 1) Set up data ####

# Specify data to work with for this example
# Define "dat" as the data to work with for the rest of the script 
#   to avoid having to change the name for new data 

data(meuse) # data on soils analysis near Meuse River in Netherlands
meuse_sf = st_as_sf(meuse, coords = c("x", "y"), crs = 28992, agr = "constant")

dat_sf = meuse_sf
dat = dat_sf |> 
  as.data.frame() |> 
  select(-geometry)
coords = st_coordinates(dat_sf)
dat$x = coords[,1]
dat$y = coords[,2]


#### 2) Fit a gam as a test run for looking at model residuals ####

# Remove NAs 
dat = dat[complete.cases(dat), ]

# Fit the gam model 
m <- gam(log(zinc) ~ s(elev) + landuse, data = dat) 

# Get the model residuals and add to data frame 
dat$resid <- residuals(m, type = "response") 
      # Or could use type = "response" to get on response scale

# Give coordinates to "dat" as a spatial object
coordinates(dat) <- ~x + y  # Note this data set uses x and y in meters E&N


#### 3) Visualize model residuals and make empirical variogram ####

# Visualize / map residuals by size and color
dat_df <- as.data.frame(dat)
dat_df$resid <- resid(m)
ggplot(dat_df, aes(x = x, y = y, size = abs(resid), color = resid > 0)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("red", "blue"),
                     labels = c("Negative", "Positive"))

# Make an empirical variogram 
vgm_emp <- variogram(resid ~ 1, data = dat)
vgm_emp
plot(vgm_emp)


#### 4) Make a Moran's I correlogram for the residuals ####

moran.correl <- correlog(
  x         = dat_df$x,
  y         = dat_df$y,
  z         = dat_df$resid,
  increment = 100,    # bin width in data coordinate units (meters for Meuse)
  resamp    = 999    # reps for permutation test
)

plot(moran.correl) # Filled points are significantly different from zero 
                    # based on permutation test. 

#### 5) Fit a Gaussian process spatial model to the residuals ####

resids_model_brm <- brm(resid ~ gp(x, y, scale = FALSE), 
                        data = dat_df, 
                        chains = 4, 
                        cores = 4)

summary(resids_model_brm)
plot(resids_model_brm)
# sdgp parameter is the spatial variance (larger means more important)
# lscale parameter is the approximate spatial range of the 
#   spatial autocorrelation patterns. 
#   I think we could use 
#     these outputs as a measure of the scale of spatial synchrony as well. 

# Visualize the spatial surface by predicting to a grid, then interpolating 

# Make the grid
grid_resolution <- 20 # Choose how fine to make the grid 
                        # in fractions of the distance across the 
                        # current dataset 
pred_grid <- expand.grid(
  x = seq(min(dat_df$x), max(dat_df$x), length.out = grid_resolution),
  y = seq(min(dat_df$y), max(dat_df$y), length.out = grid_resolution)
)
# Predict to points in the grid using posterior_epred()
epred_draws <- posterior_epred(
  resids_model_brm, 
  newdata = pred_grid, 
  ndraws = 100 # Keep this number pretty small or it gets slow
)
pred_grid$mean_estimate <- colMeans(epred_draws)
pred_grid$sd_estimate   <- apply(epred_draws, 2, sd)

# Plot an interpolated surface as a raster 
ggplot(pred_grid, aes(x=x, y=y, fill = mean_estimate)) + 
  geom_raster(interpolate = TRUE) + 
  theme_minimal() + 
  coord_equal() 
  

#### 6. For fun, include a spatial term in the model, then recheck 

# Using a GAM model 

m_spatial <- gam(log(zinc) ~ s(elev, k = 5) + landuse + s(x,y, k = 10), data = dat)

#plot(m_spatial)
# Note the wiggliness of the spatial smooth depends on k 
#   at k = 5 it's mostly a trend surface, at k = 10 it's hillier 

dat$resid_spatial <- residuals(m_spatial, type = "response") 
moran.correl_spatial <- correlog(
  x         = dat$x,
  y         = dat$y,
  z         = dat$resid_spatial,
  increment = 100,    # bin width in data coordinate units (meters for Meuse)
  resamp    = 999    # reps for permutation test
)
plot(moran.correl_spatial) # autocorrelation is gone with k = 10 
summary(m_spatial) # elev still matters a lot but is linear now 
plot(moran.correl)

# Compare the same thing using brms? 
