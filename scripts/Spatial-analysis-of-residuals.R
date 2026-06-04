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
library(patchwork)
library(tidybayes)

#### 1) Set up data ####

# Specify data to work with for this example
# Define "dat" as the data to work with for the rest of the script 
#   to avoid having to change the name for new data 

data(meuse) # data on soils analysis near Meuse River in Netherlands

head(meuse)

# Quick plots of the data 

meuse_sf = st_as_sf(meuse, coords = c("x", "y"), crs = 28992, agr = "constant")

# Soil zinc concentrations vs distance to river 
ggplot(meuse_sf, aes(size = zinc, color = elev)) + 
  geom_sf() + 
  theme_bw() + 
  scale_colour_viridis_c()

cor(meuse_sf$zinc, meuse_sf$elev)

# The rest of the script uses "dat" as the data frame to analyze, and dat_sf 
    # as the spatial points object. Set those values here. 
dat_sf = meuse_sf 

# Convert spatial object back to non-spatial data frame
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
summary(m)
plot(m)

# Diagnostic plots 
gam.check(m)

# Get the model residuals and add to data frame 
dat$resid <- residuals(m, type = "response") 

#### 3) Visualize model residuals and make empirical variogram ####
ggplot(dat, aes(x = x, y = y, size = abs(resid), color = resid > 0)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("red", "blue"),
                     labels = c("Negative", "Positive")) + 
  theme_bw()

# Make an empirical variogram 
coordinates(dat) <- c("x", "y") # need to put spatial info back in... 
vgm_emp <- variogram(resid ~ 1, data = dat)
vgm_emp
plot(vgm_emp)


#### 4) Make a Moran's I correlogram for the residuals ####

moran.correl <- correlog(
  x         = dat_df$x,
  y         = dat_df$y,
  z         = dat_df$resid,
  increment = 200,    # bin width in data coordinate units (meters for Meuse)
  resamp    = 999    # reps for permutation test
)

plot(moran.correl) # Filled points are significantly different from zero 
                    # based on permutation test. 




#### 5) Fit a Gaussian process spatial model to the residuals ####

zinc_model_brm <- brm(scale(zinc) ~ scale(elev) + gp(x, y, scale = FALSE), 
                        data = dat, 
                        chains = 4, 
                        cores = 4)

summary(zinc_model_brm)
plot(zinc_model_brm)
# sdgp parameter is the spatial variance (larger means more important)
# lscale parameter is the approximate spatial range of the 
#   spatial autocorrelation patterns. 
#   I think we could use 
#     these outputs as a measure of the scale of spatial synchrony as well. 

# Compare to gam of the same thing 

dat = as.data.frame(dat) |> 
  mutate(zinc_scaled = scale(zinc),
         elev_scaled = scale(elev))

zinc_model_gam = gam(zinc_scaled ~ elev_scaled + s(x, y, k=10), data = dat)
summary(zinc_model_gam)

gam_terms <- predict(zinc_model_gam, newdata = dat, type = "terms")
dat$gam_spatial_effect <- gam_terms[, "s(x,y)"]
gam_se_terms <- predict(zinc_model_gam, newdata = dat, type = "terms", se.fit = TRUE)
dat$gam_spatial_se <- gam_se_terms$se.fit[, "s(x,y)"]
gam_mean_plot <- ggplot(as.data.frame(dat), aes(x = x, y = y, color = gam_spatial_effect)) + 
  geom_point(size = 2.5) + 
  scale_color_distiller(palette = "RdBu", name = "GAM Smooth\nEffect") + 
  theme_bw() +
  labs(title = "Isolated GAM s(x,y) Effect")
gam_se_plot <- ggplot(as.data.frame(dat), aes(x = x, y = y, color = gam_spatial_se)) + 
  geom_point(size = 2.5) + 
  scale_color_viridis_c(option = "inferno", name = "GAM Spatial\nSE") + 
  theme_bw() +labs(title = "GAM Spatial Effect Uncertainty")
plot(gam_mean_plot + gam_se_plot)  



## Plot mean predictions and uncertainty at the point locations 
epred_draws <- posterior_epred(
  zinc_model_brm, 
  newdata = dat, 
  ndraws = 100 # Keep this number pretty small or it gets slow
)
dat$mean_estimate <- colMeans(epred_draws)
dat$sd_estimate   <- apply(epred_draws, 2, sd)

mean_plot = ggplot(as.data.frame(dat), aes(x=x, y=y, color = mean_estimate)) + 
  geom_point() + 
  scale_colour_continuous(palette = "YlOrBr") + 
  theme_bw() + 
  labs(title = "Model predictions")
sd_plot = ggplot(as.data.frame(dat), aes(x=x, y=y, color = sd_estimate)) + 
  geom_point() + 
  scale_colour_continuous(palette = "YlOrBr") + 
  theme_bw() + 
  labs(title = "Standard deviation of Model predictions")
plot(mean_plot + sd_plot)


## Let's plot JUST The spatial random effects (i.e. the gaussian process)
# Create new version of the data with fixed effect set to median (0) 
dat_gp <- as.data.frame(dat) |> 
  mutate(elev = 0)
# Extract posterior distributions for the grid
gp_epred_draws <- posterior_epred(
  zinc_model_brm, 
  newdata = dat_gp, 
  ndraws = 100
)

# Calculate the mean spatial effect and its uncertainty (SD)
dat$gp_mean <- scale(colMeans(gp_epred_draws), center = TRUE, scale= FALSE)
dat$gp_sd   <- apply(gp_epred_draws, 2, sd)

gp_mean_plot <- ggplot(as.data.frame(dat), aes(x = x, y = y, color = gp_mean)) + 
  geom_point(size = 2.5) + 
  scale_color_distiller(palette = "RdBu", name = "GP Mean\nAnomaly") + 
  theme_bw() +
  labs(title = "Isolated GP Spatial Effect")
gp_sd_plot <- ggplot(as.data.frame(dat), aes(x = x, y = y, color = gp_sd)) + 
  geom_point(size = 2.5) + 
  scale_color_viridis_c(option = "inferno", name = "GP Posterior\nSD") + 
  theme_bw() +
  labs(title = "Uncertainty of GP Spatial Effect")
plot(gp_mean_plot + gp_sd_plot)


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
summary(m_spatial) # elev still matters a lot
plot(m_spatial)

gam.check(m_spatial)

# Compare the same thing using brms? 
