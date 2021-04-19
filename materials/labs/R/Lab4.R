# LOAD REQUIRED PACKAGES AND FUNCTIONS -----------------------------------------
if (!require("pacman")) install.packages("pacman")
pkgs = c("sf", "PrevMap", "raster", "mapview", "dplyr", "tmap",
         "ggplot2") # package names
pacman::p_load(pkgs, character.only = T)

source("functions.R")

# LOAD DATA --------------------------------------------------------------------

# Load river blindness data
remo <- read.csv("data/LiberiaRemoData.csv")

# Load boundaries for Liberia
# To load vector (spatial) data we can use the st_read function from the sf
# paclage. You can find a geopackage (.gpkg file) for Liberia in the data folder.
# It contains the boundaries of the country at different administrative 
# level. We will upload layer 0 by that has the outline of the country. 
# You can check all the layers contained in a geopackage with the st_layers 
# function.
liberia <- st_read("data/gadm36_LBR.gpkg", layer = "gadm36_LBR_0")


# DATA PROCESSING --------------------------------------------------------------

# We add two extra columns, one with the prevalence and one with
# the empirical logit of prevalence
remo$prev <- remo$npos / remo$ntest
remo$elogit_prev <- log((remo$npos + 0.5) / (remo$ntest - remo$npos + 0.5))

# We now need to convert the spatial coordinates from lonlat to UTM

# Create an sf (spatial) object and use the
# st_transform function to convert the crs from 
# WGS84 long/lat (epsg = 4326) to web-mercator (epsg = 32629).
# The epsgKM function translates the coordinates
# into kilometers rather than meters (the default).
remo_sf <- remo %>% 
  st_as_sf(coords = c("long", "lat"), crs = 4326) %>% 
  st_transform(crs = epsgKM(32629))

# Extract the coordinates in UTM km from the sf object with the 
# st_coordinates function and add them to the orginal dataset
remo[c("utm_x", "utm_y")] <- st_coordinates(remo_sf)

# We also convert the boundaries for Liberia to the same CRS of the point data
liberia <- st_transform(liberia, crs = st_crs(remo_sf))

# VISUALISATION ----------------------------------------------------------------

# To have a quick look at the data you can us the mapview function from
# the mapview package

# Locations only
mapview(remo_sf)

# Color and size of the point depends on a variable
mapview(remo_sf, zcol = "prev", cex = "prev")

# RELATIONSHIP WITH ELEVATION --------------------------------------------------

# For the geostatistical binomial model with assume a linear relationship
# between prevalence on the log-odds scale and the covariates. 

# Let's see what this relationship looks like for elevation and log-elevation
ggplot(remo, aes(x = elevation, y = elogit_prev)) +
  geom_point() +
  geom_smooth(method = "gam")

ggplot(remo, aes(x = log_elevation, y = elogit_prev)) +
  geom_point() +
  geom_smooth(method = "gam")
 
# GEOSTATISTICAL BINOMIAL MODEL ------------------------------------------------

# Fit a geostatistical binomial model with log-elevation and a east-west trend
# as covariate, exponential correlation function and no nugget

# We first need to define options for the MCMC algorithm used in the Monte Carlo
# maximum likelihood method (number of samples, burnin, thinning)
cmcmc <- control.mcmc.MCML(n.sim = 10000, burnin = 2000, thin = 8)

# The binomial.logistic.MCML function from the PrevMap package allows us to 
# fit a binomial geostatistical model through Monte Carlo maximul likelihood
# estimation. Have a look at the help file of the function to see how each
# argument works help("binomial.logistic.MCML"). We are going to assume
# an exponential correlation function.

fit <- binomial.logistic.MCML(formula = npos ~ log_elevation + utm_x,
                              units.m = ~ ntest,
                              coords = ~ utm_x + utm_y, 
                              data = remo, 
                              par0 = c(-2.6, 0.3, 0, 1, 20), 
                              start.cov.pars = 20,
                              fixed.rel.nugget = 0,
                              kappa = 0.5, 
                              control.mcmc = cmcmc)

# Look at a summary of the fitted model (the argument log = F will show the )
summary(fit, log = F)

# PEDICTIONS -------------------------------------------------------------------

# Create a grid of 5km for predictions tha covers the whole country
grid_pred_sf <- st_make_grid(liberia, cellsize = 5, what = "centers")

plot(liberia$geom)
plot(grid_pred_sf, cex = .001, add = T, col = "red")

# We need to exclude the prediction locations outside of our study area.
# To do this we compare our grid with the boundaries and keep only the
# points that lies within our boundaries.
grid_pred_sf <- grid_pred_sf %>%
  st_as_sf() %>% 
  st_join(liberia, left = F)

plot(liberia$geom)
plot(grid_pred_sf, cex = .001, add = T, col = "red")

# For each prediction locations we now need to extractthe covariates used in 
# the model (log-elevation and x-coordinate).

# A raster called elevLiberia.tif is available in the data folder and can be used
# to extract elevation at the prediction locations
elev <- raster("data/elevLiberia.tif")
plot(elev)

# Extract elevation and convert to log
log_elevation <- log(extract(elev, grid_pred_sf) + 0.001)

# Create data.frame with predictors (be sure that you use the names as the ones
# specified in the dataset used to fit the model)
predictors <- data.frame(log_elevation, utm_x = st_coordinates(grid_pred_sf)[, 1])

# We extract the prediction locations from the spatial object and convert
# them to a matrix
grid_pred <- st_coordinates(grid_pred_sf)

# Predictions for a binomial geostatistical model can be computed with the
# the spatial.pred.binomial.MCML function. We can also specify for what
# threshold of prevalence we want to calculate exceedance probabilities, we
# are going to use a threshold of 20%.
predictions <- spatial.pred.binomial.MCML(object = fit, 
                                          grid.pred = grid_pred, 
                                          predictors = predictors, 
                                          control.mcmc = cmcmc, 
                                          quantiles = c(0.025, 0.975), 
                                          standard.errors = T, 
                                          thresholds = 0.2,
                                          scale.thresholds = "prevalence")

# We will now extract the predicted mean, standard errors and exceedance probs
# and put them  in a raster for plotting. To create a raster we need to feed to the 
# rasterFromXYZ functiona data.frame with the coordinates in the first two 
# columns and any other variable of  interest in the following columns.
pred_summary <- rasterFromXYZ(data.frame(predictions$grid,
                                         mean_prev = predictions$prevalence$predictions,
                                         se = predictions$prevalence$standard.errors,
                                         ex_probs = predictions$exceedance.prob), 
                              crs = crs(remo_sf))

# Plot the results
plot(pred_summary)

# EXERCISES --------------------------------------------------------------------

# 1. The predictions contain an object called "samples" with the extracted
#    samples from the marginal predictive distribution at each prediction location.
#    Using this matrix try to calculate exceedance probabilities for varying
#    threhdolds of prevalences (0.3, 0.4, 0.5...)

# 2. Fit a linear geostatistical model to the empirical logit of prevalence and
#    compares the results in terms of estimated model parameters and preditctions
#    we what you obtained with the binomial model.