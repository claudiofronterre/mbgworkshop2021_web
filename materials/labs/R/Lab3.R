# LOAD REQUIRED PACKAGES AND FUNCTIONS -----------------------------------------
if (!require("pacman")) install.packages("pacman")
pkgs = c("sf", "PrevMap", "raster", "mapview", "dplyr",
         "tmap") # package names
pacman::p_load(pkgs, character.only = T)


source("R/functions.R")

# LOAD DATA --------------------------------------------------------------------

# Load Galicia dataset
lead <- read.csv("data/GaliciaData.csv")

# DATA PROCESSING --------------------------------------------------------------

# We first need to convert the spatial coordinates from lonlat to UTM

# Create an sf (spatial) object and use the
# st_transform function to convert the crs from 
# WGS84 long/lat (epsg = 4326) to web-mercator (epsg = 3857).
# The epsgKM function translates the coordinates
# into kilometers rather than meters (the default).
lead_sf <- lead %>% 
  st_as_sf(coords = c("long", "lat"), crs = 4326) %>% 
  st_transform(crs = epsgKM(3857))

# Extract the coordinates in UTM km from the sf object with the 
# st_coordinates function and add them to the orginal dataset
lead[c("utm_x", "utm_y")] <- st_coordinates(lead_sf)

# VISUALISATION ----------------------------------------------------------------

# To have a quick look at the data you can us the mapview function from
# the mapview package

# Locations only
mapview(lead_sf)

# Color and size of the point depends on a variable
mapview(lead_sf, zcol = "loglead", cex = "loglead")

# SPATIAL MODEL - INTERCEPT ONLY -----------------------------------------------

# Fit a an intercept only linear geostatistical model with no nugget and
# exponential correlation function to the log-lead concentration
fit <- linear.model.MLE(formula = loglead ~ 1, 
                        coords = ~ utm_x + utm_y, 
                        data = lead,
                        kappa = 0.5, fixed.rel.nugget = 0, 
                        start.cov.pars = 30) 

# Look at a summary of the fitted model (the argument log = F will show the )
summary(fit, log = F)

# PEDICTIONS -------------------------------------------------------------------

# Create a grid of 1km for predictions, if we are not able to retrieve the
# boundaries for Galicia we can obtain the convex hull that contains the
# observed points and use it as our study area

# Calculate convex hull (use the st_convex_hull function)
# We also extend the convex hull by 5km such that some of the observed points
# won't lie exactly at the boundaries
boundaries <- lead_sf %>% 
  st_union() %>% 
  st_convex_hull() %>% 
  st_buffer(5)

plot(boundaries)
plot(lead_sf, add = T)

# Creat a 1km grid with the st_make_grid function
grid_pred_sf <- st_make_grid(boundaries, cellsize = 1, what = "centers")

plot(boundaries)
plot(grid_pred_sf, cex = .001, add = T, col = "red")

# We need to exclude the prediction locations outside of our study area.
# To do this we compare our grid with the boundaries and keep only the
# points that lies within our boundaries.
grid_pred_sf <- grid_pred_sf %>%
  st_as_sf() %>% 
  st_join(st_as_sf(boundaries), left = F)

plot(boundaries)
plot(grid_pred_sf, cex = .001, add = T, col = "red")

# We extract the prediction locations from the spatial object and convert
# them to a matrix
grid_pred <- st_coordinates(grid_pred_sf)

# Since we are using an intercep only model we don't need to use extract any 
# covariates and we have all the ingredients needed to obtain predictions.
# Predictions for a linear geostatistical model can be computed with the
# the spatial.pred.linear.MLE function from the PrevMap package.
predictions <- spatial.pred.linear.MLE(object = fit, 
                                       grid.pred = grid_pred, 
                                       scale.predictions = "logit", 
                                       quantiles = c(0.025, 0.975), 
                                       standard.errors = T)

# We will now extract the predicted mean and standard errors and put them
# in a raster for plotting. To create a raster we need to feed to the 
# rasterFromXYZ functiona data.frame with the coordinates in the first two 
# columns and any other variable of  interest in the following columns.
pred_mean_sd <- rasterFromXYZ(data.frame(predictions$grid.pred,
                                         mean = predictions$logit$predictions,
                                         se = predictions$logit$standard.errors), 
                              crs = crs(lead_sf))

# Plot the results
plot(pred_mean_sd)

# We can obtain nicer maps with the tmap package, an example
# is provided below
tm_shape(pred_mean_sd) +
  tm_raster(col = "mean", palette = "YlGnBu", style = "cont", midpoint = NA, 
            title = "Predicted mean\nlead concentration\n(log-scale)") +
tm_shape(boundaries) +
  tm_borders(col = "black") +
tm_compass() +
tm_scale_bar(position = c("left", "bottom")) 

# EXERCISES --------------------------------------------------------------------

# 1. Convert the predicted mean log-lead concentration back to its original scale
#    and visualise it.

# 2. Repeat this analysis by including also a set of covariates (e.g. PM10, 
#    West-East, North-South trends...). You will find a raster file for
#    PM10 in the data folder. You can load it with the raster function and 
#    then use the extract function to retrieve PM10 values for the prediction
#    locations. Make sure that the prediction locations  and the raster are in 
#    the same coordinate reference system.

# 3. Visualise side by side the predicted mean and standard errors from the 
#    model with an without covariates.
