####################################################################################


############################ Lab 1: Spatial exploratory ############################


####################################################################################


############# load packages #############
library(sf)
library(tmap)
library(dplyr)
library(tidyr)
library(geoR)
library(ggplot2)

############ source useful functions #####################
source("function.R")

######### Section 1: preliminary exploration of the dataset in R

# Load the dataset named \texttt{GaliciaData.csv} into R

data <- read.csv("data/GaliciaData.csv")


# Check the columns of the data and ensure that you understand what each column means.
View(data)

# Plot the histogram of the lead concentration

hist(data$loglead)

# summary of the data

summary(data)


####### Section 2

# Display on a map the log of the lead concentration.

## convert the R code to sf object
data_sf <- data %>% st_as_sf(., coords=c("long", "lat"), crs = 4326)
# note that the code 4326 corresponds to the default coordinate reference system for long and lat


## map the loglead
tm_shape(data_sf) + tm_symbols(col="loglead", midpoint = NA)  

## scatter plot of the log of the lead concentration and pm10
ggplot(data, aes(x = pm10, y = loglead)) + geom_point()

# put a smooth line on it 
ggplot(data, aes(x = pm10, y = loglead)) + geom_point() + geom_smooth()

## scatter plot of the log of the lead concentration and all the other variables
ggplot( data %>% gather(key, value, -loglead, -lead), 
        aes(x = value, y = loglead)) + geom_point() + 
  facet_wrap(facets = ~key, scales = "free") + geom_smooth()


########## Section 3
### construct a variogram

# First thing to note is that the variogram is constructed on the residual. 
# So fit a linear model to the loglead and extract the residual

fit <- lm(formula = loglead ~ pm10, data = data)
resi <- residuals(fit)

ggvario(coords = data[, c("long", "lat")], data = resi, bins = 15, 
        show_nbins = F, envelope = F)

# Note that we set envelop to FALSE to have just the variogram

# note that the distance is still in degrees, so we need to convert the long and lat to meters

## convert the gps coordinate to web mercator ############
utmcoords <- data_sf %>% st_transform(., crs= 3857) %>% st_coordinates()
# note that 3857 is the code for the projection to Web Mercator 

#### now construct the variogram again in kilometers
ggvario(coords = utmcoords/1000, data = resi, bins = 15, show_nbins = F)


##### Now to perform the permutation test ###########
# We will add  the envelop and nsim argument
ggvario(coords = utmcoords/1000, data = resi, bins = 15, show_nbins = F, envelop = TRUE, nsim = 999)

