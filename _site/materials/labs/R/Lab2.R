
library(PrevMap)
#library(ggplot2)

 
# Section 4 

# Question 1. Loading the data
galicia_data <- read.csv("GaliciaData.csv")

 
# Checking the structure of the data
str(galicia_data)
head(galicia_data)


# Creating new coordinates system to Web Mercator
# to a plane so that distances are in kilometers
coords.ii <- SpatialPoints(as.matrix(galicia_data[,c("long","lat")]), 
                           CRS("+init=epsg:4326"))
coords.ii.web <- spTransform(coords.ii,CRS("+init=epsg:3857"))
galicia_data$web_x <- coordinates(coords.ii.web)[,1]/1000
galicia_data$web_y <- coordinates(coords.ii.web)[,2]/1000

plot(web_y ~ web_x   , data =  galicia_data, asp=1)


# Question 1 (a). Model (1)
model1 <- linear.model.MLE(
            formula = loglead ~ 1 ,
            coords = ~ web_x + web_y , 
            data = galicia_data ,
            kappa = 0.5,
            fixed.rel.nugget = 0,
            start.cov.pars = c(50),
            method = "nlminb",
            ID.coords = NULL,
            low.rank = FALSE,
            knots = NULL,
            messages = TRUE,
            profile.llik = FALSE,
            SPDE = FALSE,
            mesh = NULL,
            SPDE.analytic.hessian = FALSE )



# Question 1 (b). Model (2)
model2 <- linear.model.MLE(
  formula = loglead ~ 1 ,
  coords = ~ web_x + web_y , 
  data = galicia_data ,
  kappa = 0.5,
  fixed.rel.nugget = NULL,
  start.cov.pars = c(50, 0.1),
  method = "nlminb",
  ID.coords = NULL,
  low.rank = FALSE,
  knots = NULL,
  messages = TRUE,
  profile.llik = FALSE,
  SPDE = FALSE,
  mesh = NULL,
  SPDE.analytic.hessian = FALSE )


summary(model2, log.cov.pars=FALSE) # What does the error message suggest?


# View the results using the following code
estimates <- data.frame(estimate=model2$estimate, 
                       std_err = sqrt(diag(model2$covariance)) )
estimates$CI_lower <- estimates$estimate - 1.96*estimates$std_err
estimates$CI_upper <- estimates$estimate + 1.96*estimates$std_err
estimates[2:4,] <- exp(estimates[2:4,])
estimates
 
 

# Likelihood ratio test between models 1 and 2
D.obs <- 2*(model2$log.lik-model1$log.lik)
1-pchisq(D.obs,1)




# Question 1 (c). Model (3)
model3 <- linear.model.MLE(
  formula = loglead ~ 1 + long,
  coords = ~ web_x + web_y , 
  data = galicia_data ,
  kappa = 0.5,
  fixed.rel.nugget = 0,
  start.cov.pars = c(50),
  method = "nlminb",
  ID.coords = NULL,
  low.rank = FALSE,
  knots = NULL,
  messages = TRUE,
  profile.llik = FALSE,
  SPDE = FALSE,
  mesh = NULL,
  SPDE.analytic.hessian = FALSE )


summary(model3, log.cov.pars=FALSE)  


# Likelihood ratio test between models 1 and 3
D.obs <- 2*(model3$log.lik-model1$log.lik)
1-pchisq(D.obs,1)
