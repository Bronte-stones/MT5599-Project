library(secr)
library(sf)
library(tidyverse)
library(raster)
library(gdistance)
library(ggplot2)
library(HelpersMG)
library(ggpubr)
library(viridis)

#write a quick aic function
aic <- function(nll, pars) {
  return (2 * nll + 2 * pars)
}

###Read the mask
load('TostMask.Rdata')
rugg <- covariates(TostMask)$stdGC


###Read the tagdata (in capture history format)
tagdata <- read.csv('tagdata_capthist.csv')


###No of detections of each individual
ndets <- colSums(tagdata)

### define encounter rate functions
##Half normal encounter rate
lambda_hner <- function(distance, lambda0, sigma){
  return(lambda0*exp(-0.5*(distance/sigma)^2))
}


#### Model 1 
#### Euclidean distance model, we assume we know the activity 
#centre which is the mean of the locations


###compute activity centres of individuals
ac_x <- sapply(1:ncol(tagdata), function(x) 
  weighted.mean(TostMask$x,tagdata[,x]))
ac_y <- sapply(1:ncol(tagdata), function(x) 
  weighted.mean(TostMask$y,tagdata[,x]))
activityCentres <- data.frame(x = ac_x, y = ac_y)

###Compute distance from activity centres to each mesh point
dist_mat <- edist(TostMask,activityCentres)


####Maximise likelihood to estimate encounter rate parameters

params = log(c(.02,6000))
euc_fixed_centre_lik <- function(params,
                                 tagdata,
                                 distance,
                                 detectfn){
  lambda0 <- exp(params[1])  ##lambda0 and sigma are positive
  sigma <- exp(params[2])
  
  lambda_mat <- detectfn(distance,lambda0,sigma)%*%
    diag(colSums(tagdata)) ## calculate the encounter rate
  
  ## calculate log likelihood
  lik = sum(sapply(1:ncol(tagdata), 
                   function(x) dpois(tagdata[,x],
                                     lambda_mat[,x],
                                     log = T))) 
  
  return(-lik)
}


mod1 = nlm(euc_fixed_centre_lik,params, tagdata = tagdata, 
           distance = dist_mat,detectfn = lambda_hner,
           print.level = 2,hessian = T)



#Non euclidean distance models

####Function to convert mesh to a raster. 
#Need raster to for the function to calculate least cost distance

maskToRaster <- function(mesh,cov){
  rast_dat <- cbind(mesh,cov)
  coordinates(rast_dat) <- ~ x + y
  gridded(rast_dat) <- TRUE
  r <- raster(rast_dat)
  return(r)
}

rugg_raster <- maskToRaster(TostMask,rugg)

##define cost functions


###Symmetric cost function, where the conductance is the mean 
#of the two values. 

costFunc1 <- function(x,a){
  return(exp(a*mean(x)))
}


### Asymmetric cost function

costFunc2 <- function(x,a){
  return(exp(a*(x[2]-x[1])))
}


#symmetric cost function with squared term
costFunc3 <- function(x, a) {
  return (exp(a * mean(x)**2))
}

#assymetric cost function with squared term
costFunc4 <- function(x, a) {
  return (exp(a* ((x[2] - x[1]) ** 2)))
}

#models with non-euclidean distance and fixed activity centres
params = c(mod1$estimate,1)

non_euc_fixed_centre_lik <- function(params,tagdata,mesh,raster,
                                     centres,detectfn,
                                     costFunction){
  lambda0 = exp(params[1])
  sigma = exp(params[2])
  alpha = params[3]
  
  transition_matrix <- transition(raster,function(x) 
    costFunction(x,alpha),directions = 16,symm  = F)
  transition_matrix_corrected <- geoCorrection(transition_matrix,
                                               type = 'c')
  non_euc_distmat <- t(costDistance(transition_matrix_corrected,
                                    as.matrix(centres),
                                    as.matrix(mesh[,1:2])))
  
  lambda_mat <- detectfn(non_euc_distmat,lambda0,sigma)%*%
    diag(colSums(tagdata)) ## calculate the encounter rate
  
  ## calculate log likelihood
  lik = sum(sapply(1:ncol(tagdata), 
                   function(x) dpois(tagdata[,x],
                                     lambda_mat[,x],
                                     log = T))) 
  
  return(-lik)
}

mod3 <- nlm(non_euc_fixed_centre_lik,params,tagdata = tagdata, 
            mesh = TostMask,raster = rugg_raster,
            centres = activityCentres,detectfn =  lambda_hner, 
            costFunction = costFunc1,print.level = 2, hessian = TRUE)

mod4 <- nlm(non_euc_fixed_centre_lik,params,tagdata = tagdata, 
            mesh = TostMask,raster = rugg_raster,
            centres = activityCentres,detectfn =  lambda_hner, 
            costFunction = costFunc2,print.level = 2, hessian = TRUE)

mod10test <- nlm(non_euc_fixed_centre_lik,params,
                 tagdata = tagdata, mesh = TostMask,
                 raster = rugg_raster,detectfn =  lambda_hner, 
                 costFunction = costFunc3,
                 centres = activityCentres,print.level = 2, hessian = TRUE)

mod12 <- nlm(non_euc_fixed_centre_lik,params,tagdata = tagdata, 
             mesh = TostMask,raster = rugg_raster,
             centres = activityCentres,detectfn =  lambda_hner, 
             costFunction = costFunc4,print.level = 2, hessian = TRUE)



#repeat these models using an exponential encounter rate
##Exponential encounter rate
lambda_exer <- function(distance, lambda0, sigma){
  return(lambda0*exp(-distance/sigma))
}


mod1exer = nlm(euc_fixed_centre_lik,mod1$estimate, 
               tagdata = tagdata, distance = dist_mat,
               detectfn = lambda_exer,print.level = 2,
               hessian = T)


mod3exer <- nlm(non_euc_fixed_centre_lik,mod3$estimate,
                tagdata = tagdata, 
                mesh = TostMask,raster = rugg_raster,
                centres = activityCentres,detectfn =  lambda_exer, 
                costFunction = costFunc1,print.level = 2, hessian = TRUE)

mod4exer <- nlm(non_euc_fixed_centre_lik, mod4$estimate,
                tagdata = tagdata, 
                mesh = TostMask,raster = rugg_raster,
                centres = activityCentres,
                detectfn =  lambda_exer, 
                costFunction = costFunc2,print.level = 2, hessian = TRUE)

mod10testexer <- nlm(non_euc_fixed_centre_lik,mod10test$estimate,
                     tagdata = tagdata, mesh = TostMask,
                     raster = rugg_raster,
                     detectfn =  lambda_exer, 
                     costFunction = costFunc3,
                     centres = activityCentres,
                     print.level = 2, hessian = TRUE)

mod12exer <- nlm(non_euc_fixed_centre_lik,mod12$estimate,
                 tagdata = tagdata, mesh = TostMask,
                 raster = rugg_raster,centres = activityCentres,
                 detectfn =  lambda_exer, 
                 costFunction = costFunc4,
                 print.level = 2, hessian = TRUE)


#calculate AICS
aic(mod1exer$minimum, length(mod1exer$estimate))
aic(mod3exer$minimum, length(mod3exer$estimate))
aic(mod4exer$minimum, length(mod4exer$estimate))
aic(mod10testexer$minimum, length(mod10testexer$estimate))
aic(mod12exer$minimum, length(mod12exer$estimate))


aic(mod1$minimum, length(mod1$estimate))
aic(mod3$minimum, length(mod3$estimate))
aic(mod4$minimum, length(mod4$estimate))
aic(mod10test$minimum, length(mod10test$estimate))
aic(mod12$minimum, length(mod12$estimate))



#calculate delta AICs
aic(mod1$minimum, length(mod1$estimate)) - 
  aic(mod4exer$minimum, length(mod4exer$estimate))
aic(mod3$minimum, length(mod3$estimate)) - 
  aic(mod4exer$minimum, length(mod4exer$estimate))
aic(mod4$minimum, length(mod4$estimate)) - 
  aic(mod4exer$minimum, length(mod4exer$estimate))
aic(mod10test$minimum, length(mod10test$estimate)) - 
  aic(mod4exer$minimum, length(mod4exer$estimate))
aic(mod12$minimum, length(mod12$estimate)) - 
  aic(mod4exer$minimum, length(mod4exer$estimate))


aic(mod1exer$minimum, length(mod1exer$estimate)) - 
  aic(mod4exer$minimum, length(mod4exer$estimate))
aic(mod3exer$minimum, length(mod3exer$estimate)) - 
  aic(mod4exer$minimum, length(mod4exer$estimate))
aic(mod4exer$minimum, length(mod4exer$estimate)) - 
  aic(mod4exer$minimum, length(mod4exer$estimate))
aic(mod10testexer$minimum, length(mod10testexer$estimate)) - 
  aic(mod4exer$minimum, length(mod4exer$estimate))
aic(mod12exer$minimum, length(mod12exer$estimate)) - 
  aic(mod4exer$minimum, length(mod4exer$estimate))




#standard errors
SEfromHessian(mod4exer$hessian)