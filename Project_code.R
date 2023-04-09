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
#input: the negative log likelihood value and number of parameters
#output: the AIC score
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
#inputs: a vector of initial parameter values, the tag data data frame, a distance matrix
#and the chosen detection function
#outputs: the negative log likelihood and parameter estimates for the euclidean distance model
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

#inputs: a vector of initial parameter values, the tag data data frame, a mask of the area,
#a raster of ruggednesses for each cell, the activity centres, 
#the chosen detection function and the chosen cost function
#outputs: the negative log likelihood and parameter estimates for the non-euclidean distance model
                   
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


########Goodness of Fit Testing#########################

#calculate the expected encounter rates for best model
#cost function 2 exponential encounter rate

params = mod4exer$estimate

lambda0 = exp(params[1])
sigma = exp(params[2])
alpha = params[3]

transition_matrix <- transition(rugg_raster,function(x) 
  costFunc2(x,alpha),directions = 16,symm  = F)
transition_matrix_corrected <- geoCorrection(transition_matrix,
                                             type = 'c')
non_euc_distmat <- t(costDistance(transition_matrix_corrected,
                                  as.matrix(activityCentres),
                                  as.matrix(TostMask[,1:2])))

lambda_mat_mod4 <- lambda_exer(non_euc_distmat,lambda0,sigma)%*%
  diag(colSums(tagdata)) 
lambda_mat_mod4_norm <- lambda_mat_mod4
#normalise the columns
for (i in 1:8) {
  lambda_mat_mod4_norm[,i] <- lambda_mat_mod4[,i] / 
    sum(lambda_mat_mod4[,i]) * sum(tagdata[,i])
}


######GOODNESS OF FIT TESTS#############

#function to calculate test statistic
#input: the observed data matrix and the expected data matrix
#output: the test statistic
good_fit_test <- function(obs,exp) {
  return (sum((sqrt(obs) - sqrt(exp)) ^ 2))
}


#the observed is the tagdata and the 
#expected is the lambda encounter rate matrix
mod4_test_stat <- good_fit_test(tagdata, lambda_mat_mod4)



#simulate data from the model using multinomial
#function to change the columns in each of the 
#lambda matrices to probabilities
column_probs <- function(matrix) {
  for (i in 1:8) {
    matrix[,i] <- matrix[,i] / sum(matrix[,i])
  }
  return (matrix)
}

#do this for all the models
prob_matrix_mod4 <- column_probs(lambda_mat_mod4)

#sample from a multinomial distribution
#do this 1000 times and check the distribution
sim_mod4_test_stat <- c()


#function to compute a simulated test statistic
sim_test_stat_func <- function(tagdat,probdat,lambdamat) {
  sim_matrix <- matrix(0,nrow(tagdat), ncol(tagdat))
  for (i in 1:8) {
    sim_matrix[,i] <- rmultinom(1,sum(tagdata[,i]), probdat[,i])
  }
  #goodness of fit
  return (good_fit_test(sim_matrix, lambdamat))
}


set.seed(49)
for (j in 1:1000) {
  
  sim_mod4_test_stat <- c(sim_mod4_test_stat, 
                          sim_test_stat_func(tagdata, 
                                             prob_matrix_mod4, 
                                             lambda_mat_mod4))
}

#function to calculate p_values
p_values <- function(simulated, observed) {
  exceed_count <- length(simulated[simulated >= observed])
  p_val <- exceed_count / length(simulated)
  return(p_val)
}

p_values(sim_mod4_test_stat, mod4_test_stat)
#p-value of 0



#do this per individual
mod4_ind_test_stats <- numeric(8)
for (i in 1:8) {
  mod4_ind_test_stats[i] <- good_fit_test(tagdata[,i],
                                          lambda_mat_mod4[,i])
}

#dataframe
obs_mod_4_ind_df <- data.frame(ind = colnames(tagdata), 
                               test_stat = mod4_ind_test_stats)

#sample from a multinomial distribution for each of the models
#do this 1000 times and check the distribution


#function to compute a simulated test statistic
sim_test_stat_func_ind_ps <- function(tagdat,probdat,lambdamat) {
  sim_matrix <- matrix(0,nrow(tagdat), ncol(tagdat))
  ind_sim_test_stat <- numeric(ncol(tagdat))
  for (i in 1:8) {
    sim_matrix[,i] <- rmultinom(1,sum(tagdata[,i]), probdat[,i])
    ind_sim_test_stat[i] <- good_fit_test(sim_matrix[,i], 
                                          lambdamat[,i])
  }
  #goodness of fit
  return (ind_sim_test_stat)
}

#simulate test statistic per individual
set.seed(49)
sim_mod4_test_stat_ind <- matrix(0,1000,8)
for (j in 1:1000) {
  sim_mod4_test_stat_ind[j,] <- sim_test_stat_func_ind_ps(
    tagdata, prob_matrix_mod4, lambda_mat_mod4)
}

#make into a dataframe
sim_mod4_test_stat_ind_vec <- c(sim_mod4_test_stat_ind)
sim_mod4_test_stat_ind_df <- data.frame(ind = rep(colnames
                                                  (tagdata), 
                                                  each = 1000), 
                                        test_stat = 
                                          sim_mod4_test_stat_ind_vec)




mod_4_ind_ps <- numeric(0)
for (i in 1:8) {
  mod_4_ind_ps[i] <- p_values(sim_mod4_test_stat_ind[,i], 
                              mod4_ind_test_stats[i])
}

mod_4_ind_ps
#all 0 apart from Individual 2 0.044

###########Now test the gof statistic actually works
#function to simulate one set of data
sim_obs_func <- function(tagdat,probdat) {
  sim_matrix <- matrix(0,nrow(tagdat), ncol(tagdat))
  for (i in 1:8) {
    sim_matrix[,i] <- rmultinom(1,sum(tagdata[,i]), probdat[,i])
  }
  return (sim_matrix)
}


#matrix for 50 x 8000 simulated
mod_4_gofs <- matrix(0,1000,50)
#vector for 100 observed
mod_4_gof_obs <- numeric(50)

set.seed(123)
#okay do this whole thing 100 times
for (k in 1:50) {
  #simulate one set of simulated for model 4 from the observed
  sim_data_mod4 <- sim_obs_func(tagdata, prob_matrix_mod4)
  
  #estimate activity centre
  ac_x_sim <- sapply(1:ncol(sim_data_mod4), function(x) 
    weighted.mean(TostMask$x,sim_data_mod4[,x]))
  ac_y_sim <- sapply(1:ncol(sim_data_mod4), function(x) 
    weighted.mean(TostMask$y,sim_data_mod4[,x]))
  activityCentres_sim <- data.frame(x = ac_x_sim, y = ac_y_sim)
  
  #fit model on simulated data
  mod4exer_sim <- nlm(non_euc_fixed_centre_lik, mod4$estimate,
                      tagdata = sim_data_mod4, 
                      mesh = TostMask,raster = rugg_raster,
                      centres = activityCentres_sim,
                      detectfn =  lambda_exer, 
                      costFunction = costFunc2,print.level = 2)
  
  #calculate lambda_matrix
  
  params = mod4exer_sim$estimate
  
  lambda0 = exp(params[1])
  sigma = exp(params[2])
  alpha = params[3]
  
  transition_matrix <- transition(rugg_raster,function(x) 
    costFunc2(x,alpha),directions = 16,symm  = F)
  transition_matrix_corrected <- geoCorrection(transition_matrix,
                                               type = 'c')
  non_euc_distmat <- t(costDistance(transition_matrix_corrected,
                                    as.matrix(activityCentres_sim),
                                    as.matrix(TostMask[,1:2])))
  
  lambda_mat_mod4_sim <- lambda_exer(non_euc_distmat,lambda0,
                                     sigma) %*%
    diag(colSums(sim_data_mod4))
  
  #normalise the columns
  for (i in 1:8) {
    lambda_mat_mod4_sim[,i] <- lambda_mat_mod4_sim[,i] / 
      sum(lambda_mat_mod4_sim[,i]) * sum(sim_data_mod4[,i])
  }
  
  #'observed' goodness of fit statistic
  mod4_test_stat_sim <- good_fit_test(sim_data_mod4, 
                                      lambda_mat_mod4_sim)
  
  #turn lambda matrix into probabilities
  prob_matrix_mod4_sim <- column_probs(lambda_mat_mod4_sim)
  sim_mod4_test_stat_sim <- c()
  
  
  for (j in 1:1000) {
    #simulate from model
    sim_mod4_test_stat_sim <- c(sim_mod4_test_stat_sim, 
                                sim_test_stat_func(
                                  sim_data_mod4, 
                                  prob_matrix_mod4_sim, 
                                  lambda_mat_mod4_sim))
  }
  mod_4_gofs[,k] <- sim_mod4_test_stat_sim
  mod_4_gof_obs[k] <- mod4_test_stat_sim
  print(k)
}



#caluculate p values
gof_p_values <- numeric(50)
for (i in 1:50) {
  gof_p_values[i] <- p_values(mod_4_gofs[,i], mod_4_gof_obs[i])
}
summary(mod_4_gofs[,2])
mod_4_gof_obs[2]
#all zeros


mod_4_gofs_2 <- matrix(0,1000,50)
#vector for 100 observed
mod_4_gof_obs_2 <- numeric(50)

set.seed(987)
#okay do this whole thing 100 times
for (k in 1:50) {
  #simulate one set of simulated for model 4 from the observed
  sim_data_mod4 <- sim_obs_func(tagdata, prob_matrix_mod4)
  
  #estimate activity centre
  ac_x_sim <- sapply(1:ncol(sim_data_mod4), function(x) 
    weighted.mean(TostMask$x,sim_data_mod4[,x]))
  ac_y_sim <- sapply(1:ncol(sim_data_mod4), function(x) 
    weighted.mean(TostMask$y,sim_data_mod4[,x]))
  activityCentres_sim <- data.frame(x = ac_x_sim, y = ac_y_sim)
  
  #fit model on simulated data
  mod4exer_sim <- nlm(non_euc_fixed_centre_lik, mod4$estimate,
                      tagdata = sim_data_mod4, mesh = TostMask,
                      raster = rugg_raster,
                      centres = activityCentres_sim,
                      detectfn =  lambda_exer, 
                      costFunction = costFunc2,
                      print.level = 2)
  
  #calculate lambda_matrix
  
  params = mod4exer_sim$estimate
  
  lambda0 = exp(params[1])
  sigma = exp(params[2])
  alpha = params[3]
  
  transition_matrix <- transition(rugg_raster,function(x) 
    costFunc2(x,alpha),directions = 16,symm  = F)
  transition_matrix_corrected <- geoCorrection(transition_matrix,
                                               type = 'c')
  non_euc_distmat <- t(costDistance(transition_matrix_corrected,
                                    as.matrix(
                                      activityCentres_sim),
                                    as.matrix(TostMask[,1:2])))
  
  lambda_mat_mod4_sim <- lambda_exer(non_euc_distmat,lambda0,
                                     sigma) %*%
    diag(colSums(sim_data_mod4)) 
  
  #normalise the columns
  for (i in 1:8) {
    lambda_mat_mod4_sim[,i] <- lambda_mat_mod4_sim[,i] / 
      sum(lambda_mat_mod4_sim[,i]) * sum(sim_data_mod4[,i])
  }
  
  #'observed' goodness of fit statistic
  mod4_test_stat_sim <- good_fit_test(sim_data_mod4, 
                                      lambda_mat_mod4_sim)
  
  #turn lambda matrix into probabilities
  prob_matrix_mod4_sim <- column_probs(lambda_mat_mod4_sim)
  sim_mod4_test_stat_sim <- c()
  
  
  for (j in 1:1000) {
    #simulate from model
    sim_mod4_test_stat_sim <- c(sim_mod4_test_stat_sim, 
                                sim_test_stat_func(sim_data_mod4, 
                                                   prob_matrix_mod4_sim, 
                                                   lambda_mat_mod4_sim))
  }
  mod_4_gofs_2[,k] <- sim_mod4_test_stat_sim
  mod_4_gof_obs_2[k] <- mod4_test_stat_sim
  print(k)
}



#caluculate p values
gof_p_values_2 <- numeric(50)
for (i in 1:50) {
  gof_p_values_2[i] <- p_values(mod_4_gofs_2[,i], 
                                mod_4_gof_obs_2[i])
}


gof_p_values <- c(gof_p_values,gof_p_values_2)
length(gof_p_values[gof_p_values <= 0.05])
#goodness of fit test statistic works




############FIT MODELS TO INDIVIDUALS########################


#fit models to different individuals


#log likelihood for fitting the models for individuals
non_euc_fixed_centre_lik_ind <- function(params,tagdata,
                                         mesh,raster,centres,
                                         detectfn,costFunction){
  lambda0 = exp(params[1])
  sigma = exp(params[2])
  alpha = params[3]
  
  transition_matrix <- transition(rugg_raster,function(x) 
    costFunction(x,alpha),directions = 16,symm  = F)
  transition_matrix_corrected <- geoCorrection(transition_matrix,
                                               type = 'c')
  non_euc_distmat <- t(costDistance(transition_matrix_corrected,
                                    as.matrix(centres),
                                    as.matrix(mesh[,1:2])))
  
  lambda_mat <- detectfn(non_euc_distmat,lambda0,sigma)%*%sum(tagdata) 
  
  lik = sum(dpois(tagdata,lambda_mat,log = T))
  
  return(-lik)
}


#a function to fit the models, compare AIC and return best model
best_model <- function(ind_tagdata, ind_ac) {
  #fit all 8 possible models
  params <- mod3$estimate
  mod3_inds <- nlm(non_euc_fixed_centre_lik_ind,params,
                   tagdata = ind_tagdata, mesh = TostMask,
                   raster = rugg_raster,centres = ind_ac,
                   detectfn =  lambda_hner, costFunction = costFunc1,
                   print.level = 2)
  params <- mod4$estimate
  mod4_inds <- nlm(non_euc_fixed_centre_lik_ind,params,
                   tagdata = ind_tagdata, mesh = TostMask,
                   raster = rugg_raster,centres = ind_ac,
                   detectfn =  lambda_hner, costFunction = costFunc2,
                   print.level = 2)
  params <- mod10test$estimate
  mod10_inds <- nlm(non_euc_fixed_centre_lik_ind,params,
                    tagdata = ind_tagdata, mesh = TostMask,
                    raster = rugg_raster,centres = ind_ac,
                    detectfn =  lambda_hner, costFunction = costFunc3,
                    print.level = 2)
  params <- mod12$estimate
  mod12_inds <- nlm(non_euc_fixed_centre_lik_ind,params,
                    tagdata = ind_tagdata, mesh = TostMask,
                    raster = rugg_raster,centres = ind_ac,
                    detectfn =  lambda_hner, costFunction = costFunc4,
                    print.level = 2)
  params <- mod3exer$estimate
  mod3exer_inds <- nlm(non_euc_fixed_centre_lik_ind,params,
                       tagdata = ind_tagdata, mesh = TostMask,
                       raster = rugg_raster,centres = ind_ac,
                       detectfn =  lambda_exer, costFunction = costFunc1,
                       print.level = 2)
  params <- mod4exer$estimate
  mod4exer_inds <- nlm(non_euc_fixed_centre_lik_ind,params,
                       tagdata = ind_tagdata, mesh = TostMask,
                       raster = rugg_raster,centres = ind_ac,
                       detectfn =  lambda_exer, costFunction = costFunc2,
                       print.level = 2)
  params <- mod10testexer$estimate
  mod10exer_inds <- nlm(non_euc_fixed_centre_lik_ind,params,
                        tagdata = ind_tagdata, mesh = TostMask,
                        raster = rugg_raster,centres = ind_ac,
                        detectfn =  lambda_exer, costFunction = costFunc3,
                        print.level = 2)
  params <- mod12exer$estimate
  mod12exer_inds <- nlm(non_euc_fixed_centre_lik_ind,params,
                        tagdata = ind_tagdata, mesh = TostMask,
                        raster = rugg_raster,centres = ind_ac,
                        detectfn =  lambda_exer, costFunction = costFunc4,
                        print.level = 2) 
  
  #calculate all AICs
  aic_mod3 <- aic(mod3_inds$minimum, length(mod3_inds$estimate))
  aic_mod4 <-aic(mod4_inds$minimum, length(mod4_inds$estimate))
  aic_mod10 <-aic(mod10_inds$minimum, length(mod10_inds$estimate))
  aic_mod12 <-aic(mod12_inds$minimum, length(mod12_inds$estimate))
  
  aic_mod3ex <-aic(mod3exer_inds$minimum, length(mod3exer_inds$estimate))
  aic_mod4ex <-aic(mod4exer_inds$minimum, length(mod4exer_inds$estimate))
  aic_mod10ex <-aic(mod10exer_inds$minimum, length(mod10exer_inds$estimate))
  aic_mod12ex <-aic(mod12exer_inds$minimum, length(mod12exer_inds$estimate))
  
  #create vector of AICs
  aic_list <- c(mod3 = aic_mod3,mod4 = aic_mod4,mod10 = aic_mod10,
                mod12 = aic_mod12, mod3ex = aic_mod3ex, mod4ex = aic_mod4ex,
                mod10ex = aic_mod10ex, mod12ex = aic_mod12ex)
  
  #find lowest and return
  return(names(aic_list)[which.min(aic_list)])
  
}


#find the best model for each individual
best_model(tagdata[,1], activityCentres[1,])
best_model(tagdata[,2], activityCentres[2,])
best_model(tagdata[,3], activityCentres[3,])
best_model(tagdata[,4], activityCentres[4,])
best_model(tagdata[,5], activityCentres[5,])
best_model(tagdata[,6], activityCentres[6,])
best_model(tagdata[,7], activityCentres[7,])
best_model(tagdata[,8], activityCentres[8,])

#run the best model for each individual
params <- mod4exer$estimate
mod4exer_ind1 <- nlm(non_euc_fixed_centre_lik_ind,params,
                     tagdata = tagdata[,1], mesh = TostMask,
                     raster = rugg_raster,centres = activityCentres[1,],
                     detectfn =  lambda_exer, costFunction = costFunc2,
                     print.level = 2, hessian = TRUE)

SEfromHessian(mod4exer_ind1$hessian)

params <- mod10testexer$estimate
mod10exer_ind2 <- nlm(non_euc_fixed_centre_lik_ind,params,
                      tagdata = tagdata[,2], mesh = TostMask,
                      raster = rugg_raster,centres = activityCentres[2,],
                      detectfn =  lambda_exer, costFunction = costFunc3,
                      print.level = 2, hessian = TRUE)
SEfromHessian(mod10exer_ind2$hessian)

params <- mod4$estimate
mod4_ind3 <- nlm(non_euc_fixed_centre_lik_ind,params,
                 tagdata = tagdata[,3], mesh = TostMask,
                 raster = rugg_raster,centres = activityCentres[3,],
                 detectfn =  lambda_hner, costFunction = costFunc2,
                 print.level = 2, hessian = TRUE)
SEfromHessian(mod4_ind3$hessian)

params <- mod4$estimate
mod4_ind4 <- nlm(non_euc_fixed_centre_lik_ind,params,
                 tagdata = tagdata[,4], mesh = TostMask,
                 raster = rugg_raster,centres = activityCentres[4,],
                 detectfn =  lambda_hner, costFunction = costFunc2,
                 print.level = 2, hessian = TRUE)
SEfromHessian(mod4_ind4$hessian)

params <- mod4$estimate
mod4_ind5 <- nlm(non_euc_fixed_centre_lik_ind,params,
                 tagdata = tagdata[,5], mesh = TostMask,
                 raster = rugg_raster,centres = activityCentres[5,],
                 detectfn =  lambda_hner, costFunction = costFunc2,
                 print.level = 2, hessian = TRUE)
SEfromHessian(mod4_ind5$hessian)

params <- mod4$estimate
mod4_ind6 <- nlm(non_euc_fixed_centre_lik_ind,params,
                 tagdata = tagdata[,6], mesh = TostMask,
                 raster = rugg_raster,centres = activityCentres[6,],
                 detectfn =  lambda_hner, costFunction = costFunc2,
                 print.level = 2, hessian = TRUE)
SEfromHessian(mod4_ind6$hessian)

params <- mod4$estimate
mod4_ind7 <- nlm(non_euc_fixed_centre_lik_ind,params,
                 tagdata = tagdata[,7], mesh = TostMask,
                 raster = rugg_raster,centres = activityCentres[7,],
                 detectfn =  lambda_hner, costFunction = costFunc2,
                 print.level = 2, hessian = TRUE)
SEfromHessian(mod4_ind7$hessian)

params <- mod4$estimate
mod4_ind8 <- nlm(non_euc_fixed_centre_lik_ind,params,
                 tagdata = tagdata[,8], mesh = TostMask,
                 raster = rugg_raster,centres = activityCentres[8,],
                 detectfn =  lambda_hner, costFunction = costFunc2,
                 print.level = 2, hessian = TRUE)
SEfromHessian(mod4_ind8$hessian)



#plot best paths for each individual
#remember to normalise the lambda_mats
#individual 1 - mod4exer
params = mod4exer_ind1$estimate
lambda0 = exp(params[1])
sigma = exp(params[2])
alpha = params[3]

transition_matrix <-
  transition(rugg_raster, function(x)
    costFunc2(x, alpha), directions = 16, symm = F)
transition_matrix_corrected <-
  geoCorrection(transition_matrix, type = 'c')
non_euc_distmat <-
  t(costDistance(
    transition_matrix_corrected,
    as.matrix(activityCentres[1, ]),
    as.matrix(TostMask[, 1:2])
  ))
lambda_mat_ind1 <-
  lambda_exer(non_euc_distmat, lambda0, sigma) %*% 
  sum(tagdata[, 1]) 
#normalise
lambda_mat_ind1 <-
  lambda_mat_ind1 / sum(lambda_mat_ind1) * sum(tagdata[, 1])


spath <- NULL
ind_locs <- ifelse(tagdata[,1] == 0, NA, tagdata[,1])

caploc <- cbind(TostMask[,1:2],cap = tagdata[,1]) %>% 
  filter(cap>0) %>% 
  dplyr::select(x,y) %>% 
  as.matrix()

spath_ind1 <- shortestPath(transition_matrix_corrected,
                           as.matrix(activityCentres)[1,],
                           caploc,output = 'SpatialLines')


#individual 2 - mod10exer
params = mod10exer_ind2$estimate
lambda0 = exp(params[1])
sigma = exp(params[2])
alpha = params[3]

transition_matrix <- transition(rugg_raster,function(x) costFunc3(x,alpha),
                                directions = 16,symm = F)
transition_matrix_corrected <- geoCorrection(transition_matrix,type = 'c')
non_euc_distmat <- t(costDistance(transition_matrix_corrected,
                                  as.matrix(activityCentres[2,]),
                                  as.matrix(TostMask[,1:2])))

lambda_mat_ind2 <- lambda_exer(non_euc_distmat,lambda0,sigma)%*%
  sum(tagdata[,2]) 
lambda_mat_ind2 <- lambda_mat_ind2 / sum(lambda_mat_ind2) * sum(tagdata[,2])

spath <- NULL
ind_locs <- ifelse(tagdata[,2] == 0, NA, tagdata[,2])

caploc <- cbind(TostMask[,1:2],cap = tagdata[,2]) %>% 
  filter(cap>0) %>% 
  dplyr::select(x,y) %>% 
  as.matrix()

spath_ind2 <- shortestPath(transition_matrix_corrected,
                           as.matrix(activityCentres)[2,],caploc,
                           output = 'SpatialLines')





#individual 3 - mod4
params = mod4_ind3$estimate
lambda0 = exp(params[1])
sigma = exp(params[2])
alpha = params[3]

transition_matrix <- transition(rugg_raster,function(x) costFunc2(x,alpha),
                                directions = 16,symm = F)
transition_matrix_corrected <- geoCorrection(transition_matrix,type = 'c')
non_euc_distmat <- t(costDistance(transition_matrix_corrected,
                                  as.matrix(activityCentres[3,]),
                                  as.matrix(TostMask[,1:2])))

lambda_mat_ind3 <- lambda_hner(non_euc_distmat,lambda0,sigma)%*%
  sum(tagdata[,3]) 
lambda_mat_ind3 <- lambda_mat_ind3 / sum(lambda_mat_ind3) * sum(tagdata[,3])


spath <- NULL
ind_locs <- ifelse(tagdata[,3] == 0, NA, tagdata[,3])

caploc <- cbind(TostMask[,1:2],cap = tagdata[,3]) %>% 
  filter(cap>0) %>% 
  dplyr::select(x,y) %>% 
  as.matrix()

spath_ind3 <- shortestPath(transition_matrix_corrected,
                           as.matrix(activityCentres)[3,],
                           caploc,output = 'SpatialLines')



#individual 4 - mod4
params = mod4_ind4$estimate
lambda0 = exp(params[1])
sigma = exp(params[2])
alpha = params[3]

transition_matrix <- transition(rugg_raster,function(x) costFunc2(x,alpha),
                                directions = 16,symm = F)
transition_matrix_corrected <- geoCorrection(transition_matrix,type = 'c')
non_euc_distmat <- t(costDistance(transition_matrix_corrected,
                                  as.matrix(activityCentres[4,]),
                                  as.matrix(TostMask[,1:2])))

lambda_mat_ind4 <- lambda_hner(non_euc_distmat,lambda0,sigma)%*%
  sum(tagdata[,4]) 
lambda_mat_ind4 <- lambda_mat_ind4 / sum(lambda_mat_ind4) * sum(tagdata[,4])


spath <- NULL
ind_locs <- ifelse(tagdata[,4] == 0, NA, tagdata[,4])

caploc <- cbind(TostMask[,1:2],cap = tagdata[,4]) %>% 
  filter(cap>0) %>% 
  dplyr::select(x,y) %>% 
  as.matrix()

spath_ind4 <- shortestPath(transition_matrix_corrected,
                           as.matrix(activityCentres)[4,],
                           caploc,output = 'SpatialLines')


#individual 5
params = mod4_ind5$estimate
lambda0 = exp(params[1])
sigma = exp(params[2])
alpha = params[3]

transition_matrix <- transition(rugg_raster,function(x) costFunc2(x,alpha),
                                directions = 16,symm = F)
transition_matrix_corrected <- geoCorrection(transition_matrix,type = 'c')
non_euc_distmat <- t(costDistance(transition_matrix_corrected,
                                  as.matrix(activityCentres[5,]),
                                  as.matrix(TostMask[,1:2])))

lambda_mat_ind5 <- lambda_hner(non_euc_distmat,lambda0,sigma)%*%sum(tagdata[,5]) 
lambda_mat_ind5 <- lambda_mat_ind5 / sum(lambda_mat_ind5) * sum(tagdata[,5])


spath <- NULL
ind_locs <- ifelse(tagdata[,5] == 0, NA, tagdata[,5])

caploc <- cbind(TostMask[,1:2],cap = tagdata[,5]) %>% 
  filter(cap>0) %>% 
  dplyr::select(x,y) %>% 
  as.matrix()

spath_ind5 <- shortestPath(transition_matrix_corrected,
                           as.matrix(activityCentres)[5,],
                           caploc,output = 'SpatialLines')



#individual 6
params = mod4_ind6$estimate
lambda0 = exp(params[1])
sigma = exp(params[2])
alpha = params[3]

transition_matrix <- transition(rugg_raster,function(x) costFunc2(x,alpha),
                                directions = 16,symm = F)
transition_matrix_corrected <- geoCorrection(transition_matrix,type = 'c')
non_euc_distmat <- t(costDistance(transition_matrix_corrected,
                                  as.matrix(activityCentres[6,]),
                                  as.matrix(TostMask[,1:2])))

lambda_mat_ind6 <- lambda_hner(non_euc_distmat,lambda0,sigma)%*%sum(tagdata[,6])
lambda_mat_ind6 <- lambda_mat_ind6 / sum(lambda_mat_ind6) * sum(tagdata[,6])


spath <- NULL
ind_locs <- ifelse(tagdata[,6] == 0, NA, tagdata[,6])

caploc <- cbind(TostMask[,1:2],cap = tagdata[,6]) %>% 
  filter(cap>0) %>% 
  dplyr::select(x,y) %>% 
  as.matrix()

spath_ind6 <- shortestPath(transition_matrix_corrected,
                           as.matrix(activityCentres)[6,],caploc,
                           output = 'SpatialLines')



#individual 7
params = mod4_ind7$estimate
lambda0 = exp(params[1])
sigma = exp(params[2])
alpha = params[3]

transition_matrix <- transition(rugg_raster,function(x) costFunc2(x,alpha),
                                directions = 16,symm = F)
transition_matrix_corrected <- geoCorrection(transition_matrix,type = 'c')
non_euc_distmat <- t(costDistance(transition_matrix_corrected,
                                  as.matrix(activityCentres[7,]),
                                  as.matrix(TostMask[,1:2])))

lambda_mat_ind7 <- lambda_hner(non_euc_distmat,lambda0,sigma)%*%sum(tagdata[,7])
lambda_mat_ind7 <- lambda_mat_ind7 / sum(lambda_mat_ind7) * sum(tagdata[,7])


spath <- NULL
ind_locs <- ifelse(tagdata[,7] == 0, NA, tagdata[,7])

caploc <- cbind(TostMask[,1:2],cap = tagdata[,7]) %>% 
  filter(cap>0) %>% 
  dplyr::select(x,y) %>% 
  as.matrix()

spath_ind7 <- shortestPath(transition_matrix_corrected,
                           as.matrix(activityCentres)[7,],
                           caploc,output = 'SpatialLines')



#individual 8 - model 4
params = mod4_ind8$estimate
lambda0 = exp(params[1])
sigma = exp(params[2])
alpha = params[3]

transition_matrix <- transition(rugg_raster,function(x) costFunc2(x,alpha),
                                directions = 16,symm = F)
transition_matrix_corrected <- geoCorrection(transition_matrix,type = 'c')
non_euc_distmat <- t(costDistance(transition_matrix_corrected,
                                  as.matrix(activityCentres[8,]),
                                  as.matrix(TostMask[,1:2])))

lambda_mat_ind8 <- lambda_hner(non_euc_distmat,lambda0,sigma)%*%sum(tagdata[,8])
lambda_mat_ind8 <- lambda_mat_ind8 / sum(lambda_mat_ind8) * sum(tagdata[,8])


spath <- NULL
ind_locs <- ifelse(tagdata[,8] == 0, NA, tagdata[,8])

caploc <- cbind(TostMask[,1:2],cap = tagdata[,8]) %>% 
  filter(cap>0) %>% 
  dplyr::select(x,y) %>% 
  as.matrix()

spath_ind8 <- shortestPath(transition_matrix_corrected,
                           as.matrix(activityCentres)[8,],
                           caploc,output = 'SpatialLines')



#plot the paths for each individual best model
path_list <- list(spath_ind1,spath_ind2,spath_ind3,spath_ind4,spath_ind5,
                  spath_ind6,spath_ind7,spath_ind8)



par(mfrow = c(2,2))

plot(rugg_raster, xlab = 'UTM X-Coordinate', ylab = 'UTM Y-Coordinate', 
     main = ('Individual 1 Best Model'))
lines(spath_ind1)
plot(rugg_raster, xlab = 'UTM X-Coordinate', ylab = 'UTM Y-Coordinate', 
     main = ('Individual 2 Best Model'))
lines(spath_ind2)
plot(rugg_raster, xlab = 'UTM X-Coordinate', ylab = 'UTM Y-Coordinate', 
     main = ('Individual 3 Best Model'))
lines(spath_ind3)
plot(rugg_raster, xlab = 'UTM X-Coordinate', ylab = 'UTM Y-Coordinate', 
     main = ('Individual 4 Best Model'))
lines(spath_ind4)
plot(rugg_raster, xlab = 'UTM X-Coordinate', ylab = 'UTM Y-Coordinate', 
     main = ('Individual 5 Best Model'))
lines(spath_ind5)
plot(rugg_raster, xlab = 'UTM X-Coordinate', ylab = 'UTM Y-Coordinate', 
     main = ('Individual 6 Best Model'))
lines(spath_ind6)
plot(rugg_raster, xlab = 'UTM X-Coordinate', ylab = 'UTM Y-Coordinate', 
     main = ('Individual 7 Best Model'))
lines(spath_ind7)
plot(rugg_raster, xlab = 'UTM X-Coordinate', ylab = 'UTM Y-Coordinate', 
     main = ('Individual 8 Best Model'))
lines(spath_ind8)




#goodness of fit testing for these
column_probs_inds <- function(lambdas) {
  vector <- lambdas / sum(lambdas)
  return (vector)
}

#function to compute a simulated test statistic
sim_test_stat_func_inds <- function(tagdat,probdat,lambdamat) {
  sim_matrix <- rmultinom(1,sum(tagdat), probdat)
  #goodness of fit
  return (good_fit_test(sim_matrix, lambdamat))
}


ind1_test_stat <- good_fit_test(tagdata[,1], lambda_mat_ind1)
ind2_test_stat <- good_fit_test(tagdata[,2], lambda_mat_ind2)
ind3_test_stat <- good_fit_test(tagdata[,3], lambda_mat_ind3)
ind4_test_stat <- good_fit_test(tagdata[,4], lambda_mat_ind4)
ind5_test_stat <- good_fit_test(tagdata[,5], lambda_mat_ind5)
ind6_test_stat <- good_fit_test(tagdata[,6], lambda_mat_ind6)
ind7_test_stat <- good_fit_test(tagdata[,7], lambda_mat_ind7)
ind8_test_stat <- good_fit_test(tagdata[,8], lambda_mat_ind8)


#probabilities
prob_matrix_ind1 <- column_probs_inds(lambda_mat_ind1)
prob_matrix_ind2 <- column_probs_inds(lambda_mat_ind2)
prob_matrix_ind3 <- column_probs_inds(lambda_mat_ind3)
prob_matrix_ind4 <- column_probs_inds(lambda_mat_ind4)
prob_matrix_ind5 <- column_probs_inds(lambda_mat_ind5)
prob_matrix_ind6 <- column_probs_inds(lambda_mat_ind6)
prob_matrix_ind7 <- column_probs_inds(lambda_mat_ind7)
prob_matrix_ind8 <- column_probs_inds(lambda_mat_ind8)


#vector for simulated test statistics
sim_test_stat_ind1 <- c()
sim_test_stat_ind2 <- c()
sim_test_stat_ind3 <- c()
sim_test_stat_ind4 <- c()
sim_test_stat_ind5 <- c()
sim_test_stat_ind6 <- c()
sim_test_stat_ind7 <- c()
sim_test_stat_ind8 <- c()




set.seed(123)
for (j in 1:1000) {
  sim_test_stat_ind1 <- c(sim_test_stat_ind1, 
                          sim_test_stat_func_inds(tagdata[,1], 
                                                  prob_matrix_ind1, 
                                                  lambda_mat_ind1))
  sim_test_stat_ind2 <- c(sim_test_stat_ind2, 
                          sim_test_stat_func_inds(tagdata[,2], 
                                                  prob_matrix_ind2, 
                                                  lambda_mat_ind2))
  sim_test_stat_ind3 <- c(sim_test_stat_ind3, 
                          sim_test_stat_func_inds(tagdata[,3], 
                                                  prob_matrix_ind3, 
                                                  lambda_mat_ind3))
  sim_test_stat_ind4 <- c(sim_test_stat_ind4, 
                          sim_test_stat_func_inds(tagdata[,4], 
                                                  prob_matrix_ind4, 
                                                  lambda_mat_ind4))
  sim_test_stat_ind5 <- c(sim_test_stat_ind5, 
                          sim_test_stat_func_inds(tagdata[,5], 
                                                  prob_matrix_ind5, 
                                                  lambda_mat_ind5))
  sim_test_stat_ind6 <- c(sim_test_stat_ind6, 
                          sim_test_stat_func_inds(tagdata[,6], 
                                                  prob_matrix_ind6, 
                                                  lambda_mat_ind6))
  sim_test_stat_ind7 <- c(sim_test_stat_ind7, 
                          sim_test_stat_func_inds(tagdata[,7], 
                                                  prob_matrix_ind7, 
                                                  lambda_mat_ind7))
  sim_test_stat_ind8 <- c(sim_test_stat_ind8, 
                          sim_test_stat_func_inds(tagdata[,8], 
                                                  prob_matrix_ind8, 
                                                  lambda_mat_ind8))
  
}


p_values(sim_test_stat_ind1, ind1_test_stat)
p_values(sim_test_stat_ind2, ind2_test_stat)
p_values(sim_test_stat_ind3, ind3_test_stat)
p_values(sim_test_stat_ind4, ind4_test_stat)
p_values(sim_test_stat_ind5, ind5_test_stat)
p_values(sim_test_stat_ind6, ind6_test_stat)
p_values(sim_test_stat_ind7, ind7_test_stat)
p_values(sim_test_stat_ind8, ind8_test_stat)


