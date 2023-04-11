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