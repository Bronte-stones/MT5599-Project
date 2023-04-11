
#goodness of fit testing for individual models
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