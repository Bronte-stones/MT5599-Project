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