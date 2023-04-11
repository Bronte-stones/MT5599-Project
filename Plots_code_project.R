#survey area
#plot the landscape
ggplot(TostMask, aes(x = x, y = y))+
  geom_tile(aes(fill = rugg))+
  coord_equal() + 
  ggtitle("Survey Area") +
  xlab("UTM x coordinate") + ylab("UTM y coordinate") + 
  labs(fill = "Ruggedness (std GC)") + 
  theme(plot.title = element_text(hjust = 0.5))



#histograms of observations
rugg_hists <- list()
for (i in (1:8)) {
  rugg_hists[[i]] <- ggplot(rugg_tag_data, aes(x = rugg_bin, y = rugg_tag_data[,i])) +
    geom_col() +
    xlab('Ruggedness') + ylab('Frequency') +
    ggtitle(paste('Individual',i)) + 
    theme(plot.title = element_text(hjust = 0.5))
}
ggarrange(plotlist = rugg_hists, nrow = 4)



#investigating alpha plots
windows()
par(mfrow = c(2,2))
#plot individual 3 for different alphas for cost function 1
params = mod3exer$estimate
lambda0 = exp(params[1])
sigma = exp(params[2])
alpha = c(0.1,0.5,1,2)

for (i in alpha) {
  transition_matrix <- transition(rugg_raster,function(x) costFunc1(x,i),directions = 16,symm = F)
  transition_matrix_corrected <- geoCorrection(transition_matrix,type = 'c')
  non_euc_distmat <- t(costDistance(transition_matrix_corrected,as.matrix(activityCentres),as.matrix(TostMask[,1:2])))
  
  lambda_mat <- lambda_exer(non_euc_distmat,lambda0,sigma)%*%diag(colSums(tagdata)) ## calculate the encounter rate
  
  spath <- NULL
  ind <- 3
  ind_locs <- ifelse(tagdata[,ind] == 0, NA, tagdata[,ind])
  
  caploc <- cbind(TostMask[,1:2],cap = tagdata[,ind]) %>% 
    filter(cap>0) %>% 
    dplyr::select(x,y) %>% 
    as.matrix()
  spath1 <- shortestPath(transition_matrix_corrected,as.matrix(activityCentres)[ind,],caploc,output = 'SpatialLines') #%>% st_as_sf()
  
  plot(rugg_raster, main = paste('Alpha =', i), xlab = 'UTM X-Coordinate', ylab = 'UTM Y-Coordinate')
  lines(spath1)
  
}


#model 4
windows()
par(mfrow = c(2,2))
params = mod4exer$estimate
lambda0 = exp(params[1])
sigma = exp(params[2])
alpha = c(0.1,0.5,1,2)

for (i in alpha) {
  transition_matrix <- transition(rugg_raster,function(x) costFunc2(x,i),directions = 16,symm = F)
  transition_matrix_corrected <- geoCorrection(transition_matrix,type = 'c')
  non_euc_distmat <- t(costDistance(transition_matrix_corrected,as.matrix(activityCentres),as.matrix(TostMask[,1:2])))
  
  lambda_mat <- lambda_exer(non_euc_distmat,lambda0,sigma)%*%diag(colSums(tagdata)) ## calculate the encounter rate
  
  
  spath <- NULL
  ind <- 3
  ind_locs <- ifelse(tagdata[,ind] == 0, NA, tagdata[,ind])
  
  caploc <- cbind(TostMask[,1:2],cap = tagdata[,ind]) %>% 
    filter(cap>0) %>% 
    dplyr::select(x,y) %>% 
    as.matrix()
  
  spath1 <- shortestPath(transition_matrix_corrected,as.matrix(activityCentres)[ind,],caploc,output = 'SpatialLines') #%>% st_as_sf()
  plot(rugg_raster, main = paste('Alpha =',i), xlab = 'UTM X-Coordinate', ylab = 'UTM Y-Coordinate')
  lines(spath1)
}



#model 10
params = mod10testexer$estimate
lambda0 = exp(params[1])
sigma = exp(params[2])
alpha = c(0.1,0.5,1,2)

for (i in alpha) {
  transition_matrix <- transition(rugg_raster,function(x) costFunc3(x,i),
                                  directions = 16,symm = F)
  transition_matrix_corrected <- geoCorrection(transition_matrix,type = 'c')
  non_euc_distmat <- t(costDistance(transition_matrix_corrected,
                                    as.matrix(activityCentres),
                                    as.matrix(TostMask[,1:2])))
  
  lambda_mat <- lambda_exer(non_euc_distmat,lambda0,sigma)%*%
    diag(colSums(tagdata))
  
  spath <- NULL
  ind <- 3
  ind_locs <- ifelse(tagdata[,ind] == 0, NA, tagdata[,ind])
  
  caploc <- cbind(TostMask[,1:2],cap = tagdata[,ind]) %>% 
    filter(cap>0) %>% 
    dplyr::select(x,y) %>% 
    as.matrix()
  
  spath1 <- shortestPath(transition_matrix_corrected,
                         as.matrix(activityCentres)[ind,],
                         caploc,output = 'SpatialLines')
  plot(rugg_raster, main = paste('Alpha = ',i), xlab = 'UTM X-Coordinate', 
       ylab = 'UTM Y-Coordinate')
  lines(spath1)
}


#model 12
params = mod12exer$estimate
lambda0 = exp(params[1])
sigma = exp(params[2])
alpha = c(0.1,0.5,1,2)

for (i in alpha) {
  transition_matrix <- transition(rugg_raster,function(x) costFunc4(x,i),
                                  directions = 16,symm = F)
  transition_matrix_corrected <- geoCorrection(transition_matrix,type = 'c')
  non_euc_distmat <- t(costDistance(transition_matrix_corrected,
                                    as.matrix(activityCentres),
                                    as.matrix(TostMask[,1:2])))
  
  lambda_mat <- lambda_exer(non_euc_distmat,lambda0,sigma)%*%
    diag(colSums(tagdata))
  
  spath <- NULL
  ind <- 3
  ind_locs <- ifelse(tagdata[,ind] == 0, NA, tagdata[,ind])
  
  caploc <- cbind(TostMask[,1:2],cap = tagdata[,ind]) %>% 
    filter(cap>0) %>% 
    dplyr::select(x,y) %>% 
    as.matrix()
  
  spath1 <- shortestPath(transition_matrix_corrected,
                         as.matrix(activityCentres)[ind,],caploc,
                         output = 'SpatialLines') 
  plot(rugg_raster, main = paste('Alpha =',i), xlab = 'UTM X-Coordinate', 
       ylab = 'UTM Y-Coordinate')
  lines(spath1)
}



#paths under best joint model
params = mod4exer$estimate
lambda0 = exp(params[1])
sigma = exp(params[2])
alpha = params[3]

transition_matrix <- transition(rugg_raster,function(x) costFunc2(x,alpha),
                                directions = 16,symm = F)
transition_matrix_corrected <- geoCorrection(transition_matrix,type = 'c')
non_euc_distmat <- t(costDistance(transition_matrix_corrected,
                                  as.matrix(activityCentres),
                                  as.matrix(TostMask[,1:2])))

lambda_mat <- lambda_exer(non_euc_distmat,lambda0,sigma)%*%
  diag(colSums(tagdata)) 

windows()
par(mfrow = c(4,2))
for (i in 1:8) {
  spath <- NULL
  ind <- i
  ind_locs <- ifelse(tagdata[,ind] == 0, NA, tagdata[,ind])
  
  caploc <- cbind(TostMask[,1:2],cap = tagdata[,ind]) %>% 
    filter(cap>0) %>% 
    dplyr::select(x,y) %>% 
    as.matrix()
  
  spath1 <- shortestPath(transition_matrix_corrected,
                         as.matrix(activityCentres)[ind,],caploc,
                         output = 'SpatialLines')
  plot(rugg_raster)
  lines(spath1, size = 0.5)
  
  
}


#goodness of fit histogram best model
par(mfrow = c(1,1))
#plot histogram of test statistics of simulated data 
#and add a line for the observed

hist(sim_mod4_test_stat, xlim = c(3100, 4300), 
     xlab = 'Goodness of Fit Test Statistic',
     main = 'Histogram of Goodness of Fit Test Statistics')
abline(v = mod4_test_stat, col = 'red')



#plot individual gof under best joint model
ggplot(sim_mod4_test_stat_ind_df) +
  geom_boxplot(aes(x=ind, y=test_stat)) +
  geom_point(data = obs_mod_4_ind_df, aes(x=ind, y=test_stat, col = 'red') ) +
  labs(x = 'Individual', y = 'Test Statistic') + 
  ggtitle('Goodness of Fit Test Statistics per Individual') +
  theme(plot.title = element_text(hjust = 0.5))



par(mfrow = c(2,2))

hist(sim_test_stat_ind1, xlim = c(0, 850), 
     xlab = 'Goodness of Fit Test Statistic',
     main = 'Individual 1')
abline(v = ind1_test_stat, col = 'red')
legend(x="topleft", 
       legend=c("Observed"),
       col=c("red"), 
       lty = 1, cex = 0.8)



hist(sim_test_stat_ind2, xlim = c(0, 850), 
     xlab = 'Goodness of Fit Test Statistic',
     main = 'Individual 2')
abline(v = ind2_test_stat, col = 'red')
legend(x="topright", 
       legend=c("Observed"),
       col=c("red"), 
       lty = 1, cex = 0.8)

hist(sim_test_stat_ind3, xlim = c(0, 850), 
     xlab = 'Goodness of Fit Test Statistic',
     main = 'Individual 3')
abline(v = ind3_test_stat, col = 'red')
legend(x="topleft", 
       legend=c("Observed"),
       col=c("red"), 
       lty = 1, cex = 0.8)

hist(sim_test_stat_ind4, xlim = c(0, 850), 
     xlab = 'Goodness of Fit Test Statistic',
     main = 'Individual 4')
abline(v = ind4_test_stat, col = 'red')
legend(x="topright", 
       legend=c("Observed"),
       col=c("red"), 
       lty = 1, cex = 0.8)

hist(sim_test_stat_ind5, xlim = c(0, 850), 
     xlab = 'Goodness of Fit Test Statistic',
     main = 'Individual 5')
abline(v = ind5_test_stat, col = 'red')
legend(x="topright", 
       legend=c("Observed"),
       col=c("red"), 
       lty = 1, cex = 0.8)

hist(sim_test_stat_ind6, xlim = c(0, 850), 
     xlab = 'Goodness of Fit Test Statistic',
     main = 'Individual 6')
abline(v = ind6_test_stat, col = 'red')
legend(x="topright", 
       legend=c("Observed"),
       col=c("red"), 
       lty = 1, cex = 0.8)

hist(sim_test_stat_ind7, xlim = c(0, 850), 
     xlab = 'Goodness of Fit Test Statistic',
     main = 'Individual 7')
abline(v = ind7_test_stat, col = 'red')
legend(x="topright", 
       legend=c("Observed"),
       col=c("red"), 
       lty = 1, cex = 0.8)

hist(sim_test_stat_ind8, xlim = c(0, 850), 
     xlab = 'Goodness of Fit Test Statistic',
     main = 'Individual 8')
abline(v = ind8_test_stat, col = 'red')
legend(x="topleft", 
       legend=c("Observed"),
       col=c("red"), 
       lty = 1, cex = 0.8)



#plot simulations for best models spatially

###plotting detections of an individual
ind_locs <- ifelse(rowSums(tagdata) == 0, NA, rowSums(tagdata))


ggplot(TostMask, aes(x = x, y = y))+
  geom_tile()+
  labs(fill = "Ruggedness (std GC)", size = 'Observations') + 
  coord_equal() + 
  geom_point(aes(size = ind_locs, alpha=0.2), color="dodgerblue", shape=16)+
  ggtitle("Observed Snow leopards") +
  xlab("UTM x coordinate") + ylab("UTM y coordinate") +
  theme(plot.title = element_text(hjust = 0.5))




#plot observations simulated from the model
#best joint model
sim_matrix_model4 <- matrix(0,nrow(tagdata), ncol(tagdata))
for (i in 1:8) {
  sim_matrix_model4[,i] <- rmultinom(1,sum(tagdata[,i]), prob_matrix_mod4[,i])
}


sim_ind_locs <- ifelse(rowSums(sim_matrix_model4) == 0, NA, rowSums(sim_matrix_model4))

ggplot(TostMask, aes(x = x, y = y))+
  geom_tile()+
  labs(fill = "Ruggedness (std GC)", size = 'Observations') + 
  coord_equal() + 
  geom_point(aes(size = sim_ind_locs, alpha=0.1), color="dodgerblue", shape=16)+
  ggtitle("Simulated Snow Leopards Joint Model") +
  xlab("UTM x coordinate") + ylab("UTM y coordinate") +
  theme(plot.title = element_text(hjust = 0.5))




#simulated observations for each individual
sim_matrix_ind1 <- rmultinom(1,sum(tagdata[,1]), prob_matrix_ind1)
sim_matrix_ind2 <- rmultinom(1,sum(tagdata[,2]), prob_matrix_ind2)
sim_matrix_ind3 <- rmultinom(1,sum(tagdata[,3]), prob_matrix_ind3)
sim_matrix_ind4 <- rmultinom(1,sum(tagdata[,4]), prob_matrix_ind4)
sim_matrix_ind5 <- rmultinom(1,sum(tagdata[,5]), prob_matrix_ind5)
sim_matrix_ind6 <- rmultinom(1,sum(tagdata[,6]), prob_matrix_ind6)
sim_matrix_ind7 <- rmultinom(1,sum(tagdata[,7]), prob_matrix_ind7)
sim_matrix_ind8 <- rmultinom(1,sum(tagdata[,8]), prob_matrix_ind8)



#dataframe of simulated individual models
sim_matrix_sep <- data.frame(sim_matrix_ind1,sim_matrix_ind2,sim_matrix_ind3,
                             sim_matrix_ind4,sim_matrix_ind5,sim_matrix_ind6,
                             sim_matrix_ind7,sim_matrix_ind8)


sim_ind_locs_sep <- ifelse(rowSums(sim_matrix_sep) == 0, NA, rowSums(sim_matrix_sep))

ggplot(TostMask, aes(x = x, y = y))+
  geom_tile()+
  labs(fill = "Ruggedness (std GC)", size = 'Observations') + 
  coord_equal() + 
  geom_point(aes(size = sim_ind_locs_sep, alpha=0.1), color="dodgerblue", shape=16)+
  ggtitle("Simulated Snow Leopards Separate Models") +
  xlab("UTM x coordinate") + ylab("UTM y coordinate") +
  theme(plot.title = element_text(hjust = 0.5))



#plot encounter rates spatially for all individuals combined

lambda_rows <- rowSums(lambda_mat_mod4)
#under joint model
ggplot(TostMask, aes(x = x, y = y))+
  geom_tile(aes(fill = log(lambda_rows)))+
  coord_equal() + 
  scale_fill_viridis(limits=c(-12, 3),discrete = FALSE, option = "G")+
  ggtitle("Combined Snow Leopards Joint Model") +
  xlab("UTM x coordinate") + ylab("UTM y coordinate") + 
  labs(fill = "Log Encounter Rate") + 
  theme(plot.title = element_text(hjust = 0.5))


lambda_inds <- data.frame(c(lambda_mat_ind1),c(lambda_mat_ind2),
                          c(lambda_mat_ind3),c(lambda_mat_ind4),
                          c(lambda_mat_ind5),c(lambda_mat_ind6),
                          c(lambda_mat_ind7),c(lambda_mat_ind8))

lambda_ind_rows <- rowSums(lambda_inds)

#under individual models
ggplot(TostMask, aes(x = x, y = y))+
  geom_tile(aes(fill = log(lambda_ind_rows)))+
  coord_equal() + 
  scale_fill_viridis(limits=c(-12, 3),discrete = FALSE, option = "G")+
  ggtitle("Combined Snow Leopards Individual Models") +
  xlab("UTM x coordinate") + ylab("UTM y coordinate") + 
  labs(fill = "Log Encounter Rate") + 
  theme(plot.title = element_text(hjust = 0.5))






#plot individual encounter rates under joint model

ind_enc_plots <- list()
for (i in (1:8)) {
  ind_enc_plots[[i]] <- ggplot(TostMask, aes(x = x, y = y))+
    geom_tile(aes(fill = lambda_mat_mod4[,i]))+
    coord_equal() + 
    scale_fill_viridis(limits=c(0, 12),discrete = FALSE,option = "H")+
    ggtitle(paste("Individual",i)) +
    xlab("UTM x coordinate") + ylab("UTM y coordinate") + 
    labs(fill = "Encounter Rate") + 
    theme(plot.title = element_text(hjust = 0.5))
}

windows()
ggarrange(plotlist = ind_enc_plots, nrow = 4)



#plot encounter rate functions for individual models
ind5_gplot <- ggplot(TostMask, aes(x = x, y = y))+
  geom_tile(aes(fill = lambda_mat_ind1))+
  coord_equal() + 
  scale_fill_viridis(limits=c(0, 3),discrete = FALSE, option = "H")+
  ggtitle("Individual 1") +
  xlab("UTM x coordinate") + ylab("UTM y coordinate") + 
  labs(fill = "Encounter Rate") + 
  theme(plot.title = element_text(hjust = 0.5))


ind6_gplot <- ggplot(TostMask, aes(x = x, y = y))+
  geom_tile(aes(fill = lambda_mat_ind2))+
  coord_equal() + 
  scale_fill_viridis(limits=c(0, 3),discrete = FALSE, option = "H")+
  ggtitle("Individual 2") +
  xlab("UTM x coordinate") + ylab("UTM y coordinate") + 
  labs(fill = "Encounter Rate") + 
  theme(plot.title = element_text(hjust = 0.5))

ind7_gplot <- ggplot(TostMask, aes(x = x, y = y))+
  geom_tile(aes(fill = lambda_mat_ind3))+
  coord_equal() + 
  scale_fill_viridis(limits=c(0, 3),discrete = FALSE, option = "H")+
  ggtitle("Individual 3") +
  xlab("UTM x coordinate") + ylab("UTM y coordinate") + 
  labs(fill = "Encounter Rate") + 
  theme(plot.title = element_text(hjust = 0.5))

ind8_gplot <-ggplot(TostMask, aes(x = x, y = y))+
  geom_tile(aes(fill = lambda_mat_ind4))+
  coord_equal() + 
  scale_fill_viridis(limits=c(0, 3),discrete = FALSE, option = "H")+
  ggtitle("Individual 4") +
  xlab("UTM x coordinate") + ylab("UTM y coordinate") + 
  labs(fill = "Encounter Rate") + 
  theme(plot.title = element_text(hjust = 0.5))

windows()
ggarrange(ind1_gplot, ind2_gplot,ind3_gplot,ind4_gplot)


ind5_gplot <- ggplot(TostMask, aes(x = x, y = y))+
  geom_tile(aes(fill = lambda_mat_ind5))+
  coord_equal() + 
  scale_fill_viridis(limits=c(0, 3),discrete = FALSE, option = "H")+
  ggtitle("Individual 5") +
  xlab("UTM x coordinate") + ylab("UTM y coordinate") + 
  labs(fill = "Encounter Rate") + 
  theme(plot.title = element_text(hjust = 0.5))


ind6_gplot <- ggplot(TostMask, aes(x = x, y = y))+
  geom_tile(aes(fill = lambda_mat_ind6))+
  coord_equal() + 
  scale_fill_viridis(limits=c(0, 3),discrete = FALSE, option = "H")+
  ggtitle("Individual 6") +
  xlab("UTM x coordinate") + ylab("UTM y coordinate") + 
  labs(fill = "Encounter Rate") + 
  theme(plot.title = element_text(hjust = 0.5))

ind7_gplot <- ggplot(TostMask, aes(x = x, y = y))+
  geom_tile(aes(fill = lambda_mat_ind7))+
  coord_equal() + 
  scale_fill_viridis(limits=c(0, 3),discrete = FALSE, option = "H")+
  ggtitle("Individual 7") +
  xlab("UTM x coordinate") + ylab("UTM y coordinate") + 
  labs(fill = "Encounter Rate") + 
  theme(plot.title = element_text(hjust = 0.5))

ind8_gplot <-ggplot(TostMask, aes(x = x, y = y))+
  geom_tile(aes(fill = lambda_mat_ind8))+
  coord_equal() + 
  scale_fill_viridis(limits=c(0, 3),discrete = FALSE, option = "H")+
  ggtitle("Individual 8") +
  xlab("UTM x coordinate") + ylab("UTM y coordinate") + 
  labs(fill = "Encounter Rate") + 
  theme(plot.title = element_text(hjust = 0.5))

windows()
ggarrange(ind5_gplot, ind6_gplot,ind7_gplot,ind8_gplot)
