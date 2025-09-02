

#### R CODE FOR AGE PREDICTIONS ####


# LIBRARIES
library(readr)
library(dplyr)
library(parallel)
library(progress)
library(ggplot2)
library(pbapply)
library(fitdistrplus)


# DATASET TO IMPORT:
# - ISO_SPOTS_ph_id_complete (it contains all the data related to the isochrones)
# - jackson_members_filt_binarie_final7000 (it contains all the stars of interest for the age predictions)



#################################
### dist_point_to_segment_2.0 ###
#################################

dist_point_to_segment_2.0 <- function(star, isocrone_p1, isocrone_p2) {
  
  px <- star[,3]
  py <- star[,2]
  x1 <- isocrone_p1[,5]
  y1 <- isocrone_p1[,4]
  x2 <- isocrone_p2[,5]
  y2 <- isocrone_p2[,4]
  
  A <- px - x1
  B <- py - y1
  C <- x2 - x1
  D <- y2 - y1
  
  m <- (y2 - y1) / (x2 - x1)
  
  dot <- A * C + B * D
  len_sq <- C * C + D * D
  param <- if (len_sq != 0) dot / len_sq else -1
  
  xx_yy <- mapply(function(p) {
    if (p < 0) {
      c(x1, y1)
    } else if (p > 1) {
      c(x2, y2)
    } else {
      c(x1 + p * C, y1 + p * D)
    }
  }, p = param)
  
  xx <- xx_yy[1, ]
  yy <- xx_yy[2, ]
  
  dx <- px - xx
  dy <- py - yy
  
  dist <- sqrt(dx * dx + dy * dy)
  
  return(dist)
}



###############################
### dist.isocrones.star_2.0 ###
###############################

dist.isocrones.star_2.0 <- function(star, isocrone) {
  cl <- makeCluster(mc <- getOption("cl.cores", 4))
  clusterExport(cl, varlist = c("dist_point_to_segment_2.0", "df_pris2022_pt1", "isocrone_beta0_ML2"))
  dist_stars_iso_i <- pbsapply(cl = cl, 1:(nrow(isocrone) - 1), function(j){
    if (isocrone[j, 7] == isocrone[j + 1, 7]) {
      dist_point_to_segment_2.0(star, isocrone[j,], isocrone[j+1,])
    } else {NA}
  })
  stopCluster(cl)
  dist_stars_iso_i
}



###############################
### dist.isocrones.star_3.0 ###
###############################

dist.isocrones.star_3.0 <- function(star_info, isocrone) {
  cl <- makeCluster(mc <- getOption("cl.cores", 10))
  
  # Carica i pacchetti necessari su ciascun nodo del cluster
  clusterEvalQ(cl, library(dplyr))
  clusterExport(cl, varlist = c("dist_point_to_segment_2.0", "isocrone"))
  
  dist_stars_iso_i <- pbsapply(cl = cl, 1:nrow(star_info), function(i) {
    iso_filtered <- isocrone %>% filter(age_yr %in% c(star_info[i, 4], star_info[i, 7]))
    sapply(1:(nrow(iso_filtered) - 1), function(j) {
      if (iso_filtered[j, 7] == iso_filtered[j + 1, 7]) {
        dist_point_to_segment_2.0(star_info[i,], iso_filtered[j,], iso_filtered[j + 1,])
      } else {
        NA
      }
    })
  })
  
  stopCluster(cl)
  dist_stars_iso_i
}



###################################
### position.point.to.curve_2.0 ###
###################################

position.point.to.curve_2.0 <- function(star, curve) {
  
  curve <- curve[order(curve$x), ]
  px <- star[,3]
  py <- star[,2]
  
  y_interpolated <- approx(curve$x, curve$y, xout = px)$y
  x_interpolated <- approx(curve$y, curve$x, xout = py)$y
  
  sign <- ifelse(!is.na(y_interpolated) & py <= y_interpolated, "DX",
                 ifelse(!is.na(y_interpolated) & py > y_interpolated, "SX",
                        ifelse(is.na(y_interpolated) & px <= x_interpolated, "DX",
                               ifelse(is.na(y_interpolated) & px > x_interpolated, "SX", "SX"))))
  return(sign)
}



###############################
### sign.isocrones.star_2.0 ###
###############################

sign.isocrones.star_2.0 <- function(star, isocrone) {
  sign_stars_iso_i <- pbsapply(unique(isocrone$isocrona), function(j){
    curve_prova <- isocrone[which(isocrone$isocrona == j),c(5,4)]
    colnames(curve_prova) <- c("x", "y")
    position.point.to.curve_2.0(star, curve_prova)
  })
  sign_stars_iso_i
}



#########################
### weighted.mean.age ###
#########################

weighted.mean.age <- function(dist1, age1, dist2, age2) {
  total_dist <- dist1 + dist2
  (age1 * 1/dist1 + age2 * 1/dist2) /(1/dist2 + 1/dist1)
}



#############################
### ISOCHRONE OF INTEREST ###
#############################

isocrone_beta0_ML2 <- ISO_SPOTS_ph_id_complete %>% filter(beta == 0 & ML == 2)
isocrone_beta0_ML2 <- isocrone_beta0_ML2[order(isocrone_beta0_ML2$age_yr),]
isocrone_beta0.2_ML2 <- ISO_SPOTS_ph_id_complete %>% filter(beta == 0.2 & ML == 2)
isocrone_beta0.2_ML2 <- isocrone_beta0.2_ML2[order(isocrone_beta0.2_ML2$age_yr),]
isocrone_beta0.4_ML2 <- ISO_SPOTS_ph_id_complete %>% filter(beta == 0.4 & ML == 2)
isocrone_beta0.4_ML2 <- isocrone_beta0.4_ML2[order(isocrone_beta0.4_ML2$age_yr),]
isocrone_beta0.6_ML2 <- ISO_SPOTS_ph_id_complete %>% filter(beta == 0.6 & ML == 2)
isocrone_beta0.6_ML2 <- isocrone_beta0.6_ML2[order(isocrone_beta0.6_ML2$age_yr),]





#### MAIN CODE FOR AGE INTERPOLATION ####


# In this code, necessary for the age prediction, we employed the predicted temperature from the Neural Network approach. In the work, we also
# used the original spectroscopic temperature from GES to compare the age predictions.

colnames(jackson_members_filt_binarie_final7000)
table(jackson_members_filt_binarie_final7000$CLUSTER)

df_jackson_for_age <- jackson_members_filt_binarie_final7000[,c(1,2,3)]

isocrone = isocrone_beta0.4_ML2
sign_str_cluster <- sign.isocrones.star_2.0(df_jackson_for_age, isocrone)
df_sign_str_cluster <- t(sign_str_cluster)

sign_df = df_sign_str_cluster

dist_iso_df <- data.frame(isocrona = isocrone$isocrona, age = isocrone$age_yr)
unique.iso <- dist_iso_df %>% distinct()

iso.sign.change <- as.data.frame(matrix(NA, ncol = 4, nrow = nrow(df_jackson_for_age)))
pb <- progress_bar$new(format = "  Interpolating [:bar] :percent in :elapsed, eta: :eta [:current/2826]", total = nrow(df_jackson_for_age))


for(i in 1:nrow(df_jackson_for_age)) {
  
  iso.sign <- as.data.frame(cbind(unique.iso, sign = sign_df[,i]))
  
  if (all(iso.sign$sign == "SX", na.rm = T)) {
    iso.sign.change[i,] <- c(200000000, "SX", "SX", 200000000)
  } else {
    
    if (all(iso.sign$sign == "DX", na.rm = T)) {
      iso.sign.change[i,] <- c(1000000, "DX", "DX", 1000000)
    } else {
      
      iso.sign <- iso.sign[which(!is.na(iso.sign$sign)),]
      
      iso.sign.change[i,] <- iso.sign[,-1] %>%
        mutate(next_sign = lead(sign), next_age = lead(age)) %>%
        filter(sign != next_sign)
    }
  }
  pb$tick()
}


colnames(iso.sign.change) <- c("age.SX", "sign.prev", "sign.post", "age.DX")
iso.sign.change$age.SX <- as.numeric(iso.sign.change$age.SX)
iso.sign.change$age.DX <- as.numeric(iso.sign.change$age.DX)

star_info <- cbind(df_jackson_for_age, iso.sign.change)

dist_prova <- dist.isocrones.star_3.0(star_info, isocrone)

iso.sign.change$dist.SX <- NA
iso.sign.change$dist.DX <- NA

pb <- progress_bar$new(format = "  Interpolating [:bar] :percent in :elapsed, eta: :eta [:current/2826]", total = nrow(df_jackson_for_age))

for(i in 1:nrow(iso.sign.change)) {
  if(iso.sign.change[i,1] == 200000000) {
    iso.sign.change[i,5] <- 1
    iso.sign.change[i,6] <- 1
  } else {
    if(iso.sign.change[i,4] == 1000000) {
      iso.sign.change[i,5] <- 1
      iso.sign.change[i,6] <- 1
    } else {
      iso.sign.change[i,5] <- min(dist_prova[[i]][1:which(is.na(dist_prova[[i]]))], na.rm = T)
      iso.sign.change[i,6] <- min(dist_prova[[i]][which(is.na(dist_prova[[i]]))+1:length(dist_prova[[i]])], na.rm = T)
    }
  }
  pb$tick()
}

age_jackson_clusters_beta0 <- weighted.mean.age(iso.sign.change$dist.SX, iso.sign.change$age.SX, iso.sign.change$dist.DX, iso.sign.change$age.DX)
age_jackson_clusters_beta0.2 <- weighted.mean.age(iso.sign.change$dist.SX, iso.sign.change$age.SX, iso.sign.change$dist.DX, iso.sign.change$age.DX)
age_jackson_clusters_beta0.4 <- weighted.mean.age(iso.sign.change$dist.SX, iso.sign.change$age.SX, iso.sign.change$dist.DX, iso.sign.change$age.DX)



df_all_age_jackson <- cbind(jackson_members_filt_binarie_final7000, age_beta0 = age_jackson_clusters_beta0,
                            age_beta0.2 = age_jackson_clusters_beta0.2, age_beta0.4 = age_jackson_clusters_beta0.4)

mean_age_per_cluster <- df_all_age_jackson %>%
  group_by(CLUSTER) %>%
  summarise(mean_age_beta0 = mean(log10(age_beta0), na.rm = TRUE),
            mean_age_beta0.2 = mean(log10(age_beta0.2), na.rm = TRUE),
            mean_age_beta0.4 = mean(log10(age_beta0.4), na.rm = TRUE),
            sd_age_beta0 = sd(log10(age_beta0), na.rm = TRUE),
            sd_age_beta0.2 = sd(log10(age_beta0.2), na.rm = TRUE),
            sd_age_beta0.4 = sd(log10(age_beta0.4), na.rm = TRUE),
            median_age_beta0 = median(log10(age_beta0), na.rm = TRUE),
            median_age_beta0.2 = median(log10(age_beta0.2), na.rm = TRUE),
            median_age_beta0.4 = median(log10(age_beta0.4), na.rm = TRUE),
            q16_age_beta0 = quantile(log10(age_beta0), 0.16, na.rm = TRUE),
            q16_age_beta0.2 = quantile(log10(age_beta0.2), 0.16, na.rm = TRUE),
            q16_age_beta0.4 = quantile(log10(age_beta0.4), 0.16, na.rm = TRUE),
            q84_age_beta0 = quantile(log10(age_beta0), 0.84, na.rm = TRUE),
            q84_age_beta0.2 = quantile(log10(age_beta0.2), 0.84, na.rm = TRUE),
            q84_age_beta0.4 = quantile(log10(age_beta0.4), 0.84, na.rm = TRUE))

as.data.frame(mean_age_per_cluster)





#### BETA SELECTION FROM STARSPOT EVOLUTIONARY MODELS ####

median_per_cluster <- df_all_age_jackson %>%
  group_by(CLUSTER) %>%
  summarise(
    median_b0_ML2_Est = median(age_beta0, na.rm = TRUE),
    median_b0.2_ML2_Est = median(age_beta0.2, na.rm = TRUE),
    median_b0.4_ML2_Est = median(age_beta0.4, na.rm = TRUE)
  )

median_per_cluster <- as.data.frame(median_per_cluster)


chisq.res.MG0 <- as.data.frame(matrix(NA, nrow = length(unique(df_LOG_age_jackson$CLUSTER)), ncol = 3))
colnames(chisq.res.MG0) <- c("beta0", "beta0.2", "beta0.4")
cluster_to_row <- setNames(seq_along(unique(df_LOG_age_jackson$CLUSTER)), unique(df_LOG_age_jackson$CLUSTER))

pb <- progress_bar$new(format = "  [:bar] :percent | Iterazione :current/:total | Trascorso :elapsed | Rimasto :eta", total = nrow(chisq.res.MG0), width = 120)
for(i in unique(df_LOG_age_jackson$CLUSTER)) {
  for(j in c(0,0.2,0.4)) {
    k <- ifelse(j == 0, 1, ifelse(j == 0.2, 2, 3))
    median.age.cluster <-  median_per_cluster[which(median_per_cluster$CLUSTER == i), k + 1]
    isochrone_set <- ISO_SPOTS_ph_id_complete %>% filter(beta == j & ML == 2)
    unique_isochrone_set_age <- sort(unique(isochrone_set$age_yr))
    diff_isochrone_median.age <- which.min(abs(unique_isochrone_set_age - median.age.cluster))
    isocrona.mediana <- isochrone_set %>% filter(age_yr == unique_isochrone_set_age[diff_isochrone_median.age]) 
    RDW_cluster <- df_LOG_age_jackson %>% filter(CLUSTER == i)
    MG0.isocrona <- approx(x = isocrona.mediana$log_Te, y = isocrona.mediana$G, xout = RDW_cluster$logTeff)$y
    res.MG0 <- RDW_cluster$MG0_ML - MG0.isocrona
    chisq.res.MG0[cluster_to_row[as.character(i)], paste0("beta", j)] <- sum((res.MG0/MG0.isocrona)^2, na.rm = T)
  }
  pb$tick()
}

table(RDW_data$Cluster)
cluster <- unique(df_LOG_age_jackson$CLUSTER)
chisq.res.MG0.cluster <- cbind(cluster, chisq.res.MG0)


beta_clusters <- NULL
for(i in 1:nrow(chisq.res.MG0.cluster)) {
  beta_selected <- which.min(chisq.res.MG0.cluster[i,c(2:4)])
  beta_clusters[i] <- ifelse(beta_selected == 1, 0, ifelse(beta_selected == 2, 0.2, 0.4))
}

# Trumpler14 was not included in the final results since its values are very uncertain due to its large distance (2.35 kpc)
beta_selected_jackson <- cbind(chisq.res.MG0.cluster, beta_clusters)





#### CODE TO CREATE THE PERTURBATED DATASET OF THE JACKSON CATALOGUE TO PREDICT THE STD DEV ASSOCIATED TO THE AGE OF THE CLUSTERS ####

# I need to simulate 100 values for each star in order to obtain a distribution of MG0_ML used in the H-R diagram for the interpolation of the age. 
# However, this magnitude depends on three components: distance, apparent G magnitude, and temperature.
# Therefore, I need to determine the sigma of each of these three components and simulate 100 values for each star.

# stars of interest for the age predictions
dim(jackson_members_filt_binarie_final7000)
head(jackson_members_filt_binarie_final7000)

# Gaia information of all the stars of the catalogue. Import the dataset expanded_GG2M.csv 
#expanded_GG2M <- read_csv("C:/Users/Marco Tarantino/OneDrive - UNIPA/Desktop/INAF WORK/new work - expanded catalogue/new work - expanded catalogue/dataset per analisi/expanded_GG2M.csv", 
#                          col_types = cols(ges_id_gaia = col_character(), 
#                                           source_id = col_character()))

expanded_G2M <- expanded_GG2M
expanded_G2M <- as.data.frame(expanded_G2M)
colnames(expanded_G2M)


### ERRORS OF THE APPARENT G MAGNITUDE OF GAIA ###

# flux formula to compute the error associated to the G apparent magnitude
expanded_G2M$phot_g_mean_mag_err <- sqrt((-2.5/log(10)*expanded_G2M$phot_g_mean_flux_error/expanded_G2M$phot_g_mean_flux)^2+(0.0027553202)^2)


### ERRORS OF PHOTOGEO DISTANCE ###
expanded_G2M$r_lo_photogeo
expanded_G2M$r_med_photogeo
expanded_G2M$r_hi_photogeo

df_error_MG0_ML <- inner_join(jackson_members_filt_binarie_final7000, expanded_G2M, 
                              c("ges_id_gaia" = "ges_id_gaia"))

colnames(df_error_MG0_ML)
df_error_MG0_ML <- df_error_MG0_ML[,c(1,2,3,4,25,24,26,16,36)]
head(df_error_MG0_ML)
summary(df_error_MG0_ML)


# selection of the cluster of interest
df_error_MG0_ML_TEFF <- df_error_MG0_ML[which(df_error_MG0_ML$CLUSTER %in% c("Cha_I", "NGC2451a", "NGC6405")),]                                            
df_error_MG0_ML_TEFF[which(df_error_MG0_ML_TEFF$CLUSTER == "NGC6709"),] # ID: 4311559274405377280
df_error_MG0_ML_TEFF <- df_error_MG0_ML_TEFF[-which(df_error_MG0_ML_TEFF$ges_id_gaia == "4311559274405377280"),]
dim(df_error_MG0_ML_TEFF)


# simulate 100 values of apparent G magnitude for each star of interest
sim100_magG <- matrix(NA, ncol = 100, nrow = nrow(df_error_MG0_ML_TEFF))
set.seed(2)
for(i in 1:nrow(df_error_MG0_ML_TEFF)) {
  sim100_magG[i,] <- rnorm(100, df_error_MG0_ML_TEFF[i,8], df_error_MG0_ML_TEFF[i,9])
}


# simulate 100 values of r_med_photogeo distances for each star of interest
hist(df_error_MG0_ML_TEFF$r_med_photogeo)
hist(df_error_MG0_ML_TEFF$r_lo_photogeo)
hist(df_error_MG0_ML_TEFF$r_hi_photogeo)

# asymmetric distribution of the distance, trying to find the best fit
fit_norm <- fitdist(df_error_MG0_ML_TEFF$r_med_photogeo, "norm")
fit_lnorm <- fitdist(df_error_MG0_ML_TEFF$r_med_photogeo, "lnorm")
fit_gamma <- fitdist(df_error_MG0_ML_TEFF$r_med_photogeo, "gamma")
fit_weibull <- fitdist(df_error_MG0_ML_TEFF$r_med_photogeo, "weibull")

summary(fit_norm)
summary(fit_lnorm)
summary(fit_gamma)
summary(fit_weibull)

# log-normal distribution selected
mu_lnorm <- (log(df_error_MG0_ML_TEFF$r_lo_photogeo) + log(df_error_MG0_ML_TEFF$r_hi_photogeo)) / 2
sigma_lnorm <- (log(df_error_MG0_ML_TEFF$r_hi_photogeo) - log(df_error_MG0_ML_TEFF$r_lo_photogeo)) / (2 * qnorm(0.84))

sim100_dist <- matrix(NA, ncol = 100, nrow = nrow(df_error_MG0_ML_TEFF))

set.seed(2)
for(i in 1:nrow(df_error_MG0_ML_TEFF)) {
  sim100_dist[i,] <- rlnorm(100, meanlog = mu_lnorm[i], sdlog = sigma_lnorm[i])
}



save(df_error_MG0_ML, sim100_magG, sim100_dist, file = "sim100_MG0_ML_error_new.RData")


# import of the files with the standard deviations related to the TEFF obtained through the bootstrap procedure on Python
#test_IC_new <- read_csv("C:/Users/Marco Tarantino/Downloads/test_IC_new.csv", 
#                        col_types = cols(ID = col_character()))
#train_IC_new <- read_csv("C:/Users/Marco Tarantino/Downloads/train_IC_new.csv", 
#                         col_types = cols(ID = col_character()))

test_IC_new <- as.data.frame(test_IC_new)
train_IC_new <- as.data.frame(train_IC_new)

length(unique(test_IC_new$ID))
length(unique(train_IC_new$ID))

df_std.dev_TEFF_train <- inner_join(df_error_MG0_ML_TEFF, train_IC_new, c("ges_id_gaia" = "ID"))
df_std.dev_TEFF_train <- df_std.dev_TEFF_train[,-c(10,11)]
colnames(df_std.dev_TEFF_train) <- c("ges_id_gaia",         "MG0_ML",              "logTeff",             "CLUSTER",             "r_lo_photogeo" ,     
                                     "r_med_photogeo",      "r_hi_photogeo" ,      "phot_g_mean_mag",     "phot_g_mean_mag_err", "TEFF",         
                                     "TEFF_DNN",            "std_dev")

df_std.dev_TEFF_test <- inner_join(df_error_MG0_ML_TEFF, test_IC_new, c("ges_id_gaia" = "ID"))
df_std.dev_TEFF_test <- df_std.dev_TEFF_test[,-c(10,11)]

colnames(df_std.dev_TEFF_test) <- c("ges_id_gaia",         "MG0_ML",              "logTeff",             "CLUSTER",             "r_lo_photogeo" ,     
                                    "r_med_photogeo",      "r_hi_photogeo" ,      "phot_g_mean_mag",     "phot_g_mean_mag_err", "TEFF",         
                                    "TEFF_DNN",            "std_dev")

df_std.dev_TEFF <- rbind(df_std.dev_TEFF_train, df_std.dev_TEFF_test)
head(df_std.dev_TEFF)
df_std.dev_TEFF <- df_std.dev_TEFF %>% group_by(CLUSTER)
df_std.dev_TEFF <- as.data.frame(df_std.dev_TEFF)


# simulations of 100 values of the TEFF for each star
sim100_TEFF_NN_new <- matrix(NA, ncol = 100, nrow = nrow(df_std.dev_TEFF))
set.seed(2)
for(i in 1:nrow(df_std.dev_TEFF)) {
  sim100_TEFF_NN_new[i,] <- rnorm(100, df_std.dev_TEFF[i,11], df_std.dev_TEFF[i,12])
}

dim(sim100_TEFF_NN_new)
hist(df_std.dev_TEFF$TEFF_DNN)
hist(sim100_TEFF_NN_new[,50])

which.max(df_std.dev_TEFF$std_dev)
df_std.dev_TEFF[1724,]
plot(df_std.dev_TEFF$TEFF_DNN, df_std.dev_TEFF$std_dev, xlab = "TEFF DNN", ylab = "std.dev TEFF DNN", pch = 20, cex = 0.7)


vec_TEFF_sim_new <- NULL
for(i in 1:nrow(sim100_TEFF_NN_new)) {
  vec_TEFF_sim.row_new <- sim100_TEFF_NN_new[i,]
  vec_TEFF_sim_new <- c(vec_TEFF_sim_new, vec_TEFF_sim.row_new)
}

length(vec_TEFF_sim_new)

# Building the dataset with the 100 pertubated values for each star. The dataset must include the ID, 
# the G magnitude, the distance r_med_photogeo, and the 100 temperature values

head(df_std.dev_TEFF)

df_stars_GG2M_complete <- read_csv("C:/Users/Marco Tarantino/OneDrive - UNIPA/Desktop/INAF WORK/new work - expanded catalogue/new work - expanded catalogue/dataset per analisi/df_stars_GG2M_complete.csv", 
                                   col_types = cols(ges_id_gaia = col_character()))

df_stars_GG2M_complete <- as.data.frame(df_stars_GG2M_complete)
dim(df_stars_GG2M_complete)
df_stars_GG2M_complete <- df_stars_GG2M_complete[,c(1,2,3,4,5, 17, 20, 26, 28, 30)]


df_Loredana_new <- inner_join(df_std.dev_TEFF, df_stars_GG2M_complete, c("ges_id_gaia"= "ges_id_gaia"))
head(df_Loredana_new)
dim(df_Loredana_new)

dataset_repeated_new <- df_Loredana_new %>% slice(rep(1:n(), each = 100))
dim(dataset_repeated_new)
head(dataset_repeated_new)

df_TEFF_sim_new <- cbind(dataset_repeated_new, vec_TEFF_sim_new)
dim(df_TEFF_sim_new)
head(df_TEFF_sim_new)
summary(df_TEFF_sim_new[1:100,])

df_TEFF_sim_new_pt2 <- df_TEFF_sim_new %>% filter(CLUSTER %in% c("Cha_I", "NGC2451a", "NGC6405"))


write.csv(df_TEFF_sim_new_pt2, "df_TEFF_sim_new_pt2.csv", row.names = FALSE)



### STD DEV DI MG0_ML ###

# The values of G magnitude, distance, and temperature have been simulated.
# Now we need to combine the three components to get 100 MG0_ML values for each star.


#merge_jackson <- read_csv("C:/Users/Marco Tarantino/OneDrive - UNIPA/Desktop/INAF WORK/new work - expanded catalogue/new work - expanded catalogue/materiale per stima eta'/std_dev_age/merge_jackson.csv", 
#                          col_types = cols(ges_id_gaia = col_character()))

merge_jackson <- as.data.frame(merge_jackson)
colnames(merge_jackson)
merge_jackson_var_sel <- merge_jackson[,c("ges_id_gaia", "e_bp_rp_pred", "AG_ML")]
head(merge_jackson_var_sel)
merge_jackson_ebprp <- merge_jackson_var_sel[which(merge_jackson_var_sel$e_bp_rp_pred < 0),]


# datasets which contain the perturbated observation for each star of interest. The second dataset contains all the info of the clusters 
# Cha_I, NGC2451a and NGC6405 which were initially excluded and then reconsidered in the analysis
#red_law_gaia_grid_NN_mont_sim <- read_csv("C:/Users/Marco Tarantino/OneDrive - UNIPA/Desktop/INAF WORK/new work - expanded catalogue/new work - expanded catalogue/materiale per stima eta'/std_dev_age/red_law_gaia_grid_NN_mont_sim.csv")
#df_AG_sim <- read_csv("red_law_gaia_grid_NN_mont_sim_chai_2451a_6405.csv", col_types = cols(ges_id_gaia = col_character()))
df_AG_sim <- as.data.frame(df_AG_sim)


# The table merge_jackson is used to select the stars with e_bp_rp_pred < 0, 
# because for these stars we cannot use the simulated AG values. 
# We will use the unperturbed Vergely values instead.

df_std_dev_ebprp_min0 <- inner_join(df_std.dev_TEFF, merge_jackson_ebprp, c("ges_id_gaia" = "ges_id_gaia"))
head(df_std_dev_ebprp_min0) # dataset con le stelle per le quali non devo perturbare i valori di AG
dim(df_std_dev_ebprp_min0)

vec_magG_sim <- NULL
for(i in 1:nrow(sim100_magG)) {
  vec_magG_sim.row <- sim100_magG[i,]
  vec_magG_sim <- c(vec_magG_sim, vec_magG_sim.row)
}

vec_dist_sim <- NULL
for(i in 1:nrow(sim100_dist)) {
  vec_dist_sim.row <- sim100_dist[i,]
  vec_dist_sim <- c(vec_dist_sim, vec_dist_sim.row)
}

length(vec_magG_sim)
length(vec_dist_sim)

dim(df_AG_sim)

df_mg_dist_Teff_100sim <- cbind(df_AG_sim, vec_magG_sim, vec_dist_sim)

head(df_mg_dist_Teff_100sim) 
df_mg_dist_Teff_100sim$MG0_ML <- df_TEFF_sim_new$MG0_ML 


# Replacing AG for stars with negative bp-rp error
head(df_std_dev_ebprp_min0) # dataset of the stars for which this replacement must be done, using Vergely's AG values


for(i in 1:nrow(df_std_dev_ebprp_min0)) {
  df_mg_dist_Teff_100sim[which(df_mg_dist_Teff_100sim[,1] == df_std_dev_ebprp_min0[i,1]),26] <- df_std_dev_ebprp_min0[i,14]
}


# Creating the dataset with the perturbed columns to generate the MG0_ML variable
df_MG0_ML_perturbated <- df_mg_dist_Teff_100sim[,c(1,2,3,4,22,32,33,26)]
head(df_MG0_ML_perturbated)
dim(df_MG0_ML_perturbated)

# Formula to calculate MG0_ML: MG0_ML = phot_g_mean_mag - AG_ML - 5 * log10(r_med_photogeo) + 5
df_MG0_ML_perturbated$vec_MG0_ML <- df_MG0_ML_perturbated$vec_magG_sim - df_MG0_ML_perturbated$AG_ML -
  5*(log10(df_MG0_ML_perturbated$vec_dist_sim)) + 5

df_jackson_perturbated <- df_MG0_ML_perturbated[,c(1,4,2,3,9,5)]
head(df_jackson_perturbated)
df_jackson_perturbated_new_pt2 <- df_jackson_perturbated
df_jackson_perturbated_new_pt2$vec_TEFF_sim <- log10(df_jackson_perturbated_new_pt2$vec_TEFF_sim)
dim(df_jackson_perturbated_new_pt2)

save(df_jackson_perturbated_new_pt2, file = "df_jackson_perturbated_new_pt2.RData")


# In df_jackson_perturbated_new_pt2, I have the 100 simulations for each TEFF and MG0_ML value. 
# To have the MG0_ML simulations centered on the initial MG0_ML value, 
# I need to calculate the standard deviation for each star and then simulate 
# 100 values with the specified mean and standard deviation.

std.dev_MG0_ML <- df_jackson_perturbated_new_pt2 %>% group_by(ges_id_gaia) %>% summarise(std.dev.MG = sd(vec_MG0_ML))
std.dev_MG0_ML <- as.data.frame(std.dev_MG0_ML)
df_jackson_perturbated_new_pt2$std.dev_MG0_ML <- NA

for(i in 1:nrow(std.dev_MG0_ML)) {
  df_jackson_perturbated_new_pt2[which(df_jackson_perturbated_new_pt2[,1] == std.dev_MG0_ML[i,1]),7] <- std.dev_MG0_ML[i,2]
}


df_jackson_perturbated_new_pt2_MG0_sim <- df_jackson_perturbated_new_pt2 %>% distinct(ges_id_gaia, .keep_all = TRUE)


sim100_MG0_ML_new <- matrix(NA, ncol = 100, nrow = nrow(df_jackson_perturbated_new_pt2_MG0_sim))
set.seed(2)
for(i in 1:nrow(df_jackson_perturbated_new_pt2_MG0_sim)) {
  sim100_MG0_ML_new[i,] <- rnorm(100, df_jackson_perturbated_new_pt2_MG0_sim[i,3], df_jackson_perturbated_new_pt2_MG0_sim[i,7])
}

vec_MG0_sim <- NULL
for(i in 1:nrow(sim100_MG0_ML_new)) {
  vec_MG0_sim.row <- sim100_MG0_ML_new[i,]
  vec_MG0_sim <- c(vec_MG0_sim, vec_MG0_sim.row)
}

df_jackson_perturbated_new_pt2 <- cbind(df_jackson_perturbated_new_pt2, vec_MG0_sim)
head(df_jackson_perturbated_new_pt2)

df_jackson_perturbated_new_pt2 <- df_jackson_perturbated_new_pt2[, c(1,2,3,4,9,6)]

save(df_jackson_perturbated_new_pt2, file = "df_jackson_perturbated_new_pt2.RData")


# by this procedure we obtain df_jackson_perturbated_new and df_jackson_perturbated_new_pt2 with 100 perturbated values for each star of 
# the 22 clusters of interest.















