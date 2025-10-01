

### RICALCOLO SPREAD DI ERRORE PER ALCUNI CLUSTER ###

# cluster da rifare: NGC6530, Cha_I, Rho_Oph, NGC2264 (avevano beta = 0), NGC2244, lam_Ori_B35, Assc50, Col197 (avevano beta = 0.6)
# utilizzo per tutti beta = 0.4, solo per NGC6530 considero beta = 0.2

# rifare anche NGC2451a (beta = 0.4) e NGC6405 (beta = 0.2) in quanto ho sbagliato il limite inferiore di età (100000 al posto di 1000000)

load("age_interp_12iso_group.RData")
load("jackson_members_filt_binarie_final7000.RData")
load("ISO_SPOTS_ph_id_complete.RData")

load("df_jackson_perturbated_new.RData")
load("df_jackson_perturbated_new_pt2.RData")

colnames(df_jackson_perturbated_new)
colnames(df_jackson_perturbated_new_pt2)
colnames(df_jackson_perturbated_new_pt2)[6] <- "vec_TEFF_sim"
df_jackson_perturbated_new_pt2$vec_TEFF_sim <- log10(df_jackson_perturbated_new_pt2$vec_TEFF_sim)

df_jackson_perturbated_final <- rbind(df_jackson_perturbated_new, df_jackson_perturbated_new_pt2)
table(df_jackson_perturbated_final$CLUSTER)

# rho oph 3600 (beta 0.4) FATTO
# NGC6530 32500 (beta 0.2) FATTO
# Cha I 7600 (beta 0.4) FATTO
# NGC2264 47100 (beta 0.4) FATTO
# NGC2244 11200 (beta 0.4)
# lam_Ori_B35 4300 (beta 0.4) FATTO
# Assc50 17500 (beta 0.4) FATTO
# Col197 9800 (beta 0.4) FATTO
# NGC2451a 3800 (beta 0.4) FATTO
# NGC6405 4900 (beta 0.2) FATTO


##############
### NGC2244 ###
##############

# considero solo l'ID, MG0_ML e Teff_NN per il calcolo dell'età
head(df_jackson_perturbated_new)

df_NGC2244_for_age <- df_jackson_perturbated_final %>% filter(CLUSTER == "NGC2244")
df_NGC2244_for_age <- df_NGC2244_for_age[,c(1,5,6)]

# codice per ritrovare tra quali isocrone cade la stella
sign_NGC2244 <- sign.isocrones.star_2.0(df_NGC2244_for_age, isocrone_beta0.4_ML2)
df_sign_NGC2244 <- t(sign_NGC2244)

isocrone = isocrone_beta0.4_ML2
sign_df = df_sign_NGC2244

dist_iso_df <- data.frame(isocrona = isocrone$isocrona, age = isocrone$age_yr)
unique.iso <- dist_iso_df %>% distinct()

iso.sign.change <- as.data.frame(matrix(NA, ncol = 4, nrow = nrow(df_NGC2244_for_age)))

pb <- progress_bar$new(format = "  Interpolating [:bar] :percent in :elapsed, eta: :eta [:current/11200]", total = nrow(df_NGC2244_for_age))

for(i in 1:nrow(df_NGC2244_for_age)) {
  
  iso.sign <- as.data.frame(cbind(unique.iso, sign = sign_df[,i]))
  
  if (all(iso.sign$sign == "SX", na.rm = T)) {
    iso.sign.change[i,] <- c(200000000, "SX", "SX", 200000000)
  } else {
    
    if (all(iso.sign$sign == "DX", na.rm = T)) {
      iso.sign.change[i,] <- c(1000000, "DX", "DX", 1000000)
    } else {
      
      iso.sign_SX <- iso.sign %>% filter(sign == "SX")
      age_SX <- max(iso.sign_SX$age)
      
      iso.sign_DX <- iso.sign %>% filter(sign == "DX")
      age_DX <- min(iso.sign_DX$age)
      
      iso.sign.change[i,] <- c(age_SX, "SX", "DX", age_DX)
      
      #iso.sign.change[i,] <- iso.sign[,-1] %>%
      #  mutate(next_sign = lead(sign), next_age = lead(age)) %>%
      #  filter(sign != next_sign)
    }
  }
  pb$tick()
}


# calcolo distanza tra le due isocrone e previsione dell'età come media pesata
colnames(iso.sign.change) <- c("age.SX", "sign.prev", "sign.post", "age.DX")
iso.sign.change$age.SX <- as.numeric(iso.sign.change$age.SX)
iso.sign.change$age.DX <- as.numeric(iso.sign.change$age.DX)
iso.sign.change$dist.SX <- NA # quinta colonna
iso.sign.change$dist.DX <- NA # sesta colonna

star_info <- cbind(df_NGC2244_for_age, iso.sign.change)


dist_beta0.4_NGC2244 <- dist.isocrones.star_3.0(star_info, isocrone_beta0.4_ML2)

df_with_dist_NGC2244 <- iso.sign.change

pb <- progress_bar$new(format = "  Interpolating [:bar] :percent in :elapsed, eta: :eta [:current/11200]", total = nrow(df_with_dist_NGC2244))

for(i in 1:nrow(df_with_dist_NGC2244)) {
  if(df_with_dist_NGC2244[i,1] == 200000000) {
    df_with_dist_NGC2244[i,5] <- 1
    df_with_dist_NGC2244[i,6] <- 1
  } else {
    if(df_with_dist_NGC2244[i,4] == 1000000) {
      df_with_dist_NGC2244[i,5] <- 1
      df_with_dist_NGC2244[i,6] <- 1
    } else {
      df_with_dist_NGC2244[i,5] <- min(dist_beta0.4_NGC2244[[i]][1:which(is.na(dist_beta0.4_NGC2244[[i]]))], na.rm = T)
      df_with_dist_NGC2244[i,6] <- min(dist_beta0.4_NGC2244[[i]][which(is.na(dist_beta0.4_NGC2244[[i]]))+1:length(dist_beta0.4_NGC2244[[i]])], na.rm = T)
    }
  }
  pb$tick()
}

age_NGC2244_beta0.4_ML2_new <- weighted.mean.age(df_with_dist_NGC2244[,5], df_with_dist_NGC2244[,1], df_with_dist_NGC2244[,6], df_with_dist_NGC2244[,4])
save(age_NGC2244_beta0.4_ML2_new, file = "age_NGC2244_beta0.4_ML2_new2.RData")




df_jackson_perturbated_final$age_pert <- NA 

df_jackson_perturbated_final[which(df_jackson_perturbated_final$CLUSTER == "NGC6530"),7] <- age_NGC6530_beta0.4_ML2_new
df_jackson_perturbated_final[which(df_jackson_perturbated_final$CLUSTER == "Cha_I"),7] <- age_Cha_I_beta0.4_ML2_new
df_jackson_perturbated_final[which(df_jackson_perturbated_final$CLUSTER == "Rho_Oph"),7] <- age_Rho_Oph_beta0.4_ML2_new
df_jackson_perturbated_final[which(df_jackson_perturbated_final$CLUSTER == "NGC2264"),7] <- age_NGC2264_beta0.4_ML2_new
df_jackson_perturbated_final[which(df_jackson_perturbated_final$CLUSTER == "NGC2244"),7] <- age_NGC2244_beta0.4_ML2_new
df_jackson_perturbated_final[which(df_jackson_perturbated_final$CLUSTER == "lam_Ori_B35"),7] <- age_lam_Ori_B35_beta0.4_ML2_new
df_jackson_perturbated_final[which(df_jackson_perturbated_final$CLUSTER == "Assc50"),7] <- age_Assc50_beta0.4_ML2_new
df_jackson_perturbated_final[which(df_jackson_perturbated_final$CLUSTER == "Col197"),7] <- age_Col197_beta0.4_ML2_new
df_jackson_perturbated_final[which(df_jackson_perturbated_final$CLUSTER == "NGC2451a"),7] <- age_NGC2451a_beta0.2_ML2_new
df_jackson_perturbated_final[which(df_jackson_perturbated_final$CLUSTER == "NGC6405"),7] <- age_NGC6405_beta0.2_ML2_new

df_jackson_perturbated_final[which(df_jackson_perturbated_final$CLUSTER == "Blanco1"),7] <- age_Blanco1_beta0.2_ML2
df_jackson_perturbated_final[which(df_jackson_perturbated_final$CLUSTER == "25_Ori"),7] <- age_25_Ori_beta0.4_ML2
df_jackson_perturbated_final[which(df_jackson_perturbated_final$CLUSTER == "lam_Ori"),7] <- age_lam_Ori_beta0.4_ML2
df_jackson_perturbated_final[which(df_jackson_perturbated_final$CLUSTER == "NGC2232"),7] <- age_NGC2232_beta0.4_ML2
df_jackson_perturbated_final[which(df_jackson_perturbated_final$CLUSTER == "NGC2451b"),7] <- age_NGC2451b_beta0.2_ML2_new
df_jackson_perturbated_final[which(df_jackson_perturbated_final$CLUSTER == "NGC2516"),7] <- age_NGC2516_beta0.2_ML2_new
df_jackson_perturbated_final[which(df_jackson_perturbated_final$CLUSTER == "gamma2_Vel"),7] <- age_gamma2_Vel_beta0.4_ML2
df_jackson_perturbated_final[which(df_jackson_perturbated_final$CLUSTER == "NGC2547"),7] <- age_NGC2547_beta0.2_ML2
df_jackson_perturbated_final[which(df_jackson_perturbated_final$CLUSTER == "IC2391"),7] <- age_IC2391_beta0.2_ML2
df_jackson_perturbated_final[which(df_jackson_perturbated_final$CLUSTER == "IC2602"),7] <- age_IC2602_beta0.2_ML2
df_jackson_perturbated_final[which(df_jackson_perturbated_final$CLUSTER == "IC4665"),7] <- age_IC4665_beta0.4_ML2
df_jackson_perturbated_final[which(df_jackson_perturbated_final$CLUSTER == "NGC6709"),7] <- age_NGC6709_beta0.4_ML2

write.csv(df_jackson_perturbated_final, "df_perturbated_final.csv", row.names = F, quote = F)





head(age_interp_12iso_group)
df_age_interp_ML2 <- age_interp_12iso_group[,c(1,2,3,4,12,21,30,39)]

cluster <- c("NGC6530", "Cha_I", "Rho_Oph", "NGC2264", "NGC2244", "lam_Ori", "lam_Ori_B35", "25_Ori", "Assc50",
             "Col197", "gamma2_Vel", "IC4665", "NGC2232", "NGC2547", "IC2602", "NGC2451b", "IC2391", "NGC2451a",
             "NGC6405", "NGC2516", "Blanco1", "NGC6709")
beta <- c(0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.4)
N <- c(325, 76, 36, 471, 112, 179, 43, 145, 175, 98, 162, 31, 44, 107, 48, 56, 31, 38, 49, 332, 83, 39)

sd.age <- c(sd(log10(df_age_interp_ML2$b0.4_ML2_Est[which(df_age_interp_ML2$CLUSTER == "NGC6530")])),
                    sd(log10(df_age_interp_ML2$b0.4_ML2_Est[which(df_age_interp_ML2$CLUSTER == "Cha_I")])),
                    sd(log10(df_age_interp_ML2$b0.4_ML2_Est[which(df_age_interp_ML2$CLUSTER == "Rho_Oph")])),
                    sd(log10(df_age_interp_ML2$b0.4_ML2_Est[which(df_age_interp_ML2$CLUSTER == "NGC2264")])),
                    sd(log10(df_age_interp_ML2$b0.4_ML2_Est[which(df_age_interp_ML2$CLUSTER == "NGC2244")])),
                    sd(log10(df_age_interp_ML2$b0.4_ML2_Est[which(df_age_interp_ML2$CLUSTER == "lam_Ori")])),
                    sd(log10(df_age_interp_ML2$b0.4_ML2_Est[which(df_age_interp_ML2$CLUSTER == "lam_Ori_B35")])),
                    sd(log10(df_age_interp_ML2$b0.4_ML2_Est[which(df_age_interp_ML2$CLUSTER == "25_Ori")])),
                    sd(log10(df_age_interp_ML2$b0.4_ML2_Est[which(df_age_interp_ML2$CLUSTER == "Assc50")])),
                    sd(log10(df_age_interp_ML2$b0.4_ML2_Est[which(df_age_interp_ML2$CLUSTER == "Col197")])),
                    sd(log10(df_age_interp_ML2$b0.4_ML2_Est[which(df_age_interp_ML2$CLUSTER == "gamma2_Vel")])),
                    sd(log10(df_age_interp_ML2$b0.4_ML2_Est[which(df_age_interp_ML2$CLUSTER == "IC4665")])),
                    sd(log10(df_age_interp_ML2$b0.4_ML2_Est[which(df_age_interp_ML2$CLUSTER == "NGC2232")])),
                    sd(log10(df_age_interp_ML2$b0.2_ML2_Est[which(df_age_interp_ML2$CLUSTER == "NGC2547")])),
                    sd(log10(df_age_interp_ML2$b0.2_ML2_Est[which(df_age_interp_ML2$CLUSTER == "IC2602")])),
                    sd(log10(df_age_interp_ML2$b0.2_ML2_Est[which(df_age_interp_ML2$CLUSTER == "NGC2451b")])),
                    sd(log10(df_age_interp_ML2$b0.2_ML2_Est[which(df_age_interp_ML2$CLUSTER == "IC2391")])),
                    sd(log10(df_age_interp_ML2$b0.2_ML2_Est[which(df_age_interp_ML2$CLUSTER == "NGC2451a")])),
                    sd(log10(df_age_interp_ML2$b0.2_ML2_Est[which(df_age_interp_ML2$CLUSTER == "NGC6405")])),
                    sd(log10(df_age_interp_ML2$b0.2_ML2_Est[which(df_age_interp_ML2$CLUSTER == "NGC2516")])),
                    sd(log10(df_age_interp_ML2$b0.2_ML2_Est[which(df_age_interp_ML2$CLUSTER == "Blanco1")])),
                    sd(log10(df_age_interp_ML2$b0.4_ML2_Est[which(df_age_interp_ML2$CLUSTER == "NGC6709")])))


# calcolo le stesse statistiche ma sui valori di eta' perturbati
df_perturbated_final <- df_jackson_perturbated_final

cluster_stats <- as.data.frame(matrix(NA, ncol = 8, nrow = 22))
colnames(cluster_stats) <- c("cluster", "mean_age_err", "sd_age_err", "Q16_age_err", "Q84_age_err", "median_age_err",
                             "MAD_age_err", "sd_sd_age_err")
cluster_stats$cluster <- cluster
i = "NGC6709"
for(i in cluster) {
  df_perturbated_cluster <- df_perturbated_final[which(df_perturbated_final$CLUSTER == i),]
  df_stats <- df_perturbated_cluster %>% group_by(ges_id_gaia) %>% summarise(mean_age_err = mean(log10(age_pert)),
                                                                             sd_age_err = sd(log10(age_pert)),
                                                                             Q16_age_err = quantile(log10(age_pert), 0.16),
                                                                             Q84_age_err = quantile(log10(age_pert), 0.84),
                                                                             median_age_err = median(log10(age_pert)),
                                                                             MAD_age_err = mad(log10(age_pert)))
  df_stats <- as.data.frame(df_stats)
  j <- which(i == cluster_stats[,1])
  cluster_stats[j,2:6] <- colMeans(df_stats[,c(2:6)])
  cluster_stats[j,7] <- mean(df_stats[which(df_stats$MAD_age_err > 0),7])
  cluster_stats[j,8] <- sd(df_stats$sd_age_err)
}





df_stats_clusters <- data.frame(cluster = cluster, N = N, beta = beta, mean_age_est = mean.age, sd_age = sd.age,
                                Q16_age = quantile16.age, Q84_age = quantile84.age,
                                median_age_est = median.age, MAD_age_est = mad.age,
                                L_bound_est_tot = median.age - mad.age, U_bound_est_tot = median.age + mad.age)

df_stats_clusters_final <- cbind(df_stats_clusters, cluster_stats[,-1])
df_stats_clusters_final$sd_age > df_stats_clusters_final$sd_age_err
df_stats_clusters_final$MAD_age_est > df_stats_clusters_final$MAD_age_err


# inserisco le MAD*k, dove k e' la costante di normalizazzione per confrontare la MAD con la std.dev
k = 1.4826
df_stats_clusters_final$MAD_age_est_k <- df_stats_clusters_final$MAD_age_est * k
df_stats_clusters_final$MAD_age_err_k <- df_stats_clusters_final$MAD_age_err * k


# eta' di letteratura (Jeffries) ed intervallo di confidenza

Lit_age <- c(6.7, 6.26, 6.91, 6.78, 6.76, 7, 6.94, 7.18, 6.85, 6.85, 7.2, 7.76, 7.45, 7.53, 7.63, 7.64,
             7.8, 7.86, 7.88, 8.14, 7.86, 8.15)

L_bound_Lit <- c(6, NA, 6, 6.26, 6, 6.56, NA, 7.17, 6, 6.32, 7.19, 7.72, 7.43, 7.52, 7.6, 7.61, 7.72, 7.8,
                 7.84, 8.13, 7.82, 8.11)

U_bound_Lit <- c(6.7, NA, 6.91, 6.83, 6.76, 7.02, NA, 7.19, 6.85, 6.92, 7.21, 7.81, 7.47, 7.55, 7.67, 7.67,
                 7.96, 7.93, 7.93, 8.16, 7.93, 8.19)


# eta' di training utilizzate da Jeffries nel suo paper

Train_age <- c(6.34, 6.26, 6.57, 6.65, 6.78, 6.94, 6.94, 7.17, 6.91, 7.07, 7.21, 7.54, 7.44, 7.55, 7.63,
               7.58, 7.62, 7.7, 7.8, 8.32, 8.04, 8.22)

Train_err <- c(0.36, NA, 0.12, 0.18, 0.29, 0.16, NA, 0.1, 0.15, 0.1, 0.11, 0.18, 0.17, 0.04, 0.06, 0.09,
               0.14, 0.13, 0.23, 0.16, 0.05, 0.05)

L_bound_Train <- Train_age - Train_err
U_bound_Train <- Train_age + Train_err

df_stats_clusters_final <- cbind(df_stats_clusters_final, Lit_age, L_bound_Lit, U_bound_Lit, Train_age, L_bound_Train, U_bound_Train)
df_stats_clusters_final_new <- df_stats_clusters_final

write.csv(df_stats_clusters_final_new, "df_stats_clusters_final_new.csv", row.names = F, quote = F)









