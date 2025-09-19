

### Custom function to predict the age of the stars ###

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
