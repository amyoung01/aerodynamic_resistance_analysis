rm(list = ls())
cat("\14")
graphics.off()

library(anytime)
require(circular)
require(phenocamr)
library(lubridate)
library(boot)

wdir <- "/Volumes/GoogleDrive/My Drive/Young_aerodynamic_resistance_analysis"

setwd(paste0(wdir,"/data/ancillary_data"))
phenoflux_metadata_table <- read.csv("pheno_flux_sites_to_use.csv")
phenoflux_metadata_table$vegtype[phenoflux_metadata_table$fluxsite == "US-Ro4"] <- "GR"


sites <- c(as.character(phenoflux_metadata_table$fluxsite),"US-Ne-corn","US-Ne-soybean")
phenos <- c(as.character(phenoflux_metadata_table$phenosite),"mead_corn","mead_soybean")
vegtype <- c(as.character(phenoflux_metadata_table$vegtype),"AG","AG")
vegtype[sites == "US-Ton"] <- "SA"
ra_vars <- c("smooth_kB_inv")

export_veg <- c()
export_sites <- c()
export_dates <- c()

bool_site_id <- logical(length = length(sites))

dates <- c()

td_sites <- c()

td_boot_results <- matrix(NA,nrow = length(sites),ncol = 3)
R <- 200

for (i in 1:length(sites)){
  
  # Don't process for alfalfa sites
  if (sites[i] == "US-Bi1" | sites[i] == "US-Tw3"){

    next

  }
  
  bool_site_id[i] <- TRUE
  
  setwd(paste0(wdir,"/results/5_resistance_values"))
  data <- read.csv(sprintf("%s_resistance_values.csv",sites[i]))
  data[data == -9999] <- NA
  
  if (sites[i] == "US-Ne-soybean"){
   
    data$smooth_kB_inv[data$smooth_kB_inv < -2] <- NA
     
  }
  
  setwd(paste0(wdir,"/results/3_processed_phenocam_data/transition_dates"))
  td <- read.csv(sprintf("%s_gcc_transition_dates.csv",phenos[i]))

  td_10 <- td$transition_10
  td_90 <- td$transition_90
  
  td_10_deg = 360*(yday(td_10)/365)
  td_10_rad = circular::rad(td_10_deg)
  td_sites <- c(td_sites,sites[i])
  
  fx <- function(td_10_rad,indices) median.circular(td_10_rad[indices],units="radians",na.rm = TRUE)
  
  td_mean_boot <- boot(td_10_rad,statistic = fx,R = R)
  td_ci_boot <- quantile.circular(td_mean_boot$t,probs=c(0.025,0.975),na.rm=TRUE)
  
  td_boot_results[i,2] <- as.numeric(td_mean_boot[1])
  td_boot_results[i,-2] <- as.numeric(td_ci_boot)
  
  for (j in 1:length(td_10)){
    
    # if (i == 24 | i == 25){
    # } else {
    # }
    
    reldays <- as.numeric(as.Date(data$date) - as.Date(td_10[j]))
    id <- reldays >= -182 & reldays <= 182
    
    reldays <- reldays[id]
    radians <- rad(360 * ((reldays + 183)/365))
    n = length(radians)
      
    export_hi_vals_ij <- matrix(NA,nrow=n,ncol=length(ra_vars))
    export_lo_vals_ij <- matrix(NA,nrow=n,ncol=length(ra_vars))
    
    dates <- c(dates,as.Date(data$date[id]))
    
    for (v in 1:length(ra_vars)){
     
      vals <- data[id,ra_vars[v]]
      
      if (sum(is.na(vals)) == n){
        
        next
        
      }
      
      prctile_vals <- quantile(vals,
                               probs = c(0.25,0.75),
                               na.rm = TRUE)
      
      lo_id <- which(vals < prctile_vals[1])
      hi_id <- which(vals > prctile_vals[2])
    
      export_hi_vals_ij[hi_id,v] <- radians[hi_id]
      export_lo_vals_ij[lo_id,v] <- radians[lo_id]
      
    }
    
    if (i == 1 & j == 1){
      
      export_hi_vals <- export_hi_vals_ij
      export_lo_vals <- export_lo_vals_ij
      
      colnames(export_hi_vals) <- paste0("hi_",ra_vars)
      colnames(export_lo_vals) <- paste0("lo_",ra_vars)
      
    } else {
      
      export_hi_vals <- rbind(export_hi_vals,export_hi_vals_ij)
      export_lo_vals <- rbind(export_lo_vals,export_lo_vals_ij)
      
    }
    
    export_sites <- c(export_sites,rep(as.character(sites[i]),n))
    export_veg <- c(export_veg,rep(as.character(vegtype[i]),n))
    
  }
  
   td_diff_num <- as.numeric(as.Date(td_90) - as.Date(td_10))
   td_diff_rad <- rad(360 * ((td_diff_num + 183)/365))
    
   if (i == 1){
     
      td_diff <- as.numeric(mean.circular(as.circular(td_diff_rad,units="radians"),na.rm = TRUE))
     
   } else {
     
     td_diff_i <- as.numeric(mean.circular(as.circular(td_diff_rad,units="radians"),na.rm = TRUE))
     td_diff <- c(td_diff,td_diff_i)
     
   }
  
}

export_table <- cbind(data.frame(site = export_sites,date = anydate(dates),vegtype = export_veg),
                      as.data.frame(export_hi_vals),
                      as.data.frame(export_lo_vals))

hi_agg <- aggregate(export_table[paste0("hi_",ra_vars)],
                    by=list(export_table$site),
                    FUN = quantile.circular,
                    na.rm = TRUE,
                    probs = c(0.1,0.5,0.9))

lo_agg <- aggregate(export_table[paste0("lo_",ra_vars)],
                    by=list(export_table$site),
                    FUN = quantile.circular,
                    na.rm = TRUE,
                    probs = c(0.1,0.5,0.9))

boot_fx <- function(x,indices){
  
  dt <- as.circular(x[indices],units="radians")
  c(as.numeric(quantile.circular(dt,probs=0.9,na.rm=TRUE)))
  
}

summary_stats <- data.frame(site = hi_agg$Group.1,
                            vegtype = vegtype[bool_site_id][sort.int(sites[bool_site_id],index.return = TRUE)$ix],
                            td_diff = td_diff,
                            hi_10 = as.data.frame(hi_agg[paste0("hi_",ra_vars)])[[1]][,1],
                            hi_50 = as.data.frame(hi_agg[paste0("hi_",ra_vars)])[[1]][,2],
                            hi_90 = as.data.frame(hi_agg[paste0("hi_",ra_vars)])[[1]][,3],
                            lo_10 = as.data.frame(lo_agg[paste0("lo_",ra_vars)])[[1]][,1],
                            lo_50 = as.data.frame(lo_agg[paste0("lo_",ra_vars)])[[1]][,2],
                            lo_90 = as.data.frame(lo_agg[paste0("lo_",ra_vars)])[[1]][,3])

transition_hi2lo <- numeric(length = nrow(summary_stats))

for (i in 1:length(transition_hi2lo)){
 
 hi2lo_i <- mean.circular(c(summary_stats$hi_90[i],summary_stats$lo_10[i]),units = "radians")
 transition_hi2lo[i] <- hi2lo_i
  
}

transition_hi2lo <- ifelse(transition_hi2lo < 0,2 * pi - abs(transition_hi2lo),transition_hi2lo)

summary_stats$hi2lo <- transition_hi2lo

B <- matrix(NA,nrow = length(sites),ncol = 3)

for (i in 1:length(sites)){
  
  # Don't process for alfalfa sites
  if (sites[i] == "US-Bi1" | sites[i] == "US-Tw3"){

    next

  }
 
  id = export_table$site == sites[i]
  
  df = export_table[id,4]
  b_i = boot(df,boot_fx,R = R)
  
  B[i,2] = median.circular(as.circular(b_i$t,units = "radians"))
  B[i,-2] = as.numeric(quantile.circular(b_i$t,probs=c(0.025,0.975)))
  
}

td_boot_results <- ifelse(td_boot_results < 0,2*pi - abs(td_boot_results),td_boot_results)
B <- ifelse(B<0,2*pi - abs(B),B)

export_boot_df <- as.data.frame(cbind(td_boot_results,B))
colnames(export_boot_df) <- c("gcc_loci","gcc_median","gcc_hici",
                              paste0(ra_vars,"_loci"),paste0(ra_vars,"_median"),paste0(ra_vars,"_hici"))

export_boot_df$site <- sites
export_boot_df$vegtype <- vegtype

export_boot_df <- export_boot_df[,c(ncol(export_boot_df)-1,ncol(export_boot_df),1:(ncol(export_boot_df)-2))]

export_boot_df[is.na(export_boot_df)] <- -9999
export_table[is.na(export_table)] <- -9999

write.csv(export_table,
          file = paste0(wdir,"/results/6_seasonality_stats/seasonality_",ra_vars,".csv"),
          row.names = FALSE)

write.csv(summary_stats,file = paste0(wdir,"/results/6_seasonality_stats/summary_stats_",ra_vars,".csv"),row.names = FALSE)

write.csv(export_boot_df,file = paste0(wdir,"/results/6_seasonality_stats/bootstrap_results_",ra_vars,".csv"),row.names = FALSE)
