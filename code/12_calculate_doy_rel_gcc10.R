rm(list = ls()) # Clear all data from workspace
graphics.off() # Close all current figures and plots
cat("\14") # Clear command console

# Load required packages
require(circular)  # Circular Statistics (version 0.4.93)
require(lubridate) # Make Dealing with Dates a Little Easier (version 1.7.10)

# Set working directory path
wdir <- "/Volumes/GoogleDrive/My Drive/Young_aerodynamic_resistance_analysis"

# Load in metadata table with site-specific information
setwd(paste0(wdir,"/data/ancillary_data"))
phenoflux_metadata_table <- read.csv("pheno_flux_sites_to_use.csv")

sites  <- c(as.character(phenoflux_metadata_table$fluxsite),"US-Ne-corn","US-Ne-soybean")
phenos <- c(as.character(phenoflux_metadata_table$phenosite),"mead_corn","mead_soybean")

for (i in 1:length(sites)){
  
   # Ignore the alfalfa sites for this analysis
   if (sites[i] == "US-Bi1" | sites[i] == "US-Tw3"){

    next

  }
  
  setwd(paste0(wdir,"/results/5_resistance_values"))
  fluxdat <- read.csv(sprintf("%s_resistance_values.csv",sites[i]))

  setwd(paste0(wdir,"/results/3_processed_phenocam_data/time_series"))
  pheno_ts <- read.csv(sprintf("%s_gcc_time_series.csv",phenos[i]))

  setwd(paste0(wdir,"/results/3_processed_phenocam_data/transition_dates"))
  pheno_td <- read.csv(sprintf("%s_gcc_transition_dates.csv",phenos[i]))
  pheno_td <- pheno_td[pheno_td$direction == "rising",]

  flux_doy <- lubridate::yday(as.Date(fluxdat$date))
  pheno_ts_doy <- lubridate::yday(as.Date(pheno_ts$date))
  
  gcc_10_doy <- lubridate::yday(as.Date(pheno_td$transition_10))
  
  flux_rad <- circular::as.circular(rad(360*(flux_doy/365)),unit="radians")
  pheno_ts_rad <- circular::as.circular(rad(360*(pheno_ts_doy/365)),unit="radians")
  gcc_10_rad <- circular::as.circular(rad(360*(gcc_10_doy/365)),unit="radians")
  
  avg_gcc_10 <- circular::mean.circular(gcc_10_rad,na.rm = TRUE)
  
  flux_rel_doy <- round(365 * (as.numeric(deg(flux_rad - avg_gcc_10))/360))
  pheno_ts_rel_doy <- round(365 * (as.numeric(deg(pheno_ts_rad - avg_gcc_10))/360))
  
  max_val = max(flux_rel_doy,na.rm = TRUE)
  min_val = min(flux_rel_doy,na.rm = TRUE)
  
  if (max_val > 182) {
    
    doy_gt_182 <- flux_rel_doy > 182
    flux_rel_doy[doy_gt_182] = min_val - (max_val - flux_rel_doy[doy_gt_182])
    
  } else if (min_val < -182) {
    
    doy_gt_182 = flux_rel_doy < -182
    flux_rel_doy[doy_gt_182] = max_val + (flux_rel_doy[doy_gt_182] - min_val)
    
  }

  max_val <- max(pheno_ts_rel_doy,na.rm = TRUE)
  min_val <- min(pheno_ts_rel_doy,na.rm = TRUE)
  
  if (max_val > 182) {
    
    doy_gt_182 <- pheno_ts_rel_doy > 182
    pheno_ts_rel_doy[doy_gt_182] = min_val - (max_val - pheno_ts_rel_doy[doy_gt_182])
    
  } else if (min_val < -182) {
    
    doy_gt_182 = pheno_ts_rel_doy < -182
    pheno_ts_rel_doy[doy_gt_182] = max_val + (pheno_ts_rel_doy[doy_gt_182] - min_val)
    
  }  
  
  fluxdat$rel_doy <- flux_rel_doy
  pheno_ts$rel_doy <- pheno_ts_rel_doy
  
  setwd(paste0(wdir,"/results/5_resistance_values"))
  write.csv(fluxdat,sprintf("%s_resistance_values.csv",sites[i]),row.names = FALSE)

  setwd(paste0(wdir,"/results/3_processed_phenocam_data/time_series"))
  write.csv(pheno_ts,sprintf("%s_gcc_time_series.csv",phenos[i]),row.names = FALSE)
    
}
