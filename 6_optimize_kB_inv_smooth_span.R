rm(list = ls())
cat("\14")
graphics.off()

require(phenocamr)
require(lubridate)

wdir <- "/Volumes/GoogleDrive/My Drive/Young_aerodynamic_resistance_analysis"

setwd(paste0(wdir,"/data/ancillary_data"))
phenoflux_metadata_table <- read.csv("pheno_flux_sites_to_use.csv")
phenoflux_metadata_table$vegtype[phenoflux_metadata_table$fluxsite == "US-Ro4"] <- "GR"

sites <- phenoflux_metadata_table$fluxsite

span_vals_all <- matrix(NA,nrow=length(sites),ncol=4)

setwd(paste0(wdir,"/results/5_resistance_values"))

for (i in 1:length(sites)){

  if (sites[i] == "US-Bi1" | sites[i] == "US-Tw3") {
    
    next
    
  }
  
  fluxdat <- read.csv(sprintf("%s_resistance_values.csv",sites[i]))
  fluxdat[fluxdat == -9999] <- NA
  
  unique_yr = unique(year(as.Date(fluxdat$date)))
  
  x <- 1:nrow(fluxdat)
  
  y <- matrix(NA,nrow=nrow(fluxdat),ncol=2)
  
  y[,1] <- fluxdat$r_h_pred
  y[,2] <- fluxdat$kB_inv
  
  span_vals <- numeric(length = 2)
  
  smooth_vals = matrix(NA,nrow=nrow(fluxdat),ncol=2)
  
  for (j in 1:length(span_vals)){
    
    span_j <- optimal_span(y[,j],x,step = 0.001,plot = FALSE)
    
    if (i == 13){
      
      span_j <- 0.2
      
    }
    
    if (span_j > 0 & span_j < 0.4){
      
      smooth_vals_j <- loess(y[,j] ~ x,span = span_j)
      smooth_vals[,j] <- round(predict(smooth_vals_j,newdata = x),digits = 2)
      span_vals_all[i,j] <- span_j
      
      
    } else {

      smooth_vals_j <- loess(y[,j] ~ x,span = 0.4/length(unique_yr))
      smooth_vals[,j] <- round(predict(smooth_vals_j,newdata = x),digits = 2)

    }

  }
  
  fluxdat$smooth_r_h_pred = smooth_vals[,1]
  fluxdat$smooth_kB_inv = smooth_vals[,2]

  fluxdat[is.na(fluxdat)] = -9999;
  
  write.csv(x = fluxdat, 
            sprintf("%s_resistance_values.csv",sites[i]),
            append = FALSE,
            row.names = FALSE)  
  
}
