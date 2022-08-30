# This script goes through and organizes and gathers PhenoCam Gcc data into what
# is needed for the aerodynamic resistance analysis. Specifically, it identifies 
# the Gcc statistic (e.g., mean, median, 75th percentile, 90th percentile) that
# minimizes the RMSE. This is the Gcc time series we use. We also recalculate
# transition dates for when the Gcc magnigutde is at 10%, 50%, and 90% percent,
# instead of the published 10%, 25%, and 50%.

rm(list=ls())
graphics.off()

library(phenocamr)
library(anytime)

wdir <- "/Volumes/GoogleDrive/My Drive/Young_aerodynamic_resistance_analysis"

setwd(paste0(wdir,"/data/ancillary_data"))
phenocam_flux_metadata_table <- read.csv("pheno_flux_sites_to_use.csv")

phenos <- phenocam_flux_metadata_table$phenosite
vegtype <- phenocam_flux_metadata_table$vegtype

rmse_stats <- c("mean","50","75","90")

for (i in 1:length(phenos)){
    
    setwd(sprintf('%s//data//raw_data//phenocam//%s',wdir,as.character(phenos[i])))
    pheno_dir <- getwd()
    versions <- files <- list.files()
    
    # Go through each version for each site. There are different versions depending
    # on if the field of view of the camera changes.
    for (v in 1:length(versions)){
        
        setwd(sprintf('%s//%s',pheno_dir,versions[v]));
        file_to_read = paste0(sprintf("%s_%s_1day_transition_dates.csv",
                                      phenos[i],versions[v]))
        
        # The next several lines reads in the RMSE data from the file 
        # header. We find which statistic minimizes this RMSE.
        fid <- file(description = file_to_read,open = "rt")
        header_lines <- readLines(con = fid,n = 15)
        close(fid)
        
        header_lines <- header_lines[12:15]
        rmse <- NA * numeric(length = length(header_lines))
        
        for (r in 1:4) {
          
          str_r <- unlist(strsplit(header_lines[r]," "))
          rmse[r] <- as.numeric(str_r[5])
          
        }
        
        min_rmse <- which(rmse == min(rmse))
        min_rmse <- ifelse(length(min_rmse) > 1,min_rmse[1],min_rmse)
        
        # Read in PhenoCam time series data for 1-day step
        pheno_td <- read_phenocam(filename = sprintf("%s_%s_1day.csv",
                                                     phenos[i],versions[v]))
        
        # identify transition dates fro 10%, 50%, and 90% of the annual Gcc
        # magnitude.
        new_dates <- transition_dates(pheno_td,
                                      lower_thresh = 0.1,
                                      middle_thresh = 0.5,
                                      upper_thresh = 0.9,
                                      percentile = rmse_stats[min_rmse],
                                      reverse = FALSE,
                                      plot = FALSE)
        
        # Store new dates
        new_dates$direction <- "rising"
        new_dates <- new_dates[,c(ncol(new_dates),1:9)]
        new_dates[,2] <- anydate(new_dates[,2])
        new_dates[,3] <- anydate(new_dates[,3])
        new_dates[,4] <- anydate(new_dates[,4])
        new_dates[,5] <- anydate(new_dates[,5])
        new_dates[,6] <- anydate(new_dates[,6])
        new_dates[,7] <- anydate(new_dates[,7])
        new_dates[,8] <- anydate(new_dates[,8])
        new_dates[,9] <- anydate(new_dates[,9])
        new_dates[,10] <- anydate(new_dates[,10])

        if (v == 1){
            
            gcc_trans_dates <- new_dates

        } else {

            gcc_trans_dates <- rbind(gcc_trans_dates, new_dates)
            
        }
        
        # Read in three day data and use only time series that minimizes RMSE
        pheno_ts_name <- paste0(sprintf("%s_%s_3day.csv",
                                        phenos[i],versions[v]))
        pheno_ts_v <- read.csv(pheno_ts_name,header = TRUE,sep = ",",skip = 24)
        
        col_id_1 <- which(colnames(pheno_ts_v) == paste0("gcc_",rmse_stats[min_rmse]))
        col_id_2 <- which(colnames(pheno_ts_v) == paste0("smooth_gcc_",rmse_stats[min_rmse]))
        
        # Make simple data frame to store time series data for given site
        new_ts <- data.frame(date = pheno_ts_v$date,
                             gcc = pheno_ts_v[,col_id_1],
                             smooth_gcc = pheno_ts_v[,col_id_2])
        
        if (v == 1){
            
            pheno_ts <- new_ts
            
        } else {
            
            if (as.Date(new_ts$date[1]) < as.Date(pheno_ts$date[nrow(pheno_ts)])){
                
                id <- which(as.Date(new_ts$date) == as.Date(pheno_ts$date[nrow(pheno_ts)]))
                new_ts <- new_ts[(id+1):nrow(new_ts),]
                
            }
            
            pheno_ts <- rbind(pheno_ts, new_ts)
            
        }
        
    }
    
    pheno_ts[is.na(pheno_ts)] <- -9999
   
    # Export time series
    pheno_ts_file_to_write <- sprintf("%s_gcc_time_series.csv",phenos[i])
    write.csv(pheno_ts,paste0(wdir,
                              "/results/3_processed_phenocam_data/time_series/",
                              pheno_ts_file_to_write),
              row.names = FALSE)
    
    pheno_td_file_to_write <- sprintf("%s_gcc_transition_dates.csv",phenos[i])
    write.csv(gcc_trans_dates,
              paste0(wdir,
                     "/results/3_processed_phenocam_data/transition_dates/",
                     pheno_td_file_to_write),
              row.names = FALSE) 
    
}
