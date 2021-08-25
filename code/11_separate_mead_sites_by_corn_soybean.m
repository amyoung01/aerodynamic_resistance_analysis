clear; clc; close all;

wdir = '/Volumes/GoogleDrive/My Drive/Young_aerodynamic_resistance_analysis';

addpath([wdir,'/code/z_functions']);

cd([wdir,'/data/ancillary_data']);
mead_crops = readtable('mead2_mead3_crop_types.csv');

sites = {'US-Ne1','US-Ne2','US-Ne3'};
phenos = {'mead1','mead2','mead3'};

for i = 1:length(sites)
    
    cd([wdir,'/results/5_resistance_values']);
    fluxdat = readtable(sprintf('%s_resistance_values.csv',char(sites(i))));
    fluxdat = standardizeMissing(fluxdat,-9999);
    
    cd([wdir,'//results/3_processed_phenocam_data/time_series']);
    pheno_ts = readtable(sprintf('%s_gcc_time_series.csv',char(phenos(i))));
    
    cd([wdir,'//results/3_processed_phenocam_data/transition_dates']);
    pheno_td = readtable(sprintf('%s_gcc_transition_dates.csv',char(phenos(i))));
    
    if strcmp(sites(i),'US-Ne1')
        
        crop = repmat({'corn'},height(fluxdat),1);
        crop_pheno = repmat({'corn'},height(pheno_ts),1);
        crop_td = repmat({'corn'},height(pheno_td),1);        
        
    else
        
        crop = cell(height(fluxdat),1);
        crop_pheno = cell(height(pheno_ts),1);
        crop_td = cell(height(pheno_td),1);
                
        yr = year(fluxdat.date);
        unique_yr = mead_crops.year;
        
        for y = 1:length(unique_yr)
            
            yr_id = yr == unique_yr(y);
            
            if strcmp(sites(i),'US-Ne2')
                
                crop(yr_id) = table2array(mead_crops(y,2));                
                crop_pheno(year(pheno_ts.date) == unique_yr(y)) = table2array(mead_crops(y,2));
                crop_td(year(pheno_td.transition_10) == unique_yr(y)) = table2array(mead_crops(y,2));
                
            else
                
                crop(yr_id) = table2array(mead_crops(y,3));
                crop_pheno(year(pheno_ts.date) == unique_yr(y)) = table2array(mead_crops(y,3));
                crop_td(year(pheno_td.transition_10) == unique_yr(y)) = table2array(mead_crops(y,3));                
                
            end
            
            
        end
        
    end
    
    fluxdat.crop = crop;
    pheno_ts.crop = crop_pheno;
    pheno_td.crop = crop_td;
    
    fluxdat = setNaN(fluxdat,-9999);
    pheno_ts = setNaN(pheno_ts,-9999);
    pheno_td = setNaN(pheno_td,-9999);
    
    cd([wdir,'/results/5_resistance_values']);
    writetable(fluxdat,sprintf('%s_resistance_values.csv',char(sites(i))));
    
    cd([wdir,'//results/3_processed_phenocam_data/time_series']);
    writetable(pheno_ts,sprintf('%s_gcc_time_series.csv',char(phenos(i))));
    
    cd([wdir,'//results/3_processed_phenocam_data/transition_dates']);
    writetable(pheno_td,sprintf('%s_gcc_transition_dates.csv',char(phenos(i))));    
    
end
    
for j = 1:3
    
    for i = 1:length(sites)
        
        if j == 1
            
            cd([wdir,'/results/5_resistance_values']);
            data_j = readtable(sprintf('%s_resistance_values.csv',char(sites(i))));
            data_j = standardizeMissing(data_j,-9999);            
            
        elseif j == 2
            
            cd([wdir,'//results/3_processed_phenocam_data/time_series']);
            data_j = readtable(sprintf('%s_gcc_time_series.csv',char(phenos(i))));
            
        else
            
            cd([wdir,'//results/3_processed_phenocam_data/transition_dates']);
            data_j = readtable(sprintf('%s_gcc_transition_dates.csv',char(phenos(i))));
            
        end
        
        if i == 1
            
            data_j.site = repmat({'US-Ne1'},height(data_j),1);
            
        elseif i == 2
            
            data_j.site = repmat({'US-Ne2'},height(data_j),1);
            
        else
            
            data_j.site = repmat({'US-Ne3'},height(data_j),1);
            
        end
        
        if i == 1
            
            data_j_corn = data_j;
            continue;
            
        end
        
        corn_id = strcmp(data_j.crop,'corn');
        soyb_id = strcmp(data_j.crop,'soybean');
        
        data_j_corn = [data_j_corn; data_j(corn_id,:)];
        
        if i == 2
            
            data_j_soyb = data_j(soyb_id,:);
            
        elseif i == 3 
            
            data_j_soyb = [data_j_soyb; data_j(soyb_id,:)];
            
        end

    end
    
    if j == 1
        
        cd([wdir,'/results/5_resistance_values']);
        data_j_corn = setNaN(data_j_corn,-9999);
        writetable(data_j_corn,'US-Ne-corn_resistance_values.csv');
        
        data_j_soyb = setNaN(data_j_soyb,-9999);
        writetable(data_j_soyb,'US-Ne-soybean_resistance_values.csv');
        
    elseif j == 2
        
        cd([wdir,'/results/3_processed_phenocam_data/time_series']);
        data_j_corn = setNaN(data_j_corn,-9999);
        writetable(data_j_corn,'mead_corn_gcc_time_series.csv');
        
        data_j_soyb = setNaN(data_j_soyb,-9999);
        writetable(data_j_soyb,'mead_soybean_gcc_time_series.csv');        
        
    else
        
        cd([wdir,'/results/3_processed_phenocam_data/transition_dates']);
        data_j_corn = setNaN(data_j_corn,-9999);
        writetable(data_j_corn,'mead_corn_gcc_transition_dates.csv');
        
        data_j_soyb = setNaN(data_j_soyb,-9999);
        writetable(data_j_soyb,'mead_soybean_gcc_transition_dates.csv');
        
    end
    
end
