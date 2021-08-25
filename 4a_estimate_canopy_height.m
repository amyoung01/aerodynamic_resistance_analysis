clear all; 
clc; 
close all;

wdir = '/Volumes/GoogleDrive/My Drive/Young_aerodynamic_resistance_analysis';

% Add path to load in custom written functions to workspace
addpath([wdir,'/code/z_functions']);

cd([wdir,'/data/ancillary_data']);
phenoflux_metadata = readtable('pheno_flux_sites_to_use.csv','Delimiter',',');

sites = phenoflux_metadata.fluxsite;
vegtype = phenoflux_metadata.vegtype;
primary_wd = phenoflux_metadata.primary_wd;
measurement_heights = phenoflux_metadata.zr;

var_names = readtable('variables_to_import_for_fluxsites.csv','Delimiter',',');

ZL_var_name = var_names.ZL;
zr_var_name = var_names.zr;

obs_canopy_heights = [phenoflux_metadata.hc_low, ...
                      phenoflux_metadata.hc_high];

obs_canopy_heights(obs_canopy_heights == -9999) = NaN;
avg_canopy_heights = nanmean(obs_canopy_heights,2);

wd_range = [primary_wd - 90,primary_wd + 90]; % [degrees]
wd_range(wd_range < 0) = wd_range(wd_range < 0) + 360;
wd_range(wd_range > 360) = wd_range(wd_range > 360) - 360;

k  = 0.41; % Von-Karman Constant

window_size = 3;
prop_missing_allowed = 0.95;

for i = 1:length(sites)

    % -----------------------------------------------  
    zr = measurement_heights(i); % Measurement height for wind speed data at site
    
    if zr == -9999, continue; end
    
    obs_hc = avg_canopy_heights(i); % BADM Reported canopy height
    
    if isnan(obs_hc), obs_hc = 0; end       
    
    % -----------------------------------------------
    cd([wdir,'/results/2_filtered_flux_data/24hr']);
    
    file_to_import = dir([char(sites(i)),'_*.csv']);
    filename = file_to_import.name;
    filename_parts = strsplit(filename,'_');         
    
    fluxdat = readtable(file_to_import.name);
    fluxdat = standardizeMissing(fluxdat,-9999);
            
    if strcmp(zr_var_name(i),'zr'), zr = fluxdat.zr; end
    
    % -----------------------------------------------  
    timestep = strsplit(char(filename_parts(2)),'.');
    timestep = timestep(1);
    
    t = 2;
    if strcmp(timestep,'HR'), t = 1; end      
    
    % -----------------------------------------------      
    n = height(fluxdat);
    n_per_window = 24 * t * window_size;
    n_samp = ceil(n / n_per_window);
    n_missing_allowed = floor(n_per_window * prop_missing_allowed);
    el = window_size - median(1:window_size);
    
    % -----------------------------------------------  
    t_air = fluxdat.t_air; % [C]
    H = fluxdat.H; % [W m-2]
    pressure = fluxdat.pressure; % [kPa]
        
    wind_speed = fluxdat.wind_speed; % [m s^-1]
    wind_direction = fluxdat.wind_direction; % [degrees]
    ustar = fluxdat.ustar; % [m s^-1]
    
    % Find and filter by neutral log-wind profile conditions ------------------------
    
    if strcmp(ZL_var_name(i),'NA')
        
        L = Monin_Obukhov_Length(t_air,ustar,H,pressure,false); % [m]
        
        if strcmp(vegtype(i),'DB') || strcmp(vegtype(i),'EN')
            
            zeta = (zr - 0.7 * obs_hc) ./ L;
            
        else
            
            zeta = zr ./ L;
            
        end
        
    else
        
        zeta = fluxdat.ZL;
        
    end
    
    zeta_neutral_bool = abs(zeta) < 0.1;
    
    % -----------------------------------------------    
    max_ustar = 0.6;
    if strcmp(vegtype(i),'DB') || strcmp(vegtype(i),'EN'), max_ustar = 1.0; end 
    ustar_bool = ustar >= 0.2 & ustar <= max_ustar;
    
    % -----------------------------------------------  
    if wd_range(i,1) < wd_range(i,2)
        
        wd_bool = wind_direction > wd_range(i,1) & wind_direction < wd_range(i,2);
        
    else
        
        wd_bool = wind_direction > wd_range(i,1) | wind_direction < wd_range(i,2);
        
    end
    
    % -----------------------------------------------  
    to_keep_id = zeta_neutral_bool & ustar_bool & wd_bool;
    
    % -----------------------------------------------  
    ws_ustar_ratio_obs = k * wind_speed ./ ustar;
    ws_ustar_ratio_obs(~to_keep_id) = NaN;
    
    % -----------------------------------------------  
    ha_mod = @(ha,zr) log((zr - 0.7 * ha) ./ (0.1 * ha));
    
    if strcmp(vegtype(i),'DB') || strcmp(vegtype(i),'EN')
        
        if zr < 1.5 * max(obs_hc)
            
            ha_mod = @(ha,zr) log((zr - 0.7 * ha) ./ (0.1 * ha)) + log(1.25);
            
        end
        
    end
    
    % -----------------------------------------------  
 
    x = fluxdat.datetime_start(t*12:t*24:height(fluxdat));
    
    export_table = table;
    export_table.date = x;
    export_table.ha = NaN(height(export_table),1);

    day_cnt = median(1:window_size);
    
    for j = 1:n_samp
        
        center_id = window_size * 24 * t * (j - 1) + 0.5 * window_size * 24 * t;
        start_id = center_id - 0.5 * window_size * 24 * t + 1;
        end_id = center_id + 0.5 * window_size * 24 * t;
                
        if end_id > n | center_id > n | start_id > n, break; end
        
        id = start_id:end_id;
        
        obs = ws_ustar_ratio_obs(id);    
                        
        if sum(isnan(obs)) > n_missing_allowed
            day_cnt = day_cnt + window_size;
            continue; 
        end      
        
        if length(zr) > 1
            
            zr_j = zr(center_id);
            
        else
            
            zr_j = zr;
            
        end        
        
        if isnan(zr_j)
            
            day_cnt = day_cnt + window_size;
            continue; 
            
        end
        
        objfun = @(x) nansum(abs(ha_mod(x,zr_j) - obs)) ./ sum(~isnan(obs));                
        ha_j = fminbnd(objfun,0.001,zr_j);
                                
        export_table.ha(day_cnt) = ha_j;
        
        day_cnt = day_cnt + window_size;

    end
    
    export_table.date.Format = 'yyyy-MM-dd';
    export_table = setNaN(export_table,-9999);
 
    cd([wdir,'/results/4_canopy_height']);
    writetable(export_table,sprintf('%s_canopy_height.csv',char(sites(i))));
    
end
% End of Script ----------------------------------------------------------------