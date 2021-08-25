clear all; 
clc; 
close all;

wdir = '/Volumes/GoogleDrive/My Drive/Young_aerodynamic_resistance_analysis';

% Add path to load in custom written functions to workspace
addpath([wdir,'/code/z_functions']);

cd([wdir,'/data/ancillary_data']);
phenoflux_metadata = readtable('profile_metadata.csv','Delimiter',',');
vegtype = phenoflux_metadata.vegtype;

zr = [phenoflux_metadata.zr_1, ...
      phenoflux_metadata.zr_2];

sites = phenoflux_metadata.fluxsite;

var_names = readtable('variables_to_import_for_fluxsites.csv', ...
    'Delimiter',',');

sites_to_keep = ismember(var_names.Site,sites);
var_names = var_names(sites_to_keep,:);

primary_wd = phenoflux_metadata.primary_wd;
wd_range = [primary_wd - 90,primary_wd + 90];
wd_range(wd_range < 0) = wd_range(wd_range < 0) + 360;
wd_range(wd_range > 360) = wd_range(wd_range > 360) - 360;

window_size = 3;
prop_missing_allowed = 0.95;

k = 0.41;

for i = 1:length(sites)

    % -----------------------------------------------  
    zr_i = zr(i,:);        
    obs_hc = phenoflux_metadata.hc(i);
        
    % -----------------------------------------------
    cd([wdir,'/results/2_filtered_flux_data/24hr_windprofiles']);
    
    file_to_import = dir([char(sites(i)),'_*.csv']);
    filename = file_to_import.name;
    filename_parts = strsplit(filename,'_');         
    
    fluxdat = readtable(file_to_import.name);
    fluxdat = standardizeMissing(fluxdat,-9999);
    
    if i == 4
        
        to_keep = fluxdat.datetime_start > '2006-01-01';
        fluxdat = fluxdat(to_keep,:);
        
    end
                
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
        
    ws_1 = fluxdat.wind_speed; % [m s^-1]
    ws_2 = fluxdat.wind_speed_2; % [m s^-1]
    
    wind_direction = fluxdat.wind_direction; % [degrees]
    ustar = fluxdat.ustar; % [m s^-1]
    
    % Find and filter by neutral log-wind profile conditions ------------------------    
    L = Monin_Obukhov_Length(t_air,ustar,H,pressure,false); % [m]    
    zeta = zr_i(1) ./ L;        
    zeta_neutral_bool = abs(zeta) < 0.10;
    
    % -----------------------------------------------    
    max_ustar = 0.6;
    if strcmp(vegtype(i),'DB'), max_ustar = 1.0; end 
    ustar_bool = ustar >= 0.2 & ustar <= max_ustar;
    
    % -----------------------------------------------  
    if wd_range(i,1) < wd_range(i,2)
        
        wd_bool = wind_direction > wd_range(i,1) & wind_direction < wd_range(i,2);
        
    else
        
        wd_bool = wind_direction > wd_range(i,1) | wind_direction < wd_range(i,2);
        
    end
    
    ws_height_bool = ws_1 > ws_2;
    
    % -----------------------------------------------  
    
    to_keep_id = zeta_neutral_bool & ws_height_bool & ustar_bool & wd_bool;
    
    d_mdl = @(d,zr) log((zr(1) - d) ./ (zr(2) - d));
    
    % -----------------------------------------------  
    ws_ustar_ratio_obs = k * (ws_1 - ws_2) ./ ustar;
    ws_ustar_ratio_obs(~to_keep_id) = NaN;    

    ws_ustar_ratio_obs_z0 = exp(k * ws_1 ./ ustar);
    ws_ustar_ratio_obs_z0(~to_keep_id) = NaN;
    
    % -----------------------------------------------       
    x = fluxdat.datetime_start(t*12:t*24:height(fluxdat));
    
    export_table = table;
    export_table.date = x;
    export_table.d = NaN(height(export_table),1);
    export_table.z0m = NaN(height(export_table),1);
    
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
                    
        objfun = @(x) nansum(abs(d_mdl(x,zr_i) - obs)) ./ sum(~isnan(obs));
        
        d_j = fminbnd(objfun,0.001,zr_i(2));
        
        export_table.d(day_cnt) = d_j;
        
        export_table.z0m(day_cnt) = (zr_i(1) - export_table.d(day_cnt)) ./ nanmedian(ws_ustar_ratio_obs_z0(id));
                
        day_cnt = day_cnt + window_size;

    end
    
    export_table.date.Format = 'yyyy-MM-dd';
    
    export_table = setNaN(export_table,-9999);
 
    cd([wdir,'/results/4_canopy_height']);
    writetable(export_table,sprintf('%s_wind_profile_z0m_d.csv',char(sites(i))));
        
end