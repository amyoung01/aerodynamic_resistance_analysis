clear; clc; close all;

wdir = '/Volumes/GoogleDrive/My Drive/Young_aerodynamic_resistance_analysis';

% Add path to load in custom written functions to workspace
addpath([wdir,'/code/z_functions']);

cd([wdir,'/data/ancillary_data']);
phenoflux_metadata = readtable('pheno_flux_sites_to_use.csv', 'Delimiter',',');

sites = phenoflux_metadata.fluxsite;
phenos = phenoflux_metadata.phenosite;
primary_veg = phenoflux_metadata.vegtype;
primary_veg(strcmp(sites,'US-Ro4')) = {'GR'};
emissivity = phenoflux_metadata.emissivity;

cp = 1004.834; % Specific heat of air for constant pressure [J K^-1 kg^-1]
k  = 0.41; % Von-Karman Constant

timeinterval = [10,14];

for i = 1:length(sites)    
    
    window_size = 1;
    
    % Change window sites for alfalfa sites
    if strcmp(sites(i),'US-Bi1') || strcmp(sites(i),'US-Tw3')
        
        window_size = 1; 
    
    end
    
    if phenoflux_metadata.zr(i) == -9999, continue; end
        
    cd([wdir,'/results/2_filtered_flux_data/midday']);
    file_to_import = dir([char(sites(i)),'_*.csv']);
    file_name_parts = strsplit(file_to_import.name,'_');
    timestep = char(file_name_parts(2)); timestep = timestep(1:2);
    
    t = 2;
    if strcmp(timestep,'HR'), t = 1; end
        
    n_NaN = 0.5 * diff(timeinterval) * t;        
    
    fluxdat = readtable(file_to_import.name);
    fluxdat = standardizeMissing(fluxdat,-9999);
    fluxdat_var_names = fluxdat.Properties.VariableNames;
    
    cd([wdir,'/results/4_canopy_height/']);
    canopy_height = readtable(sprintf('%s_canopy_height.csv',char(sites(i))));
    canopy_height = standardizeMissing(canopy_height,-9999);
    ha = canopy_height.ha;
         
    lw_in = fluxdat.lw_in;           % [W m^-2]
    lw_out = fluxdat.lw_out;         % [W m^-2]
    netrad = fluxdat.netrad;
    H = fluxdat.H;                   % [W m^-2]
    LE = fluxdat.LE;
    t_air = fluxdat.t_air + 273.15;  % Air temperature, convert from C to K
    pressure = fluxdat.pressure;     % [kPa]
    wind_speed = fluxdat.wind_speed; % [m s^-1]
    ustar = fluxdat.ustar;           % [m s^-1]
    
    % Surf Temp [K]
    t_surf = radiometric_surface_temperature(lw_in, lw_out, emissivity(i));    
    delta_t = t_surf - t_air;
    bowen = H ./ LE;
    
    to_remove_id = delta_t .* H < 0 | bowen < 0 | bowen > 10;
            
    delta_t(to_remove_id) = NaN;
    t_air(to_remove_id) = NaN;
    t_surf(to_remove_id) = NaN;
    netrad(to_remove_id) = NaN;
    H(to_remove_id) = NaN;
    LE(to_remove_id) = NaN;
    ustar(to_remove_id) = NaN;
    pressure(to_remove_id) = NaN;
    
    rho = air_density(t_air,pressure);
    
    r_m = wind_speed ./ ustar.^2;
    
    dn_hh = floor(datenum(fluxdat.datetime_start + minutes(15)));
    dn_daily = floor(datenum(canopy_height.date));
    
    kB_inv_optim = NaN(length(dn_hh),1);
    ha_hh = NaN(length(dn_hh),1);
    
    el = window_size - median(1:window_size);
        
    for j = median(1:window_size):window_size:length(dn_daily)        
        
        if dn_daily(j) + el > max(dn_daily), continue; end
        
        id = dn_hh >= dn_daily(j) - el & dn_hh <= dn_daily(j) + el;       
        
        if sum(isnan(H(id) .* ...
                     r_m(id) .* ...
                     ustar(id) .* ...
                     delta_t(id) .* ...
                     rho(id))) > 0.75 * t * 4 * window_size
            
            continue;                 
        
        end               
                
        H_pred = @(kB_inv) (cp * rho(id) .* delta_t(id)) ./ ...
                           (r_m(id) + kB_inv ./ (k * ustar(id)));
                       
        objfun = @(x) nansum(abs(H_pred(x) - H(id))) ./ sum(~isnan(H(id)));
        param = fminbnd(objfun,-5,30);
                               
        kB_inv_optim(id) = param;
        ha_hh(id) = ha(j);
                
    end
                                        
    r_b = kB_inv_optim ./ (k * ustar);
    
    complete_cases_id = isnan(sum([H,r_m,r_b],2));
    r_m(complete_cases_id) = NaN;
    r_b(complete_cases_id) = NaN;
    
    r_h_pred = r_m + r_b;
    r_h_obs = (cp * rho .* delta_t) ./ H;
            
    % Create Input Table of half-hour values to create daily summaries using
    % dailyFluxStats.m. 
    hh_input_table = table;
    hh_input_table.datetime_start = fluxdat.datetime_start;
    hh_input_table.datetime_end = fluxdat.datetime_end;
    hh_input_table.netrad = netrad;
    hh_input_table.H = H;
    hh_input_table.LE = LE;
    hh_input_table.t_air = t_air;    
    hh_input_table.t_surf = t_surf;
    hh_input_table.lw_in = lw_in;
    hh_input_table.lw_out = lw_out;
    hh_input_table.pressure = pressure;
    hh_input_table.wind_speed = wind_speed;
    hh_input_table.ustar = ustar;
    hh_input_table.r_h_obs = r_h_obs;
    hh_input_table.r_h_pred = r_h_pred;
    hh_input_table.r_m = r_m;
    hh_input_table.r_b = r_b;
    hh_input_table.kB_inv = kB_inv_optim;    
    hh_input_table.z0m = 0.1 * ha_hh;
    hh_input_table.z0h = hh_input_table.z0m ./ exp(hh_input_table.kB_inv);
    
    export_table = dailyFluxStats(hh_input_table, ...
                                  {'datetime_start','datetime_end'}, ...
                                  timeinterval, ... % Time of day interval
                                  'mean', ... % Functions to calculate
                                  n_NaN); % number of NaN's allowed
                              
    dates = export_table.date;
    
    mask = false(height(export_table),1);
    mask(median(1:window_size):window_size:height(export_table)) = true;
    
    export_table_arr = table2array(export_table(:,2:end));
    export_table_arr(~mask,:) = NaN;
    
    export_table = array2table(export_table_arr, ...
                               'VariableNames',export_table.Properties.VariableNames(2:end));
    export_table.date = dates;
    
    export_table = export_table(:,[end,1:end-1]);

    % Set NaN values as -9999
    export_table = setNaN(export_table,-9999); 
    
    % Export table for each site.
    cd([wdir,'/results/5_resistance_values/daily']);
    writetable(export_table, sprintf('%s_resistance_values.csv', ...
                                     char(sites(i))));
        
end
% End of Script ----------------------------------------------------------------