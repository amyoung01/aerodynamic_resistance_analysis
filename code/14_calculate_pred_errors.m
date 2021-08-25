clear; clc; close all;

wdir = '/Volumes/GoogleDrive/My Drive/Young_aerodynamic_resistance_analysis';

% Add path to load in custom written functions to workspace
addpath([wdir,'/code/z_functions']);

cd([wdir,'/data/ancillary_data']);
phenoflux_metadata = readtable('pheno_flux_sites_to_use.csv','Delimiter',',');
phenoflux_metadata = phenoflux_metadata(~strcmp(phenoflux_metadata.fluxsite,'US-Ne3'),:);

sites = phenoflux_metadata.fluxsite;
phenos = phenoflux_metadata.phenosite;
primary_veg = phenoflux_metadata.vegtype; 
primary_veg(strcmp(sites,'US-Ro4')) = {'GR'};

sites(strcmp(sites,'US-Ne1')) = {'US-Ne-corn'};
sites(strcmp(sites,'US-Ne2')) = {'US-Ne-soybean'};

% cd([wdir,'/results/6_prediction_errors']);
% optimal_spans = readtable('z_kB_inv_span_values.csv','Delimiter',',','TreatAsEmpty','NA');
% span = optimal_spans.kB_inv_span;

cp = 1004.834; % Specific heat of air for constant pressure [J K^-1 kg^-1]
k  = 0.41; % Von-Karman Constant

vegtype = {'DB','EN','GR','SH','AG'};    
kB_inv_obs = [0.5,1,3.25,4.25,1.88];

for i = 1:length(sites)
    
    if strcmp(sites(i),'US-Bi1') || strcmp(sites(i),'US-Tw3'), continue; end
    
    cd([wdir,'/results/5_resistance_values']);
    file_to_import = dir([char(sites(i)),'_*.csv']);
    file_name_parts = strsplit(file_to_import.name,'_');
    timestep = char(file_name_parts(2)); timestep = timestep(1:2);
    
    t = 2;
    if strcmp(timestep,'HR'), t = 1; end
    
    fluxdat = readtable(file_to_import.name);
    fluxdat = standardizeMissing(fluxdat,-9999);
    
    t_air = fluxdat.t_air;
    t_surf = fluxdat.t_surf;
    ustar = fluxdat.ustar;
    pressure = fluxdat.pressure;
    H_obs = fluxdat.H;
    
    r_m = fluxdat.r_m;
    
    kB_inv_zero = zeros(height(fluxdat),1);
    kB_inv_constant = kB_inv_obs(strcmp(vegtype,primary_veg(i))) * ones(height(fluxdat),1);
    kB_inv_smooth = fluxdat.smooth_kB_inv;
    
%     to_keep_kB_inv = false(length(kB_inv_optim),1);
%     to_keep_kB_inv(3:5:end) = true;
%     kB_inv_optim(~to_keep_kB_inv) = NaN;
    
%     x = datenum(fluxdat.date);
%     y = kB_inv_optim;
    
%     kB_inv_smooth = fluxdat.smooth_kB_inv;
            
%     kB_inv_smooth(to_keep_kB_inv) = smooth(x,y,span(i),'loess');
    
    rho = air_density(t_air,pressure,true);
    delta_t = t_surf - t_air;
    
    r_h_zero = r_m + kB_inv_zero ./ (k * ustar);
    r_h_contstant = r_m + kB_inv_constant ./ (k * ustar);
    r_h_optim = r_m + kB_inv_smooth ./ (k * ustar);
    
    H_pred_zero = (cp * rho .* delta_t) ./ r_h_zero;
    H_pred_constant = (cp * rho .* delta_t) ./ r_h_contstant;
    H_pred_optim = (cp * rho .* delta_t) ./ r_h_optim;
   
    E1 = (H_pred_zero ./ H_obs);
    E2 = (H_pred_constant ./ H_obs);
    E3 = (H_pred_optim ./ H_obs);
    
    export_table = table;
    export_table.date = fluxdat.date;
    export_table.rel_doy = fluxdat.rel_doy;
    export_table.H_obs = H_obs;    
    export_table.error_1 = E1;
    export_table.error_2 = E2;
    export_table.error_3 = E3;
    
    export_table = setNaN(export_table,-9999);
    
    cd([wdir,'/results/7_prediction_errors']);
    writetable(export_table,sprintf('%s_H_pred_errors.csv',char(sites(i))));
    
end