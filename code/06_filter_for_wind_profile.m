% This script is written to go through and filter flux data and export data tables
% for each flux site. This script requires custom-written functions and a list of
% flux site codes to process.
%
% Date Created: 2019-01-04
% Date Modified:
%
% -----------------------------------------------------------------------------------

% INITIALIZE WORKSPACE
clear all;
clc;
close all;

% SET WORKING DIRECTORY
wdir = '/Volumes/GoogleDrive/My Drive/Young_aerodynamic_resistance_analysis';

% Add path to load in custom written functions to workspace
addpath([wdir,'/code/z_functions']);

cd([wdir,'/data/ancillary_data']);

% LOAD IN METADATA TABLE CONTAINING LIST OF SITES TO PROCESS
phenocam_flux_metadata_table = ...
    readtable('profile_metadata.csv','Delimiter',',');

sites = phenocam_flux_metadata_table.fluxsite;
start_date = floor(datenum(datetime('2000-01-01','InputFormat','yyyy-MM-dd')));
end_date   = floor(datenum(datetime('2018-12-31','InputFormat','yyyy-MM-dd')));

vars_names = readtable('variables_to_import_for_fluxsites.csv','Delimiter',',');

% INITIALIZE CRITERIA FOR FILTERING FLUX DATA
vars_of_interest = vars_names.Properties.VariableNames(2:end);

time_of_day_interval = [0,24];

for i = 1:length(sites)
    
    cd([wdir,'/results/1_summarized_precipitation']);
    precip_i = readtable(sprintf('%s_precip_days_to_remove.csv',char(sites(i))));
    precip_i = standardizeMissing(precip_i,-9999);
    
    cd([wdir,'/data/raw_data/ameriflux/BASE']);
    file_to_read = dir(sprintf('AMF_%s_*',char(sites(i))));
    fluxdat = readtable(file_to_read.name,'HeaderLines',2,'Delimiter',',');
    fluxdat = standardizeMissing(fluxdat,-9999);
    
    if strcmp(sites(i),'US-Syv')
        
        PA_id = ~isnan(fluxdat.PA);
        fluxdat.PA_1_1_1(PA_id) = fluxdat.PA(PA_id);
    
    end
    
    filename_parts = strsplit(file_to_read.name,'_');
    timestep = char(filename_parts(4));
    
    flux_time_scale = minutes(30);
    if strcmp(timestep,'HR'), flux_time_scale = minutes(60); end
    
    no_time_id = isnan(fluxdat.TIMESTAMP_START);
    fluxdat = fluxdat(~no_time_id,:);
            
    raw_data_var_names = fluxdat.Properties.VariableNames;
    
    % First, find the two date columns and convert to separate 'datetime' objects.
    % These will later be merged with the reduced data table for the analsysis.
    flux_datetime_start = datetime(num2str(fluxdat.TIMESTAMP_START), ...
        'InputFormat','yyyyMMddHHmm');
    flux_datetime_end = datetime(num2str(fluxdat.TIMESTAMP_END), ...
        'InputFormat','yyyyMMddHHmm');
    
    flux_datetime_start.Format = 'yyyy-MM-dd HH:mm';
    flux_datetime_end.Format = 'yyyy-MM-dd HH:mm';    
        
    hr = hour(flux_datetime_start + flux_time_scale/2);
    
    midday_bool = hr >= time_of_day_interval(1) & hr < time_of_day_interval(2);    
            
    day_values = floor(datenum(flux_datetime_start + flux_time_scale/2));
    
    day_value_bool = day_values >= start_date & day_values <= end_date;
    
    to_keep_id = midday_bool & day_value_bool;
    
    fluxdat = fluxdat(to_keep_id,:);
    flux_datetime_start = flux_datetime_start(to_keep_id);
    flux_datetime_end = flux_datetime_end(to_keep_id);    

    fluxsite_row = find(strcmp(vars_names.Site,char(sites(i))));
    
    fluxdat_filtered = table;
    
    for k = 1:length(vars_of_interest)
        
        if strcmp(char(vars_of_interest(k)),'precip')
            
            continue;
            
        end
        
        variable_column = ...
            find(strcmp(vars_names.Properties.VariableNames, ...
            vars_of_interest(k)));
        
        varname_to_import = vars_names(fluxsite_row,variable_column);
        varname_to_import = char(table2array(varname_to_import));
        
        if i == 5 & strcmp(vars_of_interest(k),'wind_speed')
            
            varname_to_import = 'WS_1_2_1';
            
        end
        
        if i == 5 & strcmp(vars_of_interest(k),'wind_speed_2')
            
            varname_to_import = 'WS_1_3_1';
            
        end
        
        if strcmp(varname_to_import,'NA')
            
             fluxdat_filtered.(char(vars_of_interest(k))) = NaN(height(fluxdat),1);
            
        else
            
             fluxdat_filtered.(char(vars_of_interest(k))) = fluxdat.(varname_to_import);
            
        end
        
    end
    
    fluxdat =  fluxdat_filtered; clear  fluxdat_filtered;
    
    sum_netrad_component_values = ...
        fluxdat.sw_in - fluxdat.sw_out + fluxdat.lw_in - fluxdat.lw_out;
    
    both_nan_idx = isnan(fluxdat.netrad) & (isnan(sum_netrad_component_values) == false);
    fluxdat.netrad(both_nan_idx) = sum_netrad_component_values(both_nan_idx);
    
    fluxdat = setNaN(fluxdat,-9999);
    
    fluxdat.datetime_start = flux_datetime_start;
    fluxdat.datetime_end = flux_datetime_end;
    
    fluxdat = fluxdat(:,[end-1,end,1:end-2]);
    
    export_filename = ...
        sprintf('%s//results//2_filtered_flux_data//24hr_windprofiles//%s_%s.csv',...
                wdir,char(sites(i)),timestep);
    
    % Export data tables
    writetable(fluxdat,export_filename);
    
end
% End of script ----------------------------------------------------------------