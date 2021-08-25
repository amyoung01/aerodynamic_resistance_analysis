% This script is written to go through and filter flux data and export data tables
% for each flux site. This script requires custom-written functions and a list of
% flux site codes to process.
%
% Date Created: 2019-01-04
% Date Modified for publication: 2020-11-11
%
% -----------------------------------------------------------------------------------

% Initialize workspace
clear all;
clc;
close all;

% Set working directory
wdir = '/Volumes/GoogleDrive/My Drive/Young_aerodynamic_resistance_analysis';

% Add path to load in custom written functions to workspace
addpath([wdir,'/code/z_functions']);

cd([wdir,'/data/ancillary_data']);

% Load in list of sites to process
phenocam_flux_metadata_table = ...
    readtable('pheno_flux_sites_to_use.csv','Delimiter',',');

sites = phenocam_flux_metadata_table.fluxsite;
start_date = floor(datenum(phenocam_flux_metadata_table.start_date));
end_date = floor(datenum(phenocam_flux_metadata_table.end_date));

% Variable names for each site to include in exported data table
vars_names = readtable('variables_to_import_for_fluxsites.csv','Delimiter',',');

% Names of variables we are interested in keeping from each flux data file
vars_of_interest = vars_names.Properties.VariableNames(2:end);

% Filtering for only midday values
time_of_day_interval = [10,14];

% Variable filters at the half-hour (or hour) timescale
variable_filters = {{'netrad < 50'}, ...
                    {'ustar < 0.2'}, ...
                    {'H < 50'}};

for i = 1:length(sites)
    
    % Load in list of days to filter out based on precipitation
    cd([wdir,'/results/1_summarized_precipitation']);
    precip_i = readtable(sprintf('%s_precip_days_to_remove.csv',char(sites(i))));
    precip_i = standardizeMissing(precip_i,-9999);
    
    % Load in raw BASE ameriflux data file
    cd([wdir,'/data/raw_data/ameriflux/BASE']);
    file_to_read = dir(sprintf('AMF_%s_*',char(sites(i))));
    fluxdat = readtable(file_to_read.name,'HeaderLines',2,'Delimiter',',');
    fluxdat = standardizeMissing(fluxdat,-9999);
    
    filename_parts = strsplit(file_to_read.name,'_');
    timestep = char(filename_parts(4));
    
    flux_time_scale = minutes(30);
    if strcmp(timestep,'HR'), flux_time_scale = minutes(60); end
    
    % only include half hour values where time is recorded. This affected maybe
    % one site in total.
    no_time_id = isnan(fluxdat.TIMESTAMP_START);
    fluxdat = fluxdat(~no_time_id,:);
    
    % Get variable names from the downloaded ameriflux BASE file
    raw_data_var_names = fluxdat.Properties.VariableNames;
    
    % Convert TIMESTAMP integers to datetime values
    flux_datetime_start = datetime(num2str(fluxdat.TIMESTAMP_START), ...
        'InputFormat','yyyyMMddHHmm');
    flux_datetime_end = datetime(num2str(fluxdat.TIMESTAMP_END), ...
        'InputFormat','yyyyMMddHHmm');
    
    flux_datetime_start.Format = 'yyyy-MM-dd HH:mm';
    flux_datetime_end.Format = 'yyyy-MM-dd HH:mm';    
    
    % Hour values for filtering for midday
    hr = hour(flux_datetime_start + flux_time_scale/2);    
    midday_bool = hr >= time_of_day_interval(1) & hr < time_of_day_interval(2);    
    
    % Get unique integers for each day
    day_values = floor(datenum(flux_datetime_start + flux_time_scale/2));
    
    % Keep only days during specified time span in metadata file. This was
    % determined by looking at time series of different variables for each site
    % and looking for completeness and jumps.
    day_value_bool = day_values >= start_date(i) & day_values <= end_date(i);
    
    % only keep these days and midday vallues
    to_keep_id = midday_bool & day_value_bool;
    
    fluxdat = fluxdat(to_keep_id,:);
    flux_datetime_start = flux_datetime_start(to_keep_id);
    flux_datetime_end = flux_datetime_end(to_keep_id);    
    
    % Get ready to extract only the variables we are interested in. This uses
    % custom written functions found in the z_functions folder
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
        
        if strcmp(varname_to_import,'NA')
            
             fluxdat_filtered.(char(vars_of_interest(k))) = NaN(height(fluxdat),1);
            
        else
            
             fluxdat_filtered.(char(vars_of_interest(k))) = fluxdat.(varname_to_import);
            
        end
        
    end
    
    fluxdat =  fluxdat_filtered; clear  fluxdat_filtered;
    
    sum_netrad_component_values = ...
        fluxdat.sw_in - fluxdat.sw_out + fluxdat.lw_in - fluxdat.lw_out;
    
    both_nan_idx = isnan(fluxdat.netrad) & ~isnan(sum_netrad_component_values);
    fluxdat.netrad(both_nan_idx) = sum_netrad_component_values(both_nan_idx);
    
    fluxdat = filterfluxdata(fluxdat,variable_filters,{'netrad','ustar','H'});

    datenum_to_remove = datenum(precip_i.date(logical(precip_i.precip_bool)));
    datenum_hh_values = floor(datenum(flux_datetime_start + minutes(15)));
    
    Lia = ismember(datenum_hh_values,datenum_to_remove);
    
    blank_table = array2table(NaN(sum(Lia),width(fluxdat)));
    fluxdat(Lia,:) = blank_table;
    
    fluxdat = setNaN(fluxdat,-9999);
    
    fluxdat.datetime_start = flux_datetime_start;
    fluxdat.datetime_end = flux_datetime_end;
    
    fluxdat = fluxdat(:,[end-1,end,1:end-2]);
    
    export_filename = ...
        sprintf('%s//results//2_filtered_flux_data//midday//%s_%s.csv',...
                wdir,char(sites(i)),timestep);
    
    % Export data tables
    writetable(fluxdat,export_filename);
    
end
% End of script ----------------------------------------------------------------