% This script is written to go through and find days to exclude based on a
% precipitation filter.
%
% Date Created: 2019-01-04
% Date Modified for publication: 2020-11-11
% 
% -----------------------------------------------------------------------------------

% INITIALIZE WORKSPACE
clear all;
clc;
close all;

% Set working directory
wdir = '/Volumes/GoogleDrive/My Drive/Young_aerodynamic_resistance_analysis';

% Add path to load in custom written functions to workspace
addpath([wdir,'/code/z_functions']);

% Load in list of sites to process
phenocam_flux_metadata_table = readtable([wdir,'/data/ancillary_data/pheno_flux_sites_to_use.csv']);

sites = phenocam_flux_metadata_table.fluxsite;

% load in variable names for precip for each site
variables_to_import = readtable([wdir,'/data/ancillary_data/variables_to_import_for_fluxsites.csv'],'Delimiter',',');
P_labels = variables_to_import.precip;
     
variable_filters = {'x<0'};
               
for i = 1:length(sites)
    
    % Just interested in precip here
    variables_of_interest = {'datetime_start','datetime_end',char(P_labels(i))};
    
    % Read in raw ameriflux data file
    cd([wdir,'/data/raw_data/ameriflux/BASE']);
    file_to_read = dir(sprintf('AMF_%s_*',char(sites(i))));
    fluxdat_i = readtable(file_to_read.name,'HeaderLines',2);
    
    % Convert TIMESTAMP integers in raw ameriflux files to a datetime format
    flux_datetime_start = datetime(num2str(fluxdat_i.TIMESTAMP_START), ...
                                    'InputFormat','yyyyMMddHHmm');
    flux_datetime_end   = datetime(num2str(fluxdat_i.TIMESTAMP_END), ...
                                    'InputFormat','yyyyMMddHHmm');
                                
    flux_datetime_start.Format = 'yyyy-MM-dd HH:mm';
    flux_datetime_end.Format   = 'yyyy-MM-dd HH:mm';
    
    % Add converted datetime values to imported data table
    fluxdat_i.datetime_start = flux_datetime_start; clear flux_datetime_start;
    fluxdat_i.datetime_end   = flux_datetime_end; clear flux_datetime_end;
    
    % Find the variables of interest in the flux data file
    varNames = fluxdat_i.Properties.VariableNames;
    [cols_to_keep,new_var_order,vars_available] = ...
        findFluxVars(varNames,variables_of_interest);
    fluxdat_i = fluxdat_i(:,cols_to_keep); fluxdat_i = fluxdat_i(:,new_var_order);
    fluxdat_i.Properties.VariableNames = variables_of_interest(vars_available);
    fluxdat_i.Properties.VariableNames = {variables_of_interest{1:2},'precip'};   
    
    % Hour values
    hr_vals = hour(fluxdat_i.datetime_start);
    
    % Need to do this in two parts to find days where it rained from 8:00 pm
    % night before to 2:00 pm day of. Use custom function 'summarizeFluxPrecip'
    % to do this. This function can be found in z_functions folder.
    fluxdat_1 = fluxdat_i;
    fluxdat_2 = fluxdat_i;
    
    to_remove_idx = hr_vals >= 14;
    fluxdat_1 = fluxdat_1(to_remove_idx,:);
    daily_precip_1 = summarizeFluxPrecip(fluxdat_1.precip, ...
                                         fluxdat_1.datetime_start, ...
                                         0.75, ...
                                         false);    
    
    to_remove_idx = hr_vals < 20;
    fluxdat_2 = fluxdat_2(to_remove_idx,:);    

    daily_precip_2 = summarizeFluxPrecip(fluxdat_2.precip, ...
                                         fluxdat_2.datetime_start, ...
                                         0.75, ...
                                         false);

    date = daily_precip_1.date;
    
    p1 = daily_precip_1.precip;
    p2 = daily_precip_2.precip;
    
    P = zeros(length(p1),2);
    P(2:end,1) = p2(2:end);
    P(1:end-1,2) = p1(1:end-1);
    precip_bool = sum(P,2) > 0;
    
    precip_days = table;
    precip_days.date = date;
    precip_days.precip_bool = precip_bool;
                                         
    precip_days = setNaN(precip_days,-9999);
    
    % Export a table that gives a boolean value of whether the day should be
    % excluded due to precip.
    writetable(precip_days, ...
        sprintf('%s//results//1_summarized_precipitation//%s_precip_days_to_remove.csv',wdir,char(sites(i)))); 
        
end