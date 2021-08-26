% Same as script '04_filter_for_midday.m' but does not include any filtering.
% Returns all half-hour (or hour) values. Used for calculating aerodynamic canopy height
%
% Date Created: 2019-01-04
% Date Modified for publication: 2020-11-11
%
% -----------------------------------------------------------------------------------
clear all;
clc;
close all;

wdir = '/Volumes/GoogleDrive/My Drive/Young_aerodynamic_resistance_analysis';
addpath([wdir,'/code/z_functions']);

cd([wdir,'/data/ancillary_data']);
phenocam_flux_metadata_table = ...
    readtable('pheno_flux_sites_to_use.csv','Delimiter',',');

sites = phenocam_flux_metadata_table.fluxsite;
start_date = floor(datenum(phenocam_flux_metadata_table.start_date));
end_date = floor(datenum(phenocam_flux_metadata_table.end_date));

vars_names = readtable('variables_to_import_for_fluxsites.csv','Delimiter',',');
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
    
    filename_parts = strsplit(file_to_read.name,'_');
    timestep = char(filename_parts(4));
    
    flux_time_scale = minutes(30);
    if strcmp(timestep,'HR'), flux_time_scale = minutes(60); end
    
    no_time_id = isnan(fluxdat.TIMESTAMP_START);
    fluxdat = fluxdat(~no_time_id,:);
            
    raw_data_var_names = fluxdat.Properties.VariableNames;
    
    flux_datetime_start = datetime(num2str(fluxdat.TIMESTAMP_START), ...
        'InputFormat','yyyyMMddHHmm');
    flux_datetime_end = datetime(num2str(fluxdat.TIMESTAMP_END), ...
        'InputFormat','yyyyMMddHHmm');
    
    flux_datetime_start.Format = 'yyyy-MM-dd HH:mm';
    flux_datetime_end.Format = 'yyyy-MM-dd HH:mm';    
        
    hr = hour(flux_datetime_start + flux_time_scale/2);
    
    midday_bool = hr >= time_of_day_interval(1) & hr < time_of_day_interval(2);    
            
    day_values = floor(datenum(flux_datetime_start + flux_time_scale/2));
    
    day_value_bool = day_values >= start_date(i) & day_values <= end_date(i);
    
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
        sprintf('%s//results//2_filtered_flux_data//24hr//%s_%s.csv',...
                wdir,char(sites(i)),timestep);
    
    % Export data tables
    writetable(fluxdat,export_filename);
    
end
% End of script ----------------------------------------------------------------
