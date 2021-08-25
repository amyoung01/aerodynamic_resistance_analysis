% This is to reformat the raw data files for the mead sites to include how the
% measurement height (zr) changes depending on the height of vegetation.

clear;
clc;
close all;

fluxsites = {'US-Ne1','US-Ne2','US-Ne3'};

var_names = {'SW_IN_1_1_1','SW_OUT_1_1_1','LW_IN_1_1_1','LW_OUT_1_1_1', ...
             'NETRAD_1_1_1','H_1_1_1','LE_1_1_1','G_PI_F_1_1_1','TA_1_1_1','P_PI_F_1_1_1','PA_PI_F_1_1_1', ...
             'WS_1_1_1','WS_2_1_1','WS_1_2_1','WS_1_3_1','WD_1_1_1','USTAR_1_1_1','ZL_PI_F_1_1_1','RH_1_1_1'};
         

for i = 1:3
    
    cd('/Volumes/GoogleDrive/My Drive/Young_aerodynamic_resistance_analysis/data/raw_data/ameriflux/mead_sites_raw_data');

    fluxdat = readtable(sprintf('AMF_%s_BASE_HR_9-5.csv',char(fluxsites(i))));
      
    flux_datetime_start = fluxdat.TIMESTAMP_START;
    flux_datetime_end = fluxdat.TIMESTAMP_END;
    
    new_table_arr = -9999 * ones(height(fluxdat),length(var_names));
    
    for j = 1:length(var_names)
        
        id = find(strcmp(fluxdat.Properties.VariableNames,var_names(j)));
        
        if ~isempty(id)
            new_table_arr(:,j) = table2array(fluxdat(:,id));
        else
            new_table_arr(:,j) = NaN(height(fluxdat),1);            
        end
        
    end
    
    zr = -9999 * ones(height(fluxdat),1);
    
    WS = -9999 * ones(height(fluxdat),1);
    WS_1_col = strcmp(var_names,'WS_1_1_1');
    WS_2_col = strcmp(var_names,'WS_2_1_1');
    
    id_1 = new_table_arr(:,WS_1_col) > -9999 & new_table_arr(:,WS_2_col) == -9999;
    id_2 = new_table_arr(:,WS_2_col) > -9999 & new_table_arr(:,WS_1_col) == -9999;
    
    WS(id_1) = new_table_arr(id_1,WS_1_col);
    WS(id_2) = new_table_arr(id_2,WS_2_col);
    
    zr(id_1) = 6;
    zr(id_2) = 3;

    new_table_arr = new_table_arr(:,~WS_2_col);
    new_table_arr(:,WS_1_col) = WS;
    
    export_table = array2table(new_table_arr,'VariableNames',var_names(~WS_2_col));
    export_table.zr = zr;
    export_table.TIMESTAMP_START = flux_datetime_start;
    export_table.TIMESTAMP_END = flux_datetime_end;
    
    export_table = export_table(:,[end-1,end,1:end-2]);    
    
    cd /Volumes/GoogleDrive/'My Drive'/Young_phenology_ET_analysis/data/raw_data/ameriflux/BASE;
    writetable(export_table,sprintf('AMF_%s_BASE_HR_9-5.csv', ...
                                    char(fluxsites(i))));
    
    
    
end