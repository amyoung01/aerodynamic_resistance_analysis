function [fluxsummary_table,times_to_keep] = dailyFluxStats(input_table, ...
                                                          vars2exclude_names, ...
                                                          timeofday_interval, ...
                                                          function_handle, ...
                                                          n_NaN_allowed)
                                          
if nargin < 5 || isempty(n_NaN_allowed) == true
    
   n_NaN_allowed = 0; 
   
end

if nargin < 3 || isempty(timeofday_interval) == true
    
   timeofday_interval = [0,23]; 
   
end

T = input_table;
varnames = T.Properties.VariableNames;

if sum(strcmp(varnames,'datetime_start')) == 0
    
   error('No datetime_start variable in table.');
   
end

datetime_start = T.datetime_start;

if iscell(vars2exclude_names) == false
    
   vars2exclude_names = {vars2exclude_names}; 
   
end

if nargin < 2 || isempty(vars2exclude_names) == true
    
   vars2exclude_names = {'datetime_start'}; 
   
end

vars2exclude_idx = false(1,length(varnames));

for v = 1:length(vars2exclude_names)
    
    v_idx = strcmp(varnames,char(vars2exclude_names(v)));
    
    if sum(v_idx) == 0
        
        continue;
        
    end
    
    vars2exclude_idx(v_idx) = true;
    
end

T = T(:,logical(1-vars2exclude_idx));

if nargin < 4 || isempty(function_handle) == true
    
   function_handle = repmat({'mean'},size(T,2),1); 
   
end

if size(function_handle,1) ~= size(T,2) || ischar(function_handle)
    
    if ischar(function_handle) == true
        
        function_handle = {function_handle}; 
        
    end
    
    function_handle = repmat(function_handle,size(T,2),1);
    
    warning('Not the correct number of function handles.'); 
    
end

% Get datevec format for datetime and then subset data by time interval
datevec_start = datevec(datetime_start);

timeofday_idx = (datevec_start(:,4)>=timeofday_interval(1) & ...
                 datevec_start(:,4)<timeofday_interval(2));
             
% **************
times_to_keep = timeofday_idx;
                        
data_array_halfhourly = table2array(T(timeofday_idx,:));
nanidx_array_halfhourly = isnan(data_array_halfhourly);

% Get unique day identifier
individual_day_idx = datenum(datevec_start(timeofday_idx,1), ... Years
                             datevec_start(timeofday_idx,2), ... Months
                             datevec_start(timeofday_idx,3));  % Days

% min_dayidx = min(individual_day_idx); max_dayidx = max(individual_day_idx);
unique_dayidx = unique(individual_day_idx);

output_array = NaN(length(unique(individual_day_idx)),size(T,2),size(function_handle,2));

for i = 1:size(T,2)
    
    for f = 1:size(function_handle,2)
        
        function_handle_i = char(function_handle(i,f));
        function_handle_i = str2func(['nan',function_handle_i]);
        
        
        variable_daily_i = accumarray(individual_day_idx, ...
            data_array_halfhourly(:,i), ...
            [], ...
            function_handle_i);
        
        variable_daily_i = variable_daily_i(unique_dayidx);
        
        n_NaN_perday_array_i = accumarray(individual_day_idx, ...
            nanidx_array_halfhourly(:,i), ...
            [], ...
            @sum);
        
        n_NaN_perday_array_i = n_NaN_perday_array_i(unique_dayidx);
        
        bool_n_NaN_array = n_NaN_perday_array_i > n_NaN_allowed;
                
        variable_daily_i(bool_n_NaN_array) = NaN;
        output_array(:,i,f) = variable_daily_i;
        
    end
    
end

day_identifier_datevec = datevec(unique_dayidx);
date = datetime(day_identifier_datevec(:,1), ...
                day_identifier_datevec(:,2), ...
                day_identifier_datevec(:,3));
               
date.Format = 'yyyy-MM-dd';

if size(function_handle,2) == 1
    
    fluxsummary_table = array2table(output_array,'VariableNames',varnames(logical(1-vars2exclude_idx)));
    fluxsummary_table.date = date;
    fluxsummary_table = fluxsummary_table(:,[end,1:(end-1)]);
    
else
    
    fluxsummary_table = struct;
    fluxsummary_table.date = date;
    
    for f = 1:size(function_handle,2)
        
        table_f = array2table(output_array(:,:,f),'VariableNames',varnames(logical(1-vars2exclude_idx)));
        eval(['fluxsummary_table.table_',num2str(f),' = table_f;']);
        
    end
    
    fluxsummary_table.daily_stat = function_handle';
    
end

end