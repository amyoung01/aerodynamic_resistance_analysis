function output_table = setNaN(input_table,nan_value)
    
    if nargin < 2
        nan_value = -9999;
    end
    
    varnames = input_table.Properties.VariableNames;
    
    remove_datetime_idx = true(1,length(varnames));
        
    for i = 1:width(input_table)

        if ~isnumeric(input_table.(char(varnames(i))))
            remove_datetime_idx(i) = false;
        end
        
    end
    
    A = table2array(input_table(:,remove_datetime_idx));
    A(isnan(A)) = nan_value;
    
    T = array2table(A,'VariableNames',varnames(remove_datetime_idx));
        
    output_table = [input_table(:,~remove_datetime_idx),T];
    
    [~,I] = sort([find(~remove_datetime_idx),find(remove_datetime_idx)]);
    
    output_table = output_table(:,I);

end