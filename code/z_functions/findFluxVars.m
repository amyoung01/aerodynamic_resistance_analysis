function [column_linidx,new_order_of_vars,vars_available] = ...
                                              findFluxVars(varNames,variables_of_interest)

flux_variable_names = varNames;
n_flux_variables    = length(flux_variable_names);

% Find columns of fluxdat_i table that are needed for analysis. Use
% variables_of_interest to find these columns.
columns_tokeep_idx      = false(1,n_flux_variables);
reordering_columns_idx  = zeros(1,n_flux_variables);
vars_available          = false(1,length(variables_of_interest));

for j = 1:length(variables_of_interest)
    
    var_idx_1 = strcmp(flux_variable_names,char(variables_of_interest(j)));
    var_idx_2 = strcmp(flux_variable_names,[char(variables_of_interest(j)),'_1_1_1']);
    var_idx_3 = strcmp(flux_variable_names,[char(variables_of_interest(j)),'_PI']);
    var_idx_4 = strcmp(flux_variable_names,[char(variables_of_interest(j)),'_PI_1_1_1']);
    
    variable_name_mtx = [var_idx_1',var_idx_2',var_idx_3',var_idx_4'];
    
    if sum(variable_name_mtx(:)) == 1
        variable_name_idx = variable_name_mtx(:,find(sum(variable_name_mtx,1)));
    elseif sum(variable_name_mtx(:)) > 1
        variable_name_idx = variable_name_mtx(:,min(find(sum(variable_name_mtx,1))));
    elseif sum(variable_name_mtx(:)) == 0
        variable_name_idx = 0;
    end
    
    if variable_name_idx == 0
        warning(sprintf('%s is not available in data table.', ...
            char(variables_of_interest(j))));
        continue;
    end
    
    vars_available(j) = true;
    columns_tokeep_idx(variable_name_idx==true) = true;
    reordering_columns_idx(variable_name_idx==1) = j;
    
end
reordering_columns_idx(reordering_columns_idx==0) = [];
[~,new_order_of_vars] = sort(reordering_columns_idx);
column_linidx = columns_tokeep_idx;

end