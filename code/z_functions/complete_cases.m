function return_table = complete_cases(input_table,id_to_remove,cols_to_skip)

ncol = width(input_table);

if isempty(cols_to_skip)
    
   cols_to_skip = false(1,ncol);
    
end

if ~islogical(cols_to_skip)
    
   cols_to_skip_bool = false(1,ncol);
   cols_to_skip_bool(cols_to_skip) = true;
   cols_to_skip = cols_to_skip_bool;
    
end

leave_out_table = input_table(:,cols_to_skip);
df_arr = table2array(input_table(:,~cols_to_skip));

for i = 1:size(df_arr,2)
   
   df_arr(id_to_remove,i) = NaN;
    
end

return_table = array2table(df_arr, ...
    'VariableNames',input_table.Properties.VariableNames(~cols_to_skip));
return_table = [leave_out_table,return_table];

end