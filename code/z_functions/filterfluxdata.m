function output_table = filterfluxdata(input_table,filters,rmByVar)

T = input_table;

for i = 1:length(filters)
    
    filters_i = char(filters{i}); 
    var2filter = strsplit(filters_i,' '); var2filter = var2filter{1};
    
    eval(['idx = find(T.',filters_i,');']);
    T.(var2filter)(idx) = -9999;
    
end

output_array = table2array(T);
output_array(output_array == -9999) = NaN;
                        
if nargin == 3 && ~isempty(rmByVar)
    
    if numel(rmByVar) == 1
        
        if strcmp(rmByVar,'all')
            
            output_array(sum(isnan(output_array),2) > 0,:) = NaN;
            
        else
            
            var2process_idx = strcmp(T.Properties.VariableNames,rmByVar);
            
            for k = 1:size(output_array,2)
                
                output_array(isnan(output_array(:,var2process_idx)),k) = NaN;
                
            end
            
        end
        
    elseif numel(rmByVar) > 1
        
        for j = 1:numel(rmByVar)
            
            var2process_idx = strcmp(T.Properties.VariableNames,char(rmByVar(j)));
                        
            for k = 1:size(output_array,2)
                
                output_array(isnan(output_array(:,var2process_idx)),k) = NaN;

                
            end
            
        end
        
    end
    
end

output_array(isnan(output_array)) = -9999;
output_table = array2table(output_array,'VariableNames',T.Properties.VariableNames);

end