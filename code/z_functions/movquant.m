function [yq,xq] = movquant(x,y,percentile,window,middle_value_bool,end_flag)

if nargin < 5
    middle_value_bool = false; 
end

if nargin < 6
    end_flag = 'constant'; 
end

if mod(window,2) == 0
    error('window size must be odd');
end

if isdatetime(x)
    x = floor(datenum(x));
    datetime_bool = true;
else
    datetime_bool = false;
end

unique_x = unique(x); 

x_start = unique_x(median([1,window]));
x_end   = max(unique_x) - floor(window/2); 
x_to_fill = x_start:x_end;

yq = NaN(length(unique_x(1):unique_x(end)),1);

k = x_start - unique_x(1) + 1;

for i = x_to_fill
    
    idx_i = x>=(i-floor(window/2)) & x<=(i+floor(window/2));
    
    if sum(idx_i) == 0
        yq(k) = NaN; 
    end
        
    yq(k) = prctile(y(idx_i),percentile);
    
    k = k + 1;
    
end

if window > 1
    
    idx_start = 1:floor(window/2);
    idx_end   = (length(unique_x)-floor(window/2)+1):length(unique_x);
    
    if strcmp(end_flag,'constant')
        
        idx = (x>=unique_x(idx_start(1)) & x<=unique_x(idx_start(end)));
        yq(idx_start) = prctile(y(idx),percentile);
        
        idx = (x>=unique_x(idx_end(1)) & x<=unique_x(idx_end(end)));
        yq(idx_end) = prctile(y(idx),percentile);
        
    elseif strcmp(end_flags,'nan')
        
        idx = ((x>=idx_start(1) & x<=idx_start(end)) | (x>=idx_end(1) & x<=idx_end(end)));
        yq(idx) = NaN;
        
    end
end

xq = unique_x(1):unique_x(end);

if datetime_bool
    dv = datevec(xq);
    days = datetime(dv(:,1),dv(:,2),dv(:,3));
    days.Format = 'yyyy-MM-dd';
    xq = days;
end

if window > 1 && middle_value_bool
    remove_bool = false(length(unique_x),1);
    remove_bool((x_start - unique_x(1) + 1):window:length(unique_x)) = true;
    remove_bool = logical(true(length(unique_x),1) - remove_bool);
    yq(remove_bool) = NaN;
end

end