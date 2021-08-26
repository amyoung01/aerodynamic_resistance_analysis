function [smooth_vals,new_x,new_y,subs] = smooth_across_yr(x,y,window_size,span,x_extend)

%     y = peaks(365);
%     y = y(:,175);
%     y = repmat(y,4,1);
%     y = [y; y(1:100)];
%     y = y + normrnd(0,3,length(y),1);

    n = length(x);
    unique_x = unique(x);    
    subs = zeros(n,1);
    
    mv = median(1:window_size);
    el = window_size - mv;
    x_step = mv:window_size:length(unique_x);    
    
    k = 1;
    
    for j = x_step
        
        if j + el > length(unique_x)
            
            continue;
            
        end
        
        id = x >= unique_x(j-el) & x <= unique_x(j+el);
        subs(id) = k;
        
        k = k + 1;
        
    end
    
    y(subs == 0) = [];    
    subs(subs == 0) = [];
    
    y_med = accumarray(subs,y,[],@nanmedian);
    y_25 = accumarray(subs,y,[],@(x) prctile(x,25));
    y_75 = accumarray(subs,y,[],@(x) prctile(x,75));
    
    x_extended = [unique_x(1)-x_extend-1:unique_x(end)+x_extend];
    x_extended = x_extended(mv):window_size:max(x_extended);
    
    y_med = [y_med(end-x_extend/window_size:end); y_med; y_med(1:x_extend/window_size)];
    y_25  = [y_25(end-x_extend/window_size:end); y_25; y_25(1:x_extend/window_size)];
    y_75  = [y_75(end-x_extend/window_size:end); y_75; y_75(1:x_extend/window_size)];
    
    smooth_y_med = smooth(y_med,span,'loess');
    smooth_y_25 = smooth(y_25,span,'loess');
    smooth_y_75 = smooth(y_75,span,'loess');
    
    smooth_vals = [smooth_y_med, smooth_y_25, smooth_y_75];
    new_x = x_extended;
    new_y = y;
    
end