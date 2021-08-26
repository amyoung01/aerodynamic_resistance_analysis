function [precip_table] = summarizeFluxPrecip(P,datetime_vector,nan_threshold,returnMonthly)

P(P<0) = NaN;

if nargin < 4
   returnMonthly = false; 
end

if nargin < 3
   nan_threshold = 1; 
end

if isempty(nan_threshold) == true
   nan_threshold = 1; 
end

datevec_for_P = datevec(datetime_vector);
dn = datenum(datevec_for_P(:,1),datevec_for_P(:,2),datevec_for_P(:,3));

precip_daily = accumarray(dn,P,[],@nansum); precip_daily = precip_daily(min(dn):max(dn));

nanidx = isnan(P);
sum_nanidx = accumarray(dn,nanidx,[],@sum); sum_nanidx = sum_nanidx(min(dn):max(dn));
sum_nanidx_toreport = sum_nanidx;
sum_nanidx_prop = sum_nanidx./(length(dn)/length(unique(dn)));

to_remove_idx = false(length(precip_daily),1);
to_remove_idx(sum_nanidx_prop>=nan_threshold) = true;
precip_daily(to_remove_idx) = NaN;

precip_vals = precip_daily;
daily_datevec  = datevec(min(dn):max(dn));
datetime_vals = datetime(daily_datevec(:,1),daily_datevec(:,2),daily_datevec(:,3));
datetime_vals.Format = 'yyyy-MM-dd';

if returnMonthly == true
    unique_years = unique(daily_datevec(:,1));
    precip_monthly = zeros(length(unique_years)*12,1);
    monthly_datevec = zeros(length(unique_years)*12,3);
    monthly_datevec(:,3) = 15;
    sum_nanidx_toreport = zeros(length(unique_years)*12,1);
    cnt = 1;
    for i = 1:length(unique_years)
        for j = 1:12
            month_idx = find(daily_datevec(:,1)==unique_years(i) & daily_datevec(:,2)==j);
            precip_monthly(cnt) = nansum(precip_daily(month_idx));
            sum_nanidx_toreport(cnt) = sum(isnan(precip_daily(month_idx)));
            monthly_datevec(cnt,1) = unique_years(i);
            monthly_datevec(cnt,2) = j;
            cnt = cnt + 1;
        end
    end
    precip_vals = precip_monthly;
    datetime_vals = datetime(monthly_datevec(:,1),monthly_datevec(:,2),monthly_datevec(:,3));
    datetime_vals.Format = 'yyyy-MM';
end

precip_table = table;
precip_table.date = datetime_vals;
precip_table.precip = precip_vals;
precip_table.n_missing = sum_nanidx_toreport;
precip_table.Properties.VariableUnits = {'day','mm','count'};

end