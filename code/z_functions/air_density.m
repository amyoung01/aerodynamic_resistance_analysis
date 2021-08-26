function rho = calc_air_density(T,PA,K_bool)

Rd = 287.0586; % Gas constant of dry air [J kg-1 K-1]

if nargin < 3
   K_bool = true; 
end

if K_bool == false
   T = T + 273.15; 
end

PA = PA*1000; % Convert kPa to Pa

rho = PA./(Rd*T);

end