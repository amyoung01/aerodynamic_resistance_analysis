function L = Monin_Obukhov_Length(t_air,ustar,H,pressure,K_bool)

if nargin < 5
    K_bool = true;
end

if K_bool == false
    t_air = t_air + 273.15;
end

rho = air_density(t_air,pressure);

% Constants
von_karman = 0.41;
cp = 1004.834;
g = 9.81;

L = -(cp * rho .* ustar.^3 .* t_air) ./ (von_karman * g * H);

end