function v = kinematic_molecular_viscosity(t_air,pressure,K_bool)

    if nargin < 3
        K_bool = true;
    end
    
    if K_bool == false
        t_air = t_air + 273.15;
    end
    
    pressure = 1000 * pressure;  % Convert from kPa to Pa
    
    v = (1.327e-05 * (101325 ./ pressure)) .* ((t_air ./ 273.15).^1.81);

end