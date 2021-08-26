function Re = Reynolds_number(ustar,t_air,pressure,z0m,K_bool)

    if nargin < 5
       K_bool = true; 
    end

    v = kinematic_molecular_viscosity(t_air,pressure,K_bool);

    Re = z0m .* ustar ./ v;

end