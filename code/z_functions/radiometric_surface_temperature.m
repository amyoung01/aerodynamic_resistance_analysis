function TS = radiometric_surface_temperature(LW_IN,LW_OUT,emissivity,K_bool)

    stefanboltzmann = 5.67e-8;

    if nargin < 3
        
       emissivity = 0.97;
       
    end    
    
    if nargin < 4
        
        K_bool = true;
        
    end
    
    TS = nthroot((LW_OUT - (1 - emissivity) * LW_IN) ./ (emissivity * stefanboltzmann),4);
    
    if K_bool == false
        
        TS = TS - 273.15;
        
    end
    
end