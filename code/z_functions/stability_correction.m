function [psi_h,psi_m] = stability_correction(zeta)

psi_h = zeros(numel(zeta),1); 
psi_m = zeros(numel(zeta),1);
              
x = (1 - 16 * zeta).^0.25;

stable = (zeta >= 0 | isnan(zeta));
psi_h(stable) = -5 * zeta(stable);
psi_m(stable) = -5 * zeta(stable);

unstable = (zeta < 0 | isnan(zeta));
psi_h(unstable) = 2 * log((1 + x(unstable).^2)./2);
psi_m(unstable) = 2 * log((1 + x(unstable))./2) + log((1 + x(unstable).^2)./2) - 2 * atan(x(unstable)) + pi/2;

end
