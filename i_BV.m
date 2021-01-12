%BV equation
%overpot = Electrode overpotential (V)
%i0 = Exchange current density (Use same units as desired for output current density i)
%a = charge transfer coefficient
%z = # electrons transferred
%temp = temperature (Kelvin)
function [i] = i_BV(overpot, i0, a, z, temp)
	R = 8.314; %J/mol/K
    F = 96485.3329; %C/mol e-
    i = i0*(exp((a*z*F)/(R*temp)*(overpot))-exp(-((1-a)*z*F)/(R*temp)*(overpot)));
end