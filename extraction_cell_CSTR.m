%{
	Continuous stirred tank model for base metal recovery cell
	Author: Matt/Ryan
	To Add: Applied Voltage calc. Non-constant volume?
%} 
clear all
global R F temp pressure volume Q;
%universal constants
R = 8.314; %J/(mol K)
F = 96485.3329; %C/mol
% system parameters
temp = 298; %K
pressure = 1; % atm
volume = 10; %L
Q = 0.01; % L/s (flowrate)
%{
	Reactions
	Anodic: Fe2+ --> Fe3+ + e-
	Cathodic: Cu2+ + 2e- --> Cu(s)
%}
% # of electrons in rxn
z_Cu = 2;
z_Fe = 1;
% exchange current densities
i0_Cu = 5E-5; %A/cm^2
i0_Fe = 5E-5; %A/cm^2
% charge transfer coefficients
alpha_Cu = 0.5;
alpha_Fe = 0.5;
% Standard half reaction potentials vs. SHE @ 298 K, 1 atm, 1 M https://en.wikipedia.org/wiki/Standard_electrode_potential_(data_page)
Eo_Cu = 0.337; %V
Eo_Fe = 0.77; %V
%activity coefficients of ions in solution
gamma_Cu2 = 1;
gamma_Fe2 = 1;
gamma_Fe3 = 1;
%initial concentrations in mol/L
%Outlet (In CSTR)
Ci_Cu2 = 0.5;
Ci_Fe2 = 0.1;
Ci_Fe3 = 0.001; 
%inlet concentrations
C_Cu2_o = 0.7;
C_Fe2_o = 0.4;
C_Fe3_o = 0.1;

%Electrode areas
S_cat = 100; %cm^2
S_an = 100; %cm^2

%Applied Voltage (potentiostat)
Vapp = 0.8; %V

%Mass balance on generic ion in reactor
dC_dt = @(C, C_o, I, z) ((C_o-C)*Q + I*F/z)/volume;

%Solver
tfinal = 1; %total time in seconds
h = 1; %Step size = 1 s

%initializing concentrations
C_Cu2 = zeros(1, tfinal/h+1);
C_Fe2 = zeros(1, tfinal/h+1);
C_Fe3 = zeros(1, tfinal/h+1);
C_Cu2(1) = Ci_Cu2;
C_Fe2(1) = Ci_Fe2;
C_Fe3(1) = Ci_Fe3;

%Nernst potential arrays
Erev_Cu = zeros(1, tfinal/h+1);
Erev_Fe = zeros(1, tfinal/h+1);
Erev_Cu(1) = Eo_Cu - R*temp/(z_Cu*F)*log(1/(gamma_Cu2*C_Cu2(1)));
Erev_Fe(1) = Eo_Fe - R*temp/(z_Fe*F)*log(gamma_Fe2*C_Fe2(1)/gamma_Fe3*C_Fe3(1));
r_sol = ones(1, tfinal/h)*10; %ohms, add calculation for this later
r_hardware = 0; %ohms
for k = 1:1:(tfinal/h)
    %Loop to solve for current
    j = 1;
    max_iter = 10000;
    I_an = zeros(1, max_iter);
    eta_an = ones(1, max_iter);
    eta_Cu = zeros(1, max_iter);
    I_Cu = zeros(1, max_iter);
    I_cat = zeros(1, max_iter);
    error = zeros(1, max_iter);
    pct_error = ones(1, max_iter);
    I_an(j) = 0.1*(Vapp - abs(Erev_Cu(k)-Erev_Fe(k)) - 0.1)/(r_sol(j)+r_hardware); %Guess for current, in A, assuming no overpotential
    while pct_error > 0.1
        %solve for anodic overpotential using B-V equation
        i_BV_an = @(overpot) i_BV(overpot, i0_Fe, alpha_Fe, z_Fe, temp) - I_an(j)/S_an;
        eta_an(j) = fzero(i_BV_an, 0.1);
        %solve for cathodic overpotential of copper reaction
        %eta_Cu(j) = Vapp - abs(Erev_Cu(k)-Erev_Fe(k)) - I_an(j)*(r_sol+r_hardware) - abs(eta_an(j));
        %eta_Cu(j) = eta_Cu(j) * (-1); %negative (cathodic overpotential)
        eta_Cu(j) = -Vapp + Erev_Fe(k) + eta_an(j) + I_an(j)*(r_sol+r_hardware) - Erev_Cu(k);
        I_Cu(j) = i_BV(eta_Cu(j), i0_Cu, alpha_Cu, z_Cu, temp)*S_cat;
        I_cat(j) = I_Cu(j);
        error(j) = 2*(I_an(j)-(-I_cat(j)))/((abs(I_cat(j))+abs(I_an(j)))); %I_cathode should be negative
        pct_error(j) = abs(error(j))*100;
        I_an(j+1) = I_an(j)-I_cat(j)/1000;
        j = j + 1;
        if j > max_iter
            break;
        end
    end
	%Runge Kutta 4 to solve for ionic concentrations
	k1 = dC_dt(C_Cu2(k), C_Cu2_o, I_Cu(j-1), z_Cu);
	k2 = dC_dt(C_Cu2(k)+h*k1/2, C_Cu2_o, I_Cu(j-1), z_Cu);
	k3 = dC_dt(C_Cu2(k)+h*k2/2, C_Cu2_o, I_Cu(j-1), z_Cu);
	k4 = dC_dt(C_Cu2(k)+h*k3, C_Cu2_o, I_Cu(j-1), z_Cu);
	C_Cu2(k+1) = C_Cu2(k) + 1/6*h*(k1+2*k2+2*k3+k4);
	%Nernst
	Erev_Cu(k+1) = Eo_Cu - R*temp/(z_Cu*F)*log(1/(gamma_Cu2*C_Cu2(k+1)));
	Erev_Fe(k+1) = Eo_Fe - R*temp/(z_Fe*F)*log(gamma_Fe2*C_Fe2(k+1)/gamma_Fe3*C_Fe3(k+1));
end
