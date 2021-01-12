clear all
tic
%inlet concs for now
nAg0(1) = 0.10; %mol
nAu0(1) = 0.05; %mol
nPd0(1) = 0.001; %mol
nH0(1) = 1e-1; %pH of 1
%conductivity values, check on this (S m2 /mol)
gammaAg = 61.90e-4;
gammaAu = 4.1e7;
gammaPd = 0.1;
gammaH = 349.81e-4;
%divergence from ideality values (check on this shit too but low
%priority)... something about fugacity comes to mind???
gamAg = 1;%assume all ideal unless proven otherwise
gamAu = 1;
gamPd = 1;
gamH = 1;
%setting up initial conc in cell
nAg(1) = nAg0(1); % moles of silver remaining
nAu(1) = nAu0(1);
nPd(1) = nPd0(1);
nH(1) = nH0(1);
aO2 = 0.21;%atm, its always this cause atmosphere,but maybe we wanna pressurize

%Universal constants (keep em handy)
F = 96500; %C/mol
Rgas = 8.3145; %J/molK
mmAg = 196.96657; %g/mol
mmAu = 107.8682; %g/mol
mmPd = 106.42; %g/mol

%System variables
Fin = 50; %L/s
Fout = Fin; %kept same for now, may want to make into different values to account for accumulation and controls
T = 298.15; %K, we can play around with this but if we want to vary this kinetically then I oop
V = 1000; %L
cursivel = 0.1; %m, characteristic distance
Acat = 1000; %m2, area
Aan = 100; %m2, area
Vapp = 8; %V
IRhardware = 0; %set to 0 for now before i do something with it

%calculated vars (will be iterated)
K(1) = gammaAg*nAg(1)+gammaAu*nAu(1)+gammaPd*2*nPd(1)+gammaH*nH(1);
R(1) = cursivel/(Acat*K(1));
Ean(1) = -1.23 - (Rgas*T/(4*F))*log(aO2*(gamH*nH(1)/V)^4);
ErevAg(1) = 0.7989 + (Rgas*T/(F))*log(1/(gamAg*nAg(1)/V));
ErevAu(1) = 1.692 + (Rgas*T/(F))*log(1/(gamAg*nAu(1)/V));
ErevPd(1) = 0.951 + (Rgas*T/(2*F))*log(1/(((gamAg*nPd(1))/V)^2));
I(1) = abs((Vapp - abs(ErevAg(1)-Ean(1)) - IRhardware - 0.6)/R(1)); %this is the initial initial guess for I

%assume well mixed, or bulk and surface concentrations are the same
alphaAn = 0.5;
alphaAg = 0.5; %base assumption as most systems are close to this number
alphaAu = 0.5; %base assumption as most systems are close to this number
alphaPd = 0.5; %base assumption as most systems are close to this number
iAg0 = 7.3e-9; %A/m2, from https://link-springer-com.proxy.lib.uwaterloo.ca/article/10.1007/s10800-007-9434-x - verify this link applies
iAu0 = 2E-8; %A/m2, value of dissolution including Na2S in dissolution reaction from https://link-springer-com.proxy.lib.uwaterloo.ca/article/10.1134/S1023193506040021, converting A/cm2 to A/m2 
iPd0 = 4e-9; %using a mix of a few from that reference
iAn0 = 1E-8; %A/cm2

%electrochemical equations and relations
dnAgdt = @(nAg0,nAg,t,iAg) nAg0*Fin/V-nAg*Fout/V-(1/F)*iAg*Acat;
dnAudt = @(nAu0,nAu,t,iAu) nAu0*Fin/V-nAu*Fout/V-(1/F)*iAu*Acat;
dnPddt = @(nPd0,nPd,t,iPd) nPd0*Fin/V-nPd*Fout/V-(1/(2*F))*iPd*Acat;

%modelling method vars
t(1) = 0; %time initial
h = 1/60; %step size, basically does a minute of time
tfinal = 72;

for iter = 1:1:(tfinal/h)
    %solving for I at current step
    loopcheck = 0;
    numCurrentIter(iter) = 1;
    while loopcheck == 0
       etaAn(iter) = Rgas*T*log((I(iter)/Aan)/iAn0)/(4*F); %solving for anodic overpotential from estimated current
       anodefunc = @(etaAn)Aan*iAn0*(exp(alphaAn*4*F*etaAn/(Rgas*T))-exp(-(1-alphaAn)*4*F*etaAn/(Rgas*T)))-I(iter);
       
       etaAn(iter) = fzero(anodefunc,etaAn(iter));
       etaAg(iter) = abs(Vapp-abs(ErevAg(iter)-Ean(iter)) - IRhardware - I(iter)*R(iter) - abs(etaAn(iter)));
       iAg(iter) = iAg0*(exp(alphaAg*1*F*(etaAg(iter))/(Rgas*T))-exp(-(1-alphaAg)*1*F*etaAg(iter)/(Rgas*T)));

       etaAu(iter) = abs(Vapp-abs(ErevAu(iter)-Ean(iter)) - IRhardware - I(iter)*R(iter) - abs(etaAn(iter)));
       iAu(iter) = iAu0*(exp(alphaAu*1*F*(etaAu(iter))/(Rgas*T))-exp(-(1-alphaAu)*1*F*etaAu(iter)/(Rgas*T)));

       etaPd(iter) = abs(Vapp-abs(ErevPd(iter)-Ean(iter)) - IRhardware - I(iter)*R(iter) - abs(etaAn(iter)));
       iPd(iter) = iPd0*(exp(alphaPd*2*F*(etaPd(iter))/(Rgas*T))-exp(-(1-alphaPd)*2*F*etaPd(iter)/(Rgas*T)));

       Icathode(iter) = (iAu(iter)+iAg(iter)+iPd(iter))*Acat;

       error = Icathode(iter)-I(iter);
       if isinf(error) | isnan(error)
           disp("ENDED PROGRAM: Currents do not converge, change how currents are chosen or make Vapp higher")
           return
       end
       %pause%use to debug solver
       if abs(error) < 0.05
            loopcheck = 1;
       else
           I(iter) = I(iter)+0.5*error;
           numCurrentIter(iter) = numCurrentIter(iter)+1;
           loopcheck = 0;
           if numCurrentIter(iter) > 10000
               disp("ENDED PROGRAM: Too many iterations")
               return
           end
        end
    end
    t(iter+1) = h*iter;
    %Ag
    k1 = dnAgdt(nAg0(iter),nAg(iter),t(iter),iAg(iter));
    k2 = dnAgdt((nAg0(iter)+(1/2)*h*k1),(nAg(iter)+(1/2)*h*k1),(t(iter)+(1/2)*h),(iAg(iter)+(1/2)*h*k1));
    k3 = dnAgdt((nAg0(iter)+(1/2)*h*k2),(nAg(iter)+(1/2)*h*k2),(t(iter)+(1/2)*h),(iAg(iter)+(1/2)*h*k2));
    k4 = dnAgdt((nAg0(iter)+h*k3),(nAg(iter)+h*k3),(t(iter)+h),(iAg(iter)+h*k3));
    nAg(iter+1) = subplus(nAg(iter) + (1/6)*h*(k1+2*k2+2*k3+k4))+eps;
    nAg0(iter+1) = nAg0(iter);%nAg(iter+1);%for continuous recycle, use the commented out section. May add functionality for further coupled interactions.
    %Au
    k1 = dnAudt(nAu0(iter),nAu(iter),t(iter),iAu(iter));
    k2 = dnAudt((nAu0(iter)+(1/2)*h*k1),(nAu(iter)+(1/2)*h*k1),(t(iter)+(1/2)*h),(iAu(iter)+(1/2)*h*k1));
    k3 = dnAudt((nAu0(iter)+(1/2)*h*k2),(nAu(iter)+(1/2)*h*k2),(t(iter)+(1/2)*h),(iAu(iter)+(1/2)*h*k2));
    k4 = dnAudt((nAu0(iter)+h*k3),(nAu(iter)+h*k3),(t(iter)+h),(iAu(iter)+h*k3));
    nAu(iter+1) = subplus(nAu(iter) + (1/6)*h*(k1+2*k2+2*k3+k4))+eps;
    nAu0(iter+1) = nAu0(iter);%nAu(iter+1);
    %Pd
    k1 = dnPddt(nPd0(iter),nPd(iter),t(iter),iPd(iter));
    k2 = dnPddt((nPd0(iter)+(1/2)*h*k1),(nPd(iter)+(1/2)*h*k1),(t(iter)+(1/2)*h),(iPd(iter)+(1/2)*h*k1));
    k3 = dnPddt((nPd0(iter)+(1/2)*h*k2),(nPd(iter)+(1/2)*h*k2),(t(iter)+(1/2)*h),(iPd(iter)+(1/2)*h*k2));
    k4 = dnPddt((nPd0(iter)+h*k3),(nPd(iter)+h*k3),(t(iter)+h),(iPd(iter)+h*k3));
    nPd(iter+1) = subplus(nPd(iter) + (1/6)*h*(k1+2*k2+2*k3+k4))+eps;
    nPd0(iter+1) = nPd0(iter);%nPd(iter+1);
    %H
    nH(iter+1) = nH(iter); % refine this, im assuming something constantly balances pH here (or the conc is so big it dont matta)
    %other variables to calc
    K(iter+1) = gammaAg*nAg(iter+1)+gammaAu*nAu(iter+1)+gammaPd*2*nPd(iter+1)+gammaH*nH(iter+1);
    R(iter+1) = cursivel/(Acat*K(iter+1));
    Ean(iter+1) = -1.23 - (Rgas*T/(4*F))*log(aO2*(gamH*nH(iter+1)/V)^4);
    ErevAg(iter+1) = 0.7989 + (Rgas*T/(F))*log(1/(gamAg*nAg(iter+1)/V));
    ErevAu(iter+1) = 1.692 + (Rgas*T/(F))*log(1/(gamAg*nAu(iter+1)/V));
    ErevPd(iter+1) = 0.951 + (Rgas*T/(2*F))*log(1/(((gamAg*nPd(iter+1))/V)^2));
    I(iter+1) = I(iter); %this is the initial guess for the next step, to be looped on
    %additional processing equations
    wAg(iter+1) = (Acat*mmAg/F)*h*60*(sum(iAg)-iAg(end));
    wAu(iter+1) = (Acat*mmAu/F)*h*60*(sum(iAu)-iAu(end));
    wPd(iter+1) = (Acat*mmPd/(2*F))*h*60*(sum(iPd)-iPd(end));
end
CAg = nAg./V;
CAu = nAu./V;
CPd = nPd./V;
RemAg = CAg./CAg(1);
RemAu = CAu./CAu(1);
RemPd = CPd./CPd(1);
plot(t,RemAg,t,RemAu,t,RemPd)
xlabel('Time (h)')
ylabel('% consumption from initial')
title('Remaining in solution over time, step size 1 min')
legend('Silver','Gold','Palladium')
figure
plot(t,wAg,t,wAu,t,wPd)
xlabel('Time (h)')
ylabel('Total amount deposited (g)')
title('Deposited solution over time, step size 1 min')
legend('Silver','Gold','Palladium')
toc
