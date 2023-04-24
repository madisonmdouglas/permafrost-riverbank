function [y, yblock, Ebank, Eblock, Eent, Ethaw, Tbank, Tblock, Tthaw, Tent, Fbank, Fthaw] = ...
    RunPfModel_Vslump(H, U, Tw, S, zmax, dz, tmax, dt, f_ice, rhoB, Ti, Tau_crit, M, sigma_T, sigma_C, sigma_S)
%% run timeseries to calculate the rate of permafrost riverbank erosion by thaw- and entrainment-limited conditions and considering slump blocks
% inputs:
%   H = water depth (m)
%   U = average water velocity (m/s)
%   Tw = water temperature (degC)
%   S = channel slope (m/m)
%   zmax = bank height (m)
%   dz = vertical resolution for bank profile (m)
%   tmax = number of timesteps (days)
%   dt = timestep (days)
%   f_ice = mass fraction of riverbank comprised of ice (kg/kg)
%   rhoB = permafrost bank bulk density (kg/m3)
%   Ti = initial permafrost temperature (degC)
%   Tau_crit = critical Shields stress to entrain sediment (dimensionless)
%   M = coefficient for entrainment equation (kg/m2/s)
%   sigma_T = bank tensile strength (Pa)
%   sigma_C = bank compressive strength (Pa)
%   sigma_S = bank shear strength (Pa)
% outputs:
%   y = matrix of updated bank ycoordinates (m) [zmax/dz x tmax/dt] 
%   y_block = vector of slump block area (m2) [1 x tmax/dt]
%   Ebank = bank erosion rate (m/day)
%   Eblock = slump block erosion rate (m/day)
%   Eent = entrainment-limited erosion rate (m/day)
%   Ethaw = thaw-limited erosion rate (m/day)
%   Tbank = fraction of time spent eroding bank (days/days)
%   Tblock = fraction of time spent eroding slump blocks (days/days)
%   Tthaw = fraction of time with thaw-limited erosion (days/days)
%   Tent = fraction of time with entrainment-limited erosion (days/days)
%   Fbank = fraction of erosion on the bank versus on blocks (m/m)
%   Fthaw = fraction of erosion that was thaw-limited (m/m)

f_sed = 1 - f_ice;                      % mass fraction of sediment in bank (kg/kg)
if f_sed>0.8
    f_sed=0.8;                          % equivalent to 40% volumetric porosity
end
Lf = 334000 * f_ice;                    % permafrost latent heat of fusion (J/kg)
cp = f_ice*2108 + f_sed*700;            % specific heat of permafrost (J/kg/K)

y = zeros(zmax/dz+1, tmax/dt);
yblock = zeros(1, tmax/dt);
Ebank = zeros(zmax/dz+1, tmax/dt);
Eblock = zeros(1, tmax/dt);
Eent = zeros(1, tmax/dt);
Ethaw = zeros(1, tmax/dt);
thawTF = zeros(1, tmax/dt-1);
bankTF = zeros(1, tmax/dt-1);
for ts = 2:dt*tmax
    [y(:,ts), yblock(ts), Ebank(:,ts), Eblock(ts), Eent(ts), Ethaw(ts), thawTF(ts-1), bankTF(ts-1)] = ...
        RunTimestep_Vslump(y(:,ts-1),yblock(ts-1),zmax,dz,Tw(ts),U(ts),H(ts),S,dt,f_sed,rhoB,Lf,cp,Ti,Tau_crit,M,sigma_T,sigma_C,sigma_S);
end

% calculate fraction of time thaw- vs entrainment-limited
Tthaw = sum(thawTF)/length(thawTF);
Tent = 1 - Tthaw;

% calculate bank fraction of erosion thaw- vs entrainment-limited
Fthaw = sum(thawTF.*max(Ebank(:,2:end)))/sum(max(Ebank(:,2:end)));

% calculate fraction of time spent eroding bank vs slump blocks
Tbank = sum(bankTF)/length(bankTF);
Tblock = 1 - Tbank;

% calculate bank fraction of erosion on the bank versus slump blocks
Fbank = sum(max(Ebank))/(sum(max(Ebank)) + sum(Eblock));

end